//! Main BLASTN run function
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process.

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::{mpsc::channel, Arc};

use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
use crate::stats::lookup_nucl_params;
use crate::core::blast_encoding::{
    encode_iupac_to_blastna,
    COMPRESSION_RATIO,
};

// Import from parent blastn module (super::super:: because we're in blast_engine/)
use super::super::args::BlastnArgs;
use super::super::alignment::{
    build_blastna_matrix,
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_traceback_with_scratch,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    greedy_gapped_alignment_score_only,
    greedy_gapped_alignment_with_traceback,
    GreedyAlignScratch,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    GapAlignScratch,
};
use super::super::extension::{
    build_nucl_score_table,
    extend_hit_ungapped_approx_ncbi,
    extend_hit_ungapped_exact_ncbi,
    type_of_word,
};
use crate::utils::dust::MaskedInterval;
use super::super::constants::TWO_HIT_WINDOW;
use super::super::coordination::{configure_task, finalize_task_config, read_sequences, prepare_sequence_data, build_lookup_tables};
use super::super::lookup::{build_unmasked_ranges, reverse_complement};
use super::super::ncbi_cutoffs::{compute_blastn_cutoff_score, GAP_TRIGGER_BIT_SCORE_NUCL};
use super::super::blast_extend::DiagStruct;
use super::super::interval_tree::{BlastIntervalTree, TreeHsp, IndexMethod};
use super::super::filtering::{
    blast_hsp_test_identity_and_length,
    hsp_test,
    purge_hsps_with_common_endpoints,
    purge_hsps_with_common_endpoints_ex,
    reevaluate_hsp_with_ambiguities_gapped_ex,
    subject_best_hit,
    ReevalParams,
};
use super::super::hsp::{
    BlastnHitList,
    BlastnHsp,
    BlastnHspList,
    compare_hsp_lists as compare_blastn_hsp_lists,
    prune_hitlists_by_size,
    score_compare_hsps as score_compare_blastn_hsps,
    trim_by_max_hsps,
    update_best_evalue,
    write_output_blastn_hitlists,
};

// Import from this module (blast_engine)
use super::calculate_evalue;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/ncbi_math.h:160-161
// ```c
// #define NCBIMATH_LN2 0.69314718055994530941723212145818
// ```
const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:123-143
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

#[derive(Clone)]
struct SubjectChunk {
    offset: usize,
    length: usize,
    residual: usize,
    seq_ranges: Vec<(i32, i32)>,
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

        let no_chunking =
            offset == 0 && residual == 0 && self.next == self.full_right;
        let seq_ranges = if no_chunking {
            self.soft_ranges.clone()
        } else if !subject_masked {
            vec![(residual as i32, length as i32)]
        } else {
            debug_assert_eq!(residual, 0);
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
            let mut ranges = Vec::with_capacity(count);
            for j in 0..count {
                let left = self.soft_ranges[start_idx + j].0 - offset;
                let right = self.soft_ranges[start_idx + j].1 - offset;
                ranges.push((left, right));
            }
            if let Some(first) = ranges.first_mut() {
                if first.0 < 0 {
                    first.0 = 0;
                }
            }
            if let Some(last) = ranges.last_mut() {
                if last.1 > length as i32 {
                    last.1 = length as i32;
                }
            }
            ranges
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
            hit_len_array.fill(0);
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
        let bucket = ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
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
        let bucket = ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
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
            self.chain.resize(self.capacity as usize, DiagHashCell::default());
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
    let (top_shift, bottom_shift) = if (start % 2) == 0 { (2u32, 6u32) } else { (0u32, 4u32) };
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1647-1662
// ```c
// scan_range[1] = 0;  /* start pos of scan */
// scan_range[2] = subject->length - lut_word_length; /* end pos */
// if (subject->mask_type != eNoSubjMasking) {
//     scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
//     scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
// }
// ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1651-1667
    // ```c
    // if (subject->mask_type != eNoSubjMasking) {
    //     scansub = (TNaScanSubjectFunction)
    //           BlastChooseNucleotideScanSubjectAny(lookup_wrap);
    //     ...
    //     scan_range[1] = subject->seq_ranges[0].left + word_length - lut_word_length;
    //     scan_range[2] = subject->seq_ranges[0].right - lut_word_length;
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2631-2677
    // ```c
    // switch (mb_lt->lut_word_length) {
    // case 9:
    //     if (scan_step == 1)
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_9_1;
    //     if (scan_step == 2)
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_9_2;
    //     else
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
    //     break;
    // case 10:
    //     if (scan_step == 1)
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_10_1;
    //     else if (scan_step == 2)
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_10_2;
    //     else if (scan_step == 3)
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_10_3;
    //     else
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
    //     break;
    // case 11:
    //     switch (scan_step % COMPRESSION_RATIO) {
    //     case 0:
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
    //         break;
    //     case 1:
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_11_1Mod4;
    //         break;
    //     case 2:
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_11_2Mod4;
    //         break;
    //     case 3:
    //         mb_lt->scansub_callback = (void *)s_MBScanSubject_11_3Mod4;
    //         break;
    //     }
    //     break;
    // case 12:
    // case 16:
    //     mb_lt->scansub_callback = (void *)s_MBScanSubject_Any;
    //     break;
    // }
    // ```
    if lut_word_length >= 9 && lut_word_length <= 12 {
        if subject_masked || start % COMPRESSION_RATIO != 0 {
            scan_subject_kmers_range_mb_any(
                packed,
                subject_len,
                lut_word_length,
                scan_step,
                subject_masked,
                start,
                end,
                on_kmer,
            );
            return;
        }

        match lut_word_length {
            9 => {
                if scan_step == 1 {
                    scan_subject_kmers_range_mb_9_1(
                        packed,
                        subject_len,
                        scan_step,
                        start,
                        end,
                        on_kmer,
                    );
                } else if scan_step == 2 {
                    scan_subject_kmers_range_mb_9_2(
                        packed,
                        subject_len,
                        scan_step,
                        start,
                        end,
                        on_kmer,
                    );
                } else {
                    scan_subject_kmers_range_mb_any(
                        packed,
                        subject_len,
                        lut_word_length,
                        scan_step,
                        subject_masked,
                        start,
                        end,
                        on_kmer,
                    );
                }
            }
            10 => {
                if scan_step == 1 {
                    scan_subject_kmers_range_mb_10_1(
                        packed,
                        subject_len,
                        scan_step,
                        start,
                        end,
                        on_kmer,
                    );
                } else if scan_step == 2 {
                    scan_subject_kmers_range_mb_10_2(
                        packed,
                        subject_len,
                        scan_step,
                        start,
                        end,
                        on_kmer,
                    );
                } else if scan_step == 3 {
                    scan_subject_kmers_range_mb_10_3(
                        packed,
                        subject_len,
                        scan_step,
                        start,
                        end,
                        on_kmer,
                    );
                } else {
                    scan_subject_kmers_range_mb_any(
                        packed,
                        subject_len,
                        lut_word_length,
                        scan_step,
                        subject_masked,
                        start,
                        end,
                        on_kmer,
                    );
                }
            }
            11 => match scan_step % COMPRESSION_RATIO {
                0 => scan_subject_kmers_range_mb_any(
                    packed,
                    subject_len,
                    lut_word_length,
                    scan_step,
                    subject_masked,
                    start,
                    end,
                    on_kmer,
                ),
                1 => scan_subject_kmers_range_mb_11_1mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    on_kmer,
                ),
                2 => scan_subject_kmers_range_mb_11_2mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    on_kmer,
                ),
                _ => scan_subject_kmers_range_mb_11_3mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    on_kmer,
                ),
            },
            12 => scan_subject_kmers_range_mb_any(
                packed,
                subject_len,
                lut_word_length,
                scan_step,
                subject_masked,
                start,
                end,
                on_kmer,
            ),
            _ => {}
        }
        return;
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
    mut on_kmer: F,
) where
    F: FnMut(usize, u64),
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

    while determine_scanning_offsets(
        seq_ranges,
        word_length as i32,
        lut_word_length as i32,
        &mut scan_range,
    ) {
        let start = scan_range[1];
        let end = scan_range[2];
        if start >= 0 && end >= 0 {
            scan_subject_kmers_range(
                packed,
                subject_len,
                lut_word_length,
                scan_step,
                subject_masked,
                start as usize,
                end as usize,
                &mut on_kmer,
            );
        }
        scan_range[1] = scan_range[2] + 1;
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
    qs: usize,        // query start
    qe: usize,        // query end (exclusive)
    ss: usize,        // subject start in search_seq coordinates
    se: usize,        // subject end in search_seq coordinates
    score: i32,       // ungapped score
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
}

/// Internal 0-based coordinates for Phase 2 interval tree processing
/// NCBI reference: blast_traceback.c:679-692 - Phase 2 uses internal coordinates
#[derive(Clone)]
struct InternalHitData {
    q_offset_0: usize,  // 0-based query start
    q_end_0: usize,     // 0-based query end
    s_offset_0: usize,  // 0-based subject start (canonical, always < s_end_0)
    s_end_0: usize,     // 0-based subject end (canonical, always > s_offset_0)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3910-3913
    // ```c
    // tmp_hsp.query.offset = q_start;
    // tmp_hsp.query.end = q_end;
    // tmp_hsp.query.frame = query_info->contexts[context].frame;
    // ```
    query_frame: i32,
    query_context_offset: i32,
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
// ```c
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// ```
fn purge_prelim_hits_with_common_endpoints(mut hits: Vec<PrelimHit>) -> Vec<PrelimHit> {
    if hits.len() <= 1 {
        return hits;
    }

    // Pass 1: common starts
    hits.sort_by(|a, b| {
        a.context_idx
            .cmp(&b.context_idx)
            .then_with(|| a.prelim_qs.cmp(&b.prelim_qs))
            .then_with(|| a.prelim_ss.cmp(&b.prelim_ss))
            .then_with(|| b.prelim_score.cmp(&a.prelim_score))
            .then_with(|| b.prelim_qe.cmp(&a.prelim_qe))
            .then_with(|| b.prelim_se.cmp(&a.prelim_se))
    });

    let mut i = 0usize;
    while i + 1 < hits.len() {
        let j = i + 1;
        let same_context = hits[i].context_idx == hits[j].context_idx;
        let same_q_start = hits[i].prelim_qs == hits[j].prelim_qs;
        let same_s_start = hits[i].prelim_ss == hits[j].prelim_ss;
        if same_context && same_q_start && same_s_start {
            hits.remove(j);
        } else {
            i += 1;
        }
    }

    // Pass 2: common ends
    hits.sort_by(|a, b| {
        a.context_idx
            .cmp(&b.context_idx)
            .then_with(|| a.prelim_qe.cmp(&b.prelim_qe))
            .then_with(|| a.prelim_se.cmp(&b.prelim_se))
            .then_with(|| b.prelim_score.cmp(&a.prelim_score))
            .then_with(|| b.prelim_qs.cmp(&a.prelim_qs))
            .then_with(|| b.prelim_ss.cmp(&a.prelim_ss))
    });

    let mut i = 0usize;
    while i + 1 < hits.len() {
        let j = i + 1;
        let same_context = hits[i].context_idx == hits[j].context_idx;
        let same_q_end = hits[i].prelim_qe == hits[j].prelim_qe;
        let same_s_end = hits[i].prelim_se == hits[j].prelim_se;
        if same_context && same_q_end && same_s_end {
            hits.remove(j);
        } else {
            i += 1;
        }
    }

    hits
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
        && (hsp1.prelim_ss as isize - hsp2.prelim_ss as isize
            - hsp1.prelim_qs as isize
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
        let score_density = (hsp1.prelim_score as f64 + hsp2.prelim_score as f64)
            / (len1 + len2);

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2860-3055
// ```c
// Int2 Blast_HSPListsMerge(..., Int4 contexts_per_query, Int4 *split_offsets,
//                          Int4 chunk_overlap_size, Boolean allow_gap, Boolean short_reads)
// ```
fn merge_prelim_hits_subject_split(
    combined: &mut Vec<PrelimHit>,
    mut incoming: Vec<PrelimHit>,
    split_offset: usize,
    chunk_overlap_size: usize,
    allow_gap: bool,
) {
    if incoming.is_empty() {
        return;
    }
    if combined.is_empty() {
        *combined = incoming;
        return;
    }

    let mut combined_overlap = Vec::new();
    for (idx, hsp) in combined.iter().enumerate() {
        if hsp.prelim_se > split_offset {
            combined_overlap.push(idx);
        }
    }

    let mut incoming_overlap = Vec::new();
    for (idx, hsp) in incoming.iter().enumerate() {
        if hsp.prelim_ss < split_offset.saturating_add(chunk_overlap_size) {
            incoming_overlap.push(idx);
        }
    }

    let mut deleted = vec![false; incoming.len()];

    for &i in combined_overlap.iter() {
        let hsp1_context = combined[i].context_idx;
        let end_diag = combined[i].prelim_qe as i32 - combined[i].prelim_se as i32;
        for &j in incoming_overlap.iter() {
            if deleted[j] {
                continue;
            }
            if incoming[j].context_idx != hsp1_context {
                continue;
            }
            let start_diag = incoming[j].prelim_qs as i32 - incoming[j].prelim_ss as i32;
            if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE {
                if merge_two_prelim_hits(&mut combined[i], &incoming[j], allow_gap) {
                    deleted[j] = true;
                }
            }
        }
    }

    incoming = incoming
        .into_iter()
        .enumerate()
        .filter_map(|(idx, hsp)| if deleted[idx] { None } else { Some(hsp) })
        .collect();

    combined.extend(incoming);
    combined.sort_by(score_compare_prelim_hits);
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
    hits_with_internal: Vec<(BlastnHsp, InternalHitData)>,
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
    fn new(query_count: usize) -> Self {
        Self {
            hits_with_internal: Vec::new(),
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
            hit_level_array: Vec::new(),
            hit_len_array: Vec::new(),
            diag_hash: DiagHashTable::new(TWO_HIT_WINDOW as i32),
            diag_table_offset: TWO_HIT_WINDOW as i32,
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
            offset_pairs: Vec::new(),
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
    let mut approx_entries = 0usize;
    let mut max_q_off = 0usize;

    for (ctx, masks) in contexts.iter().zip(context_masks.iter()) {
        let ranges = build_unmasked_ranges(ctx.seq.len(), masks);
        let ctx_offset = ctx.query_offset.max(0) as usize;
        for (start, end) in ranges {
            if end <= start {
                continue;
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_util.c:190-203
            // ```c
            // num_entries += loc->ssr->right - loc->ssr->left;
            // curr_max = MAX(curr_max, loc->ssr->right);
            // ```
            approx_entries = approx_entries.saturating_add(end.saturating_sub(start));
            let abs_right = ctx_offset.saturating_add(end);
            if abs_right > max_q_off {
                max_q_off = abs_right;
            }
        }
    }

    (approx_entries, max_q_off)
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
fn build_subject_seq_ranges_from_masks(
    masks: &[MaskedInterval],
    subject_len: usize,
) -> Vec<(i32, i32)> {
    if subject_len == 0 {
        return Vec::new();
    }
    if masks.is_empty() {
        return vec![(0, subject_len as i32)];
    }

    let mut sorted = masks.to_vec();
    sorted.sort_by_key(|m| m.start);

    let mut merged: Vec<MaskedInterval> = Vec::new();
    for mask in sorted {
        let start = mask.start.min(subject_len);
        let end = mask.end.min(subject_len);
        if start >= end {
            continue;
        }
        match merged.last_mut() {
            Some(last) if start <= last.end => {
                last.end = last.end.max(end);
            }
            _ => merged.push(MaskedInterval::new(start, end)),
        }
    }

    let mut ranges: Vec<(i32, i32)> = Vec::new();
    let mut cursor = 0usize;
    for mask in merged {
        if mask.start > cursor {
            ranges.push((cursor as i32, mask.start as i32));
        }
        cursor = cursor.max(mask.end);
    }
    if cursor < subject_len {
        ranges.push((cursor as i32, subject_len as i32));
    }

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
fn post_process_hits_and_write(
    all_hits: Vec<BlastnHsp>,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    out_path: &Option<std::path::PathBuf>,
    verbose: bool,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    if verbose {
        eprintln!(
            "[INFO] Received {} raw hits total, starting post-processing...",
            all_hits.len()
        );
    }
    let chain_start = std::time::Instant::now();

    // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
    // Per-subject traceback already applies purge/reevaluation; skip duplicate global pass.
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-180
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;
    //    Int4 query_index;
    //    BlastHSP** hsp_array;
    //    Int4 hspcnt;
    //    double best_evalue;
    // } BlastHSPList;
    // typedef struct BlastHitList {
    //    Int4 hsplist_count;
    //    BlastHSPList** hsplist_array;
    // } BlastHitList;
    // ```
    let mut per_query: Vec<FxHashMap<u32, Vec<BlastnHsp>>> =
        vec![FxHashMap::default(); query_ids.len()];
    for hsp in all_hits {
        let q_idx = hsp.q_idx as usize;
        per_query[q_idx].entry(hsp.s_idx).or_default().push(hsp);
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:183-187
    // ```c
    // typedef struct BlastHSPResults {
    //    Int4 num_queries;
    //    BlastHitList** hitlist_array;
    // } BlastHSPResults;
    // ```
    let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(query_ids.len());
    hit_lists.resize_with(query_ids.len(), || None);
    for (q_idx, subject_map) in per_query.into_iter().enumerate() {
        if subject_map.is_empty() {
            continue;
        }

        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;
        //    Int4 query_index;
        //    BlastHSP** hsp_array;
        //    Int4 hspcnt;
        // } BlastHSPList;
        // ```
        let mut hsplist_array: Vec<BlastnHspList> = Vec::with_capacity(subject_map.len());
        for (s_idx, mut hsps) in subject_map {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
            // ```c
            // int ScoreCompareHSPs(const void* h1, const void* h2) { ... }
            // ```
            hsps.sort_by(score_compare_blastn_hsps);

            let mut hsp_list = BlastnHspList {
                oid: s_idx,
                query_index: q_idx as u32,
                hsps,
                best_evalue: f64::MAX,
            };

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
            // ```c
            // Int2 Blast_TrimHSPListByMaxHsps(...)
            // ```
            if max_hsps_per_subject > 0 {
                trim_by_max_hsps(&mut hsp_list, max_hsps_per_subject);
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1748
            // ```c
            // static double s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
            // ```
            update_best_evalue(&mut hsp_list);

            hsplist_array.push(hsp_list);
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3071-3106
        // ```c
        // static int s_EvalueCompareHSPLists(const void* v1, const void* v2) { ... }
        // ```
        hsplist_array.sort_by(compare_blastn_hsp_lists);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
        // ```c
        // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
        // ```
        prune_hitlists_by_size(&mut hsplist_array, hitlist_size);

        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:168-180
        // ```c
        // typedef struct BlastHitList {
        //    Int4 hsplist_count;
        //    BlastHSPList** hsplist_array;
        // } BlastHitList;
        // ```
        hit_lists[q_idx] = Some(BlastnHitList { hsplist_array });
    }

    let total_hits: usize = hit_lists
        .iter()
        .filter_map(|list| list.as_ref())
        .map(|list| list.hsplist_array.iter().map(|l| l.hsps.len()).sum::<usize>())
        .sum();

    if verbose {
        eprintln!(
            "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
            chain_start.elapsed().as_secs_f64(),
            total_hits
        );
    }
    let write_start = std::time::Instant::now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3071-3106
    // ```c
    // static int s_EvalueCompareHSPLists(const void* v1, const void* v2) { ... }
    // ```
    write_output_blastn_hitlists(&hit_lists, out_path.as_ref(), query_ids, subject_ids)?;

    if verbose {
        eprintln!(
            "[INFO] Output written in {:.2}s",
            write_start.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

pub fn run(args: BlastnArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
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
    let use_parallel = num_threads > 1;

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
    if use_parallel {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .context("Failed to build thread pool")?;
    }

    // Configure task-specific parameters (initial configuration)
    let mut config = configure_task(&args);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1415-1426
    // ```c
    // if (options->db_filter && options->word_size < 16) {
    //    Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //       "The limit_lookup option can only be used with word size >= 16");
    //    return BLASTERR_OPTION_VALUE_INVALID;
    // }
    // ```
    if args.limit_lookup && config.effective_word_size < 16 {
        anyhow::bail!("The limit_lookup option can only be used with word size >= 16");
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:1495-1501
    // ```c
    // arg_desc.AddDefaultKey(kArgMaxDbWordCount, ...);
    // arg_desc.SetConstraint(kArgMaxDbWordCount,
    //                        new CArgAllowValuesBetween(2, 255, true));
    // ```
    if args.limit_lookup && args.max_db_word_count < 2 {
        anyhow::bail!("The max_db_word_count option must be >= 2");
    }

    // Read sequences
    let (queries, query_ids, subjects) = read_sequences(&args)?;
    if queries.is_empty() || subjects.is_empty() {
        return Ok(());
    }

    // Prepare sequence data (including DUST masking)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1401-1405
    // ```c
    // /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
    //   sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
    // seq_arg.encoding = eBlastEncodingProtein;
    // ```
    let seq_data = prepare_sequence_data(&args, queries, query_ids, subjects);

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
            .subjects
            .iter()
            .map(|record| {
                Arc::<str>::from(
                    record
                        .id()
                        .split_whitespace()
                        .next()
                        .unwrap_or("unknown"),
                )
            })
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
    let mut query_context_records: Vec<bio::io::fasta::Record> = Vec::new();
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

        query_context_records.push(bio::io::fasta::Record::with_attrs(
            q_record.id(),
            q_record.desc(),
            seq,
        ));
        query_context_records.push(bio::io::fasta::Record::with_attrs(
            q_record.id(),
            q_record.desc(),
            rc_seq.as_slice(),
        ));
        query_context_masks.push(plus_masks);
        query_context_masks.push(minus_masks);
        query_concat_offset += q_len * 2 + 1;
    }
    let query_concat_length = query_concat_offset;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // ```c
    // query_info->contexts[0].query_offset = 0;
    // query_info->contexts[1].query_offset = kStrandLength + 1;
    // ```
    let query_context_offsets: Vec<i32> = query_contexts
        .iter()
        .map(|ctx| ctx.query_offset)
        .collect();
    let query_context_index = QueryContextIndex::new(&query_contexts);

    // Finalize configuration with query-dependent parameters (adaptive lookup table selection)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:46-47
    // ```c
    // BlastChooseNaLookupTable(const LookupTableOptions* lookup_options,
    //                          Int4 approx_table_entries, Int4 max_q_off,
    //                          Int4 *lut_width)
    // ```
    let discontig_template = args.task == "dc-megablast";
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    // ```c
    // diag_array_length = 1;
    // while (diag_array_length < (qlen+window_size))
    //     diag_array_length = diag_array_length << 1;
    // diag_table->diag_array_length = diag_array_length;
    // diag_table->diag_mask = diag_array_length-1;
    // ```
    let diag_array_lengths: Vec<usize> = seq_data
        .queries
        .iter()
        .map(|q| {
            let qlen = q.seq().len();
            let mut diag_array_length = 1usize;
            while diag_array_length < qlen + TWO_HIT_WINDOW {
                diag_array_length <<= 1;
            }
            diag_array_length
        })
        .collect();
    let diag_masks: Vec<usize> = diag_array_lengths
        .iter()
        .map(|len| len.saturating_sub(1))
        .collect();

    // Get Karlin-Altschul parameters
    // CRITICAL: NCBI uses DIFFERENT params for ungapped and gapped calculations!
    // - Ungapped params (kbp_std): Used for gap_trigger calculation
    // - Gapped params (kbp_gap): Used for cutoff_score_max and length adjustment
    // Reference: blast_parameters.c:343-344 uses kbp_std for gap_trigger
    use crate::config::NuclScoringSpec;

    // Gapped params (for length adjustment and cutoff_score_max)
    let scoring_spec_gapped = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    let params_gapped = lookup_nucl_params(&scoring_spec_gapped);

    // Ungapped params (for gap_trigger calculation)
    // NCBI uses gap_open=0, gap_extend=0 to get ungapped params
    let scoring_spec_ungapped = NuclScoringSpec {
        reward: config.reward,
        penalty: config.penalty,
        gap_open: 0,
        gap_extend: 0,
    };
    let params_ungapped = lookup_nucl_params(&scoring_spec_ungapped);

    // Use gapped params for E-value calculations (same as before)
    let params = params_gapped.clone();

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
    let query_eff_searchsp: Vec<i64> = query_contexts
        .iter()
        .map(|ctx| {
            let query_len = ctx.seq.len() as i64;
            let result = compute_length_adjustment_ncbi(
                query_len,
                db_len_total_i64,
                db_num_seqs_i64,
                &params_gapped,
            );
            let length_adjustment = result.length_adjustment;
            let effective_db_length =
                (db_len_total_i64 - db_num_seqs_i64 * length_adjustment).max(1);
            let effective_query_length = (query_len - length_adjustment).max(1);
            effective_db_length * effective_query_length
        })
        .collect();

    // Build lookup tables
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // ```c
    // pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
    // mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
    // ```
    let (lookup_tables, scan_step) = build_lookup_tables(
        &config,
        &args,
        &query_context_records,
        &query_context_masks,
        &query_context_offsets,
        &seq_data.subjects,
        approx_table_entries,
    );

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
    let progress_bar = if use_parallel {
        let bar = ProgressBar::new(seq_data.subjects.len() as u64);
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
    let max_hsps_per_subject = args.max_hsps_per_subject;

    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    // BLASTN-specific debug mode: set LOSAT_DEBUG_BLASTN=1 to enable detailed hit loss diagnostics
    let blastn_debug = std::env::var("LOSAT_DEBUG_BLASTN").is_ok();
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
    let hash_lookup_ref = lookup_tables.hash_lookup.as_ref();
    let queries_ref = &seq_data.queries;
    let query_contexts_ref = &query_contexts;
    let query_base_offsets_ref = &query_base_offsets;
    let query_eff_searchsp_ref = &query_eff_searchsp;
    let diag_array_lengths_ref = &diag_array_lengths;
    let diag_masks_ref = &diag_masks;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:812-826
    // ```c
    // if (subjects.GetMask(i).NotEmpty()) {
    //     s_SeqLoc2MaskedSubjRanges(...);
    //     BlastSeqBlkSetSeqRanges(..., eSoftSubjMasking);
    // }
    // ```
    let subject_masks_ref = &seq_data.subject_masks;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // ```c
    // BlastSeqBlkSetSequence(subj, sequence.data.release(), ...);
    // ...
    // BlastSeqBlkSetCompressedSequence(subj,
    //                                  compressed_seq.data.release());
    // ```
    let subject_encodings_ref = &seq_data.subject_encodings;
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

    // Capture config values for use in closure
    let effective_word_size = config.effective_word_size;
    let min_ungapped_score = config.min_ungapped_score;
    let use_dp = config.use_dp;
    let use_direct_lookup = config.use_direct_lookup;
    let reward = config.reward;
    let penalty = config.penalty;
    let gap_open = config.gap_open;
    let gap_extend = config.gap_extend;
    let params_for_closure = params.clone();  // Gapped params for E-value
    let params_ungapped_for_closure = params_ungapped.clone();  // Ungapped params for gap_trigger
    let params_gapped_for_closure = params_gapped.clone();  // Gapped params for cutoff_score_max
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:454-468
    // ```c
    // double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
    // params->gap_x_dropoff = (Int4)(options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    // params->gap_x_dropoff_final = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
    // ```
    let min_lambda = params_gapped_for_closure.lambda;
    let x_drop_gapped = ((config.x_drop_gapped as f64 * NCBIMATH_LN2) / min_lambda) as i32;
    let x_drop_final = ((config.x_drop_final as f64 * NCBIMATH_LN2) / min_lambda)
        .max(x_drop_gapped as f64) as i32;
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
    let x_dropoff_init = ((super::super::constants::X_DROP_UNGAPPED as f64 * NCBIMATH_LN2)
        / params_ungapped_for_closure.lambda)
        .ceil() as i32;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:866-869 (query blastna, subject ncbi2na)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
    let encoded_queries_blastna: Vec<Vec<u8>> = query_contexts_ref
        .iter()
        .map(|ctx| encode_iupac_to_blastna(ctx.seq.as_slice()))
        .collect();

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
    let process_subject = |s_idx: usize,
                           s_record: &bio::io::fasta::Record,
                           gap_scratch: &mut GapAlignScratch,
                           subject_scratch: &mut SubjectScratch,
                           subject_hits: &mut Option<Vec<BlastnHsp>>| {
            let queries = queries_ref;
            let query_contexts = query_contexts_ref;
            let query_base_offsets = query_base_offsets_ref;
            let query_eff_searchsp = query_eff_searchsp_ref;
            let diag_array_lengths = diag_array_lengths_ref;
            let diag_masks = diag_masks_ref;
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
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:812-826
            // ```c
            // if (subjects.GetMask(i).NotEmpty()) {
            //     s_SeqLoc2MaskedSubjRanges(...);
            //     BlastSeqBlkSetSeqRanges(..., eSoftSubjMasking);
            // }
            // ```
            let subject_masks = subject_masks_ref
                .get(s_idx)
                .map(|m| m.as_slice())
                .unwrap_or(&[]);
            let subject_masked = !subject_masks.is_empty();
            let soft_ranges = if subject_masked {
                build_subject_seq_ranges_from_masks(subject_masks, s_len_full)
            } else {
                vec![(0i32, s_len_full as i32)]
            };

            if s_len_full < effective_word_size {
                return;
            }

            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
            // ```c
            // BlastSeqBlkSetSequence(subj, sequence.data.release(), ...);
            // ```
            let subject_encoding = &subject_encodings_ref[s_idx];
            let s_seq_blastna_full = subject_encoding.blastna.as_slice();
            let s_seq_packed_full = subject_encoding.ncbi2na_packed.as_slice();

            let dbseq_chunk_overlap = DBSEQ_CHUNK_OVERLAP;
            let mut split_state = SubjectSplitState::new(s_len_full, soft_ranges);
            let mut subject_chunks: Vec<SubjectChunk> = Vec::new();
            loop {
                match split_state.next_chunk(subject_masked, dbseq_chunk_overlap) {
                    SubjectChunkStatus::Done => break,
                    SubjectChunkStatus::NoRange => continue,
                    SubjectChunkStatus::Ok(chunk) => subject_chunks.push(chunk),
                }
            }
            if subject_chunks.is_empty() {
                return;
            }

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
            let subject_len = s_len_full as i64;
            let mut cutoff_scores: Vec<i32> = Vec::with_capacity(queries.len());
            for q_record in queries.iter() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-843
                // ```c
                // BLAST_ComputeLengthAdjustment(..., query_length, db_length, ...);
                // effective_search_space = effective_db_length *
                //                         (query_length - length_adjustment);
                // ```
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-369
                // ```c
                // if (!gapped_calculation || sbp->matrix_only_scoring) {
                //     Int4 query_length = query_info->contexts[context].query_length;
                //     if (program_number == eBlastTypeBlastn) query_length *= 2;
                // } else {
                //     new_cutoff = gap_trigger;
                // }
                // ```
                let query_len = q_record.seq().len() as i64;
                let cutoff = compute_blastn_cutoff_score(
                    query_len,
                    subject_len,
                    evalue_threshold,
                    GAP_TRIGGER_BIT_SCORE_NUCL,
                    &params_ungapped_for_closure,  // UNGAPPED for gap_trigger (NCBI: kbp_std)
                    &params_gapped_for_closure,    // GAPPED for cutoff_score_max (NCBI: kbp_gap)
                    1.0,
                );
                cutoff_scores.push(cutoff);
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:380-383
            // ```c
            // if (curr_cutoffs->x_dropoff_init == 0)
            //    curr_cutoffs->x_dropoff = new_cutoff;
            // else
            //    curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
            // ```
            let mut x_dropoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
            for cutoff in cutoff_scores.iter() {
                let x_dropoff = if x_dropoff_init == 0 { *cutoff } else { x_dropoff_init };
                x_dropoff_scores.push(x_dropoff);
            }
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:408-412
            // ```c
            // curr_cutoffs->reduced_nucl_cutoff_score = (Int4)(0.8 * new_cutoff);
            // ```
            let mut reduced_cutoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
            for cutoff in cutoff_scores.iter() {
                reduced_cutoff_scores.push((0.8 * (*cutoff as f64)) as i32);
            }

            // BLASTN debug: Log cutoff scores for this query-subject pair
            if blastn_debug {
                eprintln!(
                    "[BLASTN_DEBUG] Subject {} (len={}): cutoff_scores={:?}, params_ungapped=(lambda={:.4}, K={:.4}), params_gapped=(lambda={:.4}, K={:.4})",
                    s_id, s_len_full, cutoff_scores,
                    params_ungapped_for_closure.lambda, params_ungapped_for_closure.k,
                    params_gapped_for_closure.lambda, params_gapped_for_closure.k
                );
            }

            let collect_prelim_hits_for_chunk = |chunk: &SubjectChunk,
                                                 gap_scratch: &mut GapAlignScratch,
                                                 subject_scratch: &mut SubjectScratch|
                                                 -> Vec<PrelimHit> {
                let hits_with_internal = &mut subject_scratch.hits_with_internal;
                hits_with_internal.clear();
                let prelim_hits = &mut subject_scratch.prelim_hits;
                prelim_hits.clear();
                let greedy_align_scratch = &mut subject_scratch.greedy_align_scratch;

                let s_len = chunk.length;
                let subject_seq_ranges = chunk.seq_ranges.as_slice();
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
                    0,                           // q_min
                    (query_concat_length + 1) as i32,  // q_max
                    0,                           // s_min
                    (s_len + 1) as i32,          // s_max
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
                const MAX_ARRAY_DIAG_SIZE: usize = 12_000_000;
                let (diag_array_length_single, diag_mask_single) = if queries.len() == 1 {
                    (diag_array_lengths[0], diag_masks[0] as isize)
                } else {
                    (0usize, 0isize)
                };
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
                let use_array_indexing = !use_diag_hash
                    && queries.len() == 1
                    && diag_array_length_single <= MAX_ARRAY_DIAG_SIZE;
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
                    let search_seq: &[u8] = s_seq_blastna;
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
                        diag_array_length_single
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
                    let init_diag = DiagStruct {
                        last_hit: -diag_window,
                        flag: 0,
                    };
                    if use_array_indexing {
                        if hit_level_array.len() != diag_array_size {
                            hit_level_array.resize(diag_array_size, init_diag);
                            hit_level_array.fill(init_diag);
                            hit_len_array.resize(diag_array_size, 0);
                            hit_len_array.fill(0);
                        } else {
                            hit_len_array.resize(diag_array_size, 0);
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

                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
                        // ```c
                        // diag = s_off + diag_table->diag_array_length - q_off;
                        // real_diag = diag & diag_table->diag_mask;
                        // ```
                        let (diag_array_length, diag_mask) = if use_array_indexing {
                            (diag_array_length_single as isize, diag_mask_single)
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
                            let window_size_i32 = TWO_HIT_WINDOW as i32;
                            let delta_calc = window_size_i32 - word_length as i32;
                            let delta_limit = if delta_calc < 0 {
                                0
                            } else {
                                (scan_range as i32).min(delta_calc)
                            };
                            window_size_i32 + delta_limit + 1
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
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
                        // ```c
                        // offset_array_size = OFFSET_ARRAY_SIZE +
                        //     ((BlastMBLookupTable*)lookup->lut)->longest_chain;
                        // ```
                        let offset_array_size = OFFSET_ARRAY_SIZE + two_stage.longest_chain();
                        let offset_pairs = &mut subject_scratch.offset_pairs;
                        offset_pairs.clear();
                        if offset_pairs.capacity() < offset_array_size {
                            offset_pairs.reserve(offset_array_size - offset_pairs.capacity());
                        }

                        let mut process_offset_pairs = |offset_pairs: &mut Vec<OffsetPair>| {
                            for pair in offset_pairs.drain(..) {
                                let q_off0 = pair.q_off;
                                let kmer_start = pair.s_off;
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:730-733
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
                                let cutoff_score = cutoff_scores[query_idx];

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
                                if q_pos_usize + lut_word_length > q_seq.len()
                                    || kmer_start + lut_word_length > s_len
                                {
                                    continue;
                                }

                                // For two-stage lookup, we use the lut_word_length match position
                                // The actual word_length matching will happen in extension phase
                                // NCBI BLAST: s_BlastnExtendInitialHit does left+right extension to verify word_length
                                let extended = 0; // Will be calculated in extension phase

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
                                // ```c
                                // diag = s_off + diag_table->diag_array_length - q_off;
                                // real_diag = diag & diag_table->diag_mask;
                                // ```
                                let diag = if use_array_indexing {
                                    (kmer_start as isize + diag_array_length - q_pos_usize as isize)
                                        & diag_mask
                                } else {
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-829
                                    // ```c
                                    // diag = s_off - q_off;
                                    // s_end = s_off + word_length;
                                    // ```
                                    kmer_start as isize - q_pos_usize as isize
                                };

                                // Check if this seed is in the debug window
                                let s_pos = kmer_start + chunk.offset;
                                let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                                    q_pos_usize >= q_start && q_pos_usize <= q_end &&
                                    s_pos >= s_start && s_pos <= s_end
                                } else {
                                    false
                                };

                                if debug_enabled && in_window {
                                    dbg_window_seeds += 1;
                                }

                                // NCBI BLAST does NOT check mask_array/mask_hash at seed level
                                // NCBI reference: na_ungapped.c:671-672
                                // NCBI only checks: if (s_off_pos < last_hit) return 0;
                                // This check happens below (line 753) using hit_level_array[real_diag].last_hit
                                // The mask_array/mask_hash check above was LOSAT-specific and caused excessive filtering
                                // REMOVED to match NCBI BLAST behavior

                                // NCBI reference: na_ungapped.c:656-672
                                // Two-hit filter: check if there's a previous hit within window_size
                                // NCBI: Boolean two_hits = (window_size > 0);
                                // NCBI: last_hit = hit_level_array[real_diag].last_hit;
                                // NCBI: hit_saved = hit_level_array[real_diag].flag;
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
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:833-836
                                    // ```c
                                    // rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);
                                    // if(!rc)  last_hit = 0;
                                    // ```
                                    let (level, _hit_len, hit_saved) = diag_hash
                                        .retrieve(diag as i32)
                                        .unwrap_or((0, 0, false));
                                    (level, hit_saved)
                                } else {
                                    (0, false)
                                };

                                // NCBI reference: na_ungapped.c:1081-1140
                                // For two-stage lookup, verify word_length match BEFORE ungapped extension
                                // This is done by left+right extension from the lut_word_length match position
                                let (q_ext_start, s_ext_start) = if word_length > lut_word_length {
                                    // NCBI BLAST: s_BlastNaExtend left/right extension uses packed subject
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1093-1136
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

                                    // LEFT extension (backwards from lut_word_length start)
                                    // CRITICAL: NCBI relies on MIN(ext_to, s_offset) limit, NO explicit boundary check
                                    // NCBI BLAST: MIN(ext_to, s_offset) - limit left extension by s_offset
                                    // s_offset (kmer_start) is the maximum number of positions we can extend left
                                    // NCBI reference: na_ungapped.c:1106-1114
                                    // for (; ext_left < MIN(ext_to, s_offset); ++ext_left) {
                                    //     s_off--;
                                    //     q--;
                                    //     if (s_off % COMPRESSION_RATIO == 3)
                                    //         s--;
                                    //     if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6) != *q)
                                    //         break;
                                    // }
                                    // CRITICAL: NCBI only limits by MIN(ext_to, s_offset) - NO query bounds check
                                    // But LOSAT must also limit by q_pos to prevent underflow
                                    // NCBI assumes sequences are long enough, but we need explicit bounds
                                    let max_ext_left = ext_to.min(kmer_start).min(q_pos_usize);

                                    // NCBI BLAST: Loop condition is ext_left < MIN(ext_to, s_offset)
                                    // NO explicit boundary check inside loop (NCBI doesn't have it)
                                    // SAFETY: We've limited by subject bounds. For query, we trust sequences are long enough
                                    // (NCBI doesn't check query bounds either - it's undefined behavior if out of bounds)
                                    let mut ext_left = 0usize;
                                    let mut q_left = q_pos_usize;
                                    let mut s_off = kmer_start;
                                    let mut s_idx = s_off / COMPRESSION_RATIO;
                                    while ext_left < max_ext_left {
                                        if debug_enabled {
                                            dbg_left_ext_iters += 1;
                                        }
                                        // NCBI BLAST: s_off--; q--;
                                        // Reference: na_ungapped.c:1107-1108
                                        q_left -= 1;
                                        s_off -= 1;
                                        if s_off % COMPRESSION_RATIO == COMPRESSION_RATIO - 1 {
                                            s_idx -= 1;
                                        }
                                        // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by max_ext_left.
                                        let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
                                        let s_base = ((s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6) as u8;
                                        // SAFETY: q_left is bounded by max_ext_left and q_pos_usize.
                                        let q_base = unsafe { *q_seq_blastna.get_unchecked(q_left) };
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
                                    if ext_left < ext_to {
                                        // NCBI BLAST: s_off = s_offset + lut_word_length;
                                        // NCBI BLAST: if (s_off + ext_to - ext_left > s_range) continue;
                                        // Reference: na_ungapped.c:1122-1124
                                        let mut s_off = kmer_start + lut_word_length;
                                        if s_off + (ext_to - ext_left) > s_len {
                                            // Not enough room in subject, skip this seed (NCBI behavior)
                                            continue;
                                        }

                                        let mut ext_right = 0usize;
                                        let q_off = q_pos_usize + lut_word_length;
                                        // LOSAT-specific: Also check query bounds (NCBI doesn't, causing UB)
                                        if q_off + (ext_to - ext_left) > q_seq_blastna.len() {
                                            // Not enough room in query, skip this seed
                                            continue;
                                        }
                                        let mut q_right = q_off;
                                        let mut s_idx = s_off / COMPRESSION_RATIO;

                                        // NCBI BLAST: for (; ext_right < ext_to - ext_left; ++ext_right)
                                        // Reference: na_ungapped.c:1128-1136
                                        // NO bounds checks in loop condition (NCBI relies on pre-loop subject check only)
                                        // NCBI doesn't check query bounds - it would access out of bounds (undefined behavior)
                                        // For Rust safety, we use unsafe indexing but trust sequences are long enough
                                        // (lookup table only finds matches within valid ranges, so sequences should be long enough)
                                        while ext_right < (ext_to - ext_left) {
                                            if debug_enabled {
                                                dbg_right_ext_iters += 1;
                                            }
                                            // NCBI BLAST: if (base mismatch) break;
                                            // Reference: na_ungapped.c:1129-1131
                                            // SAFETY: s_idx tracks s_off / COMPRESSION_RATIO within bounds by s_len check above.
                                            let s_byte = unsafe { *search_seq_packed.get_unchecked(s_idx) };
                                            let s_base = ((s_byte << (2 * (s_off % COMPRESSION_RATIO))) >> 6) as u8;
                                            // SAFETY: q_right is bounded by q_seq_blastna length check above.
                                            let q_base = unsafe { *q_seq_blastna.get_unchecked(q_right) };
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
                                    (q_pos_usize - ext_left, kmer_start - ext_left)
                                } else {
                                    // word_length == lut_word_length, no extension needed
                                    (q_pos_usize, s_pos)
                                };

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:656-672
                                // ```c
                                // last_hit = hit_level_array[real_diag].last_hit;
                                // s_off_pos = s_off + diag_table->offset;
                                // if (s_off_pos < last_hit) return 0;
                                // ```
                                // NCBI calls s_Blastn* after left/right word-length verification, so s_off is
                                // already adjusted by ext_left (s_BlastNaExtendDirect).
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
                                    let hits = two_stage.get_hits(kmer);
                                    if hits.is_empty() {
                                        return true;
                                    }
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
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:674-680
                                    // ```c
                                    // word_type = s_TypeOfWord(query, subject, &q_off, &s_off,
                                    //                          query_mask, query_info, s_range,
                                    //                          word_length, lut_word_length, lut, TRUE, &extended);
                                    // ```
                                    // s_TypeOfWord uses query_mask to skip masked seeds
                                    let (wt, ext, q_off_adj, s_off_adj) = type_of_word(
                                        q_seq,
                                        search_seq,
                                        q_off,
                                        s_off,
                                        query_mask,
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
                                                    && (off_s_end as isize - off_s_l as isize) <= s_b {
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
                                                    && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                                                    diag_hash.retrieve(diag_i32 + delta as i32) {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize - delta) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize) <= s_b {
                                                        off_found = true;
                                                        break;
                                                    }
                                                }
                                                if let Some((off_s_end, off_s_l, _)) =
                                                    diag_hash.retrieve(diag_i32 - delta as i32) {
                                                    if off_s_l > 0
                                                        && (off_s_end as isize) >= s_a
                                                        && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                                        q_seq,
                                        search_seq,
                                        q_off,
                                        s_off,
                                        query_mask,
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

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:728-766
                                // ```c
                                // if (hit_ready) {
                                //     if (word_params->ungapped_extension) {
                                //         ...
                                //     }
                                // }
                                // ```
                                if !hit_ready {
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                                    // ```c
                                    // hit_level_array[real_diag].last_hit = s_end_pos;
                                    // hit_level_array[real_diag].flag = hit_ready;
                                    // diag_table->hit_len_array[real_diag] =
                                    //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
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
                                        hit_len_array[diag_idx] = (s_end_pos - s_off_pos) as u8;
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
                                let x_dropoff = x_dropoff_scores[query_idx];
                                let reduced_cutoff = reduced_cutoff_scores[query_idx];
                                let ungapped = if word_length < 11 {
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
                                let qs = ungapped.q_start;
                                let qe = ungapped.q_start + ungapped.length;
                                let ss = ungapped.s_start;
                                let ungapped_se = ungapped.s_start + ungapped.length;
                                let ungapped_score = ungapped.score;

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
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                                    // ```c
                                    // hit_level_array[real_diag].last_hit = s_end_pos;
                                    // hit_level_array[real_diag].flag = hit_ready;
                                    // diag_table->hit_len_array[real_diag] =
                                    //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
                                    // ```
                                    if use_array_indexing && diag_idx < diag_array_size {
                                        hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                                        hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                                        hit_len_array[diag_idx] = if hit_ready {
                                            0
                                        } else {
                                            (s_end_pos - s_off_pos) as u8
                                        };
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
                                    seed_q_off: q_off,
                                    seed_s_off: s_off,
                                    qs,
                                    qe,
                                    ss,
                                    se: ungapped_se,
                                    score: ungapped_score,
                                });

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                                // ```c
                                // hit_level_array[real_diag].last_hit = s_end_pos;
                                // hit_level_array[real_diag].flag = hit_ready;
                                // diag_table->hit_len_array[real_diag] =
                                //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
                                // ```
                                // Update hit_level_array after collecting ungapped hit
                                hit_ready = true;
                                if use_array_indexing && diag_idx < diag_array_size {
                                    hit_level_array[diag_idx].last_hit = ungapped_s_end_pos as i32;
                                    hit_level_array[diag_idx].flag = 1;
                                    hit_len_array[diag_idx] = 0;
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
                        };

                        scan_subject_kmers_with_ranges(
                            search_seq_packed,
                            s_len,
                            word_length,
                            lut_word_length,
                            scan_step,
                            &subject_seq_ranges,
                            subject_masked,
                            |kmer_start, current_lut_kmer| {
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
                                let matches_slice = two_stage.get_hits(current_lut_kmer);

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
                                if (debug_mode || blastn_debug) && matches_slice.len() > 1000 {
                                    eprintln!("[WARN] Large matches_slice: len={} for kmer at position {}", matches_slice.len(), kmer_start);
                                }

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1504-1533
                                // ```c
                                // max_hits -= mb_lt->longest_chain;
                                // if (total_hits >= max_hits)
                                //     break;
                                // total_hits += s_BlastMBLookupRetrieve(mb_lt,
                                //     index, offset_pairs + total_hits, s_off);
                                // ```
                                if offset_pairs.len() >= OFFSET_ARRAY_SIZE {
                                    process_offset_pairs(offset_pairs);
                                }

                                // For each match, add to the offset_pairs buffer.
                                for &q_off_1 in matches_slice {
                                    if debug_enabled {
                                        dbg_seeds_found += 1;
                                    }
                                    if q_off_1 == 0 {
                                        continue;
                                    }
                                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
                                    // ```c
                                    // offset_pairs[i].qs_offsets.q_off   = q_off - 1;
                                    // offset_pairs[i++].qs_offsets.s_off = s_off;
                                    // ```
                                    offset_pairs.push(OffsetPair {
                                        q_off: q_off_1 as usize - 1,
                                        s_off: kmer_start,
                                    });
                                } // end of for matches_slice in two-stage lookup

                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
                                // ```c
                                // offset_array_size = OFFSET_ARRAY_SIZE +
                                //     ((BlastMBLookupTable*)lookup->lut)->longest_chain;
                                // ```
                                if offset_pairs.len() >= offset_array_size {
                                    process_offset_pairs(offset_pairs);
                                }
                            },
                        );

                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:451-497
                        // ```c
                        // hitsfound = scansub(lookup_wrap, subject, offset_pairs,
                        //                     kScanSubjectOffsetArraySize, &scan_range[1]);
                        // if (hitsfound == 0) continue;
                        // hits_extended += extend(offset_pairs, hitsfound, ...);
                        // ```
                        if !offset_pairs.is_empty() {
                            process_offset_pairs(offset_pairs);
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
                        scan_subject_kmers_with_ranges(
                            search_seq_packed,
                            s_len,
                            safe_k,
                            safe_k,
                            scan_step,
                            &subject_seq_ranges,
                            subject_masked,
                            |kmer_start, current_kmer| {
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

                            // Phase 2: Use PV-based direct lookup (O(1) with fast PV filtering) for word_size <= 13
                            // For word_size > 13, use hash-based lookup
                            let matches_slice: &[u32] = if use_direct_lookup {
                                // Use PV for fast filtering before accessing the lookup table
                                pv_direct_lookup_ref.map(|pv_dl| pv_dl.get_hits_checked(current_kmer)).unwrap_or(&[])
                            } else {
                                // Use hash-based lookup for larger word sizes
                                hash_lookup_ref.and_then(|hl| hl.get(&current_kmer).map(|v| v.as_slice())).unwrap_or(&[])
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
                                let cutoff_score = cutoff_scores[query_idx];

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:663-664
                            // ```c
                            // diag = s_off + diag_table->diag_array_length - q_off;
                            // real_diag = diag & diag_table->diag_mask;
                            // ```
                            let (diag_array_length, diag_mask) = if use_array_indexing {
                                (diag_array_length_single as isize, diag_mask_single)
                            } else {
                                (0isize, 0isize)
                            };
                            let diag = if use_array_indexing {
                                (kmer_start as isize + diag_array_length - q_pos_usize as isize)
                                    & diag_mask
                            } else {
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:828-829
                                // ```c
                                // diag = s_off - q_off;
                                // s_end = s_off + word_length;
                                // ```
                                kmer_start as isize - q_pos_usize as isize
                            };

                            // Check if this seed is in the debug window
                            let s_pos = kmer_start + chunk.offset;
                            let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                                q_pos_usize >= q_start && q_pos_usize <= q_end &&
                                s_pos >= s_start && s_pos <= s_end
                            } else {
                                false
                            };

                            if debug_enabled && in_window {
                                dbg_window_seeds += 1;
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
                                let (level, _hit_len, hit_saved) = diag_hash
                                    .retrieve(diag as i32)
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
                                // if (Delta < 0) Delta = 0;
                                // s_BlastDiagHashInsert(hash_table, diag, s_end_pos,
                                //                       (hit_ready) ? 0 : s_end_pos - s_off_pos,
                                //                       hit_ready, s_off_pos, window_size + Delta + 1);
                                // ```
                                let window_size_i32 = TWO_HIT_WINDOW as i32;
                                let delta_calc = window_size_i32 - safe_k as i32;
                                let delta_limit = if delta_calc < 0 {
                                    0
                                } else {
                                    (scan_range as i32).min(delta_calc)
                                };
                                window_size_i32 + delta_limit + 1
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
                                    hash_lookup_ref
                                        .and_then(|hl| hl.get(&kmer).map(|v| v.as_slice()))
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
                                    q_seq,
                                    search_seq,
                                    q_off,
                                    s_off,
                                    query_mask,
                                    safe_k, // word_length
                                    safe_k, // lut_word_length (same for non-two-stage)
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
                                    let s_a = s_off_pos as isize + safe_k as isize - window_size as isize;
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
                                                    let off_entry = &hit_level_array[off_diag_idx];
                                                    (off_entry.last_hit, hit_len_array[off_diag_idx])
                                                } else {
                                                    (0, 0)
                                                }
                                            };
                                            if off_s_l > 0
                                                && (off_s_end as isize - delta) >= s_a
                                                && (off_s_end as isize - off_s_l as isize) <= s_b {
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
                                                    let off_entry = &hit_level_array[off_diag_idx];
                                                    (off_entry.last_hit, hit_len_array[off_diag_idx])
                                                } else {
                                                    (0, 0)
                                                }
                                            };
                                            if off_s_l > 0
                                                && (off_s_end as isize) >= s_a
                                                && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                                                diag_hash.retrieve(diag_i32 + delta as i32) {
                                                if off_s_l > 0
                                                    && (off_s_end as isize - delta) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize) <= s_b {
                                                    off_found = true;
                                                    break;
                                                }
                                            }
                                            if let Some((off_s_end, off_s_l, _)) =
                                                diag_hash.retrieve(diag_i32 - delta as i32) {
                                                if off_s_l > 0
                                                    && (off_s_end as isize) >= s_a
                                                    && (off_s_end as isize - off_s_l as isize + delta) <= s_b {
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
                                    q_seq,
                                    search_seq,
                                    q_off,
                                    s_off,
                                    query_mask,
                                    safe_k, // word_length
                                    safe_k, // lut_word_length (same for non-two-stage)
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

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:728-766
                            // ```c
                            // if (hit_ready) {
                            //     if (word_params->ungapped_extension) {
                            //         ...
                            //     }
                            // }
                            // ```
                            if !hit_ready {
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                                // ```c
                                // hit_level_array[real_diag].last_hit = s_end_pos;
                                // hit_level_array[real_diag].flag = hit_ready;
                                // diag_table->hit_len_array[real_diag] =
                                //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
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
                                    hit_len_array[diag_idx] = (s_end_pos - s_off_pos) as u8;
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
                            let x_dropoff = x_dropoff_scores[query_idx];
                            let reduced_cutoff = reduced_cutoff_scores[query_idx];
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
                            let qs = ungapped.q_start;
                            let qe = ungapped.q_start + ungapped.length;
                            let ss = ungapped.s_start;
                            let ungapped_se = ungapped.s_start + ungapped.length;
                            let ungapped_score = ungapped.score;

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
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                                // ```c
                                // hit_level_array[real_diag].last_hit = s_end_pos;
                                // hit_level_array[real_diag].flag = hit_ready;
                                // diag_table->hit_len_array[real_diag] =
                                //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
                                // ```
                                if use_array_indexing && diag_idx < diag_array_size {
                                    hit_level_array[diag_idx].last_hit = s_end_pos as i32;
                                    hit_level_array[diag_idx].flag = if hit_ready { 1 } else { 0 };
                                    hit_len_array[diag_idx] = if hit_ready {
                                        0
                                    } else {
                                        (s_end_pos - s_off_pos) as u8
                                    };
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

                            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:768-771
                            // ```c
                            // hit_level_array[real_diag].last_hit = s_end_pos;
                            // hit_level_array[real_diag].flag = hit_ready;
                            // diag_table->hit_len_array[real_diag] =
                            //     (hit_ready) ? 0 : s_end_pos - s_off_pos;
                            // ```
                            // Update hit_level_array after collecting ungapped hit
                            hit_ready = true;
                            if use_array_indexing && diag_idx < diag_array_size {
                                hit_level_array[diag_idx].last_hit = ungapped_s_end_pos as i32;
                                hit_level_array[diag_idx].flag = 1;
                                hit_len_array[diag_idx] = 0;
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
                    } // end of if/else for two-stage vs original lookup
                } // end of strand loop

                // ============================================================
                // NCBI ARCHITECTURE: Batch process ungapped hits in score order
                // ============================================================
                // NCBI reference: blast_gapalign.c:3824 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))
                // NCBI processes init_hsp_array sorted by score DESCENDING
                // High-score HSPs are processed first, their GAPPED HSPs added to interval tree
                // Lower-score HSPs are then checked for containment against these GAPPED HSPs

                // Sort by score DESCENDING (highest first)
                ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

                // Debug counters for containment analysis
                let mut dbg_containment_skipped = 0usize;
                let total_ungapped = ungapped_hits.len();
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
                    let is_contained = interval_tree.contains_hsp(
                        &ungapped_tree_hsp,
                        uh.query_context_offset,
                        min_diag_separation,
                    );

                    if is_contained {
                        // NCBI: Skip gapped extension if ungapped HSP is contained
                        dbg_containment_skipped += 1;
                        continue;
                    }

                    // Select gapped-start seed within the ungapped HSP.
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4012-4046
                    let cutoff_score = cutoff_scores[uh.query_idx as usize];
                    let (
                        prelim_qs,
                        prelim_qe,
                        prelim_ss,
                        prelim_se,
                        prelim_score,
                        seed_qs,
                        seed_ss,
                    ) = if use_dp {
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

                        // Preliminary DP gapped extension (score-only)
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:176-183
                        let (p_qs, p_qe, p_ss, p_se, p_score, _, _, _, _, _) = extend_gapped_heuristic_with_scratch(
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
                    });
                    continue;

                }

                // DEBUG: Log gapped extension stats
                let gapped_elapsed = gapped_start_time.elapsed();
                if verbose {
                    eprintln!(
                        "[INFO] Gapped extension took {:?}: calls={}, skipped={}, total_ungapped={}",
                        gapped_elapsed,
                        dbg_gapped_calls,
                        dbg_containment_skipped,
                        total_ungapped
                    );
                }


                let mut chunk_prelim_hits: Vec<PrelimHit> = prelim_hits.drain(..).collect();
                if !chunk_prelim_hits.is_empty() {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
                    // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
                    chunk_prelim_hits = purge_prelim_hits_with_common_endpoints(chunk_prelim_hits);
                    chunk_prelim_hits.sort_by(score_compare_prelim_hits);
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
            if use_parallel && subject_chunks.len() > 1 {
                let mut chunk_results: Vec<(usize, usize, Vec<PrelimHit>)> = subject_chunks
                    .par_iter()
                    .map_init(
                        || (GapAlignScratch::new(), SubjectScratch::new(queries_ref.len())),
                        |state, chunk| {
                            let (gap_scratch, subject_scratch) = state;
                            let hits = collect_prelim_hits_for_chunk(chunk, gap_scratch, subject_scratch);
                            (chunk.offset, chunk.overlap, hits)
                        },
                    )
                    .collect();
                chunk_results.sort_by_key(|(offset, _, _)| *offset);
                for (offset, overlap, hits) in chunk_results {
                    merge_prelim_hits_subject_split(&mut combined_prelim_hits, hits, offset, overlap, true);
                }
            } else {
                for chunk in subject_chunks.iter() {
                    let hits = collect_prelim_hits_for_chunk(chunk, gap_scratch, subject_scratch);
                    merge_prelim_hits_subject_split(&mut combined_prelim_hits, hits, chunk.offset, chunk.overlap, true);
                }
            }

            if combined_prelim_hits.is_empty() {
                return;
            }

            let hits_with_internal = &mut subject_scratch.hits_with_internal;
            hits_with_internal.clear();
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

            let prelim_hits = &mut combined_prelim_hits;
            if !prelim_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:358-365
                // ```c
                // /* Make sure the HSPs in the HSP list are sorted by score, as they should be. */
                // ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
                // ```
                prelim_hits.sort_by(score_compare_prelim_hits);

                interval_tree.reset();

                for prelim in prelim_hits.iter() {
                    let ctx = &query_contexts[prelim.context_idx as usize];
                    let q_seq_blastna = encoded_queries_blastna[prelim.context_idx as usize].as_slice();

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:403-405
                    // ```c
                    // if (program_number == eBlastTypeRpsBlast ||
                    //     !BlastIntervalTreeContainsHSP(tree, hsp, query_info,
                    //                          hit_options->min_diag_separation)) {
                    // ```
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
                    if interval_tree.contains_hsp(
                        &prelim_tree_hsp,
                        prelim.query_context_offset,
                        min_diag_separation,
                    ) {
                        continue;
                    }

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
                            None => {
                                continue;
                            }
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
                    let mut start_shift: usize = 0;
                    let mut adjusted_subject = s_seq_blastna;
                    let mut adjusted_s_len = s_seq_blastna.len();
                    let mut trace_s_start_adj = trace_s_start;
                    let mut s_start_i32 = trace_s_start as i32;
                    let mut s_len_i32 = adjusted_s_len as i32;
                    let start_shift_i32 = adjust_subject_range(
                        &mut s_start_i32,
                        &mut s_len_i32,
                        trace_q_start as i32,
                        q_seq_blastna.len() as i32,
                    );
                    start_shift = start_shift_i32 as usize;
                    adjusted_s_len = s_len_i32 as usize;
                    trace_s_start_adj = s_start_i32 as usize;
                    adjusted_subject = &s_seq_blastna[start_shift..start_shift + adjusted_s_len];

                    let x_drop_trace = x_drop_final;
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
                        edit_ops,
                    ) = if use_dp {
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
                            Some(value) => value,
                            None => {
                                continue;
                            }
                        }
                    };

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:598-600
                    // ```c
                    // Blast_HSPAdjustSubjectOffset(hsp, start_shift);
                    // ```
                    if start_shift != 0 {
                        final_ss = final_ss.saturating_add(start_shift);
                        final_se = final_se.saturating_add(start_shift);
                    }

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

                    let aln_len = matches + mismatches + gap_letters;

                    if use_dp && hsp_test(matches, aln_len, percent_identity, min_hit_length) {
                        continue;
                    }

                    let (bit_score, eval) = calculate_evalue(
                        score,
                        ctx.seq.len(),
                        db_len_total,
                        db_num_seqs,
                        &params_for_closure,
                    );

                    if eval > evalue_threshold {
                        continue;
                    }

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

                    let internal = InternalHitData {
                        q_offset_0: final_qs,
                        q_end_0: final_qe,
                        s_offset_0: final_ss,
                        s_end_0: final_se,
                        query_frame: prelim.query_frame,
                        query_context_offset: prelim.query_context_offset,
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
                    hits_with_internal.push((BlastnHsp {
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
                        query_frame: prelim.query_frame,
                        query_length,
                        q_idx: prelim.query_idx,
                        s_idx: s_idx as u32,
                        raw_score: score,
                        gap_info,
                    }, internal));
                }

                prelim_hits.clear();
            }

            // =================================================================
            // NCBI POST-GAPPED PROCESSING
            // Reference: blast_traceback.c:633-692
            // =================================================================

            // Step 1: Extract hits and internal coordinates
            let (local_hits, internals): (Vec<BlastnHsp>, Vec<InternalHitData>) =
                hits_with_internal.drain(..).unzip();

            // Step 2: Endpoint purging pass 1 (trim, purge=false)
            // NCBI reference: blast_traceback.c:637-638
            // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
            let hits_step1 = local_hits.len();
            let (mut local_hits, mut extra_start) = purge_hsps_with_common_endpoints_ex(local_hits, false);
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
                let cutoff = cutoff_scores.get(hit.q_idx as usize).copied().unwrap_or(0);
                let context_idx = (hit.q_idx as usize) * 2 + if hit.query_frame < 0 { 1 } else { 0 };

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1891
                // ```c
                // hsp->evalue =
                //   BLAST_KarlinStoE_simple(score, kbp[kbp_context],
                //        query_info->contexts[hsp->context].eff_searchsp);
                // ```
                let eff_searchsp = query_eff_searchsp
                    .get(context_idx)
                    .copied()
                    .unwrap_or(0);
                let reeval_params = ReevalParams {
                    lambda: params_for_closure.lambda,
                    k: params_for_closure.k,
                    eff_searchsp,
                    db_len: db_len_total,
                    db_num_seqs,
                };

                // NCBI reference: blast_traceback.c:653-665 (reevaluate with blastna sequences)
                let q_seq_blastna = encoded_queries_blastna[context_idx].as_slice();
                let s_seq_eval = s_seq_blastna;
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
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:656-663
                // ```c
                // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
                // if (!delete_hsp)
                //     delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
                //                                                 query_nomask, subject,
                //                                                 score_options, hit_options);
                // ```
                if delete {
                    hit.raw_score = i32::MIN; // Mark for removal
                    continue;
                }

                if blast_hsp_test_identity_and_length(
                    hit,
                    q_seq_blastna,
                    s_seq_eval,
                    percent_identity,
                    min_hit_length,
                ) {
                    hit.raw_score = i32::MIN; // Mark for removal
                }
            }
            local_hits.retain(|h| h.raw_score != i32::MIN);
            let hits_step3 = local_hits.len();

            // Step 4: Endpoint purging pass 2 (delete, purge=true) - BLASTN only
            // NCBI reference: blast_traceback.c:667-668
            // if(program_number == eBlastTypeBlastn) {
            //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
            // }
            let mut local_hits = purge_hsps_with_common_endpoints(local_hits);
            let hits_step4 = local_hits.len();

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
            local_hits.sort_by(score_compare_blastn_hsps);

            // Step 6: Phase 2 - Tree reset and second containment pass
            // NCBI reference: blast_traceback.c:678
            // Blast_IntervalTreeReset(tree);
            interval_tree.reset();

            // NCBI reference: blast_traceback.c:679-692
            // Remove any HSPs that are contained within other HSPs.
            // Since the list is sorted by score already, any HSP
            // contained by a previous HSP is guaranteed to have a
            // lower score, and may be purged.
            let mut final_hits = Vec::with_capacity(local_hits.len());
            for hit in local_hits {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
                // ```c
                // if (hsp->query.frame != hsp->subject.frame) {
                //    *q_end = query_length - hsp->query.offset;
                //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
                // }
                // ```
                let (q_offset_0, q_end_0) = if hit.query_length > 0 && hit.query_frame < 0 {
                    (
                        hit.query_length.saturating_sub(hit.q_end),
                        hit.query_length
                            .saturating_sub(hit.q_start)
                            .saturating_add(1),
                    )
                } else {
                    (hit.q_start.saturating_sub(1), hit.q_end)
                };
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:537-550
                // ```c
                // query_start = query_context_offset;
                // region_start = query_start + hsp->query.offset;
                // region_end = query_start + hsp->query.end;
                // ```
                let base_offset = query_base_offsets[hit.q_idx as usize];
                let query_context_offset = if hit.query_frame < 0 {
                    (base_offset + hit.query_length + 1) as i32
                } else {
                    base_offset as i32
                };
                // Reconstruct TreeHsp from hit (convert 1-based back to 0-based)
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3908-3913
                // ```c
                // tmp_hsp.query.offset = q_start;
                // tmp_hsp.query.end = q_end;
                // tmp_hsp.query.frame = query_info->contexts[context].frame;
                // tmp_hsp.subject.offset = s_start;
                // tmp_hsp.subject.end = s_end;
                // ```
                // NCBI uses canonical coordinates: subject.offset < subject.end always
                let tree_hsp = TreeHsp {
                    query_offset: q_offset_0 as i32,
                    query_end: q_end_0 as i32,
                    subject_offset: (hit.s_start.min(hit.s_end).saturating_sub(1)) as i32,
                    subject_end: hit.s_start.max(hit.s_end) as i32,
                    score: hit.raw_score,
                    query_frame: hit.query_frame,
                    query_length: hit.query_length as i32,
                    query_context_offset,
                    subject_frame_sign: 1,
                };

                // NCBI reference: blast_traceback.c:682-688
                if !interval_tree.contains_hsp(&tree_hsp, query_context_offset, min_diag_separation) {
                    interval_tree.add_hsp(tree_hsp, query_context_offset, IndexMethod::QueryAndSubject);
                    final_hits.push(hit);
                }
                // else: HSP is contained within another, skip (implicit delete)
            }

            if subject_besthit && !final_hits.is_empty() {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:858-862
                // ```c
                // if ((hit_options->hsp_filt_opt != NULL) &&
                //     (hit_options->hsp_filt_opt->subject_besthit_opts != NULL)) {
                //    Blast_HSPListSubjectBestHit(program_number,
                //        hit_options->hsp_filt_opt->subject_besthit_opts,
                //        query_info, hsp_list);
                // }
                // ```
                let mut grouped: FxHashMap<u32, Vec<BlastnHsp>> = FxHashMap::default();
                for hit in final_hits.drain(..) {
                    grouped.entry(hit.q_idx).or_default().push(hit);
                }
                let mut filtered = Vec::new();
                for (q_idx, mut hits) in grouped {
                    let query_len = queries[q_idx as usize].seq().len();
                    subject_best_hit(&mut hits, query_len);
                    filtered.extend(hits);
                }
                final_hits = filtered;
            }

            // Phase 2 debug output
            if debug_mode || blastn_debug {
                let phase2_tree_removed = hits_step4 - final_hits.len();
                eprintln!(
                    "[DEBUG] Phase 2 breakdown: step1(extract)={}, step2(purge1)={} (-{}), step3(reeval)={} (-{}), step4(purge2)={} (-{}), step6(tree)={} (-{})",
                    hits_step1,
                    hits_step2, hits_step1 - hits_step2,
                    hits_step3, hits_step2 - hits_step3,
                    hits_step4, hits_step3 - hits_step4,
                    final_hits.len(), phase2_tree_removed
                );
            }

            // Suppress unused warning
            let _ = internals;

            if !final_hits.is_empty() {
                // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
                *subject_hits = Some(final_hits);
            }
            if let Some(bar) = progress_bar.as_ref() {
                bar.inc(1);
            }
        };

    if use_parallel {
        // Channel for sending hits
        // Use Option to signal completion: None means "all subjects processed"
        let (tx, rx) = channel::<Option<Vec<BlastnHsp>>>();
        let out_path = args.out.clone();
        // Keep a sender for the main thread to send the completion signal
        let tx_main = tx.clone();
        let verbose = args.verbose;
        let query_ids_arc = Arc::clone(&query_ids_arc);
        let subject_ids_arc = Arc::clone(&subject_ids_arc);

        let writer_handle = std::thread::spawn(move || -> Result<()> {
            if verbose {
                eprintln!("[INFO] Writer thread started, waiting for hits...");
            }
            let mut all_hits = Vec::new();
            let mut messages_received = 0usize;

            while let Ok(msg) = rx.recv() {
                match msg {
                    Some(hits) => {
                        messages_received += 1;
                        if verbose && (messages_received == 1 || messages_received % 100 == 0) {
                            eprintln!(
                                "[INFO] Received message #{}, {} hits so far",
                                messages_received,
                                all_hits.len() + hits.len()
                            );
                        }
                        all_hits.extend(hits);
                    }
                    None => {
                        // Completion signal received
                        if verbose {
                            eprintln!(
                                "[INFO] Completion signal received after {} messages",
                                messages_received
                            );
                        }
                        break;
                    }
                }
            }

            post_process_hits_and_write(
                all_hits,
                hitlist_size,
                max_hsps_per_subject,
                &out_path,
                verbose,
                query_ids_arc.as_ref(),
                subject_ids_arc.as_ref(),
            )
        });

        seq_data.subjects
            .par_iter()
            .enumerate()
            .for_each_init(
                || {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:313-319 (BLAST_GapAlignStructNew)
                    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct per thread)
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
                    (tx.clone(), GapAlignScratch::new(), SubjectScratch::new(queries_ref.len()))
                },
                |state, (s_idx, s_record)| {
                    let (tx, gap_scratch, subject_scratch) = state;
                    let mut subject_hits: Option<Vec<BlastnHsp>> = None;
                    process_subject(s_idx, s_record, gap_scratch, subject_scratch, &mut subject_hits);
                    if let Some(hits) = subject_hits {
                        // NCBI reference: blast_traceback.c:633-692 (post-gapped processing is per-subject)
                        tx.send(Some(hits)).unwrap();
                    }
                },
            );

        if let Some(bar) = progress_bar.as_ref() {
            bar.finish();
        }
        if args.verbose {
            eprintln!("[INFO] Parallel processing complete, sending completion signal...");
        }

        // Send completion signal to writer thread
        // This ensures the writer thread exits even if sender-drop semantics are delayed
        tx_main.send(None).unwrap();
        drop(tx_main); // Explicitly drop to ensure channel closes

        writer_handle.join().unwrap()?;
        return Ok(());
    }

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
    let mut all_hits = Vec::new();
    let mut gap_scratch = GapAlignScratch::new();
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
    let mut subject_scratch = SubjectScratch::new(queries_ref.len());
    for (s_idx, s_record) in seq_data.subjects.iter().enumerate() {
        let mut subject_hits: Option<Vec<BlastnHsp>> = None;
        process_subject(
            s_idx,
            s_record,
            &mut gap_scratch,
            &mut subject_scratch,
            &mut subject_hits,
        );
        if let Some(hits) = subject_hits {
            all_hits.extend(hits);
        }
    }

    if let Some(bar) = progress_bar.as_ref() {
        bar.finish();
    }
    post_process_hits_and_write(
        all_hits,
        hitlist_size,
        max_hsps_per_subject,
        &args.out,
        args.verbose,
        query_ids_arc.as_ref(),
        subject_ids_arc.as_ref(),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:689-706
    // ```c
    // p.first = (slp->GetInt().GetFrom() > offset)? slp->GetInt().GetFrom() - offset : 0;
    // p.second = MIN(slp->GetInt().GetTo() - offset, length-1);
    // if (slp->GetInt().GetTo() >= offset && p.first < length) {
    //     output.push_back(p);
    // }
    // ```
    #[test]
    fn test_build_subject_seq_ranges_from_masks() {
        let masks = vec![
            MaskedInterval::new(2, 5),
            MaskedInterval::new(8, 12),
        ];
        let ranges = build_subject_seq_ranges_from_masks(&masks, 10);
        assert_eq!(ranges, vec![(2, 4), (8, 9)]);
    }
}
