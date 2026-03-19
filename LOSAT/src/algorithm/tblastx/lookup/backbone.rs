//! Backbone lookup table for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c

use super::super::diagnostics::diagnostics_enabled;
use super::super::translation::QueryFrame;
use super::{
    compute_backbone_size, compute_mask, compute_unmasked_intervals, encode_kmer_3, ilog2, pv_set,
    QueryContext, LOOKUP_ALPHABET_SIZE, LOOKUP_WORD_LENGTH, PV_ARRAY_BTS,
};
use crate::config::ScoringMatrix;
use crate::stats::karlin_calc::{
    apply_check_ideal, compute_aa_composition, compute_karlin_params_ungapped,
    compute_score_freq_profile, compute_std_aa_composition,
};
use crate::stats::{lookup_protein_params_ungapped, KarlinParams};
use crate::utils::matrix::blosum62_score;

pub const AA_HITS_PER_CELL: usize = 3;

/// NCBI AaLookupBackboneCell - EXACT port
#[repr(C)]
#[derive(Clone, Copy)]
pub struct BackboneCell {
    pub num_used: i32,
    pub entries: [i32; AA_HITS_PER_CELL],
}

impl Default for BackboneCell {
    fn default() -> Self {
        Self {
            num_used: 0,
            entries: [0; AA_HITS_PER_CELL],
        }
    }
}

/// NCBI-style lookup table with flat backbone
pub struct BlastAaLookupTable {
    pub backbone: Vec<BackboneCell>,
    pub overflow: Vec<i32>,
    // Presence vector: PV_ARRAY_TYPE = Uint4
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:41-44
    pub pv: Vec<u32>,
    pub frame_bases: Vec<i32>,
    pub num_contexts: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:246-250
    // ```c
    // Uint4 QueryInfo_GetSeqBufLen(const BlastQueryInfo* qinfo)
    // {
    //     BlastContextInfo * cinfo = & qinfo->contexts[qinfo->last_context];
    //     return cinfo->query_offset + cinfo->query_length + (cinfo->query_length ? 2 : 1);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:112-127
    // ```c
    // (*seq_blk)->sequence_start = (Uint1 *) buffer;
    // (*seq_blk)->sequence = (*seq_blk)->sequence_start+1;
    // (*seq_blk)->length = length;
    // ```
    pub query_length: i32,
    pub word_length: i32,
    pub alphabet_size: i32,
    pub charsize: i32,
    pub mask: i32,
    pub longest_chain: i32,
    // NCBI BlastAaLookupTable has no lazy-neighbor flag; only NCBI fields are stored.
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:95-140
    pub threshold: i32,
    pub row_max: Vec<i32>,
}

impl BlastAaLookupTable {
    #[inline]
    pub fn get_context_idx(&self, concat_off: i32) -> usize {
        let bases = &self.frame_bases;
        let mut lo = 0usize;
        let mut hi = self.num_contexts;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if concat_off < bases[mid] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        lo.saturating_sub(1)
    }

    #[inline(always)]
    pub fn get_hits(&self, index: usize) -> &[i32] {
        let cell = unsafe { self.backbone.get_unchecked(index) };
        let num = cell.num_used as usize;
        if num == 0 {
            return &[];
        }
        if num <= AA_HITS_PER_CELL {
            unsafe { std::slice::from_raw_parts(cell.entries.as_ptr(), num) }
        } else {
            let cursor = cell.entries[0] as usize;
            unsafe { std::slice::from_raw_parts(self.overflow.as_ptr().add(cursor), num) }
        }
    }
}

type LookupChain = Vec<i32>;
type LookupBackboneChains = Vec<Option<LookupChain>>;

#[derive(Default)]
struct LookupBuildStats {
    total_exact_positions: usize,
    skipped_invalid_residue: usize,
    skipped_seg_mask: usize,
    exact_added_count: usize,
    neighbor_added_count: usize,
    words_with_exact_only: usize,
    words_with_neighbors: usize,
    max_neighbors_for_single_word: usize,
    max_neighbor_word_idx: usize,
    neighbor_words_generated: usize,
}

struct NeighborInfo<'a> {
    thin_backbone: &'a mut LookupBackboneChains,
    query_word: &'a [u8],
    subject_word: [u8; 32],
    alphabet_size: usize,
    wordsize: usize,
    charsize: usize,
    row_max: &'a [i32],
    offset_list: &'a [i32],
    threshold: i32,
    query_bias: i32,
    neighbor_added_count: usize,
    neighbor_words_generated: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:96-97
// ```c
// BlastAaLookupIndexQuery((BlastAaLookupTable*) lookup_wrap->lut, matrix,
//                         query, lookup_segments, 0);
// ```
struct PreparedLookupQuery {
    concat_query: Vec<u8>,
    lookup_locations: Vec<(i32, i32)>,
    frame_bases: Vec<i32>,
    contexts: Vec<QueryContext>,
    skipped_seg_mask: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1032-1114
// ```c
// start_offset = query_info->contexts[context].query_offset;
// end_offset = query_info->contexts[context].query_length + start_offset - 1;
// ...
// if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
//     tail = BlastSeqLocNew(tail ? &tail : complement_mask,
//                           start_offset, end_offset);
//     continue;
// }
// ...
// filter_start = start_offset + seq_range->left;
// filter_end = start_offset + seq_range->right;
// ...
// tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
// ```
fn build_lookup_segments_for_context(
    query_offset: i32,
    query_length: usize,
    seg_masks: &[(usize, usize)],
) -> Vec<(i32, i32)> {
    if query_length == 0 {
        return Vec::new();
    }

    let start_offset = query_offset;
    let end_offset = query_offset
        + i32::try_from(query_length).expect("NCBI BLAST query length must fit in Int4")
        - 1;

    if seg_masks.is_empty() {
        return vec![(start_offset, end_offset)];
    }

    let mut sorted_masks = seg_masks.to_vec();
    sorted_masks.sort_by_key(|&(start, _)| start);

    let mut lookup_segments = Vec::new();
    let mut first = true;
    let mut last_interval_open = true;
    let mut left = 0i32;

    for (start, end) in sorted_masks {
        let filter_start =
            start_offset + i32::try_from(start).expect("NCBI BLAST mask start must fit in Int4");
        let filter_end =
            start_offset + i32::try_from(end).expect("NCBI BLAST mask end must fit in Int4") - 1;

        if first {
            last_interval_open = true;
            first = false;

            if filter_start > start_offset {
                left = start_offset;
            } else {
                left = filter_end + 1;
                continue;
            }
        }

        let right = filter_start - 1;
        if left <= right {
            lookup_segments.push((left, right));
        }
        if filter_end >= end_offset {
            last_interval_open = false;
            break;
        }
        left = filter_end + 1;
    }

    if last_interval_open && left <= end_offset {
        lookup_segments.push((left, end_offset));
    }

    lookup_segments
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-100
// ```c
// if (word_length > to - from + 1)
//     continue;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:102-128
// ```c
// seq = query->sequence + from;
// word_target = seq + lut_word_length;
// for (offset = from; offset <= to; offset++, seq++) { ... }
// ```
#[inline]
fn count_lookup_positions(lookup_segments: &[(i32, i32)], word_length: usize) -> usize {
    lookup_segments
        .iter()
        .map(|&(left, right)| {
            let seg_len = usize::try_from(right - left + 1)
                .expect("NCBI BLAST lookup segment length must be non-negative");
            seg_len.saturating_sub(word_length - 1)
        })
        .sum()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-314
// ```c
// info->contexts[i].query_offset = new_offsets[i];
// distance = new_offsets[i+1] - new_offsets[i];
// info->contexts[i].query_length = distance ? distance-1 : 0;
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:621-636
// ```c
// if (!mask_at_hash) {
//     BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc, program_number);
// }
// if (lookup_segments) {
//     BLAST_ComplementMaskLocations(program_number, query_info,
//                                   filter_maskloc, lookup_segments);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
// ```c
// /* Increment offset by 1 extra byte for the sentinel NULLB
//    between frames. */
// offset += length + 1;
// frame_offsets[context+1] = offset;
// ```
fn prepare_lookup_query(
    queries: &[Vec<QueryFrame>],
    ideal_params: KarlinParams,
    std_comp: &[f64; LOOKUP_ALPHABET_SIZE],
    word_length: usize,
) -> PreparedLookupQuery {
    let mut prepared = PreparedLookupQuery {
        concat_query: Vec::new(),
        lookup_locations: Vec::new(),
        frame_bases: Vec::new(),
        contexts: Vec::new(),
        skipped_seg_mask: 0,
    };
    let mut query_offset = 0i32;

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            prepared.frame_bases.push(query_offset);
            prepared.concat_query.extend_from_slice(&frame.aa_seq[1..]);

            let lookup_segments =
                build_lookup_segments_for_context(query_offset, frame.aa_len, &frame.seg_masks);
            let total_aa_positions = frame.aa_len.saturating_sub(word_length - 1);
            let lookup_positions = count_lookup_positions(&lookup_segments, word_length);
            prepared.skipped_seg_mask += total_aa_positions.saturating_sub(lookup_positions);
            prepared
                .lookup_locations
                .extend_from_slice(&lookup_segments);

            let ctx_comp = compute_aa_composition(&frame.aa_seq, frame.aa_len);
            let score_min = -4;
            let score_max = 11;
            let sfp = compute_score_freq_profile(&ctx_comp, std_comp, score_min, score_max);
            let computed_params = compute_karlin_params_ungapped(&sfp).unwrap_or(ideal_params);
            let final_params = apply_check_ideal(computed_params, ideal_params);

            prepared.contexts.push(QueryContext {
                q_idx: q_idx as u32,
                f_idx: f_idx as u8,
                frame: frame.frame,
                aa_seq: frame.aa_seq.clone(),
                aa_seq_nomask: frame.aa_seq_nomask.clone(),
                aa_len: frame.aa_len,
                orig_len: frame.orig_len,
                frame_base: query_offset,
                karlin_params: final_params,
            });

            query_offset += frame.aa_seq.len() as i32 - 1;
        }
    }

    prepared
}

#[inline(always)]
fn lookup_chain_num_used(chain: &[i32]) -> usize {
    usize::try_from(chain[1]).expect("NCBI BLAST lookup chain hit count must fit in usize")
}

#[inline(always)]
fn lookup_chain_entries(chain: &[i32]) -> &[i32] {
    let num_used = lookup_chain_num_used(chain);
    &chain[2..2 + num_used]
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:33-77
// ```c
// index = ComputeTableIndex(wordsize, charsize, seq);
// if (backbone[index] == NULL) {
//     chain_size = 8;
//     hits_in_chain = 0;
//     chain = (Int4 *) malloc(chain_size * sizeof(Int4));
//     chain[0] = chain_size;
//     chain[1] = hits_in_chain;
//     backbone[index] = chain;
// }
// if ((hits_in_chain + 2) == chain_size) {
//     chain_size = chain_size * 2;
//     chain = (Int4 *) realloc(chain, chain_size * sizeof(Int4));
//     backbone[index] = chain;
//     chain[0] = chain_size;
// }
// chain[chain[1] + 2] = query_offset;
// chain[1]++;
// ```
#[inline]
fn blast_lookup_add_word_hit(
    backbone: &mut LookupBackboneChains,
    wordsize: usize,
    charsize: usize,
    seq: &[u8],
    query_offset: i32,
) {
    let mut index = 0usize;
    for &residue in seq.iter().take(wordsize) {
        index = (index << charsize) | residue as usize;
    }

    let chain = backbone[index].get_or_insert_with(|| {
        let mut chain = vec![0; 8];
        chain[0] = 8;
        chain[1] = 0;
        chain
    });

    let hits_in_chain = lookup_chain_num_used(chain);
    let chain_size = usize::try_from(chain[0]).expect("NCBI BLAST lookup chain size must fit");
    if (hits_in_chain + 2) == chain_size {
        let new_size = chain_size * 2;
        chain.resize(new_size, 0);
        chain[0] = i32::try_from(new_size).expect("NCBI BLAST lookup chain size must fit Int4");
    }

    chain[hits_in_chain + 2] = query_offset;
    chain[1] += 1;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:79-132
// ```c
// for (loc = locations; loc; loc = loc->next) {
//     Int4 from = loc->ssr->left;
//     Int4 to = loc->ssr->right;
//     if (word_length > to - from + 1)
//         continue;
//     seq = query->sequence + from;
//     word_target = seq + lut_word_length;
//     for (offset = from; offset <= to; offset++, seq++) {
//         if (seq >= word_target) {
//             BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                                   seq - lut_word_length, offset - lut_word_length);
//         }
//         if (*seq & invalid_mask)
//             word_target = seq + lut_word_length + 1;
//     }
//     if (seq >= word_target) {
//         BlastLookupAddWordHit(backbone, lut_word_length, charsize,
//                               seq - lut_word_length, offset - lut_word_length);
//     }
// }
// ```
fn blast_lookup_index_query_exact_matches(
    backbone: &mut LookupBackboneChains,
    word_length: i32,
    charsize: i32,
    lut_word_length: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    skipped_invalid_residue: &mut usize,
) -> usize {
    let invalid_mask: u8 = 0xffu8 << (charsize as u32);
    let mut total_exact_positions = 0usize;

    for &(from, to) in locations {
        if word_length > to - from + 1 {
            continue;
        }

        let mut offset = from;
        let mut seq = usize::try_from(from).expect("NCBI BLAST query offset must fit usize");
        let mut word_target = seq + usize::try_from(lut_word_length).expect("word length fits");

        while offset <= to {
            if seq >= word_target {
                let word_start = seq - usize::try_from(lut_word_length).expect("word length fits");
                blast_lookup_add_word_hit(
                    backbone,
                    usize::try_from(lut_word_length).expect("word length fits"),
                    usize::try_from(charsize).expect("charsize fits"),
                    &query[word_start..],
                    offset - lut_word_length,
                );
                total_exact_positions += 1;
            }

            if (query[seq] & invalid_mask) != 0 {
                word_target = seq + usize::try_from(lut_word_length).expect("word length fits") + 1;
                *skipped_invalid_residue += 1;
            }
            offset += 1;
            seq += 1;
        }

        if seq >= word_target {
            let word_start = seq - usize::try_from(lut_word_length).expect("word length fits");
            blast_lookup_add_word_hit(
                backbone,
                usize::try_from(lut_word_length).expect("word length fits"),
                usize::try_from(charsize).expect("charsize fits"),
                &query[word_start..],
                offset - lut_word_length,
            );
            total_exact_positions += 1;
        }
    }

    total_exact_positions
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:546-606
// ```c
// score -= info->row_max[query_word[current_pos]];
// row = info->matrix[query_word[current_pos]];
// if (current_pos == info->wordsize - 1) {
//     for (i = 0; i < alphabet_size; i++) {
//         if (score + row[i] >= threshold) {
//             subject_word[current_pos] = i;
//             for (j = 0; j < offset_list[1]; j++) {
//                 BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
//                                       charsize, subject_word,
//                                       query_bias + offset_list[j + 2]);
//             }
//         }
//     }
//     return;
// }
// for (i = 0; i < alphabet_size; i++) {
//     if (score + row[i] >= threshold) {
//         subject_word[current_pos] = i;
//         s_AddWordHitsCore(info, score + row[i], current_pos + 1);
//     }
// }
// ```
fn blast_aa_add_word_hits_core(info: &mut NeighborInfo<'_>, score: i32, current_pos: usize) {
    let query_residue = info.query_word[current_pos] as usize;
    let score = score - info.row_max[query_residue];

    if current_pos == info.wordsize - 1 {
        for residue in 0..info.alphabet_size {
            let residue_score = blosum62_score(query_residue as u8, residue as u8);
            if score + residue_score >= info.threshold {
                info.subject_word[current_pos] = residue as u8;
                for &query_offset in lookup_chain_entries(info.offset_list) {
                    blast_lookup_add_word_hit(
                        info.thin_backbone,
                        info.wordsize,
                        info.charsize,
                        &info.subject_word,
                        info.query_bias + query_offset,
                    );
                    info.neighbor_added_count += 1;
                }
                info.neighbor_words_generated += 1;
            }
        }
        return;
    }

    for residue in 0..info.alphabet_size {
        let residue_score = blosum62_score(query_residue as u8, residue as u8);
        if score + residue_score >= info.threshold {
            info.subject_word[current_pos] = residue as u8;
            blast_aa_add_word_hits_core(info, score + residue_score, current_pos + 1);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:472-544
// ```c
// w = query + offset_list[2];
// score = matrix[w[0]][w[0]];
// for (i = 1; i < lookup->word_length; i++)
//     score += matrix[w[i]][w[i]];
// if (lookup->threshold == 0 || score < lookup->threshold) {
//     for (i = 0; i < offset_list[1]; i++) {
//         BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
//                               lookup->charsize, w,
//                               query_bias + offset_list[i + 2]);
//     }
// }
// if (lookup->threshold == 0)
//     return;
// score = row_max[w[0]];
// for (i = 1; i < lookup->word_length; i++)
//     score += row_max[w[i]];
// s_AddWordHitsCore(&info, score, 0);
// ```
fn blast_aa_add_word_hits(
    thin_backbone: &mut LookupBackboneChains,
    alphabet_size: usize,
    word_length: usize,
    charsize: usize,
    threshold: i32,
    query: &[u8],
    offset_list: &[i32],
    query_bias: i32,
    row_max: &[i32],
) -> (usize, usize, usize) {
    let query_offset = usize::try_from(offset_list[2]).expect("NCBI BLAST query offset must fit");
    let query_word = &query[query_offset..];
    let mut self_score = blosum62_score(query_word[0], query_word[0]);
    for residue in query_word.iter().take(word_length).skip(1) {
        self_score += blosum62_score(*residue, *residue);
    }

    let mut exact_added_count = 0usize;
    if threshold == 0 || self_score < threshold {
        for &query_offset in lookup_chain_entries(offset_list) {
            blast_lookup_add_word_hit(
                thin_backbone,
                word_length,
                charsize,
                query_word,
                query_bias + query_offset,
            );
            exact_added_count += 1;
        }
    }
    if threshold == 0 {
        return (exact_added_count, 0, 0);
    }

    let mut info = NeighborInfo {
        thin_backbone,
        query_word,
        subject_word: [0; 32],
        alphabet_size,
        wordsize: word_length,
        charsize,
        row_max,
        offset_list,
        threshold,
        query_bias,
        neighbor_added_count: 0,
        neighbor_words_generated: 0,
    };
    let mut max_score = row_max[query_word[0] as usize];
    for residue in query_word.iter().take(word_length).skip(1) {
        max_score += row_max[*residue as usize];
    }
    blast_aa_add_word_hits_core(&mut info, max_score, 0);

    (
        exact_added_count,
        info.neighbor_added_count,
        info.neighbor_words_generated,
    )
}

/// Build NCBI-style lookup table with presence vector and overflow.
///
/// Reference: blast_aalookup.c BlastAaLookupIndexQuery + BlastAaLookupFinalize
///
/// Neighbor words are always precomputed in the lookup table (no lazy mode).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
///
/// NCBI reference (fixed alphabet and neighbor precomputation):
/// ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:436-456
/// ```c
/// ASSERT(lookup->alphabet_size <= BLASTAA_SIZE);
/// ...
/// BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
///                                   lookup->charsize, lookup->word_length,
///                                   query, location);
/// ```
pub fn build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    _karlin_params: &crate::stats::KarlinParams, // Unused - computed per context, kept for API compatibility
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let diag_enabled = diagnostics_enabled();
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE; // 28
    let charsize = ilog2(alphabet_size) + 1; // 5
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    // Compute ideal Karlin parameters (kbp_ideal) - used for check_ideal logic
    // Reference: NCBI blast_stat.c:2754 Blast_ScoreBlkKbpIdealCalc
    let ideal_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    // Compute standard amino acid composition (for database/subject)
    // Reference: NCBI blast_stat.c:2759 Blast_ResFreqStdComp
    let std_comp = compute_std_aa_composition();

    let mut build_stats = LookupBuildStats::default();
    let prepared = prepare_lookup_query(queries, ideal_params, &std_comp, word_length);
    build_stats.skipped_seg_mask = prepared.skipped_seg_mask;
    let PreparedLookupQuery {
        concat_query,
        lookup_locations,
        frame_bases,
        contexts,
        skipped_seg_mask: _,
    } = prepared;
    let query_length = i32::try_from(concat_query.len())
        .expect("NCBI BLAST concatenated query length must fit in Int4");

    // Row max for BLOSUM62 over the lookup alphabet.
    // For NCBISTDAA residues, we compute max score against any other residue.
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| {
            (0..alphabet_size)
                .map(|j| blosum62_score(i as u8, j as u8))
                .max()
                .unwrap_or(-4)
        })
        .collect();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-469
    // ```c
    // exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
    // BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                   lookup->charsize, lookup->word_length,
    //                                   query, location);
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (exact_backbone[i] != NULL) {
    //         s_AddWordHits(lookup, matrix, query->sequence,
    //                       exact_backbone[i], query_bias, row_max);
    //     }
    // }
    // ```
    let mut exact_backbone: LookupBackboneChains = vec![None; backbone_size];
    build_stats.total_exact_positions = blast_lookup_index_query_exact_matches(
        &mut exact_backbone,
        word_length as i32,
        charsize as i32,
        word_length as i32,
        &concat_query,
        &lookup_locations,
        &mut build_stats.skipped_invalid_residue,
    );

    // Phase 1a diagnostics (guarded)
    if diag_enabled {
        let unique_exact_words = exact_backbone
            .iter()
            .filter(|chain| chain.is_some())
            .count();
        let max_exact_per_word = exact_backbone
            .iter()
            .filter_map(|chain| chain.as_deref())
            .map(lookup_chain_num_used)
            .max()
            .unwrap_or(0);
        eprintln!("\n=== Phase 1a: Exact Indexing Diagnostics ===");
        eprintln!(
            "Total exact positions indexed: {}",
            build_stats.total_exact_positions
        );
        eprintln!("Unique exact words: {}", unique_exact_words);
        eprintln!("Max offsets per exact word: {}", max_exact_per_word);
        eprintln!(
            "Skipped (invalid residue): {}",
            build_stats.skipped_invalid_residue
        );
        eprintln!("Skipped (SEG mask): {}", build_stats.skipped_seg_mask);
    }

    // Phase 1b: Add neighboring words (NCBI precomputed neighbors).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
    let residue_mask: usize = (1usize << charsize) - 1;
    let mut thin_backbone: LookupBackboneChains = vec![None; backbone_size];
    let mut counts: Vec<u32> = vec![0; backbone_size];

    for idx in 0..backbone_size {
        let Some(offset_list) = exact_backbone[idx].as_deref() else {
            continue;
        };
        let (exact_added_count, neighbor_added_count, neighbor_words_generated) =
            blast_aa_add_word_hits(
                &mut thin_backbone,
                alphabet_size,
                word_length,
                charsize,
                threshold,
                &concat_query,
                offset_list,
                0,
                &row_max,
            );
        build_stats.exact_added_count += exact_added_count;
        build_stats.neighbor_added_count += neighbor_added_count;
        build_stats.neighbor_words_generated += neighbor_words_generated;
        if threshold == 0 {
            build_stats.words_with_exact_only += 1;
        }
        if neighbor_words_generated > 0 {
            build_stats.words_with_neighbors += 1;
            if neighbor_words_generated > build_stats.max_neighbors_for_single_word {
                build_stats.max_neighbors_for_single_word = neighbor_words_generated;
                build_stats.max_neighbor_word_idx = idx;
            }
        }
    }
    for idx in 0..backbone_size {
        counts[idx] = thin_backbone[idx]
            .as_deref()
            .map(lookup_chain_num_used)
            .unwrap_or(0) as u32;
    }

    // Phase 1b diagnostics (guarded)
    if diag_enabled {
        eprintln!("\n=== Phase 1b: Neighbor Generation Diagnostics ===");
        eprintln!("Threshold: {}", threshold);
        eprintln!(
            "Exact entries added (self_score < threshold): {}",
            build_stats.exact_added_count
        );
        eprintln!(
            "Neighbor entries added: {}",
            build_stats.neighbor_added_count
        );
        eprintln!(
            "Neighbor words generated (unique): {}",
            build_stats.neighbor_words_generated
        );
        eprintln!(
            "Words with exact-only additions: {}",
            build_stats.words_with_exact_only
        );
        eprintln!(
            "Words with neighbor generation: {}",
            build_stats.words_with_neighbors
        );
        eprintln!(
            "Max neighbors generated for single word: {}",
            build_stats.max_neighbors_for_single_word
        );
        if build_stats.max_neighbors_for_single_word > 0 {
            let w0 = (build_stats.max_neighbor_word_idx >> (2 * charsize)) & residue_mask;
            let w1 = (build_stats.max_neighbor_word_idx >> charsize) & residue_mask;
            let w2 = build_stats.max_neighbor_word_idx & residue_mask;
            eprintln!(
                "  Word with max neighbors: idx={} (residues {},{},{})",
                build_stats.max_neighbor_word_idx, w0, w1, w2
            );
        }
        let expansion_factor = if build_stats.total_exact_positions > 0 {
            (build_stats.exact_added_count + build_stats.neighbor_added_count) as f64
                / build_stats.total_exact_positions as f64
        } else {
            0.0
        };
        eprintln!(
            "Expansion factor (total_entries / exact_positions): {:.2}x",
            expansion_factor
        );
    }

    // Phase 2: Finalize
    let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
    // NCBI allocates pv as: (backbone_size >> PV_ARRAY_BTS) + 1.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:327-329
    let pv_size = (backbone_size >> PV_ARRAY_BTS) + 1;
    let mut pv: Vec<u32> = vec![0u32; pv_size];

    // Pre-finalize diagnostics: analyze high-frequency k-mers (guarded)
    if diag_enabled {
        eprintln!("\n=== Phase 2: High-Frequency K-mer Analysis ===");

        // Collect (count, idx) pairs and sort by count descending
        let mut count_idx_pairs: Vec<(usize, usize)> = counts
            .iter()
            .enumerate()
            .filter(|(_, &c)| c > 0)
            .map(|(i, &c)| (c as usize, i))
            .collect();
        count_idx_pairs.sort_by(|a, b| b.0.cmp(&a.0));

        // Show top 10 high-frequency k-mers
        eprintln!("Top 10 high-frequency k-mers (before suppression):");
        for (rank, &(count, idx)) in count_idx_pairs.iter().take(10).enumerate() {
            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            eprintln!(
                "  #{}: count={}, idx={} (residues {},{},{})",
                rank + 1,
                count,
                idx,
                w0,
                w1,
                w2
            );
        }

        // Distribution histogram
        let hist_buckets = [
            1,
            10,
            100,
            500,
            1000,
            5000,
            10000,
            50000,
            100000,
            usize::MAX,
        ];
        let mut hist: Vec<usize> = vec![0; hist_buckets.len()];
        let mut hist_entries: Vec<usize> = vec![0; hist_buckets.len()];
        for &(count, _) in &count_idx_pairs {
            for (i, &bucket) in hist_buckets.iter().enumerate() {
                if count <= bucket {
                    hist[i] += 1;
                    hist_entries[i] += count;
                    break;
                }
            }
        }
        eprintln!("Hit count distribution:");
        let mut prev = 0;
        for (i, &bucket) in hist_buckets.iter().enumerate() {
            if hist[i] > 0 {
                let label = if bucket == usize::MAX {
                    format!(">{}", prev)
                } else {
                    format!("{}-{}", prev + 1, bucket)
                };
                eprintln!(
                    "  {}: {} cells, {} entries",
                    label, hist[i], hist_entries[i]
                );
            }
            prev = bucket;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:282-299
    // ```c
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (lookup->thin_backbone[i]) {
    //         if (lookup->thin_backbone[i][1] > AA_HITS_PER_CELL)
    //             overflow_cells_needed += lookup->thin_backbone[i][1];
    //         if (lookup->thin_backbone[i][1] > longest_chain)
    //             longest_chain = lookup->thin_backbone[i][1];
    //     }
    // }
    // ```
    let mut overflow_size = 0usize;
    for &count in &counts {
        let count = count as usize;
        if count > AA_HITS_PER_CELL {
            overflow_size += count;
        }
    }

    let mut overflow: Vec<i32> = vec![0; overflow_size];
    let mut overflow_cursor = 0usize;

    let mut longest_chain: i32 = 0;
    for &count in &counts {
        if (count as i32) > longest_chain {
            longest_chain = count as i32;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:337-360
    // ```c
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     if (lookup->thin_backbone[i]) {
    //         PV_SET(pv, i, PV_ARRAY_BTS);
    //         bbc[i].num_used = lookup->thin_backbone[i][1];
    //         if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
    //             dest = bbc[i].payload.entries;
    //         else {
    //             bbc[i].payload.overflow_cursor = overflow_cursor;
    //             dest = (Int4 *) lookup->overflow + overflow_cursor;
    //             overflow_cursor += lookup->thin_backbone[i][1];
    //         }
    //         for (j = 0; j < lookup->thin_backbone[i][1]; j++)
    //             dest[j] = lookup->thin_backbone[i][j + 2];
    //     }
    // }
    // ```
    for idx in 0..backbone_size {
        let Some(chain) = thin_backbone[idx].as_deref() else {
            continue;
        };
        let entries = lookup_chain_entries(chain);
        let count = entries.len();
        pv_set(&mut pv, idx);
        backbone[idx].num_used = count as i32;

        if count <= AA_HITS_PER_CELL {
            for (i, &off) in entries.iter().enumerate() {
                backbone[idx].entries[i] = off;
            }
        } else {
            backbone[idx].entries[0] = overflow_cursor as i32;
            for &off in entries {
                overflow[overflow_cursor] = off;
                overflow_cursor += 1;
            }
        }
    }

    let nonempty = backbone.iter().filter(|c| c.num_used > 0).count();
    let total: usize = backbone.iter().map(|c| c.num_used.max(0) as usize).sum();
    if diag_enabled {
        eprintln!(
            "\n=== NCBI Backbone Lookup Table (BLASTAA_SIZE={}) ===",
            LOOKUP_ALPHABET_SIZE
        );
        eprintln!("Contexts: {}", contexts.len());
        eprintln!(
            "Backbone size: {} cells ({:.1} MB)",
            backbone_size,
            (backbone_size * std::mem::size_of::<BackboneCell>()) as f64 / 1e6
        );
        eprintln!("Non-empty cells: {}", nonempty);
        eprintln!("Total entries: {}", total);
        eprintln!(
            "Overflow size: {} ({:.1} MB)",
            overflow_cursor,
            (overflow_cursor * 4) as f64 / 1e6
        );
        eprintln!("Longest chain: {}", longest_chain);
        eprintln!("==========================================\n");
    }

    (
        BlastAaLookupTable {
            backbone,
            overflow,
            pv,
            frame_bases,
            num_contexts: contexts.len(),
            query_length,
            word_length: word_length as i32,
            alphabet_size: alphabet_size as i32,
            charsize: charsize as i32,
            mask: mask as i32,
            longest_chain,
            threshold,
            row_max,
        },
        contexts,
    )
}

/// Build a simple direct (exact) 3-mer lookup table for tests and diagnostics.
pub fn build_direct_lookup(queries: &[Vec<QueryFrame>]) -> Vec<Vec<(u32, u8, u32)>> {
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    let mut table: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); backbone_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            if frame.aa_len < 3 || frame.aa_seq.len() < 5 {
                continue;
            }

            let seq_ptr = &frame.aa_seq[1..frame.aa_seq.len() - 1];
            // NCBI invalid_mask = 0xff << charsize (Uint1).
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:89-120
            let invalid_mask: u8 = 0xffu8 << (charsize as u32);
            let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
            for (start, end) in unmasked {
                let from = start as i32;
                let to = end.saturating_sub(1) as i32;
                // NCBI: if (word_length > to - from + 1) continue;
                // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-100
                if word_length as i32 > to - from + 1 {
                    continue;
                }
                let mut offset = from;
                let mut word_target = from + word_length as i32;

                while offset <= to {
                    if offset >= word_target {
                        let word_start = (offset - word_length as i32) as usize;
                        let idx = encode_kmer_3(
                            seq_ptr[word_start] as usize,
                            seq_ptr[word_start + 1] as usize,
                            seq_ptr[word_start + 2] as usize,
                            charsize,
                        ) & mask;
                        table[idx].push((q_idx as u32, f_idx as u8, word_start as u32));
                    }
                    if (seq_ptr[offset as usize] & invalid_mask) != 0 {
                        // NCBI: word_target = seq + lut_word_length + 1
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                        word_target = offset + word_length as i32 + 1;
                    }
                    offset += 1;
                }

                if offset >= word_target {
                    let word_start = (offset - word_length as i32) as usize;
                    let idx = encode_kmer_3(
                        seq_ptr[word_start] as usize,
                        seq_ptr[word_start + 1] as usize,
                        seq_ptr[word_start + 2] as usize,
                        charsize,
                    ) & mask;
                    table[idx].push((q_idx as u32, f_idx as u8, word_start as u32));
                }
            }
        }
    }

    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;
    use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE};

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

    fn make_query_frame(seq: Vec<u8>, seg_masks: Vec<(usize, usize)>) -> QueryFrame {
        let mut aa_seq = Vec::with_capacity(seq.len() + 2);
        aa_seq.push(0);
        aa_seq.extend_from_slice(&seq);
        aa_seq.push(0);
        QueryFrame {
            frame: 1,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: seq.len(),
            orig_len: seq.len(),
            seg_masks,
        }
    }

    fn build_lookup_for_test(queries: Vec<Vec<QueryFrame>>, threshold: i32) -> BlastAaLookupTable {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, _) = build_ncbi_lookup(&queries, threshold, &ungapped);
        lookup
    }

    fn prepare_lookup_for_test(queries: Vec<Vec<QueryFrame>>) -> PreparedLookupQuery {
        let ideal_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let std_comp = compute_std_aa_composition();
        prepare_lookup_query(&queries, ideal_params, &std_comp, LOOKUP_WORD_LENGTH)
    }

    fn debruijn(alphabet_size: usize, order: usize) -> Vec<u8> {
        fn db(
            t: usize,
            p: usize,
            alphabet_size: usize,
            order: usize,
            a: &mut [usize],
            output: &mut Vec<u8>,
        ) {
            if t > order {
                if order % p == 0 {
                    for j in 1..=p {
                        output.push(a[j] as u8);
                    }
                }
                return;
            }

            a[t] = a[t - p];
            db(t + 1, p, alphabet_size, order, a, output);
            for residue in (a[t - p] + 1)..alphabet_size {
                a[t] = residue;
                db(t + 1, t, alphabet_size, order, a, output);
            }
        }

        let mut a = vec![0usize; alphabet_size * order + 1];
        let mut output = Vec::new();
        db(1, 1, alphabet_size, order, &mut a, &mut output);
        let wraparound = output
            .iter()
            .copied()
            .take(order.saturating_sub(1))
            .collect::<Vec<_>>();
        output.extend(wraparound);
        output
    }

    fn brute_force_lookup_count(subject_word: [u8; 3], threshold: i32) -> usize {
        let mut count = 0usize;
        for q0 in 0..BLASTAA_SIZE {
            for q1 in 0..BLASTAA_SIZE {
                for q2 in 0..BLASTAA_SIZE {
                    let score = blosum62_score(q0 as u8, subject_word[0])
                        + blosum62_score(q1 as u8, subject_word[1])
                        + blosum62_score(q2 as u8, subject_word[2]);
                    if score >= threshold
                        || (q0 as u8 == subject_word[0]
                            && q1 as u8 == subject_word[1]
                            && q2 as u8 == subject_word[2])
                    {
                        count += 1;
                    }
                }
            }
        }
        count
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-314
    // ```c
    // info->contexts[i].query_offset = new_offsets[i];
    // distance = new_offsets[i+1] - new_offsets[i];
    // info->contexts[i].query_length = distance ? distance-1 : 0;
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
    // ```c
    // offset += length + 1;
    // frame_offsets[context+1] = offset;
    // ```
    #[test]
    fn test_prepare_lookup_query_shares_single_boundary_sentinel_and_context_offsets() {
        let prepared = prepare_lookup_for_test(vec![
            vec![make_query_frame(
                vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                Vec::new(),
            )],
            vec![make_query_frame(
                vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                Vec::new(),
            )],
        ]);

        assert_eq!(
            prepared.concat_query,
            vec![
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                0,
                ncbistdaa::C,
                ncbistdaa::C,
                ncbistdaa::C,
                0,
            ]
        );
        assert_eq!(prepared.frame_bases, vec![0, 4]);
        assert_eq!(prepared.lookup_locations, vec![(0, 2), (4, 6)]);
        assert_eq!(prepared.contexts[0].frame_base, 0);
        assert_eq!(prepared.contexts[1].frame_base, 4);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1045-1056
    // ```c
    // start_offset = query_info->contexts[context].query_offset;
    // end_offset = query_info->contexts[context].query_length + start_offset - 1;
    // if (mask_loc == NULL || mask_loc->seqloc_array[context] == NULL) {
    //     tail = BlastSeqLocNew(tail ? &tail : complement_mask,
    //                           start_offset, end_offset);
    //     continue;
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1067-1114
    // ```c
    // filter_start = start_offset + seq_range->left;
    // filter_end = start_offset + seq_range->right;
    // ...
    // tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
    // ```
    #[test]
    fn test_prepare_lookup_query_converts_seg_masks_to_absolute_inclusive_segments() {
        let prepared = prepare_lookup_for_test(vec![
            vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                vec![(2, 4)],
            )],
            vec![make_query_frame(
                vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                Vec::new(),
            )],
        ]);

        assert_eq!(prepared.lookup_locations, vec![(0, 1), (4, 5), (7, 9)]);
        assert_eq!(prepared.skipped_seg_mask, 4);
    }

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
    #[test]
    fn test_build_ncbi_lookup_query_length_matches_ncbi_sequence_length() {
        let lookup = build_lookup_for_test(
            vec![
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
                vec![make_query_frame(
                    vec![ncbistdaa::C, ncbistdaa::C, ncbistdaa::C],
                    Vec::new(),
                )],
            ],
            0,
        );

        assert_eq!(lookup.query_length, 8);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:102-128
    // ```c
    // seq = query->sequence + from;
    // word_target = seq + lut_word_length;
    // ...
    // BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                       seq - lut_word_length, offset - lut_word_length);
    // ```
    #[test]
    fn test_build_ncbi_lookup_multiple_queries_do_not_cross_shared_boundary_sentinel() {
        let lookup = build_lookup_for_test(
            vec![
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
                vec![make_query_frame(
                    vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
                    Vec::new(),
                )],
            ],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
        assert_eq!(lookup.get_context_idx(0), 0);
        assert_eq!(lookup.get_context_idx(4), 1);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:340-357
    // ```c
    // if (lookup->thin_backbone[i][1] <= AA_HITS_PER_CELL)
    //     dest = bbc[i].payload.entries;
    // else {
    //     bbc[i].payload.overflow_cursor = overflow_cursor;
    //     dest = (Int4 *) lookup->overflow + overflow_cursor;
    // }
    // for (j = 0; j < lookup->thin_backbone[i][1]; j++)
    //     dest[j] = lookup->thin_backbone[i][j + 2];
    // ```
    #[test]
    fn test_build_ncbi_lookup_preserves_overflow_query_offset_append_order() {
        let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 0);
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.backbone[aaa_index].num_used, 4);
        assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
        assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:504-543
    // ```c
    // if (lookup->threshold == 0 || score < lookup->threshold) {
    //     for (i = 0; i < offset_list[1]; i++) {
    //         BlastLookupAddWordHit(lookup->thin_backbone, lookup->word_length,
    //                               lookup->charsize, w,
    //                               query_bias + offset_list[i + 2]);
    //     }
    // }
    // ...
    // for (j = 0; j < offset_list[1]; j++) {
    //     BlastLookupAddWordHit(lookup->thin_backbone, wordsize,
    //                           charsize, subject_word,
    //                           query_bias + offset_list[j + 2]);
    // }
    // ```
    #[test]
    fn test_build_ncbi_lookup_preserves_neighbor_query_offset_append_order() {
        let lookup = build_lookup_for_test(vec![vec![make_poly_a_query_frame()]], 11);
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.backbone[aaa_index].num_used, 4);
        assert_eq!(&lookup.overflow[..4], &[0, 1, 2, 3]);
        assert_eq!(lookup.get_hits(aaa_index), &[0, 1, 2, 3]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:96-128
    // ```c
    // if (word_length > to - from + 1)
    //     continue;
    // ...
    // if (*seq & invalid_mask)
    //     word_target = seq + lut_word_length + 1;
    // ```
    #[test]
    fn test_build_ncbi_lookup_skips_words_crossing_invalid_residue() {
        let lookup = build_lookup_for_test(
            vec![vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    32,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                Vec::new(),
            )]],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0, 4]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:92-100
    // ```c
    // for (loc = locations; loc; loc = loc->next) {
    //     Int4 from = loc->ssr->left;
    //     Int4 to = loc->ssr->right;
    //     if (word_length > to - from + 1)
    //         continue;
    // ```
    #[test]
    fn test_build_ncbi_lookup_respects_interval_boundaries() {
        let lookup = build_lookup_for_test(
            vec![vec![make_query_frame(
                vec![
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                    ncbistdaa::A,
                ],
                vec![(3, 4)],
            )]],
            0,
        );
        let aaa_index = encode_kmer_3(
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            ncbistdaa::A as usize,
            lookup.charsize as usize,
        ) & lookup.mask as usize;

        assert_eq!(lookup.get_hits(aaa_index), &[0]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:218-238
    // ```c
    // for(i=0;i<lookup->backbone_size;i++) {
    //   Int4 num_used = ((AaLookupSmallboneCell *)(lookup->thick_backbone))[i].num_used;
    //   if (((i & 0x1F) >= BLASTAA_SIZE) || (((i & 0x3E0) >> 5) >= BLASTAA_SIZE)
    //       || ((((i & 0x7C00) >> 10) >= BLASTAA_SIZE)))
    //       BOOST_REQUIRE_EQUAL(num_used, 0);
    //   else
    //       BOOST_REQUIRE_EQUAL(num_used, 1);
    // }
    // ```
    #[test]
    fn test_build_ncbi_lookup_debruijn_exact_words_match_ncbi_counts() {
        let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
        let lookup =
            build_lookup_for_test(vec![vec![make_query_frame(seq.clone(), Vec::new())]], 0);
        let charsize = lookup.charsize as usize;
        let residue_mask = (1usize << charsize) - 1;

        for idx in 0..lookup.backbone.len() {
            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            if w0 < BLASTAA_SIZE && w1 < BLASTAA_SIZE && w2 < BLASTAA_SIZE {
                assert_eq!(lookup.backbone[idx].num_used, 1, "valid idx={idx}");
            } else {
                assert_eq!(lookup.backbone[idx].num_used, 0, "invalid idx={idx}");
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:253-277
    // ```c
    // for each possible 3-mer:
    //   find its neighbors by brute force
    //   if (score >= 11) || exact-match
    //       count++;
    //   ensure that the number of neighbors matches the lookup table
    // ```
    #[test]
    fn test_build_ncbi_lookup_representative_neighbor_counts_match_bruteforce() {
        let seq = debruijn(BLASTAA_SIZE, LOOKUP_WORD_LENGTH);
        let lookup = build_lookup_for_test(vec![vec![make_query_frame(seq, Vec::new())]], 11);
        let representative_words = [
            [ncbistdaa::A, ncbistdaa::A, ncbistdaa::A],
            [ncbistdaa::W, ncbistdaa::W, ncbistdaa::W],
            [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N],
        ];

        for subject_word in representative_words {
            let index = encode_kmer_3(
                subject_word[0] as usize,
                subject_word[1] as usize,
                subject_word[2] as usize,
                lookup.charsize as usize,
            ) & lookup.mask as usize;
            assert_eq!(
                lookup.backbone[index].num_used as usize,
                brute_force_lookup_count(subject_word, 11),
                "subject_word={subject_word:?}"
            );
        }
    }
}
