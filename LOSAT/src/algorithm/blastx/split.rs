//! NCBI single-query BLASTX chunk metadata and score-only overlap rejoin.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    preliminary::{sort_by_score, PreliminaryHsp},
    query_setup::{prepare_queries, PreparedQueryBatch},
};
use anyhow::{ensure, Result};
use std::ops::Range;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:78-86
// ```c++
//         // multiple of 3, that way, when the nucleotide sequence(s) get(s)
//         // split, context N%6 in one chunk will have the same frame as context
//         // N%6 in the next chunk
//         case eBlastx:
//         case eTblastx:
//             // N.B.: the splitting is done on the nucleotide query sequences,
//             // then each of these chunks is translated
//             retval = 10002;
//             break;
// ```
pub const CHUNK_SIZE: usize = 10002;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:62-69
// ```c++
//     if (Blast_QueryIsTranslated(program)) {
//         // N.B.: this value must be divisible by 3 to work with translated
//         // queries, as we split them in nucleotide coordinates and then do the
//         // translation
//         retval = 297;
//     }
//     _TRACE("Using overlap chunk size " << retval);
//     return retval;
// ```
pub const OVERLAP_NT: usize = 297;
pub const OVERLAP_AA: i32 = 99;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:105-145
// ```c++
//     if ( !SplitQuery_ShouldSplit(program, *chunk_size,
//                                  concatenated_query_length, num_queries)) {
//         _TRACE("Not splitting queries");
//         return 1;
//     }
//
//     size_t overlap_size = SplitQuery_GetOverlapChunkSize(program);
//     Uint4 num_chunks = 0;
//
//     _DEBUG_ARG(size_t target_chunk_size = *chunk_size);
//
//     // For translated queries the chunk size should be divisible by CODON_LENGTH
//     if (Blast_QueryIsTranslated(program)) {
//         size_t chunk_size_delta = ((*chunk_size) % CODON_LENGTH);
//         *chunk_size -= chunk_size_delta;
//         _ASSERT((*chunk_size % CODON_LENGTH) == 0);
//     }
//
//     // Fix for small query size
//     if ((*chunk_size) > overlap_size) {
//        num_chunks = concatenated_query_length / ((*chunk_size) - overlap_size);
//     }
//
//     // Only one chunk, just return;
//     if (num_chunks <= 1) {
//        *chunk_size = concatenated_query_length;
//        return 1;
//     }
//
//     // Re-adjust the chunk_size to make load even
//     if (!Blast_QueryIsTranslated(program)) {
//        *chunk_size = (concatenated_query_length + (num_chunks - 1) * overlap_size) / num_chunks;
//        // Round up only if this will not decrease the number of chunks
//        if (num_chunks < (*chunk_size) - overlap_size ) (*chunk_size)++;
//     }
//
//     _TRACE("Number of chunks: " << num_chunks << "; "
//            "Target chunk size: " << target_chunk_size << "; "
//            "Returned chunk size: " << *chunk_size);
//
//     return num_chunks;
// ```
pub fn calculate_num_chunks(chunk_size: &mut usize, length: usize, num_queries: usize) -> usize {
    if num_queries > 1 {
        return 1;
    }
    *chunk_size -= *chunk_size % 3;
    let n = if *chunk_size > OVERLAP_NT {
        length / (*chunk_size - OVERLAP_NT)
    } else {
        0
    };
    if n <= 1 {
        *chunk_size = length;
        1
    } else {
        n
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:101-109
// ```c++
//         NCBI_THROW(CBlastException, eInvalidArgument,
//                    "Split query chunk size must be divisible by 3");
//     }
//     _TRACE("Returning query chunk size " << retval << " for " << EProgramToTaskName(program));
//     return retval;
// }
//
// CLocalBlast::CLocalBlast(CRef<IQueryFactory> qf,
//                          CRef<CBlastOptionsHandle> opts_handle,
// ```
pub fn validate_chunk_size(size: usize) -> Result<()> {
    ensure!(
        size % 3 == 0,
        "Split query chunk size must be divisible by 3"
    );
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:161-178
// ```c++
//         m_SplitBlk->SetChunkBounds(chunk_num,
//                                    TChunkRange(static_cast<unsigned int>(chunk_start), static_cast<unsigned int>(chunk_end)));
//         _TRACE("Chunk " << chunk_num << ": ranges from " << chunk_start
//                << " to " << chunk_end);
//
//         chunk_start += (m_ChunkSize - kOverlapSize);
//         if (chunk_start > m_TotalQueryLength ||
//             chunk_end == m_TotalQueryLength) {
//             break;
//         }
//     }
//
//     // For purposes of having an accurate overlap size when stitching back
//     // HSPs, save the overlap size
//     const size_t kOverlap =
//         Blast_QueryIsTranslated(m_Options->GetProgramType())
//         ? kOverlapSize / CODON_LENGTH : kOverlapSize;
//     m_SplitBlk->SetChunkOverlapSize(kOverlap);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:365-389
// ```c++
//                 for (unsigned int ctx = 0; ctx < kNumContexts; ctx++) {
//                     // handle the plus strand...
//                     if (ctx % NUM_FRAMES < CODON_LENGTH) {
//                         if (kStrand == eNa_strand_minus) {
//                             m_SplitBlk->AddContextToChunk(chunk_num,
//                                                           kInvalidContext);
//                         } else {
//                             m_SplitBlk->AddContextToChunk(chunk_num,
//                                               static_cast<Int4>(kNumContexts*queries[i]+ctx));
//                         }
//                     } else { // handle the negative strand
//                         if (kStrand == eNa_strand_plus) {
//                             m_SplitBlk->AddContextToChunk(chunk_num,
//                                                           kInvalidContext);
//                         } else {
//                             if (chunk_num == (size_t)last_query_chunk) {
//                                 // last chunk doesn't have shift
//                                 m_SplitBlk->AddContextToChunk(chunk_num,
//                                           static_cast<Int4>(kNumContexts*queries[i]+ctx));
//                             } else {
//                                 m_SplitBlk->AddContextToChunk(chunk_num,
//                                           static_cast<Int4>(kNumContexts*queries[i]+
//                                           s_AddShift(ctx, shift)));
//                             }
//                         }
// ```
#[derive(Clone, Debug)]
pub struct QueryChunk {
    pub range: Range<usize>,
    pub absolute_contexts: Vec<i32>,
    pub corrections: Vec<i32>,
    pub prepared: PreparedQueryBatch,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:50-62
// ```c++
// CQuerySplitter::CQuerySplitter(CRef<IQueryFactory> query_factory,
//                                const CBlastOptions* options)
//     : m_QueryFactory(query_factory), m_Options(options), m_NumChunks(0),
//     m_LocalQueryData(0), m_TotalQueryLength(0), m_ChunkSize(0)
// {
//     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
//     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
//     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
//     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
//         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
//     /* No split for ungapped mode JIRA SB-1082 */
//     if (!options->GetGappedMode()) m_NumChunks = 1;
//     x_ExtractCScopesAndMasks();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
// ```c++
//     size_t chunk_start = 0;
//     const size_t kOverlapSize =
//         SplitQuery_GetOverlapChunkSize(m_Options->GetProgramType());
//     for (size_t chunk_num = 0; chunk_num < m_NumChunks; chunk_num++) {
//         size_t chunk_end = chunk_start + m_ChunkSize;
//
//         // if the chunk end is larger than the sequence ...
//         if (chunk_end >= m_TotalQueryLength ||
//             // ... or this is the last chunk and it didn't make it to the end
//             // of the sequence
//             (chunk_end < m_TotalQueryLength && (chunk_num + 1) == m_NumChunks))
//         {
//             // ... assign this chunk's end to the end of the sequence
//             chunk_end = m_TotalQueryLength;
//         }
//
//         m_SplitBlk->SetChunkBounds(chunk_num,
//                                    TChunkRange(static_cast<unsigned int>(chunk_start), static_cast<unsigned int>(chunk_end)));
//         _TRACE("Chunk " << chunk_num << ": ranges from " << chunk_start
//                << " to " << chunk_end);
//
//         chunk_start += (m_ChunkSize - kOverlapSize);
//         if (chunk_start > m_TotalQueryLength ||
//             chunk_end == m_TotalQueryLength) {
//             break;
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:300-335
// ```c++
// static inline unsigned int
// s_AddShift(unsigned int context, int shift)
// {
//     _ASSERT(context == 3 || context == 4 || context == 5);
//     _ASSERT(shift == 0 || shift == 1 || shift == -1);
//
//     unsigned int retval;
//     if (shift == 0) {
//         retval = context;
//     } else if (shift == 1) {
//         retval = context == 3 ? 5 : context - shift;
//     } else if (shift == -1) {
//         retval = context == 5 ? 3 : context - shift;
//     } else {
//         abort();
//     }
//     return retval;
// }
//
// /**
//  * @brief Retrieve the shift for the negative strand
//  *
//  * @param query_length length of the query [in]
//  *
//  * @return shift (either 1, -1, or 0)
//  */
// static inline int
// s_GetShiftForTranslatedNegStrand(size_t query_length)
// {
//     int retval;
//     switch (query_length % CODON_LENGTH) {
//     case 1: retval = -1; break;
//     case 2: retval = 1; break;
//     case 0: default: retval = 0; break;
//     }
//     return retval;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:724-836
// ```c++
//             size_t correction = 0;
//             const int starting_chunk =
//                 ctx_translator.GetStartingChunk(chunk_num, ctx);
//             const int absolute_context =
//                 ctx_translator.GetAbsoluteContext(chunk_num, ctx);
//             const int last_query_chunk = qdpc.GetLastChunk(chunk_num, ctx);
//
//             if (absolute_context == kInvalidContext ||
//                 starting_chunk == kInvalidContext) {
//                 _ASSERT( !chunk_qinfo[chunk_num]->contexts[ctx].is_valid );
//                 // INT4_MAX is the sentinel value for invalid contexts
//                 m_SplitBlk->AddContextOffsetToChunk(chunk_num, INT4_MAX);
//                 continue;
//             }
//
//             // The corrections for the contexts corresponding to the negative
//             // strand in the last chunk of a query sequence are all 0
//             if (!s_IsPlusStrand(chunk_qinfo[chunk_num], ctx) &&
//                 (chunk_num == (size_t)last_query_chunk) &&
//                 (ctx % NUM_FRAMES >= 3)) {
//                 correction = 0;
//                 goto error_check;
//             }
//
//             // The corrections for the contexts corresponding to the plus
//             // strand are always the same, so only calculate the first one
//             if (s_IsPlusStrand(chunk_qinfo[chunk_num], ctx) &&
//                 (ctx % NUM_FRAMES == 1 || ctx % NUM_FRAMES == 2)) {
//                 correction = m_SplitBlk->GetContextOffsets(chunk_num).back();
//                 goto error_check;
//             }
//
//             // If the query length is divisible by CODON_LENGTH, the
//             // corrections for all contexts corresponding to a given strand are
//             // the same, so only calculate the first one
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 0) &&
//                 (ctx % NUM_FRAMES != 0) && (ctx % NUM_FRAMES != 3)) {
//                 correction = m_SplitBlk->GetContextOffsets(chunk_num).back();
//                 goto error_check;
//             }
//
//             // If the query length % CODON_LENGTH == 1, the corrections for the
//             // first two contexts of the negative strand are the same, and the
//             // correction for the last context is one more than that.
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 1) &&
//                 !s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
//
//                 if (ctx % NUM_FRAMES == 4) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back();
//                     goto error_check;
//                 } else if (ctx % NUM_FRAMES == 5) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back() + 1;
//                     goto error_check;
//                 }
//             }
//
//             // If the query length % CODON_LENGTH == 2, the corrections for the
//             // last two contexts of the negative strand are the same, which is
//             // one more that the first context on the negative strand.
//             if ((qdpc.GetQueryLength(chunk_num, ctx) % CODON_LENGTH == 2) &&
//                 !s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
//
//                 if (ctx % NUM_FRAMES == 4) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back() + 1;
//                     goto error_check;
//                 } else if (ctx % NUM_FRAMES == 5) {
//                     correction =
//                         m_SplitBlk->GetContextOffsets(chunk_num).back();
//                     goto error_check;
//                 }
//             }
//
//             if (s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
//
//                 for (int c = static_cast<int>(chunk_num); c != starting_chunk; c--) {
//                     size_t prev_len = s_GetAbsoluteContextLength(chunk_qinfo,
//                                                          c - 1,
//                                                          ctx_translator,
//                                                          absolute_context);
//                     size_t curr_len = s_GetAbsoluteContextLength(chunk_qinfo, c,
//                                                          ctx_translator,
//                                                          absolute_context);
//                     size_t overlap = min(kOverlap, curr_len);
//                     correction += prev_len - min(overlap, prev_len);
//                 }
//
//             } else {
//
//                 size_t subtrahend = 0;
//
//                 for (int c = static_cast<int>(chunk_num); c >= starting_chunk && c >= 0; c--) {
//                     size_t prev_len = s_GetAbsoluteContextLength(chunk_qinfo,
//                                                          c - 1,
//                                                          ctx_translator,
//                                                          absolute_context);
//                     size_t curr_len = s_GetAbsoluteContextLength(chunk_qinfo,
//                                                          c,
//                                                          ctx_translator,
//                                                          absolute_context);
//                     size_t overlap = min(kOverlap, curr_len);
//                     subtrahend += (curr_len - min(overlap, prev_len));
//                 }
//                 correction =
//                     global_qinfo->contexts[absolute_context].query_length -
//                     subtrahend;
//             }
//
// error_check:
//             _ASSERT((chunk_qinfo[chunk_num]->contexts[ctx].is_valid));
//             m_SplitBlk->AddContextOffsetToChunk(chunk_num, static_cast<Int4>(correction));
// ```
pub fn split_queries(
    records: &[FastaRecord],
    full: &PreparedQueryBatch,
    options: &ResolvedOptions,
) -> Result<Vec<QueryChunk>> {
    let mut size = CHUNK_SIZE;
    let length = records.iter().map(|r| r.sequence.len()).sum();
    let mut count = calculate_num_chunks(&mut size, length, records.len());
    if !options.gapped {
        count = 1;
    }
    if count == 1 {
        return Ok(vec![QueryChunk {
            range: 0..length,
            absolute_contexts: (0..full.contexts.len() as i32).collect(),
            corrections: vec![0; full.contexts.len()],
            prepared: full.clone(),
        }]);
    }
    let mut chunks = Vec::new();
    for index in 0..count {
        let start = index * (size - OVERLAP_NT);
        let end = if index + 1 == count {
            length
        } else {
            (start + size).min(length)
        };
        let mut record = records[0].clone();
        record.sequence = record.sequence[start..end].to_vec();
        record.lowercase_masks = record
            .lowercase_masks
            .iter()
            .filter_map(|&(a, b)| {
                let left = a.max(start as i32);
                let right = b.min(end as i32 - 1);
                (left <= right).then_some((left - start as i32, right - start as i32))
            })
            .collect();
        let prepared = prepare_queries(&[record], options)?;
        let shift = match length % 3 {
            1 => -1,
            2 => 1,
            _ => 0,
        };
        let absolute_contexts = (0..6)
            .map(|c| {
                if !prepared.contexts[c].is_valid {
                    -1
                } else if c < 3 || index + 1 == count || shift == 0 {
                    c as i32
                } else if shift == 1 {
                    if c == 3 {
                        5
                    } else {
                        c as i32 - 1
                    }
                } else if c == 5 {
                    3
                } else {
                    c as i32 + 1
                }
            })
            .collect();
        chunks.push(QueryChunk {
            range: start..end,
            absolute_contexts,
            corrections: vec![i32::MAX; 6],
            prepared,
        });
    }
    for index in 0..count {
        for c in 0..6 {
            let absolute = chunks[index].absolute_contexts[c];
            if absolute < 0 {
                continue;
            }
            let correction = if c == 1 || c == 2 {
                chunks[index].corrections[c - 1]
            } else if c == 0 {
                let mut sum = 0;
                for k in 1..=index {
                    let prev = absolute_length(&chunks, k as isize - 1, absolute);
                    let curr = absolute_length(&chunks, k as isize, absolute);
                    sum += prev - OVERLAP_AA.min(curr).min(prev);
                }
                sum
            } else if index + 1 == count {
                0
            } else if length % 3 == 0 && c > 3 {
                chunks[index].corrections[c - 1]
            } else if length % 3 == 1 && c == 4 {
                chunks[index].corrections[3]
            } else if length % 3 == 1 && c == 5 {
                chunks[index].corrections[4] + 1
            } else if length % 3 == 2 && c == 4 {
                chunks[index].corrections[3] + 1
            } else if length % 3 == 2 && c == 5 {
                chunks[index].corrections[4]
            } else {
                let mut sub = 0;
                for k in (0..=index).rev() {
                    let prev = absolute_length(&chunks, k as isize - 1, absolute);
                    let curr = absolute_length(&chunks, k as isize, absolute);
                    sub += curr - OVERLAP_AA.min(curr).min(prev);
                }
                full.contexts[absolute as usize].length as i32 - sub
            };
            chunks[index].corrections[c] = correction;
        }
    }
    Ok(chunks)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:477-510
// ```c++
// static string s_GetPrintableSequence(const Uint1* seq, size_t len, bool is_prot)
// {
//     string retval;
//     for (size_t i = 0; i < len; i++) {
//         retval.append(1, (is_prot
//                       ? NCBISTDAA_TO_AMINOACID[seq[i]]
//                       : BLASTNA_TO_IUPACNA[seq[i]]));
//     }
//     return retval;
// }
//
// /** Auxiliary function to validate the context offset corrections
//  * @param global global query sequence data [in]
//  * @param chunk sequence data for chunk [in]
//  * @param len length of the data to compare [in]
//  * @param is_prot whether the sequence is protein or not [in]
//  * @return true if sequence data is identical, false otherwise
//  */
// static bool cmp_sequence(const Uint1* global, const Uint1* chunk, size_t len,
//                          bool is_prot)
// {
//     bool retval = true;
//
//     for (size_t i = 0; i < len; i++) {
//         if (global[i] != chunk[i]) {
//             retval = false;
//             break;
//         }
//     }
//
//     if (retval == false) {
//         _TRACE("Comparing global: '"
//                << s_GetPrintableSequence(global, len, is_prot) << "'");
//         _TRACE("with chunk: '"
// ```
fn absolute_length(chunks: &[QueryChunk], index: isize, absolute: i32) -> i32 {
    if index < 0 {
        return 0;
    }
    let chunk = &chunks[index as usize];
    chunk
        .absolute_contexts
        .iter()
        .position(|&c| c == absolute)
        .map(|c| chunk.prepared.contexts[c].length as i32)
        .unwrap_or(0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:492-526
// ```c++
//        }
//        fprintf(stderr, "\n");
// #endif
//
//        for (j = 0; j < hitlist->hsplist_count; j++) {
//            BlastHSPList *hsplist = hitlist->hsplist_array[j];
//
//            for (k = 0; k < hsplist->hspcnt; k++) {
//                BlastHSP *hsp = hsplist->hsp_array[k];
//                Int4 local_context = hsp->context;
// #ifdef _DEBUG
//                ASSERT(local_context <= max_ctx);
//                ASSERT(local_context < num_ctx);
//                ASSERT(local_context < num_ctx_offsets);
// #endif
//
//                hsp->context = context_list[local_context];
//                hsp->query.offset += offset_list[local_context];
//                hsp->query.end += offset_list[local_context];
//                hsp->query.gapped_start += offset_list[local_context];
//                hsp->query.frame = BLAST_ContextToFrame(stream2->program,
//                                                        hsp->context);
//            }
//
//            hsplist->query_index = global_query;
//        }
//
//        Blast_HitListMerge(results1->hitlist_array + i,
//                           results2->hitlist_array + global_query,
//                           contexts_per_query, split_points,
//                           (Int4)SplitQueryBlk_GetChunkOverlapSize(squery_blk),
//                           SplitQueryBlk_AllowGap(squery_blk));
//    }
//
//    /* Sort to the canonical order, which the merge may not have done. */
// ```
pub fn adjust_chunk_hsps(hsps: &mut [PreliminaryHsp], chunk: &QueryChunk) -> Result<[i32; 6]> {
    let mut points = [-1; 6];
    for (c, &abs) in chunk.absolute_contexts.iter().enumerate() {
        if abs >= 0 {
            points[abs as usize % 6] = chunk.corrections[c];
        }
    }
    for h in hsps {
        let c = h.context;
        let abs = chunk.absolute_contexts[c];
        ensure!(abs >= 0, "HSP belongs to invalid split context");
        let offset = chunk.corrections[c];
        h.context = abs as usize;
        h.q_start += offset;
        h.q_end += offset;
        h.q_gapped_start += offset;
        h.frame = if abs % 6 < 3 {
            (abs % 6 + 1) as i8
        } else {
            -((abs % 6 - 2) as i8)
        };
    }
    Ok(points)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1488-1531
// ```c++
// s_BlastMergeTwoHSPs(BlastHSP* hsp1, BlastHSP* hsp2, Boolean allow_gap)
// {
//    ASSERT(!hsp1->gap_info || !hsp2->gap_info);
//
//    /* do not merge off-diagonal hsps for ungapped search */
//    if (!allow_gap &&
//        hsp1->subject.offset - hsp2->subject.offset -hsp1->query.offset + hsp2->query.offset)
//    {
//        return FALSE;
//    }
//
//    if(hsp1->subject.frame != hsp2->subject.frame)
// 	   return FALSE;
//
//    /* combine the boundaries of the two HSPs,
//       assuming they intersect at all */
//    if (CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end,
//                         hsp2->query.offset,
//                         hsp1->subject.offset, hsp1->subject.end,
//                         hsp2->subject.offset) ||
//        CONTAINED_IN_HSP(hsp1->query.offset, hsp1->query.end,
//                         hsp2->query.end,
//                         hsp1->subject.offset, hsp1->subject.end,
//                         hsp2->subject.end)) {
//
// 	  double score_density =  (hsp1->score + hsp2->score) *(1.0) /
// 			                  ((hsp1->query.end - hsp1->query.offset) +
// 			                   (hsp2->query.end - hsp2->query.offset));
//       hsp1->query.offset = MIN(hsp1->query.offset, hsp2->query.offset);
//       hsp1->subject.offset = MIN(hsp1->subject.offset, hsp2->subject.offset);
//       hsp1->query.end = MAX(hsp1->query.end, hsp2->query.end);
//       hsp1->subject.end = MAX(hsp1->subject.end, hsp2->subject.end);
//       if (hsp2->score > hsp1->score) {
//           hsp1->query.gapped_start = hsp2->query.gapped_start;
//           hsp1->subject.gapped_start = hsp2->subject.gapped_start;
// 	  hsp1->score = hsp2->score;
//       }
//
//       hsp1->score = MAX((int) (score_density *(hsp1->query.end - hsp1->query.offset)), hsp1->score);
//       return TRUE;
//    }
//
//    return FALSE;
// }
// ```
pub(crate) fn merge_two(a: &mut PreliminaryHsp, b: &PreliminaryHsp) -> bool {
    let inside = |q, s| a.q_start <= q && a.q_end >= q && a.s_start <= s && a.s_end >= s;
    if !inside(b.q_start, b.s_start) && !inside(b.q_end, b.s_end) {
        return false;
    }
    let density =
        (a.score + b.score) as f64 / ((a.q_end - a.q_start) + (b.q_end - b.q_start)) as f64;
    a.q_start = a.q_start.min(b.q_start);
    a.q_end = a.q_end.max(b.q_end);
    a.s_start = a.s_start.min(b.s_start);
    a.s_end = a.s_end.max(b.s_end);
    if b.score > a.score {
        a.q_gapped_start = b.q_gapped_start;
        a.s_gapped_start = b.s_gapped_start;
        a.score = b.score;
    }
    a.score = a.score.max((density * (a.q_end - a.q_start) as f64) as i32);
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2911-3003
// ```c++
//    }
//    else {            /* query seq is split */
//
//       /* An HSP can be a candidate for merging if it lies in the
//          overlap region. Whether this is true depends on whether the
//          HSP starts to the left of the split point, or ends to the
//          right of the overlap region. A complication is that 'left'
//          and 'right' have opposite meaning when the HSP is on the
//          minus strand of the query sequence */
//
//       for (index1 = 0; index1 < combined_hsp_list->hspcnt; index1++) {
//          hsp1 = combined_hsp_list->hsp_array[index1];
//          offset_idx = hsp1->context % contexts_per_query;
//          if (split_offsets[offset_idx] < 0) continue;
//          if ((hsp1->query.frame >= 0 && hsp1->query.end >
//                          split_offsets[offset_idx]) ||
//              (hsp1->query.frame < 0 && hsp1->query.offset <
//                          split_offsets[offset_idx] + chunk_overlap_size)) {
//             /* At least part of this HSP lies in the overlap strip. */
//             hsp_var = combined_hsp_list->hsp_array[hspcnt1];
//             combined_hsp_list->hsp_array[hspcnt1] = hsp1;
//             combined_hsp_list->hsp_array[index1] = hsp_var;
//             ++hspcnt1;
//          }
//       }
//       for (index2 = 0; index2 < hsp_list->hspcnt; index2++) {
//          hsp2 = hsp_list->hsp_array[index2];
//          offset_idx = hsp2->context % contexts_per_query;
//          if (split_offsets[offset_idx] < 0) continue;
//          if ((hsp2->query.frame < 0 && hsp2->query.end >
//                          split_offsets[offset_idx]) ||
//              (hsp2->query.frame >= 0 && hsp2->query.offset <
//                          split_offsets[offset_idx] + chunk_overlap_size)) {
//             /* At least part of this HSP lies in the overlap strip. */
//             hsp_var = hsp_list->hsp_array[hspcnt2];
//             hsp_list->hsp_array[hspcnt2] = hsp2;
//             hsp_list->hsp_array[index2] = hsp_var;
//             ++hspcnt2;
//          }
//       }
//    }
//
//    /* the merge process is independent of whether merging happens
//       between query chunks or subject chunks */
//
//    if (hspcnt1 > 0 && hspcnt2 > 0) {
//       hspp1 = combined_hsp_list->hsp_array;
//       hspp2 = hsp_list->hsp_array;
//
//       for (index1 = 0; index1 < hspcnt1; index1++) {
//
//          hsp1 = hspp1[index1];
//
//          for (index2 = 0; index2 < hspcnt2; index2++) {
//
//             hsp2 = hspp2[index2];
//
//             /* Skip already deleted HSPs, or HSPs from different contexts */
//             if (!hsp2 || hsp1->context != hsp2->context)
//                continue;
//
//             /* Short read qureies are shorter than the overlap region and may
//                already have a traceback */
//             if (short_reads) {
//                 hspp2[index2] = Blast_HSPFree(hsp2);
//                 continue;
//             }
//
//             /* we have to determine the starting diagonal of one HSP
//                and the ending diagonal of the other */
//
//             if (contexts_per_query < 0 || hsp1->query.frame >= 0) {
//                end_diag = s_HSPEndDiag(hsp1);
//                start_diag = s_HSPStartDiag(hsp2);
//             }
//             else {
//                end_diag = s_HSPEndDiag(hsp2);
//                start_diag = s_HSPStartDiag(hsp1);
//             }
//
//             if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
//                if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) {
//                   /* Free the second HSP. */
//                   hspp2[index2] = Blast_HSPFree(hsp2);
//                }
//             }
//          }
//       }
//
//       /* Purge the nulled out HSPs from the new HSP list */
//       Blast_HSPListPurgeNullHSPs(hsp_list);
//    }
//
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2757-2770
// ```c++
//
//    if (new_hspcnt >= hsp_list->hspcnt + combined_hsp_list->hspcnt) {
//       /* All HSPs from both arrays are saved */
//       for (index=combined_hsp_list->hspcnt, index1=0;
//            index1<hsp_list->hspcnt; index1++) {
//          if (hsp_list->hsp_array[index1] != NULL)
//             combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
//       }
//       combined_hsp_list->hspcnt = new_hspcnt;
//       Blast_HSPListSortByScore(combined_hsp_list);
//    } else {
//       /* Not all HSPs are be saved; sort both arrays by score and save only
//          the new_hspcnt best ones.
//          For the merged set of HSPs, allocate array the same size as in the
// ```
pub fn merge_chunk_hsps(
    combined: &mut Vec<PreliminaryHsp>,
    mut incoming: Vec<PreliminaryHsp>,
    points: &[i32; 6],
    max_hsps: usize,
) {
    if incoming.is_empty() {
        return;
    }
    if combined.is_empty() {
        *combined = incoming;
        return;
    }
    let mut nold = 0;
    for i in 0..combined.len() {
        let h = &combined[i];
        let point = points[h.context % 6];
        if point >= 0
            && (if h.frame >= 0 {
                h.q_end > point
            } else {
                h.q_start < point + OVERLAP_AA
            })
        {
            combined.swap(nold, i);
            nold += 1;
        }
    }
    let mut nnew = 0;
    for i in 0..incoming.len() {
        let h = &incoming[i];
        let point = points[h.context % 6];
        if point >= 0
            && (if h.frame < 0 {
                h.q_end > point
            } else {
                h.q_start < point + OVERLAP_AA
            })
        {
            incoming.swap(nnew, i);
            nnew += 1;
        }
    }
    let mut removed = vec![false; incoming.len()];
    for a in &mut combined[..nold] {
        for (i, b) in incoming[..nnew].iter().enumerate() {
            if removed[i] || a.context != b.context {
                continue;
            }
            let (end_diag, start_diag) = if a.frame >= 0 {
                (a.q_end - a.s_end, b.q_start - b.s_start)
            } else {
                (b.q_end - b.s_end, a.q_start - a.s_start)
            };
            if (end_diag - start_diag).abs() < 10 && merge_two(a, b) {
                removed[i] = true;
            }
        }
    }
    combined.extend(
        incoming
            .into_iter()
            .enumerate()
            .filter_map(|(i, h)| (!removed[i]).then_some(h)),
    );
    sort_by_score(combined);
    combined.truncate(max_hsps);
}
