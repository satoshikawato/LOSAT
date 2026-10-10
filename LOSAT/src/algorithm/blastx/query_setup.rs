//! BLASTX context/translation/mask preparation before score/statistics setup.
use super::{args::ResolvedOptions, input::FastaRecord};
use crate::algorithm::tblastx::translation::generate_frames;
use crate::core::blast_seg::SegMasker;
use crate::utils::genetic_code::GeneticCode;
use anyhow::{bail, Result};

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
// ```c++
//         for (i = 0; i < retval->last_context + 1; i++) {
//             retval->contexts[i].query_index =
//                 Blast_GetQueryIndexFromContext(i, program);
//             ASSERT(retval->contexts[i].query_index != -1);
//
//             retval->contexts[i].frame = BLAST_ContextToFrame(program,  i);
//             ASSERT(retval->contexts[i].frame != INT1_MAX);
//
//             retval->contexts[i].is_valid = TRUE;
// ```
pub struct QueryContext {
    pub query_index: usize,
    pub frame: i8,
    pub offset: usize,
    pub length: usize,
    pub is_valid: bool,
    pub lowercase_masks: Vec<(i32, i32)>,
    pub masks: Vec<(i32, i32)>,
    pub dna_masks: Vec<(i32, i32)>,
}
#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:264-266
// ```c++
//             int seg_flags = queries.GetSegmentInfo(j);
//             query_info->contexts[ctx_index].segment_flags = seg_flags;
//             query_info->contexts[ctx_index + 1].segment_flags = seg_flags;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:499-501
// ```c++
//
//     int buflen = QueryInfo_GetSeqBufLen(qinfo);
//     TAutoUint1Ptr buf((Uint1*) calloc(buflen+1, sizeof(Uint1)));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:654-659
// ```c++
//     BlastSeqBlkSetSequence(*seqblk, buf.release(), buflen - 2);
//
//     (*seqblk)->lcase_mask = mask.Release();
//     (*seqblk)->lcase_mask_allocated = TRUE;
// }
//
// ```
pub struct PreparedQueryBatch {
    pub original_lengths: Vec<usize>,
    pub contexts: Vec<QueryContext>,
    pub first_context: usize,
    pub max_length: usize,
    pub min_length: usize,
    pub sequence_start: Vec<u8>,
    pub sequence_start_nomask: Vec<u8>,
    pub lookup_segments: Vec<(i32, i32)>,
    pub split_eligible: bool,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:123-136
// ```c++
//     // in the next chunk
//     case eBlastx:
//         if (task_name == "blastx-fast" && mt_mode == true)
//         {
// 		retval = 20004;
// 		break;
// 	}
//     case eTblastx:
//         // N.B.: the splitting is done on the nucleotide query sequences, then
//         // each of these chunks is translated
//         retval = 10002;
//         break;
//     case eBlastp:
//     default:
// ```
pub const BATCH_SIZE: usize = 10002;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:139-165
// ```c++
//
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
//
//         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
//
//         if (loc->IsInt()) {
//             size_read += sequence::GetLength(loc->GetInt().GetId(),
//                                              q->GetScope());
//         } else if (loc->IsWhole()) {
//             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
//         } else {
//             // programmer error, CBlastInputSource should only return Seq-locs
// ```
pub fn batch_ranges(records: &[FastaRecord]) -> Vec<std::ops::Range<usize>> {
    let mut result = Vec::new();
    let mut start = 0;
    while start < records.len() {
        let mut end = start;
        let mut length = 0;
        while end < records.len() && length < BATCH_SIZE {
            length += records[end].sequence.len();
            end += 1;
        }
        result.push(start..end);
        start = end;
    }
    result
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-88
// ```c++
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
//         // First context
//         qinfo->contexts[0].query_offset = 0;
//         qinfo->contexts[0].query_length = length;
//         if (length == 0)
//            qinfo->contexts[0].is_valid = false;
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:191-219
// ```c++
//
//         if (translate) {
//             for (unsigned int i = 0; i < kNumContexts; i++) {
//                 unsigned int prot_length =
//                     static_cast<unsigned int>(BLAST_GetTranslatedProteinLength(length, i));
//                 max_length = MAX(max_length, prot_length);
//                 min_length = MIN(min_length, prot_length);
//
//                 Uint4 ctx_len(0);
//
//                 switch (strand) {
//                 case eNa_strand_plus:
//                     ctx_len = (i<3) ? prot_length : 0;
//                     s_QueryInfo_SetContext(query_info, ctx_index + i, ctx_len);
//                     // the missing frame is present in query_info as
//                     // zero-lenghth context
//                     min_length = 0;
//                     break;
//
//                 case eNa_strand_minus:
//                     ctx_len = (i<3) ? 0 : prot_length;
//                     s_QueryInfo_SetContext(query_info, ctx_index + i, ctx_len);
//                     min_length = 0;
//                     break;
//
//                 case eNa_strand_both:
//                 case eNa_strand_unknown:
//                     s_QueryInfo_SetContext(query_info, ctx_index + i,
//                                            prot_length);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:577-587
// ```c++
//                 for (unsigned int i = 0; i < kNumContexts; i++) {
//                     if (qinfo->contexts[ctx_index + i].query_length <= 0) {
//                         continue;
//                     }
//
//                     int offset = qinfo->contexts[ctx_index + i].query_offset;
//                     BLAST_GetTranslation(sequence.data.get() + 1,
//                                          seqbuf_rev,
//                                          na_length,
//                                          qinfo->contexts[ctx_index + i].frame,
//                                          & buf.get()[offset], gc);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:620-641
// ```c++
//
//     if (!mask_at_hash) {
//         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
//                              program_number);
//     }
//
//     if (program_number == eBlastTypeBlastx && scoring_options->is_ooframe) {
//         BLAST_CreateMixedFrameDNATranslation(query_blk, query_info);
//     }
//
//     /* Find complement of the mask locations, for which lookup table will be
//      * created. This should only be done if we do want to create a lookup table,
//      * i.e. if it is a full search, not a traceback-only search.
//      */
//     if (lookup_segments) {
//         BLAST_ComplementMaskLocations(program_number, query_info,
//                                       filter_maskloc, lookup_segments);
//     }
//
//     if (mask)
//     {
//         if (Blast_QueryIsTranslated(program_number)) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1003
// ```c++
// /// Auxiliary class to validate the genetic code input
// class CArgAllowGeneticCodeInteger : public CArgAllow
// {
// protected:
//      /// Overloaded method from CArgAllow
//      virtual bool Verify(const string& value) const {
//          static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,33};
// ```
pub fn prepare_queries(
    records: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<PreparedQueryBatch> {
    if records.is_empty() {
        bail!("cannot prepare an empty BLASTX query batch");
    }
    crate::blastinput::value_parsers::genetic_code(&options.query_gencode.to_string())
        .map_err(anyhow::Error::msg)?;
    let code = GeneticCode::try_from_id(options.query_gencode).map_err(anyhow::Error::msg)?;
    let (mut contexts, mut max_length, mut min_length) = (Vec::new(), 0, usize::MAX);
    let mut offset = 0;
    for (q, record) in records.iter().enumerate() {
        for i in 0..6 {
            let frame = if i < 3 {
                (i + 1) as i8
            } else {
                -((i - 2) as i8)
            };
            let len = record.sequence.len().saturating_sub(i % 3) / 3;
            max_length = max_length.max(len);
            min_length = min_length.min(len);
            let enabled = match options.strand.as_str() {
                "plus" => i < 3,
                "minus" => i >= 3,
                _ => true,
            };
            let length = if enabled { len } else { 0 };
            let lowercase_masks = if enabled {
                dna_to_protein_masks(&record.lowercase_masks, record.sequence.len(), frame)
            } else {
                Vec::new()
            };
            contexts.push(QueryContext {
                query_index: q,
                frame,
                offset,
                length,
                is_valid: length > 0,
                lowercase_masks,
                masks: Vec::new(),
                dna_masks: Vec::new(),
            });
            if length > 0 {
                offset += length + 1;
            }
        }
    }
    if options.strand != "both" {
        min_length = 0;
    }
    let last = contexts.last().unwrap();
    let buffer_length = last.offset + last.length + if last.length > 0 { 2 } else { 1 };
    let mut buffer = vec![0; buffer_length];
    for (q, record) in records.iter().enumerate() {
        for frame in generate_frames(&record.sequence, &code) {
            let i = if frame.frame > 0 {
                frame.frame as usize - 1
            } else {
                3 + (-frame.frame) as usize - 1
            };
            let context = &contexts[q * 6 + i];
            if context.length > 0 {
                buffer[context.offset..context.offset + context.length + 2]
                    .copy_from_slice(&frame.aa_seq);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup.hpp:192-197
    // ```c++
    //     TSeqPos size() const {
    //         TSeqPos retval = x_Size();
    //         if (retval == 0) {
    //             NCBI_THROW(CBlastException, eInvalidArgument,
    //                        "Sequence contains no data");
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:634-639
    // ```c++
    //             // to determine whether the message should contain a warning or
    //             // error?
    //             CRef<CSearchMessage> m
    //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
    //             messages[index].push_back(m);
    //             s_InvalidateQueryContexts(qinfo, index);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:650-652
    // ```c++
    //     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
    //         NCBI_THROW(CBlastException, eSetup, messages.ToString());
    //     }
    // ```
    if !contexts.iter().any(|c| c.is_valid) && records.iter().any(|q| q.sequence.is_empty()) {
        anyhow::bail!(
            "{}",
            "Warning: Sequence contains no data "
                .repeat(records.iter().filter(|q| q.sequence.is_empty()).count())
        );
    }
    let nomask = buffer.clone();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1144-1158
    // ```c++
    //         SSegOptions* seg_options = filter_options->segOptions;
    //         SegParameters* sparamsp=NULL;
    //
    //         sparamsp = SegParametersNewAa();
    //         sparamsp->overlaps = TRUE;
    //         if (seg_options->window > 0)
    //             sparamsp->window = seg_options->window;
    //         if (seg_options->locut > 0.0)
    //             sparamsp->locut = seg_options->locut;
    //         if (seg_options->hicut > 0.0)
    //             sparamsp->hicut = seg_options->hicut;
    //
    // 		status = SeqBufferSeg(sequence, length, offset, sparamsp,
    //                               seqloc_retval);
    // 		SegParametersFree(sparamsp);
    // ```
    // Same SEG parameters as the block above; the masker is built once for all contexts instead of
    // once per context (placement only).
    let seg_masker = options.seg.enabled.then(|| {
        let s = &options.seg;
        SegMasker::new(
            if s.window > 0 { s.window as usize } else { 12 },
            if s.locut > 0.0 { s.locut } else { 2.2 },
            if s.hicut > 0.0 { s.hicut } else { 2.5 },
        )
        // BLASTX keeps LOSAT's former SEG until SX (plan DW-10).
        .keeping_all_left_segments()
    });
    // EXPERIMENT (LOSAT_X_BXPAR): SEG reads only the residues of its own context
    // and the hard masking below writes only there, so the contexts of a large
    // batch (the six frames of the full query) are filtered on the search pool
    // before the loop consumes the intervals in context order.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
    // ```c
    //     for (context = query_info->first_context;
    //          context <= query_info->last_context; ++context) {
    // ...
    //         BlastSeqLoc *filter_per_context = NULL;
    //         status = s_GetFilteringLocationsForOneContext(query_blk,
    //                                                       query_info,
    //                                                       context,
    //                                                       program_number,
    //                                                       filter_options,
    //                                                       &filter_per_context,
    //                                                       blast_message);
    // ```
    // NCBI filters the contexts one after another in this loop. SEG of one context reads only that
    // context, so with LOSAT_X_BXPAR the contexts of a large batch are filtered in parallel and the
    // intervals are kept in seg_ready; the loop below still consumes them in context order.
    #[allow(unused_mut)]
    let mut seg_ready: Vec<Option<Vec<(i32, i32)>>> = Vec::new();
    #[cfg(feature = "parallel")]
    if let Some(masker) = seg_masker.as_ref() {
        if super::runtime::x_bx_parallel()
            && buffer.len() >= 60_000
            && !super::runtime::x_inner_serial()
            && rayon::current_thread_index().is_some()
            && rayon::current_num_threads() > 1
        {
            use rayon::prelude::*;
            let buffer = &buffer;
            seg_ready = contexts
                .par_iter()
                .map(|context| {
                    context.is_valid.then(|| {
                        masker
                            .mask_sequence(
                                &buffer[context.offset + 1..context.offset + 1 + context.length],
                            )
                            .into_iter()
                            .map(|r| (r.start as i32, r.end as i32 - 1))
                            .collect()
                    })
                })
                .collect();
        }
    }
    for (context_index, context) in contexts.iter_mut().enumerate() {
        if !context.is_valid {
            continue;
        }
        let mut masks = context.lowercase_masks.clone();
        if let Some(masker) = seg_masker.as_ref() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
            // ```c
            //     for (context = query_info->first_context;
            //          context <= query_info->last_context; ++context) {
            // ...
            //         status = s_GetFilteringLocationsForOneContext(query_blk,
            //                                                       query_info,
            //                                                       context,
            //                                                       program_number,
            //                                                       filter_options,
            //                                                       &filter_per_context,
            //                                                       blast_message);
            // ```
            // Dispatch point: a context filtered ahead of the loop uses its stored intervals; any
            // other context is filtered here exactly as before.
            if let Some(ready) = seg_ready.get_mut(context_index).and_then(Option::take) {
                masks.extend(ready);
            } else {
                masks.extend(
                    masker
                        .mask_sequence(
                            &buffer[context.offset + 1..context.offset + 1 + context.length],
                        )
                        .into_iter()
                        .map(|r| (r.start as i32, r.end as i32 - 1)),
                );
            }
        }
        context.masks = combine_masks(masks);
        if !options.soft_masking {
            for &(from, to) in &context.masks {
                for index in from..=to {
                    buffer[context.offset + 1 + index as usize] = 21;
                }
            }
        }
        context.dna_masks = protein_to_dna_masks(
            &context.masks,
            core_dna_length(&contexts_lengths(
                records[context.query_index].sequence.len(),
            )),
            context.frame,
        );
    }
    let lookup_segments = complement_masks(&contexts);
    Ok(PreparedQueryBatch {
        original_lengths: records.iter().map(|r| r.sequence.len()).collect(),
        contexts,
        first_context: if options.strand == "minus" { 3 } else { 0 },
        max_length,
        min_length,
        sequence_start: buffer,
        sequence_start_nomask: nomask,
        lookup_segments,
        split_eligible: records.len() == 1 && options.gapped,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:923-929
// ```c++
// BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
// {
//     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
//         return 0;
//     }
//     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
// }
// ```
fn contexts_lengths(length: usize) -> [usize; 3] {
    [
        length / 3,
        length.saturating_sub(1) / 3,
        length.saturating_sub(2) / 3,
    ]
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:143-160
// ```c++
//     Int4 start_context = NUM_FRAMES*query_index;
//     Int4 dna_length = 2;
//     Int4 index;
//
//     /* Make sure that query index is within appropriate range, and that this is
//        really a translated search */
//     ASSERT(query_index < query_info->num_queries);
//     ASSERT(start_context < query_info->last_context);
//
//     /* If only reverse strand is searched, then forward strand contexts don't
//        have lengths information */
//     if (query_info->contexts[start_context].query_length == 0)
//         start_context += 3;
//
//     for (index = start_context; index < start_context + 3; ++index)
//         dna_length += query_info->contexts[index].query_length;
//
//     return dna_length;
// ```
fn core_dna_length(lengths: &[usize; 3]) -> usize {
    lengths.iter().sum::<usize>() + 2
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1604-1635
// ```c++
//             short frame = iter->first;
//             BlastSeqLoc * bsl = iter->second;
//
//             for (BlastSeqLoc* itr = bsl; itr; itr = itr->next) {
//                 int to(0), from(0);
//
//                 if (frame < 0) {
//                     from = ((int) dna_length + frame - itr->ssr->right) / CODON_LENGTH;
//                     to = ((int) dna_length + frame - itr->ssr->left) / CODON_LENGTH;
//                 } else {
//                     from = (itr->ssr->left - frame + 1) / CODON_LENGTH;
//                     to = (itr->ssr->right - frame + 1) / CODON_LENGTH;
//                 }
//                 if (from < 0)
//                     from = 0;
//                 if (to < 0)
//                     to = 0;
//                 const int kFrameLength = frame_lengths[(CSeqLocInfo::ETranslationFrame)frame];
//                 if (from >= kFrameLength)
//                     from = kFrameLength - 1;
//                 if (to >= kFrameLength)
//                     to = kFrameLength - 1;
//
//                 _ASSERT(from >= 0 && to >= 0);
//                 _ASSERT(from < kFrameLength && to < kFrameLength);
//                 itr->ssr->left  = from;
//                 itr->ssr->right = to;
//             }
//         }
//     }
// }
//
// ```
pub fn dna_to_protein_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    let aa_length = (length - (frame.abs() - 1)).max(0) / 3;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                ((length + frame - right) / 3, (length + frame - left) / 3)
            } else {
                ((left - frame + 1) / 3, (right - frame + 1) / 3)
            };
            (from.max(0).min(aa_length - 1), to.max(0).min(aa_length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:928-954
// ```c++
//                } else {
//                    from = CODON_LENGTH*seq_range->left + frame - 1;
//                    to = CODON_LENGTH*seq_range->right + frame - 1;
//                }
//
//                if (from < 0)
//                    from = 0;
//                if (to   < 0)
//                    to   = 0;
//                if (from >= dna_length)
//                    from = dna_length - 1;
//                if (to   >= dna_length)
//                    to   = dna_length - 1;
//
//                ASSERT(from >= 0);
//                ASSERT(to   >= 0);
//                ASSERT(from < dna_length);
//                ASSERT(to   < dna_length);
//
//                seq_range->left = from;
//                seq_range->right = to;
//            }
//        }
//    }
//    return status;
// }
//
// ```
pub fn protein_to_dna_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                (length - 3 * right + frame + 1, length - 3 * left + frame)
            } else {
                (3 * left + frame - 1, 3 * right + frame - 1)
            };
            (from.max(0).min(length - 1), to.max(0).min(length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:983-1004
// ```c++
//     qsort(ptrs, (size_t)num_elems, sizeof(*ptrs),
//           s_SeqRangeSortByStartPosition);
//
//     /* Merge the overlapping elements */
//     {
//         BlastSeqLoc* curr_tail = *mask_loc = ptrs[0];
//         for (i = 0; i < num_elems - 1; i++) {
//             const SSeqRange* next_ssr = ptrs[i+1]->ssr;
//             const Int4 stop = curr_tail->ssr->right;
//
//             if ((stop + link_value) > next_ssr->left) {
//                 curr_tail->ssr->right = MAX(stop, next_ssr->right);
//                 ptrs[i+1] = BlastSeqLocNodeFree(ptrs[i+1]);
//             } else {
//                 curr_tail = ptrs[i+1];
//             }
//         }
//     }
//
//     /* Rebuild the linked list */
//     {
//         BlastSeqLoc* tail = *mask_loc;
// ```
pub fn combine_masks(mut masks: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    masks.sort_by_key(|r| r.0);
    let mut result: Vec<(i32, i32)> = Vec::new();
    for mask in masks {
        if let Some(last) = result.last_mut().filter(|r| r.1 > mask.0) {
            last.1 = last.1.max(mask.1);
        } else {
            result.push(mask);
        }
    }
    result
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1073-1115
// ```c++
//             filter_end = start_offset + seq_range->right;
//          }
//          /* The canonical "state" at the top of this
//             while loop is that both "left" and "right" have
//             been initialized to their correct values.
//             The first time this loop is entered in a call to
//             the function this is not true and the following "if"
//             statement moves everything to the canonical state. */
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
//
//          right = filter_start - 1;
//
//          /* Cache the tail of the list to avoid the overhead of traversing the
//           * list when appending to it */
//          tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
//          if (filter_end >= end_offset) {
//             /* last masked region at end of sequence */
//             last_interval_open = FALSE;
//             break;
//          } else {
//             left = filter_end + 1;
//          }
//       }
//
//       if (last_interval_open) {
//          /* Need to finish SSeqRange* for last interval. */
//          right = end_offset;
//          /* Cache the tail of the list to avoid the overhead of traversing the
//           * list when appending to it */
//          tail = BlastSeqLocNew((tail ? &tail : complement_mask), left, right);
//       }
// ```
fn complement_masks(contexts: &[QueryContext]) -> Vec<(i32, i32)> {
    let mut result = Vec::new();
    for c in contexts.iter().filter(|c| c.is_valid) {
        let start = c.offset as i32;
        let end = start + c.length as i32 - 1;
        if c.masks.is_empty() {
            result.push((start, end));
            continue;
        }
        let mut first = true;
        let mut last_interval_open = true;
        let mut left = 0;
        for &(a, b) in &c.masks {
            let filter_start = start + a;
            let filter_end = start + b;
            if first {
                last_interval_open = true;
                first = false;
                if filter_start > start {
                    left = start;
                } else {
                    left = filter_end + 1;
                    continue;
                }
            }
            result.push((left, filter_start - 1));
            if filter_end >= end {
                last_interval_open = false;
                break;
            } else {
                left = filter_end + 1;
            }
        }
        if last_interval_open {
            result.push((left, end));
        }
    }
    result
}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
// ```c++
//         for (i = 0; i < retval->last_context + 1; i++) {
//             retval->contexts[i].query_index =
//                 Blast_GetQueryIndexFromContext(i, program);
//             ASSERT(retval->contexts[i].query_index != -1);
//
//             retval->contexts[i].frame = BLAST_ContextToFrame(program,  i);
//             ASSERT(retval->contexts[i].frame != INT1_MAX);
//
//             retval->contexts[i].is_valid = TRUE;
// ```
mod tests {
    use super::*;
    #[test]
    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blastfilter_unit_test.cpp:1003-1027 BlastxLowerCaseMaskProteinLocations
    // ```c++
    //     bqff.UseProteinCoords(9180); // 9180 is length of GI|1945388
    //
    //     BlastSeqLoc* bsl = *bqff[CSeqLocInfo::eFramePlus1];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 0);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 25);
    //
    //     bsl = *bqff[CSeqLocInfo::eFramePlus2];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 0);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 24);
    //
    //     bsl = *bqff[CSeqLocInfo::eFramePlus3];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 0);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 24);
    //
    //     bsl = *bqff[CSeqLocInfo::eFrameMinus1];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 3034);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 3059);
    //
    //     bsl = *bqff[CSeqLocInfo::eFrameMinus2];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 3034);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 3058);
    //
    //     bsl = *bqff[CSeqLocInfo::eFrameMinus3];
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->left, 3034);
    //     BOOST_REQUIRE_EQUAL(bsl->ssr->right, 3058);
    // ```
    fn ncbi_active_lowercase_frame_assertions() {
        let expected = [
            (0, 25),
            (0, 24),
            (0, 24),
            (3034, 3059),
            (3034, 3058),
            (3034, 3058),
        ];
        for (frame, expected) in [1, 2, 3, -1, -2, -3].into_iter().zip(expected) {
            assert_eq!(
                dna_to_protein_masks(&[(0, 75)], 9180, frame),
                vec![expected]
            );
        }
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:995-1003
    // ```c++
    //                 ptrs[i+1] = BlastSeqLocNodeFree(ptrs[i+1]);
    //             } else {
    //                 curr_tail = ptrs[i+1];
    //             }
    //         }
    //     }
    //
    //     /* Rebuild the linked list */
    //     {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1089-1096
    // ```c++
    //                /* beginning of sequence filtered */
    //                left = filter_end + 1;
    //                continue;
    //             }
    //          }
    //
    //          right = filter_start - 1;
    //
    // ```
    fn touching_and_duplicate_masks_preserve_lookup_intervals() {
        assert_eq!(
            combine_masks(vec![(0, 0), (0, 0), (1, 2), (2, 3)]),
            vec![(0, 0), (0, 0), (1, 2), (2, 3)]
        );
        let c = QueryContext {
            query_index: 0,
            frame: 1,
            offset: 0,
            length: 2,
            is_valid: true,
            lowercase_masks: vec![],
            masks: vec![(0, 0), (0, 0)],
            dna_masks: vec![],
        };
        assert_eq!(complement_masks(&[c]), vec![(1, -1), (1, 1)]);
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:139-165
    // ```c++
    //
    //     while (size_read < GetBatchSize()) {
    //
    //         if (End())
    //             break;
    //
    //         CRef<CBlastSearchQuery> q;
    //         try { q.Reset(m_Source->GetNextSequence(scope)); }
    //         catch (const CObjReaderParseException& e) {
    //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
    //                 break;
    //             }
    //             throw;
    //         }
    //         catch (const exception&) {
    //             continue; //SB-2307. ignore well formed, not found accession
    //         }
    //
    //         CConstRef<CSeq_loc> loc = q->GetQuerySeqLoc();
    //
    //         if (loc->IsInt()) {
    //             size_read += sequence::GetLength(loc->GetInt().GetId(),
    //                                              q->GetScope());
    //         } else if (loc->IsWhole()) {
    //             size_read += sequence::GetLength(loc->GetWhole(), q->GetScope());
    //         } else {
    //             // programmer error, CBlastInputSource should only return Seq-locs
    // ```
    fn ncbi_batch_crossing_record_is_retained() {
        let record = |n| FastaRecord {
            internal_id: String::new(),
            title: String::new(),
            sequence: vec![b'A'; n],
            lowercase_masks: vec![],
            warnings: vec![],
        };
        assert_eq!(
            batch_ranges(&[record(10001), record(1), record(3)]),
            vec![0..2, 2..3]
        );
        assert_eq!(batch_ranges(&[record(10002), record(3)]), vec![0..1, 1..2]);
        assert_eq!(batch_ranges(&[record(10003), record(3)]), vec![0..1, 1..2]);
    }
}
