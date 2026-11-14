//! BLASTX context/translation/mask preparation before score/statistics setup.
use super::{args::ResolvedOptions, input::FastaRecord};
use crate::algorithm::common::masks::{combine_masks, dna_to_protein_masks, protein_to_dna_masks};
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
    /// The values NCBI takes of the query batch (reader batch) that this unit is a part of, when
    /// the batch is searched in several units (`reader_batch.rs`); `None` when the unit is the
    /// whole batch.
    pub(crate) reader_shape: Option<std::sync::Arc<super::reader_batch::ReaderBatchShape>>,
    /// The bytes of NCBI's query buffer after this unit's own when the unit is a part of a reader
    /// batch: the letters and NULLBs of the contexts of the units after it
    /// (`reader_batch::attach`, `batch_buffer_length`); 0 when the unit is the whole batch, and
    /// for a unit prepared again after phase 1 failed (`native::XUnit::take_prep`).
    pub(crate) reader_buffer_after: usize,
    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): whether `sequence_start` is on the path of the seed
    /// loops' sentinel walks (`blastx/x_ungap_sent.rs`, `query_fits`, another matrix than BLOSUM62
    /// `query_fits_of`; false with the switch off).
    pub(crate) x_ungap_sent: bool,
}
impl PreparedQueryBatch {
    /// EXPERIMENT (LOSAT_X_BXENDMEM): what the results of a searched batch keep: the batch without
    /// its residues and lookup segments, which the report does not read.
    pub(crate) fn x_for_results(&self) -> Self {
        if !super::kappa::x_end_mem() {
            return self.clone();
        }
        Self {
            original_lengths: self.original_lengths.clone(),
            contexts: self.contexts.clone(),
            first_context: self.first_context,
            max_length: self.max_length,
            min_length: self.min_length,
            sequence_start: Vec::new(),
            sequence_start_nomask: Vec::new(),
            lookup_segments: Vec::new(),
            reader_shape: self.reader_shape.clone(),
            reader_buffer_after: self.reader_buffer_after,
            x_ungap_sent: self.x_ungap_sent,
        }
    }
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
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:136-139
// ```c++
//     case eBlastp:
//     default:
//         retval = 10000;
//         break;
// ```
/// The letters of RPSTBLASTN's query batch (`GetQueryBatchSize(eRPSTblastn)`, the `default:`
/// case; DECISIONS 38 D10): its batches are read with it, never split (`rps_search::translated`).
pub const RPSTBLASTN_BATCH_SIZE: usize = 10000;
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
pub fn batch_ranges(records: &[FastaRecord], batch_size: usize) -> Vec<std::ops::Range<usize>> {
    let mut result = Vec::new();
    let mut start = 0;
    while start < records.len() {
        let mut end = start;
        let mut length = 0;
        while end < records.len() && length < batch_size {
            length += records[end].sequence.len();
            end += 1;
        }
        result.push(start..end);
        start = end;
    }
    result
}
// No NCBI counterpart of its own: context `i` (frames 1, 2, 3, -1, -2, -3) of a query of `letters`
// nucleotides as `prepare_queries` lays the contexts out (blast_setup_cxx.cpp:191-219, cited
// there): the frame's protein length, and whether the context's strand is searched (a context
// of a strand not searched has length 0).
pub(crate) fn context_frame(letters: usize, i: usize, strand: &str) -> (usize, bool) {
    let len = letters.saturating_sub(i % 3) / 3;
    let enabled = match strand {
        "plus" => i < 3,
        "minus" => i >= 3,
        _ => true,
    };
    (len, enabled)
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
/// The error of a query batch with a `-matrix` that NCBI's tables do not have
/// (`ResolvedOptions::unknown_matrix`): NCBI's (`unknown_matrix_message`), or LOSAT's refusal when
/// NCBI reads a file for the name (`ResolvedOptions::matrix_file`): NCBI 2.17.0 reads it where it
/// sets up the batch, and aborts after a matrix or fails on another file (`matrix_file_refusal`).
/// A run without a batch reads no file.
pub(crate) fn matrix_error_message(records: &[FastaRecord], options: &ResolvedOptions) -> String {
    match &options.matrix_file {
        Some(file) => crate::blastinput::ncbi_environment::matrix_file_refusal(
            &options.matrix,
            file,
            "BLASTX",
        ),
        None => unknown_matrix_message(records),
    }
}

/// NCBI's error of a query batch whose score block cannot load the `-matrix` of an ungapped search
/// (a name its matrix loader does not know, `args.rs` `unknown_matrix`): the message of each query
/// of the batch, a query without letters with its own after it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 (the score block of
/// the batch: a name that the matrix loader does not know fails with -1, blast_stat.c:1679-1680,
/// posted for all queries, blast_setup.c:498-501)
/// ```c
///     Int2 status = BLAST_MainSetUp(opts_memento->m_ProgramType,
/// ...
///     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, search_messages);
///     if (status != 0 &&
///         (!(blast_msg.Get()) || (blast_msg.Get() && blast_msg.Get()->severity == eBlastSevError)))
///     {
///         	string msg;
///         	if (search_messages.HasMessages()) {
///             		msg = search_messages.ToString();
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:246 (the text of -1)
/// ```c
///             snprintf(buf, sizeof(buf) - 1, "Unknown error code %d", error_code);
/// ```
/// One message per query (`TSearchMessages::ToString`), `Error: ` first (`GetMessage`); a query
/// without letters keeps its `Sequence contains no data` (blast_setup_cxx.cpp:631-638) after the
/// error, which `RemoveDuplicates` sorts first (`CSearchMessage::operator<`,
/// blast_types.hpp:295-304), as TBLASTX's (`tblastx/blast_engine/run_impl.rs`).
pub(crate) fn unknown_matrix_message(records: &[FastaRecord]) -> String {
    let mut message = String::new();
    for record in records {
        message.push_str("Error: Unknown error code -1 ");
        if record.sequence.is_empty() {
            message.push_str("Warning: Sequence contains no data ");
        }
    }
    message
}

pub fn prepare_queries(
    records: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<PreparedQueryBatch> {
    prepare_unit_queries(records, None, options)
}

/// `prepare_queries` for `records`, one unit of the reader batch of `batch` when it is given
/// (`reader_batch.rs`): whether the setup fails for want of data is the reader batch's.
pub(crate) fn prepare_unit_queries(
    records: &[FastaRecord],
    batch: Option<&super::reader_batch::ReaderBatchShape>,
    options: &ResolvedOptions,
) -> Result<PreparedQueryBatch> {
    if records.is_empty() {
        bail!("cannot prepare an empty BLASTX query batch");
    }
    crate::blastinput::value_parsers::genetic_code(&options.query_gencode.to_string())
        .map_err(anyhow::Error::msg)?;
    let code = GeneticCode::try_from_id(options.query_gencode).map_err(anyhow::Error::msg)?;
    // Diagnostic (LOSAT_TIMING): the contexts' layout and the six-frame translation for the
    // common line (utils/xtiming.rs). No NCBI counterpart.
    let t_translate = crate::utils::xtiming::clock();
    let (mut contexts, mut max_length, mut min_length) = (Vec::new(), 0, usize::MAX);
    let mut offset = 0;
    for (q, record) in records.iter().enumerate() {
        for i in 0..6 {
            let frame = if i < 3 {
                (i + 1) as i8
            } else {
                -((i - 2) as i8)
            };
            let (len, enabled) = context_frame(record.sequence.len(), i, &options.strand);
            max_length = max_length.max(len);
            min_length = min_length.min(len);
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
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Translate, t_translate);
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
    // `qinfo` is the query batch's: a unit of a reader batch (`batch`) fails only when no
    // context of the reader batch is valid, with the messages of all its empty queries.
    let (fails, empty) = match batch {
        Some(batch) => (batch.set_up_fails(), batch.empty_queries),
        None => {
            let empty = records.iter().filter(|q| q.sequence.is_empty()).count();
            (!contexts.iter().any(|c| c.is_valid) && empty > 0, empty)
        }
    };
    if fails {
        anyhow::bail!("{}", "Warning: Sequence contains no data ".repeat(empty));
    }
    // The matrix error comes after the set-up's errors above, with the messages of the whole
    // query batch: a unit of a reader batch raises its batch's (`ReaderBatchShape::matrix_error`).
    if options.unknown_matrix {
        let message = match batch {
            Some(batch) => batch.matrix_error.as_deref().map(str::to_string),
            None => None,
        };
        anyhow::bail!(
            "{}",
            message.unwrap_or_else(|| matrix_error_message(records, options))
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
    let seg_masker = options.seg.filter_options().map(|(window, locut, hicut)| {
        SegMasker::new(
            if window > 0 { window as usize } else { 12 },
            if locut > 0.0 { locut } else { 2.2 },
            if hicut > 0.0 { hicut } else { 2.5 },
        )
        // LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): the chunks of a split query (split.rs) are masked
        // again after the whole query, mostly over the same segments, so this masker takes the
        // trim of a segment trimmed before on this thread from its memo.
        .x_with_trim_memo()
    });
    // SEG reads only the residues of its own context and the hard masking below writes only there,
    // so the contexts of a large batch (the six frames of the full query) are filtered on the
    // search pool before the loop consumes the intervals in context order.
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
    // context, so with a pool of 2+ threads the contexts of a large batch are filtered in parallel
    // and the intervals are kept in seg_ready; the loop below still consumes them in context order.
    #[allow(unused_mut)]
    let mut seg_ready: Vec<Option<Vec<(i32, i32)>>> = Vec::new();
    #[cfg(feature = "parallel")]
    if let Some(masker) = seg_masker.as_ref() {
        if super::runtime::x_bx_par_rule(|| {
            buffer.len() >= 60_000
                && !super::runtime::x_inner_serial()
                && rayon::current_thread_index().is_some()
                && rayon::current_num_threads() > 1
        }) {
            use rayon::prelude::*;
            let buffer = &buffer;
            // Diagnostic (LOSAT_TIMING): the parallel SEG's wall time for the common line
            // (utils/xtiming.rs). No NCBI counterpart.
            let t_seg = crate::utils::xtiming::clock();
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
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
        }
    }
    for (context_index, context) in contexts.iter_mut().enumerate() {
        if !context.is_valid {
            continue;
        }
        let mut masks = Vec::new();
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
                // Diagnostic (LOSAT_TIMING): this context's SEG for the common line
                // (utils/xtiming.rs). No NCBI counterpart.
                let t_seg = crate::utils::xtiming::clock();
                masks.extend(
                    masker
                        .mask_sequence(
                            &buffer[context.offset + 1..context.offset + 1 + context.length],
                        )
                        .into_iter()
                        .map(|r| (r.start as i32, r.end as i32 - 1)),
                );
                crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
            }
        }
        // SEG's intervals first, as `BlastSeqLocAppend(filter_out, lcase_mask_slp)`: the sort of
        // `BlastSeqLocCombine` compares only the left ends, and a one-letter run on the first letter
        // of a SEG interval is merged only in this order (the results' masks).
        masks.extend_from_slice(&context.lowercase_masks);
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
    // The sentinel walks take BLOSUM62's keys, or another matrix's where it has them (item SC6b,
    // `x_ungap_sent.rs` `query_fits_of`).
    let x_ungap_sent = if options.scoring == crate::config::ScoringMatrix::Blosum62 {
        super::seed::x_ungap_sent::query_fits(&buffer, options.word_size, options.window_size)
    } else {
        super::seed::x_ungap_sent::query_fits_of(
            options.scoring,
            &buffer,
            options.word_size,
            options.window_size,
        )
    };
    Ok(PreparedQueryBatch {
        original_lengths: records.iter().map(|r| r.sequence.len()).collect(),
        contexts,
        first_context: if options.strand == "minus" { 3 } else { 0 },
        max_length,
        min_length,
        sequence_start: buffer,
        sequence_start_nomask: nomask,
        lookup_segments,
        reader_shape: None,
        reader_buffer_after: 0,
        x_ungap_sent,
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
            source: None,
        };
        assert_eq!(
            batch_ranges(&[record(10001), record(1), record(3)], BATCH_SIZE),
            vec![0..2, 2..3]
        );
        assert_eq!(
            batch_ranges(&[record(10002), record(3)], BATCH_SIZE),
            vec![0..1, 1..2]
        );
        assert_eq!(
            batch_ranges(&[record(10003), record(3)], BATCH_SIZE),
            vec![0..1, 1..2]
        );
    }

    // `BATCH_SIZE` sets the limit (`GetQueryBatchSize`): the record that crosses it ends its
    // batch, and a limit above the default keeps the longer batch whole.
    #[test]
    fn the_batch_size_is_the_limit_of_a_batch() {
        let record = |n| FastaRecord {
            internal_id: String::new(),
            title: String::new(),
            sequence: vec![b'A'; n],
            lowercase_masks: vec![],
            warnings: vec![],
            source: None,
        };
        let records = [
            record(100),
            record(100),
            record(100),
            record(100),
            record(100),
        ];
        assert_eq!(batch_ranges(&records, 200), vec![0..2, 2..4, 4..5]);
        assert_eq!(batch_ranges(&records, 201), vec![0..3, 3..5]);
        assert_eq!(
            batch_ranges(&records, 1),
            vec![0..1, 1..2, 2..3, 3..4, 4..5]
        );
        assert_eq!(batch_ranges(&records, 30000), vec![0..5]);
    }
}
