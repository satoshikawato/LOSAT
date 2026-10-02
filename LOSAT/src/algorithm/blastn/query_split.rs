//! NCBI's query splitting for BLASTN (`CQuerySplitter`): NCBI searches a query batch of at
//! least two chunks chunk by chunk in its preliminary stage, maps the HSPs of each chunk onto
//! the contexts of the batch and merges them (`run.rs`), and then runs the traceback on the
//! whole batch.

use crate::utils::dust::MaskedInterval;

/// The query chunk size of the task when the environment variable `CHUNK_SIZE` is not set
/// or blank (`SplitQuery_GetChunkSize`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:65-73
/// ```c
///         switch (program) {
///         case eBlastn:
///             retval = 1000000;
///             break;
///         case eMegablast:
///         case eDiscMegablast:
///         case eMapper:
///             retval = 5000000;
///             break;
/// ```
pub fn query_chunk_size(megablast: bool) -> usize {
    if megablast {
        5_000_000
    } else {
        1_000_000
    }
}

/// The overlap of two query chunks for a nucleotide query when the environment variable
/// `OVERLAP_CHUNK_SIZE` is not set or blank.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:51-53
/// ```c
/// SplitQuery_GetOverlapChunkSize(EBlastProgramType program)
/// {
///     size_t retval = 100;
/// ```
pub const QUERY_CHUNK_OVERLAP: usize = 100;

/// NCBI's query chunk size and chunk overlap (`size_t` values), from the environment
/// variables `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` when they are set and not blank.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:54-62
/// ```c
/// size_t
/// SplitQuery_GetChunkSize(EProgram program)
/// {
///     size_t retval = 0;
///
///     // used for experimentation purposes
///     char* chunk_sz_str = getenv("CHUNK_SIZE");
///     if (chunk_sz_str && !NStr::IsBlank(chunk_sz_str)) {
///         retval = NStr::StringToInt(chunk_sz_str);
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:53-60
/// ```c
///     size_t retval = 100;
///     // used for experimentation purposes
///     char* overlap_sz_str = getenv("OVERLAP_CHUNK_SIZE");
///     if (overlap_sz_str && !NStr::IsBlank(overlap_sz_str)) {
///         retval = NStr::StringToInt(overlap_sz_str);
///         _TRACE("Using overlap chunk size from environment " << retval);
///         return retval;
///     }
/// ```
/// The `int` converts to the 64-bit `size_t` of the oracle: a negative value becomes
/// 2^64 plus the value. `u64` keeps that arithmetic on every target.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SplitSizes {
    pub chunk_size: u64,
    pub overlap: u64,
}

impl SplitSizes {
    /// The sizes of a task, with the values of the environment variables (already read
    /// with `NStr::StringToInt`; `None` when unset or blank).
    pub fn new(megablast: bool, chunk_size: Option<i32>, overlap: Option<i32>) -> Self {
        Self {
            chunk_size: chunk_size.map_or(query_chunk_size(megablast) as u64, |value| {
                i64::from(value) as u64
            }),
            overlap: overlap.map_or(QUERY_CHUNK_OVERLAP as u64, |value| i64::from(value) as u64),
        }
    }

    /// The default sizes of a task (neither variable set).
    pub fn default_for(megablast: bool) -> Self {
        Self::new(megablast, None, None)
    }

    /// Whether `CHUNK_SIZE` was negative (`SplitQuery_GetChunkSize` returns it as a
    /// `size_t`, near 2^64).
    pub fn negative_chunk_size(self) -> bool {
        (self.chunk_size as i64) < 0
    }

    /// NCBI's `CBatchSizeMixer` maximum: the chunk size less 1000 (`size_t`), converted to
    /// its `Int4` parameter.
    ///
    /// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261
    /// ```c
    ///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
    /// ```
    pub fn mixer_max_batch_size(&self) -> i32 {
        self.chunk_size.wrapping_sub(1000) as u32 as i32
    }
}

/// One query's part in a query chunk: the query of the batch and its residues `from..to`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChunkQuery {
    pub query: usize,
    pub from: usize,
    pub to: usize,
}

impl ChunkQuery {
    fn len(&self) -> usize {
        self.to - self.from
    }
}

/// A query chunk of a split batch (`CSplitQueryBlk`): its query parts in query order, and
/// for each of its contexts (two per part, the plus strand first) the context of the batch
/// and the offset that NCBI adds to the query coordinates of the chunk's HSPs in it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryChunk {
    pub queries: Vec<ChunkQuery>,
    pub contexts: Vec<usize>,
    pub context_offsets: Vec<i32>,
}

/// The number of query chunks of a batch of `concatenated_query_length` residues, and the
/// chunk size NCBI then uses.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
/// ```c
///     // Fix for small query size
///     if ((*chunk_size) > overlap_size) {
///        num_chunks = concatenated_query_length / ((*chunk_size) - overlap_size);
///     }
///
///     // Only one chunk, just return;
///     if (num_chunks <= 1) {
///        *chunk_size = concatenated_query_length;
///        return 1;
///     }
///
///     // Re-adjust the chunk_size to make load even
///     if (!Blast_QueryIsTranslated(program)) {
///        *chunk_size = (concatenated_query_length + (num_chunks - 1) * overlap_size) / num_chunks;
///        // Round up only if this will not decrease the number of chunks
///        if (num_chunks < (*chunk_size) - overlap_size ) (*chunk_size)++;
/// ```
pub fn calculate_num_chunks(sizes: SplitSizes, concatenated_query_length: usize) -> (usize, usize) {
    let overlap_size = sizes.overlap;
    let length = concatenated_query_length as u64;
    // `Uint4 num_chunks` (split_query_aux_priv.cpp:112).
    let mut num_chunks: u32 = 0;
    if sizes.chunk_size > overlap_size {
        num_chunks = (length / (sizes.chunk_size - overlap_size)) as u32;
    }
    if num_chunks <= 1 {
        return (1, concatenated_query_length);
    }
    let chunks = u64::from(num_chunks);
    let mut chunk_size = length.wrapping_add((chunks - 1).wrapping_mul(overlap_size)) / chunks;
    if chunks < chunk_size.wrapping_sub(overlap_size) {
        chunk_size += 1;
    }
    (num_chunks as usize, chunk_size as usize)
}

/// Whether NCBI's setup of some query chunk would split that chunk again: each chunk's
/// search builds its own query splitter with the same chunk size and overlap, and only a
/// debug build asserts that it does not split. When it does (an overlap close to the chunk
/// size), the release build skips the chunk's lookup table and stops with a null reference.
/// LOSAT searches such a chunk once (approved exception 1 of PD-LOSAT-NCBI-DEFECTS), so
/// this only tells where the exception applies.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
/// ```c
///     BlastSeqSrc* seqsrc =
///         BlastSeqSrcCopy(full_data->m_SeqSrc->GetPointer());
///     CRef<SBlastSetupData> setup_data =
///         BlastSetupPreliminarySearchEx(
///                 qf, options,
///                 CRef<objects::CPssmWithParameters>(),
///                 seqsrc, num_threads);
///     BlastSeqSrcResetChunkIterator(seqsrc);
///     setup_data->m_InternalData->m_SeqSrc.Reset(new TBlastSeqSrc(seqsrc,
///                                                BlastSeqSrcFree));
///
///     _ASSERT(setup_data->m_QuerySplitter->IsQuerySplit() == false);
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux_priv.cpp:206-208
/// ```c
///     // 5. Create the lookup table
///     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
///         LookupTableWrap* lut =
/// ```
pub fn chunk_would_be_split(chunks: &[QueryChunk], sizes: SplitSizes) -> bool {
    chunks.iter().any(|chunk| {
        let length = chunk.queries.iter().map(ChunkQuery::len).sum();
        calculate_num_chunks(sizes, length).0 > 1
    })
}

/// The query chunks of a batch with queries of `lengths` residues, or `None` when NCBI
/// does not split it.
///
/// The chunk ranges on the concatenated queries:
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:145-171
/// ```c
///     size_t chunk_start = 0;
///     const size_t kOverlapSize =
///         SplitQuery_GetOverlapChunkSize(m_Options->GetProgramType());
///     for (size_t chunk_num = 0; chunk_num < m_NumChunks; chunk_num++) {
///         size_t chunk_end = chunk_start + m_ChunkSize;
///
///         // if the chunk end is larger than the sequence ...
///         if (chunk_end >= m_TotalQueryLength ||
///             // ... or this is the last chunk and it didn't make it to the end
///             // of the sequence
///             (chunk_end < m_TotalQueryLength && (chunk_num + 1) == m_NumChunks))
///         {
///             // ... assign this chunk's end to the end of the sequence
///             chunk_end = m_TotalQueryLength;
///         }
///
///         m_SplitBlk->SetChunkBounds(chunk_num,
///                                    TChunkRange(static_cast<unsigned int>(chunk_start), static_cast<unsigned int>(chunk_end)));
///         _TRACE("Chunk " << chunk_num << ": ranges from " << chunk_start
///                << " to " << chunk_end);
///
///         chunk_start += (m_ChunkSize - kOverlapSize);
///         if (chunk_start > m_TotalQueryLength ||
///             chunk_end == m_TotalQueryLength) {
///             break;
///         }
///     }
/// ```
/// A chunk the loop leaves unset keeps the empty range of `SSplitQueryBlkNew`, and has no
/// query (NCBI then fails to make its query factory, "Empty CBlastQueryVector"); the chunk
/// size of `calculate_num_chunks` leaves none.
///
/// The queries of a chunk, and their parts in it:
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:223-247
/// ```c
///     query_ranges.push_back(TChunkRange(0, static_cast<unsigned int>(m_LocalQueryData->GetSeqLength(0)))); // FIXED
///     ...
///     for (int i = 1; i < kNumQueries; i++) {
///         TSeqPos query_start = query_ranges[i-1].GetTo() + 1;
///         TSeqPos query_end = query_start + static_cast<TSeqPos>(m_LocalQueryData->GetSeqLength(i));
///         query_ranges.push_back(TChunkRange(query_start, query_end));
///     ...
///     for (size_t chunk_num = 0; chunk_num < m_NumChunks; chunk_num++) {
///         const TChunkRange chunk = m_SplitBlk->GetChunkBounds(chunk_num);
///         ...
///         for (size_t qindex = 0; qindex < query_ranges.size(); qindex++) {
///             const TChunkRange& query_range = query_ranges[qindex];
///             if ( !chunk.IntersectingWith(query_range) ) {
///                 continue;
///             }
///             m_SplitBlk->AddQueryToChunk(chunk_num, static_cast<Int4>(qindex));
/// ```
/// (`TChunkRange` is a `COpenRange`, so `GetTo() + 1` is the end of the previous query.)
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:196-210
/// ```c
///     CSeq_interval& interval = split_query_loc->SetInt();
///     const int qstart = chunk.GetFrom() - query_range.GetFrom();
///     const int qend = chunk.GetToOpen() - query_range.GetToOpen();
///
///     interval.SetFrom(max(0, qstart) + query_offset);
///
///     if (qend >= 0) {
///         interval.SetTo(query_range.GetToOpen() - query_range.GetFrom() + query_offset);
///     } else {
///         interval.SetTo(chunk.GetToOpen() - query_range.GetFrom() + query_offset);
///     }
///
///     // Note subtraction, as Seq-intervals are assumed to be
///     // open/inclusive
///     interval.SetTo() -= 1;
/// ```
/// The query location of a FASTA query starts at 0, so `query_offset` is 0.
///
/// The contexts of a chunk (both strands, as LOSAT has no `-strand`):
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:392-410
/// ```c
///             } else if (Blast_QueryIsNucleotide(kProgram)) {
///
///                 for (unsigned int ctx = 0; ctx < kNumContexts; ctx++) {
///                     // handle the plus strand...
///                     if (ctx % NUM_STRANDS == 0) {
///                         ...
///                             m_SplitBlk->AddContextToChunk(chunk_num,
///                                               static_cast<Int4>(kNumContexts*queries[i]+ctx));
///                         }
///                     } else { // handle the negative strand
///                         ...
///                             m_SplitBlk->AddContextToChunk(chunk_num,
///                                               static_cast<Int4>(kNumContexts*queries[i]+ctx));
/// ```
/// The offsets are `context_offsets`.
pub fn split_query_batch(lengths: &[usize], sizes: SplitSizes) -> Option<Vec<QueryChunk>> {
    let total_length: usize = lengths.iter().sum();
    let (num_chunks, chunk_size) = calculate_num_chunks(sizes, total_length);
    if num_chunks <= 1 {
        return None;
    }
    // A split needs `chunk_size > overlap`, so the overlap is an ordinary `int` value here.
    let overlap_size = sizes.overlap as usize;

    let mut chunk_ranges = vec![(0usize, 0usize); num_chunks];
    let mut chunk_start = 0usize;
    for (chunk_num, range) in chunk_ranges.iter_mut().enumerate() {
        let mut chunk_end = chunk_start + chunk_size;
        if chunk_end >= total_length || (chunk_end < total_length && chunk_num + 1 == num_chunks) {
            chunk_end = total_length;
        }
        *range = (chunk_start, chunk_end);
        chunk_start += chunk_size - overlap_size;
        if chunk_start > total_length || chunk_end == total_length {
            break;
        }
    }

    let mut query_ranges = Vec::with_capacity(lengths.len());
    let mut query_start = 0usize;
    for &length in lengths {
        query_ranges.push((query_start, query_start + length));
        query_start += length;
    }

    let mut chunks: Vec<QueryChunk> = chunk_ranges
        .iter()
        .map(|&(chunk_from, chunk_to)| {
            let mut queries = Vec::new();
            let mut contexts = Vec::new();
            for (query, &(query_from, query_to)) in query_ranges.iter().enumerate() {
                if chunk_from.max(query_from) >= chunk_to.min(query_to) {
                    continue;
                }
                let qstart = chunk_from as i64 - query_from as i64;
                let qend = chunk_to as i64 - query_to as i64;
                let from = qstart.max(0) as usize;
                let to = if qend >= 0 {
                    query_to - query_from
                } else {
                    chunk_to - query_from
                };
                queries.push(ChunkQuery { query, from, to });
                contexts.push(2 * query);
                contexts.push(2 * query + 1);
            }
            QueryChunk {
                queries,
                contexts,
                context_offsets: Vec::new(),
            }
        })
        .collect();

    let offsets: Vec<Vec<i32>> = (0..chunks.len())
        .map(|chunk_num| context_offsets(&chunks, chunk_num, lengths, overlap_size))
        .collect();
    for (chunk, offsets) in chunks.iter_mut().zip(offsets) {
        chunk.context_offsets = offsets;
    }
    Some(chunks)
}

/// The position of the batch context `absolute_context` among the contexts of a chunk.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:251-264
/// ```c
/// CContextTranslator::GetContextInChunk(size_t chunk_num,
///                                       int absolute_context) const
/// {
///     _ASSERT(chunk_num < m_ContextsPerChunk.size());
///     const vector<int>& context_indices = m_ContextsPerChunk[chunk_num];
///     vector<int>::const_iterator itr = find(context_indices.begin(),
///                                            context_indices.end(),
///                                            absolute_context);
///     if (itr == context_indices.end()) {
///         return kInvalidContext;
///     }
///     return static_cast<int>(itr - context_indices.begin());  // FIXED
/// }
/// ```
fn context_in_chunk(chunk: &QueryChunk, absolute_context: usize) -> Option<usize> {
    chunk
        .contexts
        .iter()
        .position(|&context| context == absolute_context)
}

/// The length of a batch context in a chunk (0 when the chunk does not have it).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:452-466
/// ```c
/// s_GetAbsoluteContextLength(const vector<const BlastQueryInfo*>& chunk_qinfo,
///                            int chunk_num,
///                            const CContextTranslator& ctx_translator,
///                            int absolute_context)
/// {
///     if (chunk_num < 0) {
///         return 0;
///     }
///
///     int pos = ctx_translator.GetContextInChunk((size_t)chunk_num,
///                                                absolute_context);
///     if (pos != kInvalidContext) {
///         return chunk_qinfo[chunk_num]->contexts[pos].query_length;
///     }
///     return 0;
/// ```
fn absolute_context_length(
    chunks: &[QueryChunk],
    chunk_num: isize,
    absolute_context: usize,
) -> usize {
    if chunk_num < 0 {
        return 0;
    }
    let chunk = &chunks[chunk_num as usize];
    context_in_chunk(chunk, absolute_context)
        .map(|position| chunk.queries[position / 2].len())
        .unwrap_or(0)
}

/// The first chunk of the run of chunks, ending with `curr_chunk`, that have the batch
/// context `absolute_context`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:266-284
/// ```c
/// CContextTranslator::GetStartingChunk(size_t curr_chunk,
///                                      Int4 context_in_chunk) const
/// {
///     int absolute_context = GetAbsoluteContext(curr_chunk, context_in_chunk);
///     if (absolute_context == kInvalidContext) {
///         return kInvalidContext;
///     }
///
///     size_t retval = curr_chunk;
///
///     for (--curr_chunk; static_cast<int>(curr_chunk) >= 0; --curr_chunk) {
///         if (GetContextInChunk(curr_chunk, absolute_context) ==
///             kInvalidContext) {
///             break;
///         }
///         retval = curr_chunk;
///     }
///     return static_cast<int>(retval);
/// }
/// ```
fn starting_chunk(chunks: &[QueryChunk], curr_chunk: usize, absolute_context: usize) -> usize {
    let mut retval = curr_chunk;
    for chunk_num in (0..curr_chunk).rev() {
        if context_in_chunk(&chunks[chunk_num], absolute_context).is_none() {
            break;
        }
        retval = chunk_num;
    }
    retval
}

/// The offsets of the contexts of a chunk (`size_t` arithmetic, stored as `Int4`): where the
/// chunk's part of a context starts in the batch's context (on the minus strand, from the
/// query's end).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:612-662
/// ```c
///         for (Int4 ctx = chunk_qinfo[chunk_num]->first_context;
///              ctx <= chunk_qinfo[chunk_num]->last_context;
///              ctx++) {
///
///             size_t correction = 0;
///             const int starting_chunk =
///                 ctx_translator.GetStartingChunk(chunk_num, ctx);
///             const int absolute_context =
///                 ctx_translator.GetAbsoluteContext(chunk_num, ctx);
///             ...
///             if (s_IsPlusStrand(chunk_qinfo[chunk_num], ctx)) {
///
///                 for (int c = static_cast<int>(chunk_num); c != starting_chunk; c--) {
///                     size_t prev_len = s_GetAbsoluteContextLength(chunk_qinfo,
///                                                          c - 1,
///                                                          ctx_translator,
///                                                          absolute_context);
///                     size_t curr_len = s_GetAbsoluteContextLength(chunk_qinfo, c,
///                                                          ctx_translator,
///                                                          absolute_context);
///                     size_t overlap = min(kOverlap, curr_len);
///                     correction += prev_len - min(overlap, prev_len);
///                 }
///
///             } else {
///
///                 size_t subtrahend = 0;
///
///                 for (int c = static_cast<int>(chunk_num); c >= starting_chunk && c >= 0; c--) {
///                     size_t prev_len = s_GetAbsoluteContextLength(chunk_qinfo,
///                                                          c - 1,
///                                                          ctx_translator,
///                                                          absolute_context);
///                     size_t curr_len = s_GetAbsoluteContextLength(chunk_qinfo,
///                                                          c,
///                                                          ctx_translator,
///                                                          absolute_context);
///                     size_t overlap = min(kOverlap, curr_len);
///                     subtrahend += (curr_len - min(overlap, prev_len));
///                 }
///                 correction =
///                     global_qinfo->contexts[absolute_context].query_length -
///                     subtrahend;
/// ```
fn context_offsets(
    chunks: &[QueryChunk],
    chunk_num: usize,
    lengths: &[usize],
    overlap_size: usize,
) -> Vec<i32> {
    chunks[chunk_num]
        .contexts
        .iter()
        .map(|&absolute_context| {
            let starting = starting_chunk(chunks, chunk_num, absolute_context);
            let mut correction = 0usize;
            if absolute_context % 2 == 0 {
                let mut c = chunk_num as isize;
                while c != starting as isize {
                    let prev_len = absolute_context_length(chunks, c - 1, absolute_context);
                    let curr_len = absolute_context_length(chunks, c, absolute_context);
                    let overlap = overlap_size.min(curr_len);
                    correction = correction.wrapping_add(prev_len - overlap.min(prev_len));
                    c -= 1;
                }
            } else {
                let mut subtrahend = 0usize;
                let mut c = chunk_num as isize;
                while c >= starting as isize && c >= 0 {
                    let prev_len = absolute_context_length(chunks, c - 1, absolute_context);
                    let curr_len = absolute_context_length(chunks, c, absolute_context);
                    let overlap = overlap_size.min(curr_len);
                    subtrahend = subtrahend.wrapping_add(curr_len - overlap.min(prev_len));
                    c -= 1;
                }
                correction = lengths[absolute_context / 2].wrapping_sub(subtrahend);
            }
            correction as i32
        })
        .collect()
}

/// The masks of a query part: NCBI restricts the query's masks (its DUST and lower-case
/// masks, `m_UserSpecifiedMasks`) to the part, and adds the DUST masks of the part itself
/// (`run.rs`). The restriction drops the last residue of the part from its range, so a mask
/// that starts there is dropped, and extends every other mask by one residue (to the end of
/// the part at most).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:278-281
/// ```c
///             CRef<CSeq_loc> mask_query_loc(new CSeq_loc);
///             s_SetSplitQuerySeqInterval(chunk, query_range, 0, mask_query_loc);
///             TMaskedQueryRegions split_mask =
///                 m_UserSpecifiedMasks[qindex].RestrictToSeqInt(mask_query_loc->GetInt());
/// ```
/// NCBI reference: ncbi-blast/c++/src/objects/seq/seqlocinfo.cpp:102-116
/// ```c
///     TSeqRange loc(location.GetFrom(), 0);
///     loc.SetToOpen(location.GetTo());
///
///     ITERATE(TMaskedQueryRegions, maskinfo, *this) {
///         const CSeq_interval& intv = (*maskinfo)->GetInterval();
///         TSeqRange mask (intv.GetFrom(), intv.GetTo());
///         TSeqRange range = loc.IntersectionWith(mask);
///         if (range.NotEmpty()) {
///             ...
///             CRef<CSeq_interval> si
///                 (new CSeq_interval(const_cast<CSeq_id&>(intv.GetId()),
///                                    range.GetFrom(),
///                                    range.GetToOpen(),
///                                    kStrand));
/// ```
/// (`TSeqRange` is a `CRange`, whose constructor takes the last position.) The chunk's query
/// block then puts the masks in the part's coordinates:
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:1041-1047
/// ```c
///    for (seqloc = *mask; seqloc; seqloc = next_loc) {
///       next_loc = seqloc->next;
///       seqloc->ssr->left = MAX(0, seqloc->ssr->left - from);
///       seqloc->ssr->right = MIN(seqloc->ssr->right, to) - from;
///       /* If this mask location does not intersect the [from,to] interval,
///          do not add it to the newly constructed list and free its contents. */
///       if (seqloc->ssr->left > seqloc->ssr->right) {
/// ```
pub fn restrict_masks(masks: &[MaskedInterval], part: &ChunkQuery) -> Vec<MaskedInterval> {
    // The part is the Seq-interval [from, to - 1]; `loc` is [from, to - 1).
    let last = part.to - 1;
    masks
        .iter()
        .filter_map(|mask| {
            // `TSeqRange mask` is [start, end) (the mask's last position is end - 1).
            let left = mask.start.max(part.from);
            let right_open = mask.end.min(last);
            // The new Seq-interval is [left, right_open], within [from, to - 1].
            (left < right_open)
                .then(|| MaskedInterval::new(left - part.from, right_open - part.from + 1))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
    // ```c
    //     if ((*chunk_size) > overlap_size) {
    //        num_chunks = concatenated_query_length / ((*chunk_size) - overlap_size);
    //     }
    //        *chunk_size = (concatenated_query_length + (num_chunks - 1) * overlap_size) / num_chunks;
    //        if (num_chunks < (*chunk_size) - overlap_size ) (*chunk_size)++;
    // ```
    #[test]
    fn the_chunk_count_and_size_are_ncbis() {
        assert_eq!(
            calculate_num_chunks(SplitSizes::default_for(false), 1_999_799),
            (1, 1_999_799)
        );
        assert_eq!(
            calculate_num_chunks(SplitSizes::default_for(false), 1_999_800),
            (2, 999_951)
        );
        // EDL933 with -task blastn: five chunks of 1,105,770 (E2f §D).
        assert_eq!(
            calculate_num_chunks(SplitSizes::default_for(false), 5_528_445),
            (5, 1_105_770)
        );
        assert_eq!(
            calculate_num_chunks(SplitSizes::default_for(true), 5_528_445),
            (1, 5_528_445)
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:612-662
    // ```c
    //                     correction += prev_len - min(overlap, prev_len);
    //                     subtrahend += (curr_len - min(overlap, prev_len));
    //                 correction =
    //                     global_qinfo->contexts[absolute_context].query_length -
    //                     subtrahend;
    // ```
    #[test]
    fn a_query_over_two_chunks_is_mapped_onto_its_contexts() {
        let chunks =
            split_query_batch(&[1_999_800], SplitSizes::default_for(false)).expect("split");
        assert_eq!(chunks.len(), 2);
        assert_eq!(
            chunks[0].queries,
            vec![ChunkQuery {
                query: 0,
                from: 0,
                to: 999_951
            }]
        );
        assert_eq!(
            chunks[1].queries,
            vec![ChunkQuery {
                query: 0,
                from: 999_851,
                to: 1_999_800
            }]
        );
        assert_eq!(chunks[0].contexts, vec![0, 1]);
        assert_eq!(chunks[0].context_offsets, vec![0, 1_999_800 - 999_951]);
        assert_eq!(chunks[1].context_offsets, vec![999_851, 0]);
    }

    #[test]
    fn parts_that_end_or_start_in_an_overlap_keep_their_starts() {
        // Chunks of 1,000,059 overlapping by 100: the first query ends 58 residues into the
        // overlap of the second and third chunks, and the second query starts there.
        let chunks = split_query_batch(&[1_999_976, 1_000_000], SplitSizes::default_for(false))
            .expect("split");
        assert_eq!(chunks.len(), 3);
        assert_eq!(
            chunks[1].queries,
            vec![
                ChunkQuery {
                    query: 0,
                    from: 999_959,
                    to: 1_999_976
                },
                ChunkQuery {
                    query: 1,
                    from: 0,
                    to: 42
                },
            ]
        );
        assert_eq!(
            chunks[2].queries,
            vec![
                ChunkQuery {
                    query: 0,
                    from: 1_999_918,
                    to: 1_999_976
                },
                ChunkQuery {
                    query: 1,
                    from: 0,
                    to: 1_000_000
                },
            ]
        );
        assert_eq!(chunks[1].contexts, vec![0, 1, 2, 3]);
        assert_eq!(
            chunks[1].context_offsets,
            vec![999_959, 0, 0, 1_000_000 - 42]
        );
        assert_eq!(chunks[2].context_offsets, vec![1_999_918, 0, 0, 0]);
    }

    // NCBI reference: ncbi-blast/c++/src/objects/seq/seqlocinfo.cpp:102-116
    // ```c
    //     TSeqRange loc(location.GetFrom(), 0);
    //     loc.SetToOpen(location.GetTo());
    //         TSeqRange range = loc.IntersectionWith(mask);
    //                                    range.GetToOpen(),
    // ```
    #[test]
    fn restricted_masks_grow_by_one_and_drop_one_at_the_last_residue() {
        let part = ChunkQuery {
            query: 0,
            from: 100,
            to: 200,
        };
        let masks = [
            MaskedInterval::new(0, 50),
            MaskedInterval::new(90, 110),
            MaskedInterval::new(150, 160),
            MaskedInterval::new(190, 300),
            MaskedInterval::new(199, 250),
        ];
        assert_eq!(
            restrict_masks(&masks, &part),
            vec![
                MaskedInterval::new(0, 11),
                MaskedInterval::new(50, 61),
                MaskedInterval::new(90, 100),
            ]
        );
    }

    // CHUNK_SIZE and OVERLAP_CHUNK_SIZE as NCBI reads them (split_query_aux_priv.cpp:100-138).
    #[test]
    fn environment_sizes_follow_ncbis_size_t_arithmetic() {
        let sizes = SplitSizes::new(false, Some(500), None);
        assert_eq!(calculate_num_chunks(sizes, 100_000), (250, 500));
        let sizes = SplitSizes::new(false, Some(40_000), Some(0));
        assert_eq!(calculate_num_chunks(sizes, 100_000), (2, 50_001));
        // An overlap at least the chunk size, a negative overlap (2^64 - 1) and a negative
        // chunk size (2^64 - 5) never split.
        assert_eq!(
            calculate_num_chunks(SplitSizes::new(false, Some(150), Some(150)), 100_000),
            (1, 100_000)
        );
        assert_eq!(
            calculate_num_chunks(SplitSizes::new(false, Some(40_000), Some(-1)), 100_000),
            (1, 100_000)
        );
        assert_eq!(
            calculate_num_chunks(SplitSizes::new(false, Some(-5), None), 100_000),
            (1, 100_000)
        );
        // The CBatchSizeMixer maximum: size_t chunk size - 1000 as Int4.
        assert_eq!(
            SplitSizes::default_for(false).mixer_max_batch_size(),
            999_000
        );
        assert_eq!(
            SplitSizes::default_for(true).mixer_max_batch_size(),
            4_999_000
        );
        assert_eq!(
            SplitSizes::new(false, Some(1001), None).mixer_max_batch_size(),
            1
        );
        assert_eq!(
            SplitSizes::new(false, Some(500), None).mixer_max_batch_size(),
            -500
        );
        assert_eq!(
            SplitSizes::new(false, Some(-5), None).mixer_max_batch_size(),
            -1005
        );
        assert_eq!(
            SplitSizes::new(false, Some(i32::MIN), None).mixer_max_batch_size(),
            i32::MAX - 999
        );
    }

    #[test]
    fn small_chunks_cover_the_queries_with_the_overlap() {
        let sizes = SplitSizes::new(false, Some(1200), Some(50));
        let chunks = split_query_batch(&[1000, 1500, 700], sizes).expect("split");
        assert!(chunks.len() > 1);
        for chunk in &chunks {
            assert!(!chunk.queries.is_empty());
        }
    }

    // Oracle (E2g audits b and round 2): a negative chunk size above a negative overlap
    // splits a batch only when the wrapped difference fits in it at least twice; -1 over
    // INT_MIN does not split 12000 residues (NCBI searches as without the variables), -5
    // over -10 does.
    #[test]
    fn a_negative_chunk_size_splits_only_a_long_enough_batch() {
        let wide = SplitSizes::new(false, Some(-1), Some(i32::MIN));
        assert!(wide.negative_chunk_size());
        assert_eq!(calculate_num_chunks(wide, 12_000).0, 1);
        let narrow = SplitSizes::new(false, Some(-5), Some(-10));
        assert!(calculate_num_chunks(narrow, 12_000).0 > 1);
        assert!(!SplitSizes::new(false, Some(3000), Some(-10)).negative_chunk_size());
    }

    // Oracle (E2g audit a): a 30000-residue query with CHUNK_SIZE 3000 fails from
    // OVERLAP_CHUNK_SIZE 1526 (chunks of 2950 over 1474) and runs at 1525.
    #[test]
    fn a_chunk_that_ncbi_would_split_again() {
        let fails = SplitSizes::new(false, Some(3000), Some(1526));
        let chunks = split_query_batch(&[30_000], fails).expect("split");
        assert!(chunk_would_be_split(&chunks, fails));
        let runs = SplitSizes::new(false, Some(3000), Some(1525));
        let chunks = split_query_batch(&[30_000], runs).expect("split");
        assert!(!chunk_would_be_split(&chunks, runs));
        let default = SplitSizes::default_for(false);
        let chunks = split_query_batch(&[1_999_800], default).expect("split");
        assert!(!chunk_would_be_split(&chunks, default));
    }
}
