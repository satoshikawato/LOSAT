//! NCBI's query splitting for BLASTP and TBLASTN (`CQuerySplitter`): NCBI searches a protein
//! query batch of at least two chunks chunk by chunk in its preliminary stage, maps the HSPs
//! of each chunk onto the queries of the batch, merges them, and then runs the traceback on
//! the whole batch. A protein query has one context, so the chunk of a query part is that
//! query's context (BLASTN's two strands are `blastn/query_split.rs`).

use crate::algorithm::blastn::query_split::{
    calculate_num_chunks, restrict_masks, ChunkQuery, SplitSizes, QUERY_CHUNK_OVERLAP,
};
use crate::utils::dust::MaskedInterval;

/// The query chunk size of blastp (10000) or tblastn (20000) when the environment variable
/// `CHUNK_SIZE` is not set or blank, and the overlap (100) when `OVERLAP_CHUNK_SIZE` is not.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:74-76,91-94
/// ```c
///         case eTblastn:
///             retval = 20000;
///             break;
///         ...
///         case eBlastp:
///         default:
///             retval = 10000;
///             break;
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:51-60
/// ```c
/// SplitQuery_GetOverlapChunkSize(EBlastProgramType program)
/// {
///     size_t retval = 100;
///     // used for experimentation purposes
///     char* overlap_sz_str = getenv("OVERLAP_CHUNK_SIZE");
///     if (overlap_sz_str && !NStr::IsBlank(overlap_sz_str)) {
///         retval = NStr::StringToInt(overlap_sz_str);
/// ```
/// `chunk_size` and `overlap` are the values of the environment variables, already read
/// with `NStr::StringToInt` (`None` when unset or blank); the `int` becomes NCBI's 64-bit
/// `size_t`.
pub fn protein_split_sizes(
    default_chunk_size: u32,
    chunk_size: Option<i32>,
    overlap: Option<i32>,
) -> SplitSizes {
    SplitSizes {
        chunk_size: chunk_size.map_or(u64::from(default_chunk_size), |value| {
            i64::from(value) as u64
        }),
        overlap: overlap.map_or(QUERY_CHUNK_OVERLAP as u64, |value| i64::from(value) as u64),
    }
}

/// A query chunk of a split protein batch (`CSplitQueryBlk`): its query parts in query order
/// (the part of query `queries[i].query` is the chunk's context `i`, whose context in the
/// batch is that query) and, for each, the offset that NCBI adds to the query coordinates of
/// the chunk's HSPs in it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProteinQueryChunk {
    pub queries: Vec<ChunkQuery>,
    pub context_offsets: Vec<i32>,
}

/// The query chunks of a protein batch with queries of `lengths` residues, or `None` when
/// NCBI does not split it.
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
/// The queries of a chunk and their parts are those of `blastn/query_split.rs`
/// `split_query_batch` (split_query_cxx.cpp:196-247). A protein query has one context:
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:415-418
/// ```c
///             } else if (Blast_QueryIsProtein(kProgram)) {
///                 m_SplitBlk->AddContextToChunk(chunk_num,
///                                               static_cast<Int4>(kNumContexts*queries[i]));
///             } else {
/// ```
/// The offsets are `context_offsets`. A chunk that the loop leaves unset keeps the empty
/// range of `SSplitQueryBlkNew` and has no query part.
pub fn split_protein_batch(lengths: &[usize], sizes: SplitSizes) -> Option<Vec<ProteinQueryChunk>> {
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

    let parts: Vec<Vec<ChunkQuery>> = chunk_ranges
        .iter()
        .map(|&(chunk_from, chunk_to)| {
            query_ranges
                .iter()
                .enumerate()
                .filter(|(_, &(query_from, query_to))| {
                    chunk_from.max(query_from) < chunk_to.min(query_to)
                })
                .map(|(query, &(query_from, query_to))| {
                    let qstart = chunk_from as i64 - query_from as i64;
                    let qend = chunk_to as i64 - query_to as i64;
                    let from = qstart.max(0) as usize;
                    let to = if qend >= 0 {
                        query_to - query_from
                    } else {
                        chunk_to - query_from
                    };
                    ChunkQuery { query, from, to }
                })
                .collect()
        })
        .collect();
    Some(
        (0..parts.len())
            .map(|chunk_num| ProteinQueryChunk {
                queries: parts[chunk_num].clone(),
                context_offsets: context_offsets(&parts, chunk_num, overlap_size),
            })
            .collect(),
    )
}

/// The length of a query's part in a chunk (0 when the chunk does not have the query or the
/// chunk number is negative).
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
fn part_length(parts: &[Vec<ChunkQuery>], chunk_num: isize, query: usize) -> usize {
    if chunk_num < 0 {
        return 0;
    }
    parts[chunk_num as usize]
        .iter()
        .find(|part| part.query == query)
        .map_or(0, |part| part.to - part.from)
}

/// The offsets of the contexts of a chunk (`size_t` arithmetic, stored as `Int4`): where the
/// chunk's part of a query starts in the query. A protein context has frame 0, so NCBI takes
/// its plus-strand branch.
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
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:630-643
/// ```c
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
/// ```
fn context_offsets(parts: &[Vec<ChunkQuery>], chunk_num: usize, overlap_size: usize) -> Vec<i32> {
    parts[chunk_num]
        .iter()
        .map(|part| {
            let query = part.query;
            let mut starting = chunk_num;
            for previous in (0..chunk_num).rev() {
                if !parts[previous].iter().any(|other| other.query == query) {
                    break;
                }
                starting = previous;
            }
            let mut correction = 0usize;
            let mut c = chunk_num as isize;
            while c != starting as isize {
                let prev_len = part_length(parts, c - 1, query);
                let curr_len = part_length(parts, c, query);
                let overlap = overlap_size.min(curr_len);
                correction = correction.wrapping_add(prev_len - overlap.min(prev_len));
                c -= 1;
            }
            correction as i32
        })
        .collect()
}

/// Whether NCBI's setup of some query chunk would split that chunk again (an overlap close
/// to the chunk size): each chunk's search builds its own query splitter, and only a debug
/// build asserts that it does not split; the release build then skips the chunk's lookup
/// table and stops with a null reference.
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
pub fn protein_chunk_would_be_split(chunks: &[ProteinQueryChunk], sizes: SplitSizes) -> bool {
    chunks.iter().any(|chunk| {
        let length = chunk.queries.iter().map(|part| part.to - part.from).sum();
        calculate_num_chunks(sizes, length).0 > 1
    })
}

/// The residues of a query part as the chunk's query reads them. With lower-case masking,
/// NCBI restricts the query's lower-case masks (`m_UserSpecifiedMasks`) to the part
/// (`blastn/query_split.rs` `restrict_masks`, which extends a mask by one residue), so the
/// part's letters are upper case except inside the restricted masks; without it the case of
/// the letters does not mask.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:278-281
/// ```c
///             CRef<CSeq_loc> mask_query_loc(new CSeq_loc);
///             s_SetSplitQuerySeqInterval(chunk, query_range, 0, mask_query_loc);
///             TMaskedQueryRegions split_mask =
///                 m_UserSpecifiedMasks[qindex].RestrictToSeqInt(mask_query_loc->GetInt());
/// ```
pub fn protein_chunk_part(query: &[u8], part: &ChunkQuery, mask_lowercase: bool) -> Vec<u8> {
    if !mask_lowercase {
        return query[part.from..part.to].to_vec();
    }
    let mut masks = Vec::new();
    let mut start = None;
    for (index, residue) in query.iter().enumerate() {
        match (residue.is_ascii_lowercase(), start) {
            (true, None) => start = Some(index),
            (false, Some(left)) => {
                masks.push(MaskedInterval::new(left, index));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(left) = start {
        masks.push(MaskedInterval::new(left, query.len()));
    }
    let mut residues: Vec<u8> = query[part.from..part.to]
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect();
    for mask in restrict_masks(&masks, part) {
        for residue in &mut residues[mask.start..mask.end] {
            *residue = residue.to_ascii_lowercase();
        }
    }
    residues
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
    // ```c
    //     if ((*chunk_size) > overlap_size) {
    //        num_chunks = concatenated_query_length / ((*chunk_size) - overlap_size);
    //     }
    // ```
    // A 60000-residue tblastn query is 3 chunks of 20067 residues that overlap by 100.
    #[test]
    fn a_long_tblastn_query_is_three_overlapping_chunks() {
        let sizes = protein_split_sizes(20000, None, None);
        let chunks = split_protein_batch(&[60_000], sizes).unwrap();
        let ranges: Vec<_> = chunks
            .iter()
            .map(|chunk| (chunk.queries[0].from, chunk.queries[0].to))
            .collect();
        assert_eq!(ranges, vec![(0, 20067), (19967, 40034), (39934, 60000)]);
        let offsets: Vec<_> = chunks
            .iter()
            .map(|chunk| chunk.context_offsets[0])
            .collect();
        assert_eq!(offsets, vec![0, 19967, 39934]);
        assert!(split_protein_batch(&[39_799], sizes).is_none());
        assert!(!protein_chunk_would_be_split(&chunks, sizes));
    }

    // Two queries share the first chunk (10051 residues, `calculate_num_chunks`); the
    // second query's later part starts where its first part's overlap begins.
    #[test]
    fn a_query_over_two_chunks_has_the_offset_of_its_first_part() {
        let sizes = protein_split_sizes(10000, None, None);
        let chunks = split_protein_batch(&[5_000, 15_000], sizes).unwrap();
        assert_eq!(chunks.len(), 2);
        assert_eq!(
            chunks[0].queries,
            vec![
                ChunkQuery {
                    query: 0,
                    from: 0,
                    to: 5_000
                },
                ChunkQuery {
                    query: 1,
                    from: 0,
                    to: 5_051
                },
            ]
        );
        assert_eq!(chunks[0].context_offsets, vec![0, 0]);
        assert_eq!(
            chunks[1].queries,
            vec![ChunkQuery {
                query: 1,
                from: 4_951,
                to: 15_000
            }]
        );
        assert_eq!(chunks[1].context_offsets, vec![4_951]);
    }

    // The lower-case mask [3, 6) restricted to the part [4, 8) is [4, 7) of the query (one
    // residue longer); the mask [10, 12) is outside the part.
    #[test]
    fn a_part_keeps_ncbis_restricted_lower_case_masks() {
        let query = b"ACDefgHIKLmn";
        let part = ChunkQuery {
            query: 0,
            from: 4,
            to: 8,
        };
        assert_eq!(protein_chunk_part(query, &part, true), b"fghI".to_vec());
        assert_eq!(protein_chunk_part(query, &part, false), b"fgHI".to_vec());
    }
}
