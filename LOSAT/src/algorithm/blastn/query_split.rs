//! NCBI's query splitting for BLASTN (`CQuerySplitter`): NCBI searches a query batch of at
//! least two chunks chunk by chunk in its preliminary stage, maps the HSPs of each chunk onto
//! the contexts of the batch and merges them (`run.rs`), and then runs the traceback on the
//! whole batch.

use super::args::QueryStrand;
use crate::utils::dust::MaskedInterval;

/// The query chunk size of the task when the environment variable `CHUNK_SIZE` is not set
/// or blank (`SplitQuery_GetChunkSize`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:65-73
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:51-53
pub const QUERY_CHUNK_OVERLAP: usize = 100;

/// NCBI's query chunk size and chunk overlap (`size_t` values), from the environment
/// variables `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` when they are set and not blank.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:53-62
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:53-60
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
    #[cfg(test)]
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
    /// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261
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
/// (`None` for a strand that is not searched: `kInvalidContext`) and the offset that NCBI adds
/// to the query coordinates of the chunk's HSPs in it (`i32::MAX` for `None`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryChunk {
    pub queries: Vec<ChunkQuery>,
    pub contexts: Vec<Option<usize>>,
    pub context_offsets: Vec<i32>,
}

/// The number of query chunks of a batch of `concatenated_query_length` residues, and the
/// chunk size NCBI then uses.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
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

/// Whether NCBI splits a query batch whose queries have `concatenated_query_length` letters into
/// chunks (`calculate_num_chunks` gives more than one; `split_query_batch` makes them).
pub fn splits(sizes: SplitSizes, concatenated_query_length: usize) -> bool {
    calculate_num_chunks(sizes, concatenated_query_length).0 > 1
}

/// The query chunks of a batch with queries of `lengths` residues, or `None` when NCBI
/// does not split it.
///
/// The chunk ranges on the concatenated queries:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
/// A chunk the loop leaves unset keeps the empty range of `SSplitQueryBlkNew`, and has no
/// query (NCBI then fails to make its query factory, "Empty CBlastQueryVector"); the chunk
/// size of `calculate_num_chunks` leaves none.
///
/// The queries of a chunk, and their parts in it:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:223-247
/// (`TChunkRange` is a `COpenRange`, so `GetTo() + 1` is the end of the previous query.)
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:196-210
/// The query location of a FASTA query starts at 0, so `query_offset` is 0.
///
/// The contexts of a chunk (the strand that is not searched is `None`):
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:392-410
/// The offsets are `context_offsets`.
pub fn split_query_batch(
    lengths: &[usize],
    sizes: SplitSizes,
    strand: QueryStrand,
) -> Option<Vec<QueryChunk>> {
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
        if chunk_end >= total_length || chunk_num + 1 == num_chunks {
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
                let (plus, minus) = strand.searched();
                contexts.push(plus.then_some(2 * query));
                contexts.push(minus.then_some(2 * query + 1));
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:251-264
fn context_in_chunk(chunk: &QueryChunk, absolute_context: usize) -> Option<usize> {
    chunk
        .contexts
        .iter()
        .position(|&context| context == Some(absolute_context))
}

/// The length of a batch context in a chunk (0 when the chunk does not have it).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:452-466
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:266-284
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:612-662
fn context_offsets(
    chunks: &[QueryChunk],
    chunk_num: usize,
    lengths: &[usize],
    overlap_size: usize,
) -> Vec<i32> {
    chunks[chunk_num]
        .contexts
        .iter()
        .map(|&context| {
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:604-627
            let Some(absolute_context) = context else {
                return i32::MAX;
            };
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqlocinfo.cpp:102-116
/// (`TSeqRange` is a `CRange`, whose constructor takes the last position.) The chunk's query
/// block then puts the masks in the part's coordinates, from the chunk's Seq-interval in
/// the query's record:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1047
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:256-262
/// With `-query_loc` the query's masks are in record coordinates (the searched letters'
/// masks moved by `offset`, the start of the query's interval) and so is the chunk's
/// Seq-interval, but the restriction above takes the part in the searched letters
/// (`query_offset` 0): a mask survives only in `[from + offset, to - 1]`. NCBI does this
/// for every split query with a range; LOSAT reproduces it. `offset` is 0 without a range.
pub fn restrict_masks(
    masks: &[MaskedInterval],
    part: &ChunkQuery,
    offset: usize,
) -> Vec<MaskedInterval> {
    // `loc` is the part [from, to - 1) in the searched letters' coordinates.
    let last = part.to - 1;
    // The chunk's Seq-interval: [from + offset, to - 1 + offset] in the record.
    let (interval_from, interval_to) = (part.from + offset, last + offset);
    masks
        .iter()
        .filter_map(|mask| {
            // `TSeqRange mask` is [start, end) in the record (the mask's last position is
            // end - 1).
            let left = (mask.start + offset).max(part.from);
            let right_open = (mask.end + offset).min(last);
            // The new Seq-interval is [left, right_open].
            if left >= right_open {
                return None;
            }
            // `BlastSeqLoc_RestrictToInterval(mask, interval_from, interval_to)`.
            let restricted_left = left.saturating_sub(interval_from);
            let restricted_right = right_open.min(interval_to) as isize - interval_from as isize;
            (restricted_left as isize <= restricted_right)
                .then(|| MaskedInterval::new(restricted_left, restricted_right as usize + 1))
        })
        .collect()
}

/// `restrict_masks` for masks given in the record's coordinates (`[start, end)` in the record),
/// as NCBI holds a query's lower-case masks: the runs of the whole record, also those before a
/// `-query_loc` interval. A run that ends on the letter before the interval then reaches the
/// first letter of the query's first part through the restriction's extra residue, which a mask
/// in the searched letters' coordinates cannot express. The steps are `restrict_masks`'s.
pub fn restrict_record_masks(
    masks: &[MaskedInterval],
    part: &ChunkQuery,
    offset: usize,
) -> Vec<MaskedInterval> {
    let last = part.to - 1;
    let (interval_from, interval_to) = (part.from + offset, last + offset);
    masks
        .iter()
        .filter_map(|mask| {
            let left = mask.start.max(part.from);
            let right_open = mask.end.min(last);
            if left >= right_open {
                return None;
            }
            let restricted_left = left.saturating_sub(interval_from);
            let restricted_right = right_open.min(interval_to) as isize - interval_from as isize;
            (restricted_left as isize <= restricted_right)
                .then(|| MaskedInterval::new(restricted_left, restricted_right as usize + 1))
        })
        .collect()
}

/// `restrict_masks` of a query whose letter before its `-query_loc` interval (record letter
/// `offset - 1`) is lower case, with `-lcase_masking`. NCBI restricts the lower-case runs of the
/// whole record (`restrict_record_masks`); of the letters outside the interval only that one can
/// reach a part, the first letter of the query's first part, whether its run ends there or goes
/// on into the interval. The masks in the searched letters' coordinates are moved into the
/// record's with that letter added. It masks the first letter of a part that starts the query and
/// is longer than `offset` letters; `restrict_masks` misses that letter when the run ends on the
/// letter before the interval, or goes on into the interval and the part is `offset + 1` letters
/// long, and masks the same letters in every other part.
pub fn restrict_masks_with_letter_before(
    masks: &[MaskedInterval],
    part: &ChunkQuery,
    offset: usize,
) -> Vec<MaskedInterval> {
    let mut record = Vec::with_capacity(masks.len() + 1);
    record.push(MaskedInterval::new(offset - 1, offset));
    record.extend(
        masks
            .iter()
            .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset)),
    );
    restrict_record_masks(&record, part, offset)
}

/// With `-filtering_db`, where the masks that NCBI's chunk keeps before a part's letters end:
/// the first restriction of `restrict_masks` keeps masks of the record that end before the
/// part's Seq-interval when the part's query has a `-query_loc` interval, and the chunk combines
/// them with the part's repeat hits (`repeat_filter::RepeatFilter::add_masks`) before its search
/// takes the part's letters. `masks` are the query's masks in the searched letters'
/// coordinates, `before` the masked letters before its interval (`repeat_filter::MasksOutside`,
/// bit `k`: record letter `offset - 1 - k`). The letters from the last of them (with the
/// restriction's extra residue) to the part's first letter; `None` without one.
pub fn masked_before_part(
    masks: &[MaskedInterval],
    before: u8,
    part: &ChunkQuery,
    offset: usize,
) -> Option<usize> {
    let last = part.to - 1;
    let interval_from = part.from + offset;
    let letters_before = (0..u8::BITS as usize)
        .filter(|&k| before & (1 << k) != 0)
        .filter_map(|k| offset.checked_sub(k + 1))
        .map(|letter| (letter, letter + 1));
    masks
        .iter()
        .map(|mask| (mask.start + offset, mask.end + offset))
        .chain(letters_before)
        .filter_map(|(start, end)| {
            let (left, right_open) = (start.max(part.from), end.min(last));
            (left < right_open && right_open < interval_from).then_some(right_open)
        })
        .max()
        .map(|end| interval_from - end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masks_before_a_part_of_a_ranged_query() {
        let part = |from, to| ChunkQuery { query: 0, from, to };
        // The query's interval starts at record letter 100: a mask of the searched letters
        // 296..298 ends (with the extra residue) 2 letters before the part 300..600's interval.
        let masks = [MaskedInterval::new(296, 298)];
        assert_eq!(masked_before_part(&masks, 0, &part(300, 600), 100), Some(2));
        // The first part keeps the lower case before the interval: record letter 97 (bit 2)
        // ends 2 letters before the interval with its extra residue; letter 99 (bit 0) reaches
        // the part's first letter, which `restrict_masks_with_letter_before` masks.
        assert_eq!(masked_before_part(&[], 0b100, &part(0, 600), 100), Some(2));
        assert_eq!(masked_before_part(&[], 0b1, &part(0, 600), 100), None);
        // Without a range nothing lies before a part's interval.
        assert_eq!(masked_before_part(&masks, 0, &part(300, 600), 0), None);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:612-662
    #[test]
    fn a_query_over_two_chunks_is_mapped_onto_its_contexts() {
        let chunks = split_query_batch(
            &[1_999_800],
            SplitSizes::default_for(false),
            QueryStrand::Both,
        )
        .expect("split");
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
        assert_eq!(chunks[0].contexts, vec![Some(0), Some(1)]);
        assert_eq!(chunks[0].context_offsets, vec![0, 1_999_800 - 999_951]);
        assert_eq!(chunks[1].context_offsets, vec![999_851, 0]);
    }

    #[test]
    fn parts_that_end_or_start_in_an_overlap_keep_their_starts() {
        // Chunks of 1,000,059 overlapping by 100: the first query ends 58 residues into the
        // overlap of the second and third chunks, and the second query starts there.
        let chunks = split_query_batch(
            &[1_999_976, 1_000_000],
            SplitSizes::default_for(false),
            QueryStrand::Both,
        )
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
        assert_eq!(chunks[1].contexts, vec![Some(0), Some(1), Some(2), Some(3)]);
        assert_eq!(
            chunks[1].context_offsets,
            vec![999_959, 0, 0, 1_000_000 - 42]
        );
        assert_eq!(chunks[2].context_offsets, vec![1_999_918, 0, 0, 0]);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:256-281
    // With `-query_loc 1000001-3200000` a mask survives in a chunk only from the chunk's
    // start plus the interval's start on (oracle: `-task blastn -lcase_masking` loses the
    // lower-case runs at record 1.2 Mb and 2.3 Mb, keeps those at 1.05 Mb and 2.12 Mb).
    #[test]
    fn ranged_chunk_masks_are_restricted_in_the_searched_coordinates() {
        let offset = 1_000_000;
        let masks = [
            MaskedInterval::new(50_000, 55_000),
            MaskedInterval::new(200_000, 205_000),
            MaskedInterval::new(1_120_000, 1_125_000),
            MaskedInterval::new(1_300_000, 1_305_000),
        ];
        let first = ChunkQuery {
            query: 0,
            from: 0,
            to: 1_100_050,
        };
        assert_eq!(
            restrict_masks(&masks, &first, offset),
            vec![MaskedInterval::new(50_000, 55_001)]
        );
        let second = ChunkQuery {
            query: 0,
            from: 1_099_950,
            to: 2_200_000,
        };
        assert_eq!(
            restrict_masks(&masks, &second, offset),
            vec![MaskedInterval::new(20_050, 25_051)]
        );
        // Without a range the same masks are restricted to each part as they are.
        assert_eq!(
            restrict_masks(&masks, &second, 0),
            vec![
                MaskedInterval::new(20_050, 25_051),
                MaskedInterval::new(200_050, 205_051)
            ]
        );
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seqlocinfo.cpp:102-116
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
            restrict_masks(&masks, &part, 0),
            vec![
                MaskedInterval::new(0, 11),
                MaskedInterval::new(50, 61),
                MaskedInterval::new(90, 100),
            ]
        );
    }

    // Masks in the record's coordinates are restricted as `restrict_masks` restricts the same
    // masks in the searched letters' coordinates; a run that ends on the letter before the
    // interval (record 9 with offset 10) reaches the first part's first letter, only while that
    // part's last letter (`to - 1`, which bounds the run in record coordinates) is not before the
    // offset.
    #[test]
    fn record_masks_restrict_as_searched_masks_and_reach_the_first_letter() {
        let offset = 10;
        let searched = [
            MaskedInterval::new(0, 3),
            MaskedInterval::new(95, 120),
            MaskedInterval::new(185, 190),
        ];
        let record: Vec<_> = searched
            .iter()
            .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset))
            .collect();
        for (from, to) in [(0, 100), (90, 200), (180, 200), (0, 11), (0, 12)] {
            let part = ChunkQuery { query: 0, from, to };
            assert_eq!(
                restrict_record_masks(&record, &part, offset),
                restrict_masks(&searched, &part, offset),
                "{from}..{to}"
            );
        }
        let before = [MaskedInterval::new(4, 10)];
        let part = |to| ChunkQuery {
            query: 0,
            from: 0,
            to,
        };
        assert_eq!(
            restrict_record_masks(&before, &part(100), offset),
            vec![MaskedInterval::new(0, 1)]
        );
        assert_eq!(
            restrict_record_masks(&before, &part(11), offset),
            vec![MaskedInterval::new(0, 1)]
        );
        assert!(restrict_record_masks(&before, &part(10), offset).is_empty());
        let later = ChunkQuery {
            query: 0,
            from: 90,
            to: 200,
        };
        assert!(restrict_record_masks(&before, &later, offset).is_empty());
        let two_before = [MaskedInterval::new(4, 9)];
        assert!(restrict_record_masks(&two_before, &part(100), offset).is_empty());
    }

    // The masks of a part from the runs of the searched letters, with the letter before the
    // interval when it is lower case, mask the residues that the runs of the whole record mask
    // (`restrict_record_masks`, NCBI's), for every lower-case pattern of a 10-letter record, four
    // intervals and every part; the searched letters' runs alone miss some.
    #[test]
    fn the_letter_before_the_interval_gives_the_whole_records_masks() {
        use super::super::coordination::collect_lowercase_masks;
        let residues = |masks: &[MaskedInterval]| -> std::collections::BTreeSet<usize> {
            masks.iter().flat_map(|mask| mask.start..mask.end).collect()
        };
        let mut missed = 0;
        for pattern in 0u32..1 << 10 {
            let record: Vec<u8> = (0..10)
                .map(|letter| {
                    if pattern >> letter & 1 == 1 {
                        b'a'
                    } else {
                        b'A'
                    }
                })
                .collect();
            let runs = collect_lowercase_masks(&record);
            for (offset, end) in [(1, 10), (2, 8), (3, 10), (4, 9)] {
                let searched = collect_lowercase_masks(&record[offset..end]);
                let before = record[offset - 1].is_ascii_lowercase();
                for from in 0..end - offset {
                    for to in from + 1..=end - offset {
                        let part = ChunkQuery { query: 0, from, to };
                        let whole = residues(&restrict_record_masks(&runs, &part, offset));
                        let masks = if before {
                            restrict_masks_with_letter_before(&searched, &part, offset)
                        } else {
                            restrict_masks(&searched, &part, offset)
                        };
                        assert_eq!(
                            residues(&masks),
                            whole,
                            "{} {offset}..{end} {from}..{to}",
                            String::from_utf8_lossy(&record)
                        );
                        missed += usize::from(
                            residues(&restrict_masks(&searched, &part, offset)) != whole,
                        );
                    }
                }
            }
        }
        assert!(missed > 0);
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
        let chunks =
            split_query_batch(&[1000, 1500, 700], sizes, QueryStrand::Both).expect("split");
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:392-410,604-626
    // A strand that is not searched has no chunk context and the offset INT4_MAX; the other
    // strand's offsets are those of the split with both strands.
    #[test]
    fn a_single_strand_split_keeps_the_offsets_of_the_strand_it_searches() {
        let lengths = [1_999_976, 1_000_000];
        let sizes = SplitSizes::default_for(false);
        let both = split_query_batch(&lengths, sizes, QueryStrand::Both).expect("split");
        for (strand, keep) in [(QueryStrand::Plus, 0usize), (QueryStrand::Minus, 1)] {
            let one = split_query_batch(&lengths, sizes, strand).expect("split");
            assert_eq!(one.len(), both.len());
            for (one, both) in one.iter().zip(&both) {
                assert_eq!(one.queries, both.queries);
                assert_eq!(one.contexts.len(), both.contexts.len());
                for (position, (&context, &offset)) in
                    one.contexts.iter().zip(&one.context_offsets).enumerate()
                {
                    if position % 2 == keep {
                        assert_eq!(context, both.contexts[position]);
                        assert_eq!(offset, both.context_offsets[position]);
                    } else {
                        assert_eq!(context, None);
                        assert_eq!(offset, i32::MAX);
                    }
                }
            }
        }
        let plus = split_query_batch(&lengths, sizes, QueryStrand::Plus).expect("split");
        assert_eq!(plus[1].contexts, vec![Some(0), None, Some(2), None]);
        let minus = split_query_batch(&lengths, sizes, QueryStrand::Minus).expect("split");
        assert_eq!(minus[1].contexts, vec![None, Some(1), None, Some(3)]);
        assert_eq!(
            minus[1].context_offsets,
            vec![i32::MAX, 0, i32::MAX, 1_000_000 - 42]
        );
    }
}
