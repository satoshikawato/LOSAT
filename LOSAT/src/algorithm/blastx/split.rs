//! NCBI single-query BLASTX chunk metadata and score-only overlap rejoin.
use super::{
    args::ResolvedOptions,
    input::FastaRecord,
    preliminary::PreliminaryHsp,
    query_setup::{prepare_queries, PreparedQueryBatch},
};
use crate::algorithm::blastn::query_split::{
    restrict_masks, restrict_record_masks, ChunkQuery, SplitSizes,
};
use crate::utils::dust::MaskedInterval;
use anyhow::{ensure, Result};
use std::ops::Range;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:78-86
pub const CHUNK_SIZE: usize = 10002;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:62-69
pub const OVERLAP_NT: usize = 297;

/// The sizes when neither `CHUNK_SIZE` nor `OVERLAP_CHUNK_SIZE` is set.
#[cfg(test)]
pub const DEFAULT_SIZES: SplitSizes = SplitSizes {
    chunk_size: CHUNK_SIZE as u64,
    overlap: OVERLAP_NT as u64,
};

/// LOSAT's refusal of a query split (`split_sizes`, `split_queries`), not an error of NCBI's
/// engine: it gets no `BLAST engine error:` (`native.rs` `search_error`).
#[derive(Debug)]
pub struct QuerySplitRefusal(pub String);

impl std::fmt::Display for QuerySplitRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for QuerySplitRefusal {}

fn refusal(message: impl Into<String>) -> anyhow::Error {
    QuerySplitRefusal(message.into()).into()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:54-62
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:50-60
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:105-111
/// The chunk size and overlap, in nucleotides, of a query batch of `query_count` queries, as
/// NCBI's query splitter reads them for each batch before the batch's queries are set up:
/// `CHUNK_SIZE` always, `OVERLAP_CHUNK_SIZE` only when the batch has one query
/// (`SplitQuery_ShouldSplit`), also for an ungapped search, which is not split. The `int`
/// becomes NCBI's 64-bit `size_t`, so -1 passes the check of 3 and -2 does not. A chunk size
/// not divisible by 3 is NCBI's engine error; a value NCBI cannot convert is its
/// `CStringException`, which ends the program (`app::read_query_split_variable`; the search's
/// error frames pass it on as it is, `native::database_search_error`).
pub fn split_sizes(query_count: usize) -> Result<SplitSizes> {
    let read = crate::blastinput::app::read_query_split_variable;
    let chunk_size = read("CHUNK_SIZE")?.map_or(CHUNK_SIZE as u64, |value| i64::from(value) as u64);
    validate_chunk_size(chunk_size)?;
    let overlap = match query_count > 1 {
        true => OVERLAP_NT as u64,
        false => {
            read("OVERLAP_CHUNK_SIZE")?.map_or(OVERLAP_NT as u64, |value| i64::from(value) as u64)
        }
    };
    Ok(SplitSizes {
        chunk_size,
        overlap,
    })
}

/// The queries of the query batch whose records (or reader-batch unit) are `records`: with a
/// reader batch of several units it is the reader batch's (`reader_batch.rs`).
fn batch_query_count(records: &[FastaRecord], full: &PreparedQueryBatch) -> usize {
    full.reader_shape
        .as_deref()
        .map_or(records.len(), |shape| shape.query_count)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:105-145
pub fn calculate_num_chunks(
    chunk_size: &mut u64,
    overlap: u64,
    length: usize,
    num_queries: usize,
) -> usize {
    if num_queries > 1 {
        return 1;
    }
    *chunk_size -= *chunk_size % 3;
    // `num_chunks` is a `Uint4` (split_query_aux_priv.cpp:112).
    let n = if *chunk_size > overlap {
        (length as u64 / (*chunk_size - overlap)) as u32
    } else {
        0
    };
    if n <= 1 {
        *chunk_size = length as u64;
        1
    } else {
        n as usize
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:101-109
pub fn validate_chunk_size(size: u64) -> Result<()> {
    ensure!(
        size.is_multiple_of(3),
        "Split query chunk size must be divisible by 3"
    );
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:161-178
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:365-389
#[derive(Clone, Debug)]
pub struct QueryChunk {
    /// Read by the tests that check the split against NCBI's chunk bounds.
    #[cfg(test)]
    pub range: Range<usize>,
    pub absolute_contexts: Vec<i32>,
    pub corrections: Vec<i32>,
    pub prepared: PreparedQueryBatch,
    /// The overlap of the batch's chunks in residues (`kOverlap`, `SetChunkOverlapSize`; 0 when
    /// the batch is not split), which the merge of two chunks' HSPs reads.
    pub overlap: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:50-62
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:300-335
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:724-836
pub fn split_queries(
    records: &[FastaRecord],
    full: &PreparedQueryBatch,
    options: &ResolvedOptions,
    sizes: SplitSizes,
) -> Result<Vec<QueryChunk>> {
    let mut size = sizes.chunk_size;
    let length = records.iter().map(|r| r.sequence.len()).sum();
    // A batch of several queries is not split.
    let num_queries = batch_query_count(records, full);
    let mut count = calculate_num_chunks(&mut size, sizes.overlap, length, num_queries);
    if !options.gapped {
        count = 1;
    }
    if count == 1 {
        return Ok(vec![QueryChunk {
            #[cfg(test)]
            range: 0..length,
            absolute_contexts: (0..full.contexts.len() as i32).collect(),
            corrections: vec![0; full.contexts.len()],
            prepared: full.clone(),
            overlap: 0,
        }]);
    }
    // A negative CHUNK_SIZE splits a batch only with a negative OVERLAP_CHUNK_SIZE below it;
    // NCBI's `size_t` chunk ranges then wrap (it stops with a CCoreException).
    if (size as i64) < 0 {
        return Err(refusal(
            crate::algorithm::blastn::query_split::negative_split_refusal("BLASTX"),
        ));
    }
    // NCBI splits the nucleotide query and translates each chunk, so an overlap that is not a
    // multiple of 3 starts the chunks after the first out of the query's reading frames
    // (split_query_aux_priv.cpp:62-64: "this value must be divisible by 3 to work with
    // translated queries"), and NCBI then moves their HSPs onto contexts of other frames.
    if !sizes.overlap.is_multiple_of(3) {
        return Err(refusal(format!(
            "OVERLAP_CHUNK_SIZE={}, which is not a multiple of 3, splits this query in chunks out of its reading frames (NCBI BLAST+ requires a multiple of 3 for a translated query); this is not supported by LOSAT's BLASTX",
            sizes.overlap
        )));
    }
    let step = (size - sizes.overlap) as usize;
    let size = size as usize;
    let overlap = (sizes.overlap / 3) as i32;
    let range = |index: usize| {
        let start = index * step;
        let end = if index + 1 == count {
            length
        } else {
            (start + size).min(length)
        };
        start..end
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-208
    // Each chunk's setup builds its own query splitter with the same sizes; when it would split
    // the chunk again (an overlap of half the chunk size or more), the release build skips the
    // chunk's lookup table and stops with a null reference. This also covers a loop that ends
    // before `count` chunks (a chunk then holds twice the step).
    if (0..count).any(|index| {
        let mut chunk_size = size as u64;
        calculate_num_chunks(&mut chunk_size, sizes.overlap, range(index).len(), 1) > 1
    }) {
        return Err(refusal(
            "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would split a query chunk again (it then stops with a null reference) is not supported by LOSAT's BLASTX",
        ));
    }
    // A chunk's query setup (translation, SEG) reads only its own slice of the query, so the chunks
    // can be prepared on the search pool; the corrections below still run in chunk order.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-249
    // NCBI prepares each query chunk in turn (query factory for the chunk, then its local query
    // data, which holds the translation and the filtering). make_chunk is the same preparation for
    // one chunk, as a closure so that several chunks can be prepared at once.
    let make_chunk = |index: usize| -> Result<QueryChunk> {
        let Range { start, end } = range(index);
        let record = if x_chunk_slice() {
            x_chunk_record(&records[0], start, end)
        } else {
            x_chunk_record_from_copy(&records[0], start, end)
        };
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
        Ok(QueryChunk {
            #[cfg(test)]
            range: start..end,
            absolute_contexts,
            corrections: vec![i32::MAX; 6],
            prepared,
            overlap,
        })
    };
    let mut chunks: Vec<QueryChunk> = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-249
    // Dispatch point: inside the search pool of 2+ threads, this thread prepares the chunks in
    // index order and the pool's idle threads prepare the next ones (`utils::xitems::x_map`: this
    // thread runs no other job of the pool meanwhile); otherwise the loop below prepares them one
    // after another as NCBI does. The chunks are collected in index order and the corrections that
    // follow still run in chunk order.
    #[cfg(feature = "parallel")]
    if super::runtime::x_bx_par_rule(|| {
        !super::runtime::x_inner_serial()
            && rayon::current_thread_index().is_some()
            && rayon::current_num_threads() > 1
    }) {
        chunks = crate::utils::xitems::x_map(count, &make_chunk)
            .into_iter()
            .collect::<Result<Vec<_>>>()?;
    }
    if chunks.is_empty() {
        for index in 0..count {
            chunks.push(make_chunk(index)?);
        }
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
                    sum += prev - overlap.min(curr).min(prev);
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
                    sub += curr - overlap.min(curr).min(prev);
                }
                full.contexts[absolute as usize].length as i32 - sub
            };
            chunks[index].corrections[c] = correction;
        }
    }
    Ok(chunks)
}

/// EXPERIMENT (LOSAT_X_BXCHUNKSLICE): whether a chunk's record is made without a copy of the whole
/// query. A copy for each chunk (569 of 5.3 MiB for a 5.5 Mb query, 3 GiB in all) was the largest
/// block each thread that prepared chunks ever held in its glibc arena, which kept that much
/// resident through the search after the block was freed.
fn x_chunk_slice() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXCHUNKSLICE"))
}

/// The lowercase masks of `whole` inside `start..end`, moved to the chunk's coordinates.
fn x_chunk_masks(whole: &FastaRecord, start: usize, end: usize) -> Vec<(i32, i32)> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281
    // Each run is one letter longer inside a chunk, and a run that starts on the chunk's last
    // letter is dropped (`restrict_masks`). With `-query_loc` NCBI restricts the record's runs by
    // the chunk in the searched letters' coordinates, so a run survives only in
    // `[start + offset, end - 1]` of the record (E2d AUTHORITY §E.1, as BLASTN's and TBLASTN's
    // split queries): the preliminary search of a chunk sees the runs of its last `offset` letters
    // unmasked, which the traceback of the whole query masks again.
    let part = ChunkQuery {
        query: 0,
        from: start,
        to: end,
    };
    let offset = whole.placement().offset;
    let masks = match whole
        .source
        .as_deref()
        .and_then(|source| source.uncut.as_deref())
    {
        // NCBI's masks are the runs of the whole record (`m_UserSpecifiedMasks`). The runs cut to
        // the interval (`lowercase_masks`) restrict to the same masks but for one letter: the run
        // that holds the letter before the interval masks the first chunk's first letter, when it
        // ends there and the chunk is longer than `offset`, or when it goes on into the interval
        // and the chunk is `offset + 1` letters long. A record cut to its interval keeps its
        // whole letters when that letter is lower case (`FastaRecord::lower_case_before`), or for
        // the archive. The restriction keeps a run only where it reaches `[start + offset,
        // end - 1]` of the record, so the runs are read from the chunk's letters in the record and
        // the one before them: a run cut at either end of these is restricted to the same letters.
        Some(uncut) => {
            let from = (start + offset).saturating_sub(1);
            let to = (end + offset).min(uncut.len());
            let runs: Vec<MaskedInterval> =
                crate::algorithm::common::masks::lowercase_runs(&uncut[from..to])
                    .into_iter()
                    .map(|(a, b)| MaskedInterval::new(from + a as usize, from + b as usize + 1))
                    .collect();
            restrict_record_masks(&runs, &part, offset)
        }
        None => {
            let runs: Vec<MaskedInterval> = whole
                .lowercase_masks
                .iter()
                .map(|&(a, b)| MaskedInterval::new(a as usize, b as usize + 1))
                .collect();
            restrict_masks(&runs, &part, offset)
        }
    };
    masks
        .into_iter()
        .map(|mask| (mask.start as i32, mask.end as i32 - 1))
        .collect()
}

// No NCBI counterpart: NCBI's chunk factory reads the chunk's range of the query it holds; here the
// record is what `prepare_queries` takes.
/// EXPERIMENT (LOSAT_X_BXCHUNKSLICE): the record of the chunk `start..end` of `whole`, with the
/// chunk's letters and masks and without a copy of the others.
fn x_chunk_record(whole: &FastaRecord, start: usize, end: usize) -> FastaRecord {
    FastaRecord {
        internal_id: whole.internal_id.clone(),
        title: whole.title.clone(),
        sequence: whole.sequence[start..end].to_vec(),
        lowercase_masks: x_chunk_masks(whole, start, end),
        warnings: whole.warnings.clone(),
        // Not the whole letters of a record cut to `-query_loc` for outfmt 11: one copy per chunk.
        source: whole
            .source
            .as_deref()
            .map(crate::blastinput::fasta_reader::RecordSource::without_uncut),
        around: None,
        skipped: None,
    }
}

/// `x_chunk_record` as `LOSAT_X_BXCHUNKSLICE=0` makes it: from a copy of the whole record.
fn x_chunk_record_from_copy(whole: &FastaRecord, start: usize, end: usize) -> FastaRecord {
    let mut record = whole.clone();
    record.sequence = record.sequence[start..end].to_vec();
    record.lowercase_masks = x_chunk_masks(whole, start, end);
    record.source = whole
        .source
        .as_deref()
        .map(crate::blastinput::fasta_reader::RecordSource::without_uncut);
    record
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:477-510
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1488-1531
// The divisor is never 0 here, unlike in TBLASTX's merge across subject chunks
// (`ChunkMergeUndefined`): only the one-hit word finder's ungapped HSPs can have a length of 0 or
// less, and an ungapped search is never split (`split_queries`, split_query_cxx.cpp:60-61), so
// the HSPs of a split query are gapped ones, each of positive length.
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

#[cfg(test)]
#[path = "split_x_chunk_record_tests.rs"]
mod x_chunk_record_tests;
