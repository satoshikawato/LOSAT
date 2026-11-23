//! NCBI's query splitting for BLASTP and TBLASTN (`CQuerySplitter`): NCBI searches a protein
//! query batch of at least two chunks chunk by chunk in its preliminary stage, maps the HSPs
//! of each chunk onto the queries of the batch, merges them, and then runs the traceback on
//! the whole batch. A protein query has one context, so the chunk of a query part is that
//! query's context (BLASTN's two strands are `blastn/query_split.rs`).

use crate::algorithm::blastn::query_split::{
    calculate_num_chunks, ChunkQuery, SplitSizes, QUERY_CHUNK_OVERLAP,
};
use crate::algorithm::common::masks::lowercase_runs;
use crate::blastinput::app::read_query_split_variable;

/// The query chunk size of blastp (10000) or tblastn (20000) when the environment variable
/// `CHUNK_SIZE` is not set or blank, and the overlap (100) when `OVERLAP_CHUNK_SIZE` is not.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:74-76,91-94
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:51-60
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
/// The queries of a chunk and their parts are those of `blastn/query_split.rs`
/// `split_query_batch` (split_query_cxx.cpp:196-247). A protein query has one context:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:415-418
/// The offsets are `context_offsets`. A chunk that the loop leaves unset keeps the empty
/// range of `SSplitQueryBlkNew` and has no query part.
pub fn split_protein_batch(lengths: &[usize], sizes: SplitSizes) -> Option<Vec<ProteinQueryChunk>> {
    let total_length: usize = lengths.iter().sum();
    let (num_chunks, chunk_size) = calculate_num_chunks(sizes, total_length);
    if num_chunks <= 1 {
        return None;
    }
    // NCBI's `size_t` overlap: the variable's value, near 2^64 for a negative OVERLAP_CHUNK_SIZE
    // (which splits a batch only below a negative CHUNK_SIZE).
    let overlap_size = sizes.overlap as usize;

    let mut chunk_ranges = vec![(0usize, 0usize); num_chunks];
    let mut chunk_start = 0usize;
    for (chunk_num, range) in chunk_ranges.iter_mut().enumerate() {
        let mut chunk_end = chunk_start + chunk_size;
        if chunk_end >= total_length || chunk_num + 1 == num_chunks {
            chunk_end = total_length;
        }
        *range = (chunk_start, chunk_end);
        // NCBI's `size_t` arithmetic (c++/src/algo/blast/api/split_query_cxx.cpp:166): an evened chunk size below a
        // positive overlap wraps the start past the queries' length at the first chunk (the loop
        // ends, the chunks after it keep no query); a negative overlap moves the start by the
        // chunk size plus the overlap's magnitude (gaps between the chunks), or past the length.
        chunk_start = chunk_start.wrapping_add(chunk_size.wrapping_sub(overlap_size));
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:452-466
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:266-284
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:630-643
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
pub fn protein_chunk_would_be_split(chunks: &[ProteinQueryChunk], sizes: SplitSizes) -> bool {
    chunks.iter().any(|chunk| {
        let length = chunk.queries.iter().map(|part| part.to - part.from).sum();
        calculate_num_chunks(sizes, length).0 > 1
    })
}

/// The residues of a query part as the chunk's query reads them with lower-case masking: NCBI
/// restricts the query's lower-case masks (`m_UserSpecifiedMasks`) to the part
/// (`blastn/query_split.rs` `restrict_masks`, which extends a mask by one residue), so the
/// part's letters are upper case except inside the restricted masks. The query's searched
/// letters start at `offset` in `record`, its whole record (`-query_loc`; `offset` 0 and the
/// query itself without a range). NCBI's masks are the lower-case runs of the whole record, so a
/// run that ends on the letter before the interval reaches the interval's first letter through
/// the restriction's extra residue (`restrict_masks`, whose steps this takes in record
/// coordinates).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:383-384
pub fn protein_chunk_part_in_record(record: &[u8], offset: usize, part: &ChunkQuery) -> Vec<u8> {
    let mut residues: Vec<u8> = record[offset + part.from..offset + part.to]
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect();
    let last = part.to - 1;
    let (interval_from, interval_to) = (part.from + offset, last + offset);
    for (start, end) in lowercase_runs(record) {
        // The run is [start, end] in the record; `loc` is [from, to - 1) in the searched letters.
        // NCBI compares the two as numbers (`RestrictToSeqInt`), so with an offset a run is cut by
        // the part's searched coordinates, not by where the part lies in the record.
        let left = (start as usize).max(part.from);
        let right_open = (end as usize + 1).min(last);
        if left >= right_open {
            continue;
        }
        let restricted_left = left.saturating_sub(interval_from);
        let restricted_right = right_open.min(interval_to) as isize - interval_from as isize;
        if restricted_left as isize <= restricted_right {
            for residue in &mut residues[restricted_left..=restricted_right as usize] {
                *residue = residue.to_ascii_lowercase();
            }
        }
    }
    residues
}

/// NCBI's query chunk size and overlap of a blastp or tblastn search
/// (`protein_split_sizes`), with the environment variables
/// `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` read as NCBI reads them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:55-58
/// A value that NCBI cannot convert is the `Err` (its exception), which stops the first query
/// batch in its query splitter before any batch is split (`query_batch::read_query_batches`).
pub fn protein_query_split_sizes(default_chunk_size: u32) -> anyhow::Result<SplitSizes> {
    Ok(protein_split_sizes(
        default_chunk_size,
        read_query_split_variable("CHUNK_SIZE")?,
        read_query_split_variable("OVERLAP_CHUNK_SIZE")?,
    ))
}

/// Where the preliminary search of a query chunk first fails to read a `-db` database
/// (`SplitBatchSteps`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadFails {
    /// In its setup, before its Karlin-Altschul check: the length of OID 0, which the setup reads
    /// where the statistics' total is not positive, of volumes without an OID (blast_setup.c:
    /// 925-937, cited at `SplitBatchSteps`). Any chunk that NCBI does not split again meets it,
    /// before the setup checks of its batch, as does a batch searched whole.
    Setup,
    /// In the search: the first chunk with a valid context that NCBI searches meets it, at the
    /// first read of a sequence (a sequence file that cannot be opened). TBLASTN also gives this
    /// kind to the subject-length error of its stores check (not a read), which the chunk's own
    /// setup checks meet first with the same text.
    Search,
}

/// Where a chunk of a split batch goes after `SplitBatchSteps::before_read`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChunkStep {
    /// The chunk's search reads the subjects.
    Read,
    /// The chunk finds nothing: none of its contexts is valid.
    Skip,
    /// The stage ends: NCBI's read of the database fails here (`ReadFails`).
    End,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:233-290
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171,234-249,876-880
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_search_runner.hpp:71-77
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-982
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1832-1836
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback_mt_priv.c:202-240
// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:233-252
/// What NCBI's search of a split query batch meets before the search of each chunk reads the
/// subjects, in NCBI's order (2.17.0 probes of line RG items X-short-bat, fix round 3, and
/// X-bpgap): one owner for BLASTP's and TBLASTN's searches of a split batch, for the set-up of a
/// batch without its search, and for where the failures of a call's setup come when a split batch
/// meets them first. Each stop is LOSAT's refusal in the program's form (`stop`: TBLASTN's
/// `blastinput::app::ErrorAfterBatchReports`, BLASTP's batch failure, psiblast's plain error); each
/// error of NCBI's setup is that program's error there.
/// - `new`: the query splitter makes the query factory of every chunk first and dereferences the
///   missing query vector of a chunk that no query reaches (a null pointer), before any chunk is
///   set up. This is where a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE (NCBI's
///   `size_t` chunk ranges wrap) whose first range covers the batch stops, and a positive pair
///   whose evened chunk size falls below the overlap (the start wraps after chunk 0).
/// - `before_read`, chunk by chunk after the chunk's setup: the setup of a chunk that would be split
///   again makes no lookup table, which the chunk's search dereferences whether or not its contexts
///   are valid, after the chunks before it were searched; then the setup of the chunk's preliminary
///   search reads the length of OID 0 where the statistics' total is not positive
///   (`ReadFails::Setup`), checks the Karlin-Altschul blocks (a chunk without a valid context ends
///   there: NCBI ignores its error), then the subject length and the hit-saving parameters (the
///   batch's setup checks, once: they are the same for every chunk); then the chunk's search reads
///   the subjects, where a read of the database that fails ends the stage (`ReadFails::Search`).
/// - `before_search`, where the chunk's search would read the subjects: the wrapped ranges of a
///   negative pair that leave no chunk without a query leave gaps between the chunks (a chunk is
///   then never split again), which NCBI searches on with the hits in the gaps lost, a setting the
///   maintainer kept refused (PD-LOSAT-NCBI-DEFECTS version 1.1, D13). It is refused at the first
///   chunk with a valid context with subjects in hand; without one nothing is lost, and the batch
///   is reported as NCBI's.
/// - `unmet_checks`, after the last chunk: the setup checks that no chunk with a valid context met
///   (each chunk without one), which the traceback's own setup of the whole batch meets, NCBI's
///   `Traceback failed` (exit status 3) after the batch's preliminary search found nothing.
pub(crate) struct SplitBatchSteps<'s> {
    program: &'static str,
    sizes: SplitSizes,
    read_fails: Option<ReadFails>,
    /// The batch's setup checks, until the first chunk with a valid context meets them.
    setup_checks: Option<anyhow::Result<()>>,
    stop: &'s dyn Fn(anyhow::Error) -> anyhow::Error,
}

impl<'s> SplitBatchSteps<'s> {
    /// The steps of `program`'s batch split into `chunks` by `sizes`, with the database's read
    /// failure `read_fails`, the batch's setup checks `setup_checks` and the program's form of a
    /// stop `stop`; a chunk without a query stops the batch here.
    pub(crate) fn new(
        program: &'static str,
        chunks: &[ProteinQueryChunk],
        sizes: SplitSizes,
        read_fails: Option<ReadFails>,
        setup_checks: anyhow::Result<()>,
        stop: &'s dyn Fn(anyhow::Error) -> anyhow::Error,
    ) -> anyhow::Result<Self> {
        // `split_protein_batch` makes no part for a chunk that no query reaches.
        if chunks.iter().any(|chunk| chunk.queries.is_empty()) {
            return Err(stop(anyhow::anyhow!(
                "a query chunk without a query is not supported by LOSAT's {program}"
            )));
        }
        Ok(Self {
            program,
            sizes,
            read_fails,
            setup_checks: Some(setup_checks),
            stop,
        })
    }

    /// The steps of `chunk` after its setup (`valid`: whether a context of it is valid) and before
    /// its search reads the subjects.
    pub(crate) fn before_read(
        &mut self,
        chunk: &ProteinQueryChunk,
        valid: bool,
    ) -> anyhow::Result<ChunkStep> {
        if protein_chunk_would_be_split(std::slice::from_ref(chunk), self.sizes) {
            return Err((self.stop)(anyhow::anyhow!(
                "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would split a query chunk again (it then stops with a null reference) is not supported by LOSAT's {}",
                self.program
            )));
        }
        if self.read_fails == Some(ReadFails::Setup) {
            return Ok(ChunkStep::End);
        }
        if !valid {
            return Ok(ChunkStep::Skip);
        }
        if let Some(checks) = self.setup_checks.take() {
            checks?;
        }
        if self.read_fails.is_some() {
            return Ok(ChunkStep::End);
        }
        Ok(ChunkStep::Read)
    }

    /// Whether `before_search` can refuse with subjects in hand (`subjects_in_hand`): only then
    /// does a caller whose chunks' searches find their own validity need it before the search.
    pub(crate) fn gap_step_pending(&self, subjects_in_hand: bool) -> bool {
        subjects_in_hand && self.sizes.negative_chunk_size()
    }

    /// The step of a chunk that `before_read` reads, where its search would read the subjects
    /// (`subjects_in_hand`: whether there is one): D13's refusal of the gaps.
    pub(crate) fn before_search(&self, subjects_in_hand: bool) -> anyhow::Result<()> {
        if self.gap_step_pending(subjects_in_hand) {
            return Err((self.stop)(anyhow::anyhow!(
                "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch with NCBI BLAST+'s size_t chunk ranges wrapped, which is not supported by LOSAT's {}",
                self.program
            )));
        }
        Ok(())
    }

    /// After the last chunk: the error of the batch's setup checks that no chunk with a valid
    /// context met, which the traceback's own setup of the batch meets; the caller makes NCBI's
    /// `Traceback failed` of its error of the setup, and a refusal (a count of 0) stays one.
    pub(crate) fn unmet_checks(&mut self) -> Option<anyhow::Error> {
        self.setup_checks.take().and_then(Result::err)
    }
}

/// Whether NCBI's search of a split batch gets to a read of the subjects (`ChunkStep::Read` or
/// `End`) with the steps `steps` of its chunks `chunks`, `valid` saying whether a context of a
/// chunk is valid: the set-up of a batch without its search, and where a failure of a call's setup
/// comes. `false`: every chunk finds nothing (`SplitBatchSteps::unmet_checks` then holds the
/// traceback's error).
pub(crate) fn split_batch_reaches_read(
    steps: &mut SplitBatchSteps<'_>,
    chunks: &[ProteinQueryChunk],
    mut valid: impl FnMut(&ProteinQueryChunk) -> anyhow::Result<bool>,
) -> anyhow::Result<bool> {
    for chunk in chunks {
        let valid = valid(chunk)?;
        match steps.before_read(chunk, valid)? {
            ChunkStep::Skip => continue,
            ChunkStep::Read | ChunkStep::End => return Ok(true),
        }
    }
    Ok(false)
}

#[cfg(test)]
#[path = "protein_query_split_tests.rs"]
mod tests;
