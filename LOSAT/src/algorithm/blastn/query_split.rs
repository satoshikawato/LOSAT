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

/// NCBI's chunk count (`Uint4`) and evened chunk size (`size_t`, 64-bit on every target).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:111-139
fn num_chunks_and_size(sizes: SplitSizes, concatenated_query_length: usize) -> (u32, u64) {
    let overlap_size = sizes.overlap;
    let length = concatenated_query_length as u64;
    // `num_chunks` is a `Uint4` (split_query_aux_priv.cpp:112).
    let mut num_chunks: u32 = 0;
    if sizes.chunk_size > overlap_size {
        num_chunks = (length / (sizes.chunk_size - overlap_size)) as u32;
    }
    if num_chunks <= 1 {
        return (1, length);
    }
    let chunks = u64::from(num_chunks);
    let mut chunk_size = length.wrapping_add((chunks - 1).wrapping_mul(overlap_size)) / chunks;
    if chunks < chunk_size.wrapping_sub(overlap_size) {
        chunk_size = chunk_size.wrapping_add(1);
    }
    (num_chunks, chunk_size)
}

/// The number of query chunks of a batch of `concatenated_query_length` residues, and the
/// chunk size NCBI then uses.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:123-138
pub fn calculate_num_chunks(sizes: SplitSizes, concatenated_query_length: usize) -> (usize, usize) {
    let (num_chunks, chunk_size) = num_chunks_and_size(sizes, concatenated_query_length);
    (num_chunks as usize, chunk_size as usize)
}

/// Whether NCBI splits a query batch whose queries have `concatenated_query_length` letters into
/// chunks (`calculate_num_chunks` gives more than one; `split_query_batch` makes them).
pub fn splits(sizes: SplitSizes, concatenated_query_length: usize) -> bool {
    num_chunks_and_size(sizes, concatenated_query_length).0 > 1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
/// LOSAT's refusal of a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE whose wrapped
/// chunk ranges leave gaps between chunks that all hold a query (decision D13 of
/// PD-LOSAT-NCBI-DEFECTS): NCBI searches on and loses the hits in the gaps. Where such a pair
/// leaves a chunk without a query, `CQuerySplitter::Split` stops first and LOSAT refuses that
/// empty chunk instead. The one text of BLASTN's, BLASTX's, BLASTP's, psiblast's and TBLASTN's
/// gap refusals.
pub fn negative_split_refusal(program: &str) -> String {
    format!(
        "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch with NCBI BLAST+'s chunk ranges wrapped, which is not supported by LOSAT's {program}"
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:234-249,876-880
/// LOSAT's refusal of a query chunk that holds no query: `CQuerySplitter::Split` dereferences its
/// missing query vector (null-pointer `CCoreException`) before any chunk is masked, set up, or
/// searched. The one text of BLASTN's, BLASTP's, psiblast's and TBLASTN's empty-chunk refusals.
pub fn empty_chunk_refusal(program: &str) -> String {
    format!("a query chunk without a query is not supported by LOSAT's {program}")
}

/// A query chunk of a split batch (`QueryChunk` for BLASTN, `ProteinQueryChunk` for protein
/// searches): its query parts in query order (`SplitBatchSteps`).
pub(crate) trait SplitChunk {
    fn chunk_queries(&self) -> &[ChunkQuery];
}

impl SplitChunk for QueryChunk {
    fn chunk_queries(&self) -> &[ChunkQuery] {
        &self.queries
    }
}

/// Whether NCBI's setup of some query chunk would split that chunk again (an overlap close to the
/// chunk size): each chunk's search builds its own query splitter, and only a debug build asserts
/// that it does not split; the release build then skips the chunk's lookup table and stops with a
/// null reference.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
pub(crate) fn chunk_would_be_split<C: SplitChunk>(chunks: &[C], sizes: SplitSizes) -> bool {
    chunks.iter().any(|chunk| {
        let length = chunk
            .chunk_queries()
            .iter()
            .map(|part| part.to - part.from)
            .sum();
        num_chunks_and_size(sizes, length).0 > 1
    })
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

/// A failed setup check of a split batch (`SplitBatchSteps::new`), as its program makes it: NCBI's
/// engine error of the setup (the subject length, the hit-saving parameters), which the
/// traceback's own setup of a batch whose chunks found nothing meets as `Traceback failed`
/// (`UnmetChecks::TracebackFailed`), or LOSAT's refusal there (a count of 0,
/// `common::database::check_setup_count`, or BLASTN's greedy setup refusal), which stays the
/// refusal.
#[derive(Debug)]
pub(crate) enum SetupCheckFailure {
    Engine(anyhow::Error),
    Refusal(anyhow::Error),
}

impl SetupCheckFailure {
    /// The error that the first chunk with a valid context meets (`SplitBatchSteps::before_read`).
    pub(crate) fn into_error(self) -> anyhow::Error {
        match self {
            SetupCheckFailure::Engine(error) | SetupCheckFailure::Refusal(error) => error,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:233-252 (cited at
// `SplitBatchSteps`)
/// The setup checks that no chunk with a valid context of a split batch met, as the traceback's
/// own setup of the batch meets them after the last chunk (`SplitBatchSteps::unmet_checks`).
#[derive(Debug)]
pub(crate) enum UnmetChecks {
    /// NCBI's engine error `Traceback failed` (`blastinput::app::TRACEBACK_FAILED`, exit status 3),
    /// which the caller makes in its program's form of an engine error of the setup.
    TracebackFailed,
    /// LOSAT's refusal of the checks (a count of 0, or BLASTN's greedy setup refusal).
    Refusal(anyhow::Error),
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
/// subjects, in NCBI's order (2.17.0 probes of line RG items X-short-bat, fix round 3, X-bpgap
/// and X-chunk): one owner for BLASTN's, BLASTP's, psiblast's and TBLASTN's searches of a split
/// batch, for the set-up of a batch without its search, and for where the failures of a call's
/// setup come when a split batch meets them first. Each stop is LOSAT's refusal in the program's
/// form (`stop`: BLASTN's and TBLASTN's `blastinput::app::ErrorAfterBatchReports`, BLASTP's batch
/// failure, psiblast's plain error); each error of NCBI's setup is that program's error there.
/// - `new`: the query splitter makes the query factory of every chunk first and dereferences the
///   missing query vector of a chunk that no query reaches (a null pointer), before any chunk is
///   set up. This is where a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE (NCBI's
///   `size_t` chunk ranges wrap) whose first range covers the batch or whose start reaches the
///   batch end stops, and a positive pair whose evened chunk size falls below or equals the
///   overlap at batch length at most the overlap (the start wraps or stops after chunk 0).
/// - `before_read`, chunk by chunk after the chunk's setup: unless `without_resplit_stop` skips
///   it for BLASTN (approved exception 1 of PD-LOSAT-NCBI-DEFECTS), the setup of a chunk that
///   would be split again makes no lookup table, which the chunk's search dereferences whether or
///   not its contexts are valid, after the chunks before it were searched; then the setup of the
///   chunk's preliminary search reads the length of OID 0 where the statistics' total is not
///   positive (`ReadFails::Setup`), checks the Karlin-Altschul blocks (a chunk without a valid
///   context ends there: NCBI ignores its error), then the subject length, hit-saving parameters
///   and greedy setup (the batch's setup checks, once: they are the same for every chunk); then
///   the chunk's search reads the subjects, where a read of the database that fails ends the stage
///   (`ReadFails::Search`).
/// - `before_search`, where the chunk's search would read the subjects: the wrapped ranges of a
///   negative pair that leave no chunk without a query leave gaps between the chunks (a chunk is
///   then never split again), which NCBI searches on with the hits in the gaps lost, a setting the
///   maintainer kept refused (PD-LOSAT-NCBI-DEFECTS version 1.1, D13). It is refused at the first
///   chunk with a valid context with subjects in hand; without one nothing is lost, and the batch
///   is reported as NCBI's.
/// - `unmet_checks`, after the last chunk: the setup checks that no chunk with a valid context met
///   (each chunk without one), which the traceback's own setup of the whole batch meets, NCBI's
///   `Traceback failed` (exit status 3) after the batch's preliminary search found nothing
///   (`UnmetChecks`: the one owner of that error of BLASTN's, BLASTP's, psiblast's and TBLASTN's).
pub(crate) struct SplitBatchSteps<'s> {
    program: &'static str,
    sizes: SplitSizes,
    /// Whether a chunk that NCBI's setup would split again stops the batch (`before_read`): true
    /// for BLASTP, psiblast and TBLASTN; false for BLASTN (approved exception 1 of
    /// PD-LOSAT-NCBI-DEFECTS).
    check_resplit: bool,
    read_fails: Option<ReadFails>,
    /// The batch's setup checks, until the first chunk with a valid context meets them.
    setup_checks: Option<Result<(), SetupCheckFailure>>,
    stop: &'s dyn Fn(anyhow::Error) -> anyhow::Error,
}

impl<'s> SplitBatchSteps<'s> {
    /// The steps of `program`'s batch split into `chunks` by `sizes`, with the database's read
    /// failure `read_fails`, the batch's setup checks `setup_checks` and the program's form of a
    /// stop `stop`; a chunk without a query stops the batch here.
    pub(crate) fn new<C: SplitChunk>(
        program: &'static str,
        chunks: &[C],
        sizes: SplitSizes,
        read_fails: Option<ReadFails>,
        setup_checks: Result<(), SetupCheckFailure>,
        stop: &'s dyn Fn(anyhow::Error) -> anyhow::Error,
    ) -> anyhow::Result<Self> {
        // `split_query_batch` / `split_protein_batch` makes no part for a chunk that no query
        // reaches.
        if chunks.iter().any(|chunk| chunk.chunk_queries().is_empty()) {
            return Err(stop(anyhow::Error::msg(empty_chunk_refusal(program))));
        }
        Ok(Self {
            program,
            sizes,
            check_resplit: true,
            read_fails,
            setup_checks: Some(setup_checks),
            stop,
        })
    }

    /// Skips the re-split stop in `before_read`: BLASTN searches a chunk that NCBI's setup would
    /// split again once (approved exception 1 of PD-LOSAT-NCBI-DEFECTS).
    pub(crate) fn without_resplit_stop(mut self) -> Self {
        self.check_resplit = false;
        self
    }

    /// The steps of `chunk` after its setup (`valid`: whether a context of it is valid) and before
    /// its search reads the subjects.
    pub(crate) fn before_read<C: SplitChunk>(
        &mut self,
        chunk: &C,
        valid: bool,
    ) -> anyhow::Result<ChunkStep> {
        if self.check_resplit && chunk_would_be_split(std::slice::from_ref(chunk), self.sizes) {
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
            checks.map_err(SetupCheckFailure::into_error)?;
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
            return Err((self.stop)(anyhow::Error::msg(negative_split_refusal(
                self.program,
            ))));
        }
        Ok(())
    }

    /// After the last chunk: the batch's setup checks that no chunk with a valid context met, which
    /// the traceback's own setup of the batch meets: NCBI's `Traceback failed` for an engine error
    /// of the setup, and a refusal (a count of 0, or BLASTN's greedy setup refusal) stays one.
    pub(crate) fn unmet_checks(&mut self) -> Option<UnmetChecks> {
        match self.setup_checks.take()? {
            Ok(()) => None,
            Err(SetupCheckFailure::Engine(_)) => Some(UnmetChecks::TracebackFailed),
            Err(SetupCheckFailure::Refusal(refusal)) => Some(UnmetChecks::Refusal(refusal)),
        }
    }
}

/// Whether NCBI's search of a split batch gets to a read of the subjects (`ChunkStep::Read` or
/// `End`) with the steps `steps` of its chunks `chunks`, `valid` saying whether a context of a
/// chunk is valid: the set-up of a batch without its search, and where a failure of a call's setup
/// comes. `false`: every chunk finds nothing (`SplitBatchSteps::unmet_checks` then holds the
/// traceback's error).
pub(crate) fn split_batch_reaches_read<C: SplitChunk>(
    steps: &mut SplitBatchSteps<'_>,
    chunks: &[C],
    mut valid: impl FnMut(&C) -> anyhow::Result<bool>,
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

/// The query chunks of a batch with queries of `lengths` residues, or `None` when NCBI
/// does not split it.
///
/// The chunk ranges on the concatenated queries:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:145-171
/// A chunk the loop leaves unset (or gives an empty range at the batch's end) keeps no query,
/// and NCBI's `CQuerySplitter::Split` then fails to make its query factory with a null-pointer
/// `CCoreException` before any chunk is set up (an evened chunk size below or equal to a positive
/// overlap when the batch is no longer than the overlap, or a negative pair whose wrapped range
/// covers or reaches the batch's end before the last chunk).
///
/// The queries of a chunk, and their parts in it:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:223-247
/// (`TChunkRange` is a `COpenRange`, so `GetTo()` plus 1 is the end of the previous query.)
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
    let (num_chunks, chunk_size) = num_chunks_and_size(sizes, total_length);
    if num_chunks <= 1 {
        return None;
    }
    let total_length_u64 = total_length as u64;
    let num_chunks_usize = num_chunks as usize;
    // NCBI's 64-bit `size_t` overlap: near 2^64 for a negative OVERLAP_CHUNK_SIZE (which splits a
    // batch only below a negative CHUNK_SIZE). Keeping the loop in `u64` preserves that wrap on
    // 32-bit targets too.
    let overlap_size = sizes.overlap;

    let mut chunk_ranges = vec![(0u64, 0u64); num_chunks_usize];
    let mut chunk_start = 0u64;
    for (chunk_num, range) in chunk_ranges.iter_mut().enumerate() {
        let mut chunk_end = chunk_start.wrapping_add(chunk_size);
        if chunk_end >= total_length_u64 || chunk_num + 1 == num_chunks_usize {
            chunk_end = total_length_u64;
        }
        *range = (chunk_start, chunk_end);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:166-170
        // An evened chunk size below a positive overlap wraps the start past the queries' length at
        // the first chunk (the loop ends, the chunks after it keep no query); an evened chunk size
        // equal to the overlap with total length equal to the overlap hits `chunk_end == total` and
        // ends after chunk 0; a negative overlap moves the start by the chunk size plus the
        // overlap's magnitude (gaps between chunks, or past/at the batch's end leaving empty chunks).
        chunk_start = chunk_start.wrapping_add(chunk_size.wrapping_sub(overlap_size));
        if chunk_start > total_length_u64 || chunk_end == total_length_u64 {
            break;
        }
    }

    let mut query_ranges = Vec::with_capacity(lengths.len());
    let mut query_start = 0u64;
    for &length in lengths {
        let len = length as u64;
        query_ranges.push((query_start, query_start + len));
        query_start += len;
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
                    (query_to - query_from) as usize
                } else {
                    (chunk_to - query_from) as usize
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

    let overlap_usize = usize::try_from(overlap_size).unwrap_or(usize::MAX);
    let offsets: Vec<Vec<i32>> = (0..chunks.len())
        .map(|chunk_num| context_offsets(&chunks, chunk_num, lengths, overlap_usize))
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
            // `BlastSeqLoc_RestrictToInterval` of `mask` to `interval_from` and `interval_to`.
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
/// `offset - 1`) is in its masks: lower case, with `-lcase_masking`, or joined by the masks of
/// `-filtering_db` (`repeat_filter::MasksOutside::before`). NCBI restricts the lower-case runs of the
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
#[path = "query_split_tests.rs"]
mod tests;
