//! The typed results of a search, one of NCBI's query batches at a time.
//!
//! NCBI's applications read a batch of queries, search it, print its results and read the next
//! batch, so a run holds the results of one batch at a time. `run_batches` is that loop for a
//! caller: each step searches one batch and hands over its records, and the caller stops the run
//! by not asking for the next batch, which is also where it reports progress.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-288
#![warn(missing_docs, unused)]

use std::io;
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use super::error::Error;
use super::input::Inputs;
use super::local_blast::{
    FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs, ReportedSubject,
};
use super::parts::Parts;
use super::run::{
    check_inputs, query_range_of, run_local, searched_queries, CountingWriter, HspRecord,
    QueryRounds, QueryStats, RunEnd,
};
use super::search::Search;
use super::start::{CallStart, RunLog};
use super::view::{self, QueryHits};
use crate::algorithm::blastn::blast_engine::QueryBatchSizes;
use crate::blastinput::app::query_batch_size;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::{next_ranged_batch_end, record_interval, RecordInterval};
use crate::cli::Commands;
use crate::report::blast_xml::XmlStats;
use crate::report::PairwiseHit;

/// The results of one query batch of a search (`run_batches`).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Batch {
    /// The query records of the batch: a range of indices in the inputs.
    pub query_indices: Range<usize>,
    /// The HSPs of the batch's queries. `q_idx` counts in the inputs and `index` over the
    /// batches of the search, as they do in the whole search's `Outcome::hsps`.
    pub hsps: Vec<HspRecord>,
    /// One entry per query of `query_indices`, as in `Outcome::query_stats`.
    pub query_stats: Vec<Option<QueryStats>>,
    /// One entry per query of `query_indices`, as in `Outcome::rounds`.
    pub rounds: Vec<QueryRounds>,
    /// What the command line writes to its standard error for this batch: the warnings of its
    /// queries, and with the first batch those of the run.
    pub diagnostics: Vec<u8>,
}

impl Batch {
    /// The batch `range` in which the run ended without a post and with exit status 0
    /// (`run::RunEnd::WithoutPost`: NCBI's end at a header record that ends before its read
    /// enters the set's elements): no HSPs, the statistics and rounds of queries that were not
    /// searched, and what the command line writes to its standard error for the batch,
    /// `diagnostics`.
    fn ended(range: Range<usize>, diagnostics: Vec<u8>) -> Self {
        let queries = range.len();
        Batch {
            query_indices: range,
            hsps: Vec::new(),
            query_stats: vec![None; queries],
            rounds: vec![
                QueryRounds {
                    stats: vec![None],
                    converged: false,
                };
                queries
            ],
            diagnostics,
        }
    }

    /// The results of each query of the batch, in the order of the inputs.
    pub fn queries(&self) -> impl Iterator<Item = QueryHits<'_>> + '_ {
        let first = self.query_indices.start;
        self.query_stats
            .iter()
            .enumerate()
            .map(move |(offset, stats)| {
                view::query_hits(
                    &self.hsps,
                    first + offset,
                    stats.as_ref(),
                    self.rounds.get(offset),
                )
            })
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-272
/// The size of NCBI's query batches: BLASTN's follows the extensions of the batch before.
enum BatchSizes {
    Fixed(u32),
    Blastn(QueryBatchSizes),
}

impl BatchSizes {
    fn current(&self) -> u32 {
        match self {
            BatchSizes::Fixed(size) => *size,
            BatchSizes::Blastn(sizes) => sizes.batch_size(),
        }
    }
}

/// The batch sizes of `command`, as each program's application reads them (`GetQueryBatchSize`);
/// BLASTN's of a database reads the call's settings (`start`, its messages on `diagnostics`).
fn batch_sizes(
    command: &Commands,
    inputs: Inputs<'_>,
    start: &mut CallStart,
    diagnostics: &mut dyn io::Write,
) -> anyhow::Result<BatchSizes> {
    Ok(match command {
        Commands::Blastn(args) => BatchSizes::Blastn(match args.database_name() {
            Some(name) => {
                let settings = start.settings(diagnostics)?;
                QueryBatchSizes::of_database(args, name, &settings)?
            }
            None => QueryBatchSizes::new(args, inputs.subject)?,
        }),
        Commands::Blastp(_) => BatchSizes::Fixed(query_batch_size(10000)?),
        Commands::Tblastn(_) => BatchSizes::Fixed(query_batch_size(20000)?),
        Commands::Tblastx(_) => BatchSizes::Fixed(query_batch_size(10002)?),
        Commands::Blastx(_) => BatchSizes::Fixed(query_batch_size(
            crate::algorithm::blastx::QUERY_BATCH_SIZE as u32,
        )?),
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:468-468
        Commands::Psiblast(_) => BatchSizes::Fixed(1),
        Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("Program::parse admits the search programs only")
        }
    })
}

/// What a batch's search tells beside its hit list.
struct BatchObserver {
    diagnostics: Arc<AtomicUsize>,
    /// Where the run's own messages end in the diagnostics (`queries_begin`).
    queries_begin: Option<usize>,
    extensions: Option<i32>,
    /// The statistics, database ranks and rounds of the batch's search calls.
    parts: Parts,
}

impl FormatObserver for BatchObserver {
    fn hsp_begin(&mut self, _format: usize, _hsp: HspIndex) {}

    fn hsp_end(&mut self, _format: usize, _hsp: HspIndex) {}

    fn query_stats(&mut self, q_idx: usize, stats: XmlStats) {
        self.parts.query_stats(q_idx, stats);
    }

    fn queries_begin(&mut self) {
        self.queries_begin = Some(self.diagnostics.load(Ordering::SeqCst));
    }

    fn batch_extensions(&mut self, extensions: i32) {
        self.extensions = Some(extensions);
    }

    fn database_subjects(&mut self, subjects: Vec<ReportedSubject>) {
        self.parts.database_subjects(subjects);
    }

    fn round_begin(&mut self, query: usize, round: u32) {
        self.parts.round_begin(query, round);
    }

    fn rounds_converged(&mut self, query: usize) {
        self.parts.rounds_converged(query);
    }

    fn pssm_query(&mut self, query: &FastaRecord) {
        self.parts.pssm_query(query);
    }
}

/// The search of one query batch at a time (`run_batches`). Each `next` searches one batch; a
/// caller that drops the iterator stops the search there.
pub struct Batches<'a> {
    search: &'a Search,
    inputs: Inputs<'a>,
    sizes: BatchSizes,
    /// The letters of each query record, which size the batches (one query of no letters for a
    /// PSSM query, one batch).
    lengths: Vec<usize>,
    /// The query records that NCBI's reader skips (`-query_loc`).
    skipped: Vec<bool>,
    next_query: usize,
    hsps_before: usize,
    searched_one: bool,
    finished: bool,
    /// The call's application start, which the batches share.
    start: CallStart,
    /// The `-logfile` target of the batches' run messages (tblastx's).
    log: RunLog,
    /// What the start wrote, which comes with the first batch.
    start_messages: Vec<u8>,
    /// Whether the call is a dry run: one batch of no query, with the start's messages.
    dry_run: bool,
    /// The run's numbers of psiblast's `-save_each_pssm` files, which go on over the batches as
    /// over the command line's (`PsiblastArgs::pssm_versions`).
    pssm_versions: crate::algorithm::psiblast::pssm::PssmVersions,
}

/// Searches the queries one of NCBI's query batches at a time and returns the batches as an
/// iterator; the search happens in the iterator's `next`. The batches are those of the command
/// line (by default 10,000 letters for BLASTP, 20,000 for TBLASTN, 10,002 for BLASTX and TBLASTX,
/// and for BLASTN a size that follows the extensions of the batch before), so the records of the
/// batches, one after the other, are those of `run_collect` with a tabular format: the same HSPs,
/// statistics and diagnostics.
///
/// A caller reports progress between batches (`Batch::query_indices`), and stops the search by not
/// asking for the next batch. A batch that fails gives its error (`Error::Search`) and ends the
/// iterator: the batches before it stand, as the reports that the command line has written by
/// then do. A batch in which the run ends with exit status 0 and no message (a database header
/// record that ends before NCBI's read enters its elements, which NCBI's search programs take as
/// the end of the run) is the last: it has no HSPs and carries its diagnostics (with the first
/// batch, the run's own), what the command line writes to its standard error up to that end. A batch is searched by one `run_local` call, each with its own thread pool of
/// `-num_threads` threads; of a database, each opens the database again, so what NCBI writes once
/// per run about a taxonomy database comes with every batch that reads one.
pub fn run_batches<'a>(search: &'a Search, inputs: Inputs<'a>) -> Result<Batches<'a>, Error> {
    let mut start_messages = Vec::new();
    let mut start = CallStart::new(search, inputs, &mut start_messages)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:534-543 (the PSSM's query
    // searched without a `CBlastInput`, cited at `blastp::blast_engine::search_pssm_query`)
    // A PSSM or alignment query (`Inputs::pssm`, `Inputs::msa`) is one batch of its own: NCBI
    // reads no query batches then.
    let lengths: Vec<usize> = match inputs.query_file {
        Some(_) => vec![0],
        None => inputs
            .query
            .iter()
            .map(|record| record.seq().len())
            .collect(),
    };
    if start.dry_run() {
        let log = start.take_log();
        return Ok(Batches {
            search,
            inputs,
            sizes: BatchSizes::Fixed(0),
            skipped: vec![false; lengths.len()],
            lengths,
            next_query: 0,
            hsps_before: 0,
            searched_one: false,
            finished: false,
            start,
            log,
            start_messages,
            dry_run: true,
            pssm_versions: Default::default(),
        });
    }
    check_inputs(search, inputs)?;
    // Under `-parse_deflines` with `-subject`, NCBI's scope keeps every batch's queries
    // (`query_batch::QueryScopes`): two queries of the run with one id are refused before the
    // first batch, as the command line refuses them; each batch's call holds its own records alone.
    if search.parse_deflines() && !inputs.database {
        let program = match search.program() {
            crate::api::Program::Blastn => "BLASTN",
            crate::api::Program::Blastp => "BLASTP",
            crate::api::Program::Blastx => "BLASTX",
            crate::api::Program::Tblastn => "TBLASTN",
            crate::api::Program::Tblastx => "TBLASTX",
            crate::api::Program::Psiblast => "PSIBLAST",
        };
        crate::blastinput::query_batch::QueryScopes::of_records(
            crate::blastinput::query_batch::QueryScopeSpan::Run,
            program,
            inputs.query,
        )
        .map_err(|error| Error::search(&error))?;
    }
    let command = search.command(inputs)?;
    // A PSSM query is one batch, whose size NCBI does not read (`BATCH_SIZE`; tblastn_app.cpp:
    // 208-221 reads it for a FASTA query only).
    let sizes = match inputs.query_file {
        Some(_) => BatchSizes::Fixed(1),
        None => batch_sizes(&command, inputs, &mut start, &mut start_messages)
            .map_err(|error| Error::search(&error))?,
    };
    let query_range = query_range_of(&command);
    let skipped = lengths
        .iter()
        .map(|&length| {
            query_range
                .is_some_and(|range| record_interval(&range, length) == RecordInterval::PastEnd)
        })
        .collect();
    Ok(Batches {
        search,
        inputs,
        sizes,
        lengths,
        skipped,
        next_query: 0,
        hsps_before: 0,
        searched_one: false,
        finished: false,
        log: start.take_log(),
        start,
        start_messages,
        dry_run: false,
        pssm_versions: Default::default(),
    })
}

impl Batches<'_> {
    /// How many queries the inputs hold (the records, or the one sequence of a PSSM or alignment
    /// file, `Inputs::pssm`, `Inputs::msa`): with `Batch::query_indices`, how far the search is.
    pub fn query_count(&self) -> usize {
        self.lengths.len()
    }

    /// The records of the next batch. A batch size of 0 reads no record: the rest goes to the
    /// program, which fails as the command line fails (no batch follows either way).
    fn next_range(&self) -> Range<usize> {
        let start = self.next_query;
        let end = next_ranged_batch_end(&self.lengths, &self.skipped, start, self.sizes.current());
        start..if end > start { end } else { self.lengths.len() }
    }

    /// Searches the batch `range`, as the command line's loop searches it, and whether the run
    /// ended in it without a post and with exit status 0 (`run::RunEnd::WithoutPost`), where
    /// NCBI's loop ends: no batch is searched after it, and the batch itself has no HSPs (NCBI
    /// wrote no report of it) and its diagnostics (`Batch::ended`).
    fn search_batch(&mut self, range: Range<usize>) -> Result<(Batch, bool), Error> {
        let inputs = match self.inputs.query_file {
            Some(_) => self.inputs,
            None => Inputs {
                query: &self.inputs.query[range.clone()],
                ..self.inputs
            },
        };
        let mut command = self.search.command(inputs)?;
        let query_range = query_range_of(&command);
        match &mut command {
            Commands::Blastn(args) => {
                args.one_query_batch = true;
                args.first_query_ordinal = range.start;
            }
            Commands::Blastp(args) => args.first_query_ordinal = range.start,
            Commands::Tblastn(args) => args.first_query_ordinal = range.start,
            Commands::Tblastx(args) => args.first_query_ordinal = range.start,
            Commands::Psiblast(args) => {
                use crate::algorithm::psiblast::pssm::PssmCheckScope;
                args.first_query_ordinal = range.start;
                args.pssm_versions = self.pssm_versions.clone();
                // The command line checks the PSSM files of all its queries before its first
                // byte: the first batch checks the run's, where the command line does.
                args.pssm_check = if self.searched_one {
                    PssmCheckScope::Checked
                } else {
                    PssmCheckScope::Run {
                        queries: self.lengths.len(),
                    }
                };
            }
            // BLASTX's messages show the record's own id; its prefetch's stop asks whether the
            // call comes after others (`local::DatabaseReports::unrecorded_calls`).
            Commands::Blastx(args) => args.first_query_ordinal = range.start,
            // rpsblast and rpstblastn are refused by `Program::parse`.
            Commands::Rpsblast(_)
            | Commands::Rpstblastn(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => {}
        }
        let ordinals = searched_queries(query_range, inputs.query);
        let mut diagnostics = CountingWriter::new();
        let mut observer = BatchObserver {
            diagnostics: diagnostics.position.clone(),
            queries_begin: None,
            extensions: None,
            parts: Parts::new(),
        };
        let mut collected: Vec<Arc<Vec<PairwiseHit>>> = Vec::new();
        let mut report = io::sink();
        // The batch's messages are gathered, so that a later batch's repeat of the run's own
        // messages is cut before they go where the run's messages go (`RunLog::route`: tblastx's
        // `-logfile` file), on the way of a failure too.
        let result = {
            let mut on_hits = |hits| collected.push(hits);
            let mut sink = self.log.sink(&mut diagnostics, false);
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut report),
                }],
                diagnostics: &mut sink,
                later_warnings: None,
                hits: Some(&mut on_hits),
                observer: Some(&mut observer),
            };
            run_local(command, inputs, &mut outputs, &mut self.start)
        };
        // The run's own messages come once, with the first batch, after the start's.
        let mut bytes = diagnostics.bytes;
        if self.searched_one {
            bytes.drain(..observer.queries_begin.unwrap_or(0));
        }
        let mut bytes = self.log.route(bytes);
        if !self.searched_one {
            bytes.splice(0..0, std::mem::take(&mut self.start_messages));
        }
        if result.map_err(|error| self.log.failed(&error))? == RunEnd::WithoutPost {
            return Ok((Batch::ended(range, bytes), true));
        }
        let finished = observer
            .parts
            .finish(
                self.search.program(),
                inputs,
                &collected,
                query_range,
                ordinals.as_deref(),
                self.hsps_before,
            )
            .map_err(|fault| {
                Error::Search(format!("internal error in the search's parts: {fault}").into())
            })?;
        let mut hsps = finished.hsps;
        for hsp in &mut hsps {
            hsp.q_idx += range.start;
        }
        self.hsps_before += hsps.len();
        if let (BatchSizes::Blastn(sizes), Some(extensions)) =
            (&mut self.sizes, observer.extensions)
        {
            sizes.searched(extensions);
        }
        Ok((
            Batch {
                query_indices: range,
                hsps,
                query_stats: finished.query_stats,
                rounds: finished.rounds,
                diagnostics: bytes,
            },
            false,
        ))
    }
}

impl Iterator for Batches<'_> {
    type Item = Result<Batch, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        // A dry run searches nothing (`Search::dry_run`).
        if self.dry_run {
            if self.finished {
                return None;
            }
            self.finished = true;
            return Some(Ok(Batch {
                query_indices: 0..0,
                hsps: Vec::new(),
                query_stats: Vec::new(),
                rounds: Vec::new(),
                diagnostics: std::mem::take(&mut self.start_messages),
            }));
        }
        // Queries without a record are one call: NCBI's `Query is Empty!`.
        let done = self.next_query >= self.lengths.len() && self.searched_one;
        if self.finished || done {
            return None;
        }
        let range = self.next_range();
        self.next_query = range.end;
        let batch = self.search_batch(range);
        self.searched_one = true;
        self.finished = !matches!(batch, Ok((_, false)));
        Some(batch.map(|(batch, _)| batch))
    }
}

#[cfg(test)]
#[path = "batches_tests.rs"]
mod tests;
