//! The typed results of a search, one of NCBI's query batches at a time.
//!
//! NCBI's applications read a batch of queries, search it, print its results and read the next
//! batch, so a run holds the results of one batch at a time. `run_batches` is that loop for a
//! caller: each step searches one batch and hands over its records, and the caller stops the run
//! by not asking for the next batch, which is also where it reports progress.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-288
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!         ...
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//!             ...
//!                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//!                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//!                 results = lcl_blast.Run();
//!             ...
//!                 ITERATE(CSearchResultSet, result, *results) {
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!                 }
//! ```
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
use super::run::{
    check_inputs, query_ordinal, query_range_of, records, run_local, searched_queries,
    CountingWriter, HspRecord, QueryStats,
};
use super::search::Search;
use super::start::{CallStart, RunLog};
use super::view::{self, QueryHits};
use crate::algorithm::blastn::blast_engine::QueryBatchSizes;
use crate::blastinput::app::query_batch_size;
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
    /// What the command line writes to its standard error for this batch: the warnings of its
    /// queries, and with the first batch those of the run.
    pub diagnostics: Vec<u8>,
}

impl Batch {
    /// The results of each query of the batch, in the order of the inputs.
    pub fn queries(&self) -> impl Iterator<Item = QueryHits<'_>> + '_ {
        let first = self.query_indices.start;
        self.query_stats
            .iter()
            .enumerate()
            .map(move |(offset, stats)| {
                view::query_hits(&self.hsps, first + offset, stats.as_ref())
            })
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-272
// ```c++
//         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
//         int batch_size = m_CmdLineArgs->GetQueryBatchSize();
//         if (batch_size) {
//             input.SetBatchSize(batch_size);
//         ...
//         } else {
//             Int8 total_len = formatter.GetDbTotalLength();
//             if (total_len > 0) {
//                 /* the optimal hits per batch scales with total db size */
//                 mixer.SetTargetHits(total_len / 3000);
//             }
//             input.SetBatchSize(mixer.GetBatchSize());
// ```
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
        // ```c++
        //             input.Reset(new CBlastInput(&*fasta, 1));
        // ```
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
    query_stats: Vec<Option<QueryStats>>,
    extensions: Option<i32>,
    /// The sequences that a database search's hit list names.
    database_subjects: Option<Vec<ReportedSubject>>,
}

impl FormatObserver for BatchObserver {
    fn hsp_begin(&mut self, _format: usize, _hsp: HspIndex) {}

    fn hsp_end(&mut self, _format: usize, _hsp: HspIndex) {}

    fn query_stats(&mut self, q_idx: usize, stats: XmlStats) {
        if self.query_stats.len() <= q_idx {
            self.query_stats.resize(q_idx + 1, None);
        }
        self.query_stats[q_idx] = QueryStats::of_context(stats);
    }

    fn queries_begin(&mut self) {
        self.queries_begin = Some(self.diagnostics.load(Ordering::SeqCst));
    }

    fn batch_extensions(&mut self, extensions: i32) {
        self.extensions = Some(extensions);
    }

    fn database_subjects(&mut self, subjects: Vec<ReportedSubject>) {
        self.database_subjects = Some(subjects);
    }
}

/// The search of one query batch at a time (`run_batches`). Each `next` searches one batch; a
/// caller that drops the iterator stops the search there.
pub struct Batches<'a> {
    search: &'a Search,
    inputs: Inputs<'a>,
    sizes: BatchSizes,
    /// The letters of each query record, which size the batches.
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
/// then do. A batch is searched by one `run_local` call, each with its own thread pool of
/// `-num_threads` threads; of a database, each opens the database again, so what NCBI writes once
/// per run about a taxonomy database comes with every batch that reads one.
pub fn run_batches<'a>(search: &'a Search, inputs: Inputs<'a>) -> Result<Batches<'a>, Error> {
    let mut start_messages = Vec::new();
    let mut start = CallStart::new(search, inputs, &mut start_messages)?;
    let lengths: Vec<usize> = inputs
        .query
        .iter()
        .map(|record| record.seq().len())
        .collect();
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
    let sizes = batch_sizes(&command, inputs, &mut start, &mut start_messages)
        .map_err(|error| Error::search(&error))?;
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
    })
}

impl Batches<'_> {
    /// How many query records the inputs hold: with `Batch::query_indices`, how far the search
    /// is.
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

    /// Searches the batch `range`, as the command line's loop searches it.
    fn search_batch(&mut self, range: Range<usize>) -> Result<Batch, Error> {
        let inputs = Inputs {
            query: &self.inputs.query[range.clone()],
            ..self.inputs
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
            Commands::Psiblast(args) => args.first_query_ordinal = range.start,
            // BLASTX's messages show the record's own id; rpsblast and rpstblastn are refused
            // by `Program::parse`.
            Commands::Blastx(_)
            | Commands::Rpsblast(_)
            | Commands::Rpstblastn(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => {}
        }
        let ordinals = searched_queries(query_range, inputs.query);
        let mut diagnostics = CountingWriter::new();
        let mut observer = BatchObserver {
            diagnostics: diagnostics.position.clone(),
            queries_begin: None,
            query_stats: Vec::new(),
            extensions: None,
            database_subjects: None,
        };
        let mut collected: Option<Arc<Vec<PairwiseHit>>> = None;
        let mut report = io::sink();
        // The batch's messages are gathered, so that a later batch's repeat of the run's own
        // messages is cut before they go where the run's messages go (`RunLog::route`: tblastx's
        // `-logfile` file), on the way of a failure too.
        let result = {
            let mut on_hits = |hits| collected = Some(hits);
            let mut sink = self.log.sink(&mut diagnostics, false);
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut report),
                }],
                diagnostics: &mut sink,
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
        result.map_err(|error| self.log.failed(&error))?;
        let ordinals = ordinals.as_deref();
        let mut hsps = match &collected {
            Some(hits) => records(
                self.search.program(),
                inputs,
                hits,
                query_range,
                ordinals,
                observer.database_subjects.as_deref(),
            ),
            None => Vec::new(),
        };
        for hsp in &mut hsps {
            hsp.index += self.hsps_before;
            hsp.q_idx += range.start;
        }
        self.hsps_before += hsps.len();
        let mut query_stats = vec![None; inputs.query.len()];
        for (q_idx, stats) in observer.query_stats.iter().enumerate() {
            if let Some(slot) = query_stats.get_mut(query_ordinal(ordinals, q_idx)) {
                *slot = *stats;
            }
        }
        if let (BatchSizes::Blastn(sizes), Some(extensions)) =
            (&mut self.sizes, observer.extensions)
        {
            sizes.searched(extensions);
        }
        Ok(Batch {
            query_indices: range,
            hsps,
            query_stats,
            diagnostics: bytes,
        })
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
        self.finished = batch.is_err();
        Some(batch)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::{run_collect, validate, FastaRecord, Program};

    const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

    #[test]
    fn a_batch_ends_where_it_reaches_its_size_and_a_query_past_the_interval_adds_nothing() {
        let query: Vec<FastaRecord> = (0..4)
            .map(|index| FastaRecord::new(format!("Query_{}", index + 1), b"q", PROTEIN))
            .collect();
        let subject = vec![FastaRecord::new("Subject_1", b"s", PROTEIN)];
        let search = validate(Program::Blastp, &[]).unwrap();
        let batches = run_batches(&search, Inputs::new(&query, &subject)).unwrap();
        assert_eq!(batches.query_count(), 4);
        let mut batches = batches;
        // The default size is 10,000 letters: one batch.
        let only = batches.next().unwrap().unwrap();
        assert_eq!(only.query_indices, 0..4);
        assert_eq!(only.hsps.len(), 4);
        assert!(batches.next().is_none());
        assert!(batches.next().is_none());
        let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
        assert_eq!(only.hsps, outcome.hsps);
        assert_eq!(only.query_stats, outcome.query_stats);
    }

    #[test]
    fn a_failure_ends_the_batches() {
        let query = vec![FastaRecord::new("Query_1", b"q", PROTEIN)];
        let search = validate(Program::Blastp, &[]).unwrap();
        let mut batches = run_batches(&search, Inputs::new(&query, &[])).unwrap();
        assert!(matches!(batches.next(), Some(Err(Error::Search(_)))));
        assert!(batches.next().is_none());
    }

    #[test]
    fn queries_without_a_record_are_one_call() {
        let subject = vec![FastaRecord::new("Subject_1", b"s", PROTEIN)];
        let search = validate(Program::Blastp, &[]).unwrap();
        let mut batches = run_batches(&search, Inputs::new(&[], &subject)).unwrap();
        let only = batches.next().unwrap().unwrap();
        assert_eq!(only.query_indices, 0..0);
        assert_eq!(
            String::from_utf8_lossy(&only.diagnostics),
            "Warning: [blastp] Query is Empty!\n"
        );
        assert!(batches.next().is_none());
    }
}
