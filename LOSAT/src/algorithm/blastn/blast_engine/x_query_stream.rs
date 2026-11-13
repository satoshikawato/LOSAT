//! EXPERIMENT (LOSAT_X_QSTREAM, in the strict set since 2026-10-18: on unless it is `0`): blastn
//! reads the records of each query batch from the query file when the batch loop of
//! `run_in_pool` takes the batch, as NCBI reads its query batches one at a time, instead of
//! holding the records of the whole query file for the whole search. It saves memory (the records of the batches not in flight) and costs
//! one more parse of the query file; where it wins and where it loses on the workstation is in its
//! entry of `docs/agents/diagnostics.md`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//! 	    BLAST_PROF_START( APP.LOOP.PRE );
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-146
//! ```c++
//! CRef<CBlastQueryVector>
//! CBlastInput::GetNextSeqBatch(CScope& scope)
//! {
//!     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
//!     TSeqPos size_read = 0;
//!
//!     while (size_read < GetBatchSize()) {
//!
//!         if (End())
//!             break;
//!
//!         CRef<CBlastSearchQuery> q;
//!         try { q.Reset(m_Source->GetNextSequence(scope)); }
//! ```
//!
//! Without the switch `search_cli` reads every record of the query file before the search
//! (`read_queries`) and the loop cuts its batches from them. LOSAT needs, before its first batch,
//! what NCBI learns only while it reads: the length of every record (the batch planner and size
//! mixer, LOSAT_X_RPTSTREAM's checks, `-query_loc`'s skipped records) and how the reading ends (an `eEOF`, a reader error, LOSAT's rejection of a Seq-id
//! line), whose message comes in the batch that NCBI reads it in. With the switch the file is read
//! twice:
//! - the first reading is today's reading, on the same source and at the same place (after
//!   `IsIStreamEmpty`, before the search's checks and the outfmt 0 prolog): every record is read
//!   with the same reader, by the loop of `read_queries` (`scan`), and dropped at once; its length
//!   is kept, and so is how the reading ended, with the reader's error and warnings as
//!   `read_queries` returns them. The query input (`QueryInput`: the placements, input positions
//!   and skipped records of `-query_loc`, the range length) and the searched lengths are made from
//!   the lengths alone (`query_input_of`, what `cut_queries` and `QueryInput::whole` make of the
//!   records). So everything that comes before the search, and every refusal, warning and error
//!   that comes before the first output today, comes from the same values at the same place;
//! - the second reading is a new source over a second handle of the same open file (`try_clone`),
//!   at the position where the first started, read forward only as far as the batches need: the
//!   loop takes its batch's records when it reaches the batch (`Queries::batch`). The records read
//!   wait in a window of chunks, each a run of consecutive records read at once (`read_to`), which
//!   the loop takes when it reaches their batch. The loop's batch is held until the batch is done,
//!   then dropped. Where the first reading was read in parts (`x_query_scan_par.rs`), the second
//!   reads the same parts again (`State::read_parts`): when the loop takes its batch, the parts up
//!   to the one that holds the batch's last record, and at least one per thread, each by a reader
//!   of its own on the search's threads (`x_query_scan_par::read_again`), one chunk per part; the
//!   records of those parts after the batch wait in the window for the next batches. Without the
//!   parts, the second reading runs on the loop's thread between the batches, while the search's
//!   threads wait: for megablast at 16 threads on the workstation that cost 1.10 to 1.15 of the
//!   time (NCMEM2 p1).
//!
//! Why the records are the same: the reader is a function of the bytes it reads and of its
//! configuration (the same `ReaderConfig`; the local-ID counter starts at 1 in a new reader and
//! moves one per record read, so the second reading numbers the records as the first). The second
//! reading starts at the same byte, with the same calls per record (`End()`, then
//! `GetNextSequence`), so its line reader refills its buffers at the same bytes (the in-avail
//! window of a regular file depends only on the file's size and the position) and reads the same
//! lines. The first reading's `IsIStreamEmpty` leaves no state but the position (it seeks back and
//! drops its buffered bytes), which the second source starts at. It reads only the records the
//! first reading read without an error, so it never meets the error or the `eEOF`. In parts, each
//! part's reader starts at a record start with the line number and the local-ID counter that the
//! one reader has there (`FastaInputSource::x_from_reader_at`), the only state a reader carries
//! from record to record in a file without `\r` (`x_query_scan_par.rs`, "Why the lengths are the
//! first reading's"), so it reads the same records, the line numbers of their messages included.
//! Each record's length is compared with the first reading's: a record whose length changed, a file
//! that ends early or a reader error (the file changed between the two readings) ends the search
//! with LOSAT's error, and every later reading of the search returns it (`State::failed`); a change
//! that keeps every record's length is read as the new content (as NCBI, which reads per batch,
//! would read it). In parts the change must also keep where each part starts and ends and its
//! lines, or the reading of that part is the error; the parts are read whole, so an error found
//! ahead of the loop's batch waits for the batch that needs its record (`State::pending`), as the
//! one reader finds it only there. The batches are cut from the first reading's lengths by the same
//! code as before (`next_ranged_batch`), so they have the same boundaries, the same records in the
//! same order, the same cut letters (`FastaRecord::cut` at the placement of `query_input_of`), the
//! same lower case letters (the `-lcase_masking` masks), the same local ids and titles (every
//! outfmt's ids and titles are made from them), and the same reader warnings per record (written
//! per batch, as before). DUST and the rest of the search read only the batch's records.
//!
//! Where the switch does not act (`Off`; the records are read before the search as without it):
//! where the reports are written after the loop (LOSAT_X_RPTSTREAM off, or one of its reasons that
//! the search shows before the reading: several formats or a file sink, verbose or debug output,
//! BLAST XML, `x_rpt_stream::refused_before_reading`; every record is
//! then read by a writer after the search, so the records would all be held anyway and the second
//! reading would only cost time), standard input (`-query -`, also with `-subject -`), a query
//! that is not a regular file (a named pipe, a character device such as `/dev/stdin`, a socket:
//! `metadata().is_file()`, which follows symbolic links), and a file whose position or second
//! handle cannot be had (`Io`, for example `try_clone` on WASI). LOSAT reads no compressed query
//! file. `run_local` (the Web adapter and ABI v1) gets its records from its caller.
//!
//! The rule (`Rule::off`): the switch does not act on a regular file shorter than twice the largest
//! query batch (`Batch`, `Rule::batch_letters`: `BATCH_SIZE`, or the `CBatchSizeMixer` maximum,
//! 999,000 letters for blastn and 4,999,000 for megablast and dc-megablast): the batch in flight
//! holds most of the file anyway. The file's bytes are compared with the batch's letters (titles
//! and line ends make the file a little longer than its letters).
//!
//! Until 2026-10-28 the rule had a third reason, `slots`, for the query-batch slots of
//! LOSAT_X_QBATCHPAR, which were removed then. Until 2026-11-03 it had `megablast`: megablast
//! searches are fast for their query bytes, so a reading on one thread while the search's threads
//! wait is several % of their time (v4-D sessions 4 and 5: 31,168 cDNAs against 3 Mb at 4 and 16
//! threads 1.052 and 1.037 of the time); it went when the second reading was read in parts.
//!
//! Elsewhere the second reading saved 3-77 MiB at 0.994-1.006 of the time (blastn and
//! dc-megablast).
//!
//! What is kept for every query (`query_input_of`): the length of each input record, whether it is
//! skipped, the input position and placement of each searched query and its searched length
//! (about 33 bytes per record without `-query_loc`). The records of every batch are kept as well
//! (cloned when the loop takes them) where a writer reads them after their batch: LOSAT_X_RPTSTREAM
//! on but not acting for the search for a reason it finds after the reading (one of the refusals
//! of `x_rpt_stream::check` that read the queries' lengths or the subjects); there the memory is
//! that of the records read before the search. With LOSAT_X_RPTSTREAM acting they are not kept.
//!
//! With `LOSAT_X_STATS` and the switch, `[X_STATS] QSTREAM records= read= batches= kept= scan_s=
//! scan_par= read_s= peak_window= off=` (records: the input records; read: the records read again;
//! off: why the switch did not act, `none` when it did; kept: the searched records the loop kept for
//! a writer after their batch, 0 where LOSAT_X_RPTSTREAM writes each batch's reports; scan_s: the
//! first reading; scan_par: 1 where it was read in parts on the search's threads,
//! `x_query_scan_par.rs`, and the second reading reads those parts, else 0; read_s: the second, on
//! the loop's clock; peak_window: the most records read ahead of the loop at once). The line is
//! printed for searches of the command line only.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::{Duration, Instant};

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::{
    FastaInputSource, FastaRecord, ParseErrorCode, QueryEnd, ReadError, ReaderConfig,
};
use crate::blastinput::seq_range::{
    record_interval, Placements, QueryInput, RecordInterval, RecordPlacement, SequenceRange,
};

// No NCBI counterpart: reads the switch once (LOSAT_X_QSTREAM is in the strict set). It chooses
// when LOSAT reads the same records.
/// Whether LOSAT_X_QSTREAM is on.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_QSTREAM"))
}

/// Why the records of a search are read before the search (`Opening::new`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The switch is off (`LOSAT_X_QSTREAM=0`).
    NotSet,
    /// The reports are written after the loop (LOSAT_X_RPTSTREAM off, or refusing the search before
    /// the reading: `x_rpt_stream::refused_before_reading`), which reads every record.
    Reports,
    /// `-query -` (standard input), which cannot be read twice.
    StandardInput,
    /// A query that is not a regular file (a named pipe, a device, a socket).
    NotAFile,
    /// The file's position or a second handle on it could not be had.
    Io,
    /// The file is shorter than twice the largest query batch (`Rule::off`).
    Batch,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Reports => "reports",
            Off::StandardInput => "standard_input",
            Off::NotAFile => "not_a_file",
            Off::Io => "io",
            Off::Batch => "batch",
        }
    }
}

/// Whether the query file of a search can be read a second time: a second handle on the open
/// file and the position its reading starts at.
pub(super) struct Opening {
    on: bool,
    off: Option<Off>,
    again: Option<(File, u64)>,
    /// Where the switch does not act on a regular file: a second handle on it and the position its
    /// reading starts at, for LOSAT_X_BNQPARTS's reading in parts (`store_parts`).
    before: Option<(File, u64)>,
}

impl Opening {
    /// `path` is `-query` and `file` the file `run` opened for it, before its reader is made;
    /// `shared_standard_input` is `-query -` with `-subject -`; `reports_streamed` is whether
    /// LOSAT_X_RPTSTREAM may write the reports batch by batch (on, and no reason of
    /// `x_rpt_stream::refused_before_reading`); `rule` is what the rule reads of the search
    /// (`Rule::off`).
    pub(super) fn new(
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
        rule: Rule,
    ) -> Opening {
        Opening::with_mode(
            on(),
            path,
            file,
            shared_standard_input,
            reports_streamed,
            rule,
        )
    }

    /// `new` with the switch's value `on` (`on()`; the unit tests set it).
    fn with_mode(
        on: bool,
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
        rule: Rule,
    ) -> Opening {
        let off_with = |off, before| {
            if on {
                crate::utils::xreach::x_reach!("LOSAT_X_QSTREAM", false);
            }
            Opening {
                on,
                off: Some(off),
                again: None,
                before,
            }
        };
        let off = |off| off_with(off, None);
        let standard_input = path.as_os_str() == "-" || shared_standard_input;
        let regular = || match file.metadata() {
            Ok(metadata)
                if metadata.is_file()
                    && crate::blastinput::input_files::has_stream_position(file) =>
            {
                Ok(metadata.len())
            }
            Ok(_) => Err(Off::NotAFile),
            Err(_) => Err(Off::Io),
        };
        // A second handle on the file and the position its reading starts at.
        let again = || {
            let mut handle = file;
            let start = handle.stream_position().ok()?;
            Some((file.try_clone().ok()?, start))
        };
        if !on || !reports_streamed {
            let before = (!standard_input && regular().is_ok()).then(again).flatten();
            return off_with(if on { Off::Reports } else { Off::NotSet }, before);
        }
        if standard_input {
            return off(Off::StandardInput);
        }
        let bytes = match regular() {
            Ok(bytes) => bytes,
            Err(why) => return off(why),
        };
        if let Some(why) = rule.off(bytes) {
            return off_with(why, again());
        }
        match again() {
            Some(again) => {
                crate::utils::xreach::x_reach!("LOSAT_X_QSTREAM", true);
                Opening {
                    on,
                    off: None,
                    again: Some(again),
                    before: None,
                }
            }
            None => off(Off::Io),
        }
    }

    /// EXPERIMENT (LOSAT_X_BNQPARTS): the parts of the store of a search whose records are read
    /// before it (the switch does not act on its regular file), read on the search's `threads`
    /// threads (`x_query_scan_par::store_parts`), or `None` where the search's reader reads them
    /// (one thread, standard input or a file that is not regular, the BLAST archive,
    /// `-parse_deflines`, whose reader numbers the records without ids across the file and adds
    /// every record to the search's scope of ids, and the cases of that module).
    pub(super) fn store_parts(
        &mut self,
        config: ReaderConfig,
        threads: usize,
    ) -> Option<Vec<super::x_query_drop::Part>> {
        if config.keep_defline || config.parse_deflines {
            return None;
        }
        let (file, start) = self.before.take()?;
        super::x_query_scan_par::store_parts(&file, start, config, threads)
    }

    /// Whether the records are read batch by batch (`Store::read`).
    pub(super) fn acts(&self) -> bool {
        self.again.is_some()
    }
}

/// No NCBI counterpart (when LOSAT reads the records only): what the rule of LOSAT_X_QSTREAM
/// reads of a search (`Rule::off`).
#[derive(Clone, Copy, Debug)]
pub(super) struct Rule {
    /// The letters of the largest query batch of the search: `BATCH_SIZE` when it is set, else
    /// the `CBatchSizeMixer` maximum of the task (`run.rs`, `x_query_batch_letters`).
    pub(super) batch_letters: u64,
}

impl Rule {
    /// Why the switch does not act on a regular query file of `file_bytes` bytes (the module's
    /// rule), or `None` where it acts: a file shorter than twice the largest query batch.
    fn off(self, file_bytes: u64) -> Option<Off> {
        if file_bytes < self.batch_letters.saturating_mul(2) {
            return Some(Off::Batch);
        }
        None
    }
}

/// The first reading: the loop of `read_queries` (`GetNextSeqBatch` until `End()`), each record's
/// length kept and the record dropped, and how the reading ended (with the warnings of the record
/// whose reading failed, as `read_queries` keeps them). With `titles` (the encoding of the BLAST
/// archive's or the search strategy's requests, `report::archive::request_title_encoding`), also
/// whether a record has a title that they cannot write (`check_request_titles`), which the batch loop
/// refuses in the batch NCBI completes with it (`run.rs` `run_in_pool`) and which keeps the
/// reports after the loop (`x_rpt_stream::Off::ArchiveTitle`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-152
/// ```c++
///     while (size_read < GetBatchSize()) {
///
///         if (End())
///             break;
///
///         CRef<CBlastSearchQuery> q;
///         try { q.Reset(m_Source->GetNextSequence(scope)); }
///         catch (const CObjReaderParseException& e) {
///             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
///                 break;
///             }
///             throw;
///         }
/// ```
fn scan<R: Read>(
    source: &mut FastaInputSource<R>,
    titles: Option<crate::report::archive::ArchiveEncoding>,
) -> (Vec<usize>, QueryEnd, bool) {
    let mut lengths = Vec::new();
    let mut archive_title = false;
    loop {
        if source.end() {
            return (lengths, QueryEnd::Input, archive_title);
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => {
                archive_title |=
                    titles.is_some_and(|encoding| !encoding.writes_text(&record.title));
                lengths.push(record.sequence.len());
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => return (lengths, QueryEnd::EofError, archive_title),
            Err(error) => return (lengths, QueryEnd::Error { error, warnings }, archive_title),
        }
    }
}

/// The records of a search's queries read batch by batch, or why they are not (its drop prints
/// the summary lines).
pub(super) struct Store {
    on: bool,
    off: Option<Off>,
    acting: Option<Acting>,
}

struct Acting {
    query_input: QueryInput,
    /// The searched length of each searched query.
    lengths: Vec<usize>,
    /// The BLAST archive (outfmt 11) or the search strategy: a record has a title that their
    /// requests cannot write (`scan`).
    archive_title: bool,
    ranged: bool,
    scan: Duration,
    /// The first reading took the parts' lengths (`x_query_scan_par.rs`).
    scan_par: bool,
    state: Mutex<State>,
}

/// The second reading.
struct State {
    second: Second,
    /// The input records read again.
    next: usize,
    /// The input records `base..next`, read and not taken by the loop, in chunks of consecutive
    /// records (first input record, records).
    base: usize,
    chunks: VecDeque<(usize, Arc<Vec<FastaRecord>>)>,
    /// The loop's batch: its first input record and its input records (up to `base`).
    current: Option<(usize, Arc<Vec<FastaRecord>>)>,
    /// The error of a failed reading: every later reading of the search returns it (the records
    /// of the failed reading are lost, and the reader is past the record that failed).
    failed: Option<String>,
    /// A change of the file that the reading in parts found ahead of the loop's batch (it reads
    /// whole parts): the first input record it cannot give and the error, which the reading of a
    /// batch that needs that record returns (then `failed`).
    pending: Option<(usize, String)>,
    stats: Stats,
}

/// How the second reading reads the records.
#[cfg_attr(
    all(unix, feature = "parallel"),
    expect(
        clippy::large_enum_variant,
        reason = "boxing a variant changes the layout and the memory; that needs an RSS and time measurement"
    )
)]
enum Second {
    /// One reader over the file, record after record (`State::read_one`), where the first reading
    /// was not read in parts.
    Serial(FastaInputSource<File>),
    /// The parts of the first reading (`x_query_scan_par::Span`), each read again by a reader of its
    /// own on the search's threads (`State::read_parts`): `next` is the first part not read.
    #[cfg(all(unix, feature = "parallel"))]
    Parts {
        file: File,
        spans: Vec<super::x_query_scan_par::Span>,
        next: usize,
        config: ReaderConfig,
    },
}

#[derive(Default)]
struct Stats {
    batches: u64,
    /// The searched records the loop kept for a writer after their batch (`Queries::note_kept`).
    kept: u64,
    read: Duration,
    peak_window: usize,
}

impl Store {
    /// The store of a search whose records are read before it (`opening` does not act).
    pub(super) fn off(opening: Opening) -> Store {
        debug_assert!(!opening.acts());
        Store {
            on: opening.on,
            off: opening.off,
            acting: None,
        }
    }

    /// The first reading, on `source` (the reader `search_cli` made of the file, after its
    /// `IsIStreamEmpty`), where `read_queries` reads without the switch: each record's length is
    /// kept and the record dropped. Returns how the reading ended. `config` is the source's
    /// configuration, `range` `-query_loc` and `threads` the search's threads (the reading in parts
    /// on them, `x_query_scan_par.rs`).
    pub(super) fn read(
        opening: Opening,
        source: &mut FastaInputSource<File>,
        config: ReaderConfig,
        range: Option<&SequenceRange>,
        threads: usize,
        titles: Option<crate::report::archive::ArchiveEncoding>,
    ) -> Result<(Store, QueryEnd)> {
        Store::read_with(opening, source, config, range, titles, |file, start| {
            super::x_query_scan_par::lengths(file, start, config, threads)
        })
    }

    /// `read` with the first reading in parts `in_parts` (from the second handle and the position
    /// the source starts from; the tests cut small parts).
    fn read_with(
        opening: Opening,
        source: &mut FastaInputSource<File>,
        config: ReaderConfig,
        range: Option<&SequenceRange>,
        titles: Option<crate::report::archive::ArchiveEncoding>,
        in_parts: impl FnOnce(&File, u64) -> Option<super::x_query_scan_par::Scan>,
    ) -> Result<(Store, QueryEnd)> {
        let (again, start) = opening
            .again
            .expect("LOSAT_X_QSTREAM reads a file it can read again");
        let started = Instant::now();
        // The BLAST archive (outfmt 11): this reading sees every record's title (`scan`).
        // The lengths read in parts on the search's threads, from the second handle at the position
        // the source starts from; where the parts are not taken, the source reads them (its
        // position is where it was: the parts read at offsets).
        // A part's reader numbers its records from the part's first (`-parse_deflines`, whose
        // records with ids take no number and whose reader adds every record to the search's
        // scope of ids, is never read batch by batch: `run.rs`).
        let parts = if !config.keep_defline && !config.parse_deflines {
            in_parts(&again, start)
        } else {
            None
        };
        let scan_par = parts.is_some();
        // (The parts are not read for the archive, which keeps the deflines.)
        let (input_lengths, spans, end, archive_title) = match parts {
            Some(parts) => (parts.lengths, Some(parts.spans), QueryEnd::Input, false),
            None => {
                let (lengths, end, archive_title) = scan(source, titles);
                (lengths, None, end, archive_title)
            }
        };
        let scan = started.elapsed();
        let (query_input, lengths) = query_input_of(input_lengths, range);
        let reopen = |handle: &File| -> Result<File> {
            let mut handle = handle
                .try_clone()
                .context("LOSAT_X_QSTREAM: the query file cannot be read again")?;
            handle
                .seek(SeekFrom::Start(start))
                .context("LOSAT_X_QSTREAM: the query file cannot be read again")?;
            Ok(handle)
        };
        // The second reading reads the parts of the first where it read in parts, else the file with
        // one reader.
        let second = match spans {
            #[cfg(all(unix, feature = "parallel"))]
            Some(spans) => Second::Parts {
                file: again,
                spans,
                next: 0,
                config,
            },
            _ => Second::Serial(FastaInputSource::from_file(reopen(&again)?, config)),
        };
        let store = Store {
            on: opening.on,
            off: None,
            acting: Some(Acting {
                ranged: range.is_some(),
                query_input,
                lengths,
                archive_title,
                scan,
                scan_par,
                state: Mutex::new(State {
                    second,
                    next: 0,
                    base: 0,
                    chunks: VecDeque::new(),
                    current: None,
                    failed: None,
                    pending: None,
                    stats: Stats::default(),
                }),
            }),
        };
        Ok((store, end))
    }

    /// Whether the records are read batch by batch.
    pub(super) fn acts(&self) -> bool {
        self.acting.is_some()
    }

    fn acting(&self) -> &Acting {
        self.acting
            .as_ref()
            .expect("LOSAT_X_QSTREAM reads the records of this search batch by batch")
    }

    /// The query input of the search (`query_input_of`).
    pub(super) fn query_input(&self) -> &QueryInput {
        &self.acting().query_input
    }
}

// No NCBI counterpart: the summary line of the switch (with LOSAT_X_STATS); it prints counters and
// changes no value NCBI computes.
impl Drop for Store {
    fn drop(&mut self) {
        // After a panic the line would claim what the run did not finish.
        if !self.on || std::thread::panicking() {
            return;
        }
        let off = self.off.map_or("none", Off::name);
        let (records, scan, scan_par, stats) = match &self.acting {
            Some(acting) => {
                let mut state = acting.lock();
                (
                    acting.query_input.input_lengths.len(),
                    acting.scan,
                    acting.scan_par,
                    (state.next, std::mem::take(&mut state.stats)),
                )
            }
            None => (0, Duration::ZERO, false, (0, Stats::default())),
        };
        let (read, stats) = stats;
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!(
                "[X_STATS] QSTREAM records={records} read={read} batches={} kept={} scan_s={:.4} scan_par={} read_s={:.4} peak_window={} off={off}",
                stats.batches,
                stats.kept,
                scan.as_secs_f64(),
                u8::from(scan_par),
                stats.read.as_secs_f64(),
                stats.peak_window,
            );
        }
    }
}

impl Acting {
    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Searched query `index` as the search reads it, from its input record: cut to its
    /// placement with `-query_loc` (as `cut_queries` cuts it), else the record itself.
    fn searched(&self, record: &FastaRecord, index: usize) -> FastaRecord {
        if self.ranged {
            let from = self.query_input.placements.offset(index);
            record.cut(from, from + self.lengths[index])
        } else {
            record.clone()
        }
    }
}

/// The error of a reading that found the file changed at input record `index`.
fn changed(index: usize, what: &str) -> String {
    format!(
        "LOSAT_X_QSTREAM: the query file changed while LOSAT read it (record {}: {what})",
        index + 1
    )
}

impl State {
    /// Reads the input records up to `end` at least (`input_lengths`: the first reading's lengths):
    /// with one reader exactly up to `end`, into one new chunk; in parts, the parts up to the one
    /// that holds input record `end - 1`, one chunk each, on `pool`'s threads.
    fn read_to(
        &mut self,
        end: usize,
        input_lengths: &[usize],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<()> {
        if let Some(failed) = &self.failed {
            return Err(anyhow::anyhow!("{failed}"));
        }
        if let Some((index, message)) = &self.pending {
            if end > *index {
                self.failed = Some(message.clone());
                return Err(anyhow::anyhow!("{message}"));
            }
        }
        if self.next >= end {
            return Ok(());
        }
        let started = Instant::now();
        if matches!(self.second, Second::Serial(_)) {
            let first = self.next;
            let mut records = Vec::with_capacity(end - first);
            while self.next < end {
                records.push(self.read_one(input_lengths)?);
                self.next += 1;
            }
            self.chunks.push_back((first, Arc::new(records)));
        } else {
            self.read_parts(end, input_lengths, pool);
        }
        self.stats.read += started.elapsed();
        self.stats.peak_window = self.stats.peak_window.max(self.next - self.base);
        if let Some((index, message)) = &self.pending {
            if end > *index {
                self.failed = Some(message.clone());
                return Err(anyhow::anyhow!("{message}"));
            }
        }
        Ok(())
    }

    /// `read_to` in parts: the parts from the first not read up to the one that holds input record
    /// `end - 1`, and with a pool at least one part per thread, each read by a reader of its own
    /// (`x_query_scan_par::read_again`) on `pool`'s threads, their records moved into the chunks in
    /// order up to the first record a part cannot give (the file changed), which is kept `pending`.
    /// No NCBI counterpart: NCBI reads its query batches with one reader (blast_input.cpp:134-171,
    /// the module doc); the parts give the same records (scheduling only).
    #[cfg(all(unix, feature = "parallel"))]
    fn read_parts(
        &mut self,
        end: usize,
        input_lengths: &[usize],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) {
        let Second::Parts {
            file,
            spans,
            next,
            config,
        } = &mut self.second
        else {
            unreachable!("LOSAT_X_QSTREAM reads its parts")
        };
        let from = *next;
        let mut to = from;
        while to < spans.len() && spans[to].first < end {
            to += 1;
        }
        // A batch smaller than a part per thread (the size mixer's first batches, a search with
        // many hits per letter) would leave threads idle, or read one part at a time on this
        // thread: the parts read ahead wait in the window, at most a part per thread.
        let parallel = pool.enabled() && pool.threads() > 1;
        if parallel {
            to = to.max(from + pool.threads()).min(spans.len());
        }
        *next = to;
        let (file, spans, config) = (&*file, &spans[..], *config);
        let read = |index: usize| {
            super::x_query_scan_par::read_again(
                file,
                &spans[index],
                index + 1 == spans.len(),
                config,
                input_lengths,
            )
        };
        let parts: Vec<_> = if parallel && to - from > 1 {
            use rayon::prelude::*;
            pool.install(|| {
                (from..to)
                    .into_par_iter()
                    .with_max_len(1)
                    .map(read)
                    .collect()
            })
        } else {
            (from..to).map(read).collect()
        };
        for (records, failure) in parts {
            if !records.is_empty() {
                let first = self.next;
                self.next += records.len();
                self.chunks.push_back((first, Arc::new(records)));
            }
            if let Some((index, what)) = failure {
                debug_assert_eq!(index, self.next);
                self.pending = Some((index, changed(index, &what)));
                return;
            }
        }
    }

    #[cfg(not(all(unix, feature = "parallel")))]
    fn read_parts(
        &mut self,
        _end: usize,
        _input_lengths: &[usize],
        _pool: &crate::utils::threading::SearchPool<'_>,
    ) {
        unreachable!("LOSAT_X_QSTREAM reads in parts on unix with the parallel feature only")
    }

    /// A chunk boundary at input record `at` (`base <= at <= next`): the chunk across it is cut in
    /// two, its records moved (copied while shared).
    fn split(&mut self, at: usize) {
        let Some(position) = self
            .chunks
            .iter()
            .position(|(start, records)| *start < at && at < start + records.len())
        else {
            return;
        };
        let (start, records) = self.chunks.remove(position).expect("the chunk");
        let (front, back) = match Arc::try_unwrap(records) {
            Ok(mut owned) => {
                let back = owned.split_off(at - start);
                (owned, back)
            }
            Err(shared) => (shared[..at - start].to_vec(), shared[at - start..].to_vec()),
        };
        self.chunks.insert(position, (at, Arc::new(back)));
        self.chunks.insert(position, (start, Arc::new(front)));
    }

    /// The loop's batch of the input records `input` (`base == input.start`), taken out of the
    /// chunks: one chunk as it is, or the records of several (moved, copied while shared).
    fn take(
        &mut self,
        input: Range<usize>,
        input_lengths: &[usize],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Arc<Vec<FastaRecord>>> {
        self.read_to(input.end, input_lengths, pool)?;
        self.split(input.end);
        let mut taken = Vec::new();
        while self
            .chunks
            .front()
            .is_some_and(|(start, _)| *start < input.end)
        {
            taken.push(self.chunks.pop_front().expect("a chunk").1);
        }
        self.base = input.end;
        let records = if taken.len() == 1 {
            taken.pop().expect("a chunk")
        } else {
            let mut records = Vec::with_capacity(input.len());
            for chunk in taken {
                match Arc::try_unwrap(chunk) {
                    Ok(owned) => records.extend(owned),
                    Err(shared) => records.extend_from_slice(&shared),
                }
            }
            Arc::new(records)
        };
        anyhow::ensure!(
            records.len() == input.len(),
            "LOSAT_X_QSTREAM: the batch of the input records {input:?} holds {} records",
            records.len()
        );
        Ok(records)
    }

    /// The next input record, read as the first reading read it (`End()`, then
    /// `GetNextSequence`), with the length the first reading found. A record of another length (or
    /// a file that ends before it, or a reader error) fails this reading and every later one
    /// (`failed`).
    fn read_one(&mut self, input_lengths: &[usize]) -> Result<FastaRecord> {
        if let Some(failed) = &self.failed {
            return Err(anyhow::anyhow!("{failed}"));
        }
        let index = self.next;
        #[cfg_attr(
            not(all(unix, feature = "parallel")),
            expect(
                clippy::infallible_destructuring_match,
                reason = "Second has a second variant (Parts) on unix with the parallel feature"
            )
        )]
        let source = match &mut self.second {
            Second::Serial(source) => source,
            #[cfg(all(unix, feature = "parallel"))]
            Second::Parts { .. } => {
                unreachable!("LOSAT_X_QSTREAM reads record after record with one reader")
            }
        };
        let what = if source.end() {
            "the file ends before it".to_string()
        } else {
            match source.next_sequence(&mut |_| Ok(())) {
                Ok(record) if record.sequence.len() == input_lengths[index] => return Ok(record),
                Ok(record) => format!(
                    "{} letters where the first reading read {}",
                    record.sequence.len(),
                    input_lengths[index]
                ),
                Err(error) => format!("{error}"),
            }
        };
        let message = changed(index, &what);
        self.failed = Some(message.clone());
        Err(anyhow::anyhow!(message))
    }

    /// Writes the reader's warnings of every input record not read yet, in input order (no batch
    /// taken): record after record with one reader, part after part on this thread in parts.
    fn write_warnings(
        &mut self,
        input_lengths: &[usize],
        out: &mut (dyn Write + Send),
    ) -> Result<()> {
        debug_assert!(self.chunks.is_empty() && self.current.is_none());
        let serial = crate::utils::threading::SearchPool::x_serial();
        while self.next < input_lengths.len() {
            if matches!(self.second, Second::Serial(_)) {
                let record = self.read_one(input_lengths)?;
                self.next += 1;
                out.write_all(&record.warnings)?;
            } else {
                let before = self.next;
                self.read_to(self.next + 1, input_lengths, &serial)?;
                debug_assert!(self.next > before, "the parts cover every record");
                for (_, records) in self.chunks.drain(..) {
                    for record in records.iter() {
                        out.write_all(&record.warnings)?;
                    }
                }
            }
            self.base = self.next;
        }
        Ok(())
    }
}

/// The query records of a search: read before it (the off path), or read batch by batch from
/// the query file (`Store`).
#[derive(Clone, Copy)]
pub(super) enum Queries<'a> {
    /// The input records and the searched ones (cut with `-query_loc`, else the input records).
    Records {
        input: &'a [FastaRecord],
        searched: &'a [FastaRecord],
    },
    Streamed(&'a Store),
    /// EXPERIMENT (LOSAT_X_BNQPACK): the records read before the search, which the loop takes out
    /// batch by batch (`x_query_drop.rs`).
    Taken(&'a super::x_query_drop::Store),
}

/// The records of one batch of the loop: its input records (skipped ones included) and its
/// searched records.
pub(super) enum BatchRecords<'a> {
    Borrowed {
        input: &'a [FastaRecord],
        searched: &'a [FastaRecord],
    },
    Read {
        input: Arc<Vec<FastaRecord>>,
        /// The searched records cut with `-query_loc` (`None`: the input records).
        cut: Option<Vec<FastaRecord>>,
    },
    /// EXPERIMENT (LOSAT_X_BNQPACK): the records taken out of the `Store`, held by the loop alone
    /// (the input and the searched records: no `-query_loc`).
    Taken(Vec<FastaRecord>),
}

impl BatchRecords<'_> {
    pub(super) fn input(&self) -> &[FastaRecord] {
        match self {
            BatchRecords::Borrowed { input, .. } => input,
            BatchRecords::Read { input, .. } => input,
            BatchRecords::Taken(records) => records,
        }
    }

    pub(super) fn searched(&self) -> &[FastaRecord] {
        match self {
            BatchRecords::Borrowed { searched, .. } => searched,
            BatchRecords::Read { input, cut } => cut.as_deref().unwrap_or(input),
            BatchRecords::Taken(records) => records,
        }
    }

    /// The searched records for the batch's search: EXPERIMENT (LOSAT_X_BNQPACK) held by the loop
    /// alone when taken out of the `Store`, else shared.
    pub(super) fn batch_queries(&mut self) -> super::x_query_drop::BatchQueries<'_> {
        match self {
            BatchRecords::Taken(records) => super::x_query_drop::BatchQueries::Held(records),
            other => super::x_query_drop::BatchQueries::Shared(other.searched()),
        }
    }

    /// The searched records for a writer after the batch: those taken out of the `Store`, else a
    /// copy.
    pub(super) fn into_searched(self) -> Vec<FastaRecord> {
        match self {
            BatchRecords::Taken(records) => records,
            other => other.searched().to_vec(),
        }
    }
}

impl<'a> Queries<'a> {
    /// Whether a batch's records leave the search with the batch (read batch by batch, or taken
    /// out of a LOSAT_X_BNQPACK `Store`): a writer after the batch reads the loop's kept records.
    pub(super) fn per_batch(&self) -> bool {
        matches!(self, Queries::Streamed(_) | Queries::Taken(_))
    }

    /// The number of input records.
    pub(super) fn input_count(&self) -> usize {
        match self {
            Queries::Records { input, .. } => input.len(),
            Queries::Streamed(store) => store.acting().query_input.input_lengths.len(),
            Queries::Taken(store) => store.lengths().len(),
        }
    }

    /// The searched length of each searched query.
    pub(super) fn searched_lengths(&self) -> Vec<usize> {
        match self {
            Queries::Records { searched, .. } => {
                searched.iter().map(|record| record.seq().len()).collect()
            }
            Queries::Streamed(store) => store.acting().lengths.clone(),
            Queries::Taken(store) => store.lengths().to_vec(),
        }
    }

    /// Whether an input record has a title that the BLAST archive (outfmt 11, `archive`) or the
    /// search strategy (`strategy`) cannot write (`check_request_titles`; the records of a
    /// LOSAT_X_BNQPACK `Store` are never theirs, whose queries' letters the requests read:
    /// `x_query_drop::Rule`).
    pub(super) fn archive_title(&self, archive: bool, strategy: bool) -> bool {
        match self {
            Queries::Records { input, .. } => crate::report::archive::check_request_titles(
                input, "query", "BLASTN", archive, strategy,
            )
            .is_err(),
            Queries::Streamed(store) => store.acting().archive_title,
            Queries::Taken(_) => false,
        }
    }

    /// Writes the reader's warnings of every input record, in input order (a search that reads
    /// every batch and searches none).
    pub(super) fn write_input_warnings(&self, out: &mut (dyn Write + Send)) -> Result<()> {
        match self {
            Queries::Records { input, .. } => {
                for record in *input {
                    out.write_all(&record.warnings)?;
                }
            }
            Queries::Taken(store) => store.write_warnings(out)?,
            Queries::Streamed(store) => {
                let acting = store.acting();
                acting
                    .lock()
                    .write_warnings(&acting.query_input.input_lengths, out)?;
            }
        }
        Ok(())
    }

    /// The records of the loop's batch of input records `input` and searched queries `searched`.
    /// The loop takes its batches in order, each starting where the one before ended.
    #[cfg(test)]
    pub(super) fn batch(
        &self,
        input: Range<usize>,
        searched: Range<usize>,
    ) -> Result<BatchRecords<'a>> {
        self.batch_in(
            input,
            searched,
            &crate::utils::threading::SearchPool::x_serial(),
        )
    }

    /// `batch`, the packed letters of a batch taken out of a LOSAT_X_BNQPACK `Store` made again by
    /// `pool`.
    pub(super) fn batch_in(
        &self,
        input: Range<usize>,
        searched: Range<usize>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<BatchRecords<'a>> {
        match *self {
            Queries::Records {
                input: input_records,
                searched: searched_records,
            } => Ok(BatchRecords::Borrowed {
                input: &input_records[input],
                searched: &searched_records[searched],
            }),
            Queries::Taken(store) => {
                debug_assert_eq!(input, searched, "no -query_loc");
                Ok(BatchRecords::Taken(store.take(input, pool)?))
            }
            Queries::Streamed(store) => {
                let acting = store.acting();
                let mut state = acting.lock();
                anyhow::ensure!(
                    input.start == state.base,
                    "LOSAT_X_QSTREAM: the batch of the input records {input:?} after those up to {}",
                    state.base
                );
                let records = state.take(input.clone(), &acting.query_input.input_lengths, pool)?;
                state.current = Some((input.start, Arc::clone(&records)));
                state.stats.batches += 1;
                let cut = acting.ranged.then(|| {
                    searched
                        .clone()
                        .map(|index| {
                            let ordinal = acting.query_input.ordinal(index);
                            acting.searched(&records[ordinal - input.start], index)
                        })
                        .collect::<Vec<_>>()
                });
                Ok(BatchRecords::Read {
                    input: records,
                    cut,
                })
            }
        }
    }

    /// The loop keeps `count` searched records of its batch for a writer after the batch (the
    /// `kept=` count of the `[X_STATS]` line).
    pub(super) fn note_kept(&self, count: usize) {
        if let Queries::Streamed(store) = self {
            store.acting().lock().stats.kept += count as u64;
        }
    }

    /// The loop is done with its batch that ends at input record `input_end` (its reports are
    /// written, or kept): its records are dropped.
    pub(super) fn done(&self, input_end: usize) {
        if let Queries::Streamed(store) = self {
            let mut state = store.acting().lock();
            if state
                .current
                .as_ref()
                .is_some_and(|(start, records)| start + records.len() == input_end)
            {
                state.current = None;
            }
        }
    }

    /// The searched records of the batches before searched query `start`, for a writer that
    /// writes their reports: the records read before the search, or `kept` (the records of the
    /// batches taken, which the loop keeps where such a writer runs: `run_in_pool`).
    pub(super) fn before<'s>(&'s self, kept: &'s [FastaRecord], start: usize) -> &'s [FastaRecord] {
        match self {
            Queries::Records { searched, .. } => &searched[..start],
            Queries::Streamed(_) | Queries::Taken(_) => {
                assert!(
                    start <= kept.len(),
                    "LOSAT_X_QSTREAM or LOSAT_X_BNQPACK: a writer reads the records of earlier batches, which were not kept"
                );
                &kept[..start]
            }
        }
    }

    /// Every searched record, for the reports after the loop (`before` with every query).
    pub(super) fn all<'s>(&'s self, kept: &'s [FastaRecord]) -> &'s [FastaRecord] {
        match self {
            Queries::Records { searched, .. } => searched,
            Queries::Streamed(_) | Queries::Taken(_) => kept,
        }
    }
}

/// The query input and the searched lengths of a search from its input records' lengths: what
/// `cut_queries` makes of the records with `-query_loc` (`range`), else `QueryInput::whole`,
/// without their letters (the tests below compare them).
fn query_input_of(
    input_lengths: Vec<usize>,
    range: Option<&SequenceRange>,
) -> (QueryInput, Vec<usize>) {
    let count = input_lengths.len();
    let Some(range) = range else {
        let lengths = input_lengths.clone();
        return (
            QueryInput {
                placements: Placements::default(),
                ordinals: (0..count).collect(),
                input_lengths,
                skipped: vec![false; count],
                range_length: None,
            },
            lengths,
        );
    };
    let mut input = QueryInput {
        placements: Placements(Vec::with_capacity(count)),
        ordinals: Vec::with_capacity(count),
        input_lengths: Vec::with_capacity(count),
        skipped: Vec::with_capacity(count),
        range_length: Some(range.length()),
    };
    let mut lengths = Vec::with_capacity(count);
    for (ordinal, length) in input_lengths.into_iter().enumerate() {
        input.input_lengths.push(length);
        let (from, to_exclusive) = match record_interval(range, length) {
            RecordInterval::Letters { from, to } => (from, to + 1),
            RecordInterval::Empty => (length, length),
            RecordInterval::PastEnd => {
                input.skipped.push(true);
                continue;
            }
        };
        input.skipped.push(false);
        lengths.push(to_exclusive - from);
        input.placements.0.push(RecordPlacement {
            offset: from,
            length,
        });
        input.ordinals.push(ordinal);
    }
    (input, lengths)
}

/// The first field in which two query inputs differ.
#[cfg(test)]
fn same_query_input(got: &QueryInput, old: &QueryInput) -> std::result::Result<(), String> {
    let field = if got.placements.0 != old.placements.0 {
        "placements"
    } else if got.ordinals != old.ordinals {
        "ordinals"
    } else if got.input_lengths != old.input_lengths {
        "input lengths"
    } else if got.skipped != old.skipped {
        "skipped records"
    } else if got.range_length != old.range_length {
        "range length"
    } else {
        return Ok(());
    };
    Err(format!(
        "the query input differs from the old one in its {field}"
    ))
}

/// How a reading ended, as text: its kind, and for an error its kind, text, line and the
/// warnings written before it.
#[cfg(test)]
fn end_text(end: &QueryEnd) -> String {
    match end {
        QueryEnd::Input => "the end of the input".to_string(),
        QueryEnd::EofError => "an eEOF".to_string(),
        QueryEnd::Error { error, warnings } => {
            let kind = match error {
                ReadError::Parse { code, line, .. } => format!("parse {code:?} line {line}"),
                ReadError::Unsupported(_) => "unsupported".to_string(),
                ReadError::Write(_) => "write".to_string(),
            };
            format!(
                "the error ({kind}) {error} after the warnings {:?}",
                String::from_utf8_lossy(warnings)
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{read_queries, QueryRecords};
    use crate::blastinput::seq_range::{cut_queries, RangeRole};

    /// The query input made from the lengths is `cut_queries`'s (with ranges that start inside,
    /// at the end of, and past some records) and `QueryInput::whole`'s.
    #[test]
    fn query_input_from_lengths_is_that_of_the_records() {
        let lengths = [0usize, 1, 5, 99, 100, 101, 250, 3, 100, 7, 1000];
        let records: Vec<FastaRecord> = lengths
            .iter()
            .enumerate()
            .map(|(index, &length)| {
                let letters: Vec<u8> = (0..length)
                    .map(|at| b"ACGTacgtN"[(at + index) % 9])
                    .collect();
                FastaRecord::new(format!("Query_{}", index + 1), b"t", &letters)
            })
            .collect();
        let (whole, whole_lengths) = query_input_of(lengths.to_vec(), None);
        same_query_input(&whole, &QueryInput::whole(&records)).unwrap();
        assert_eq!(whole_lengths, lengths);
        for text in [
            "1-5",
            "3-200",
            "100-150",
            "101-5000",
            "6-7",
            "2-3",
            "1000-2000",
        ] {
            let range = crate::blastinput::seq_range::parse_sequence_range(
                text,
                RangeRole::Query,
                "BLASTN",
            )
            .unwrap();
            let cut = cut_queries(&records, &range);
            let (input, searched) = query_input_of(lengths.to_vec(), Some(&range));
            same_query_input(&input, &cut.input).unwrap_or_else(|error| panic!("{text}: {error}"));
            let cut_lengths: Vec<usize> = cut.records.iter().map(|r| r.seq().len()).collect();
            assert_eq!(searched, cut_lengths, "{text}");
            // Each searched record is the input record cut at its placement and length.
            for (index, record) in cut.records.iter().enumerate() {
                let ordinal = input.ordinal(index);
                let from = input.placements.offset(index);
                assert_eq!(
                    &records[ordinal].cut(from, from + searched[index]),
                    record,
                    "{text}"
                );
            }
        }
    }

    fn config() -> ReaderConfig {
        ReaderConfig::query("BLASTN", false, false)
    }

    /// A file of this test in the temporary directory, removed when dropped.
    struct TempFile(std::path::PathBuf);

    impl TempFile {
        fn new(tag: &str, bytes: &[u8]) -> Self {
            static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let path = std::env::temp_dir().join(format!(
                "losat_qstream_{}_{}_{tag}.fa",
                std::process::id(),
                NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
            ));
            std::fs::write(&path, bytes).expect("write the test file");
            TempFile(path)
        }
    }

    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    /// The store of `path` as `search_cli` makes it (`IsIStreamEmpty`, then the first reading)
    /// with the switch on, and how the first reading ended (one thread: the reading in parts,
    /// `x_query_scan_par.rs`, has its own tests).
    fn store(path: &Path, range: Option<&SequenceRange>) -> (Store, QueryEnd) {
        let file = File::open(path).expect("open the test file");
        let opening = Opening::with_mode(true, path, &file, false, true, NO_RULE);
        assert!(opening.acts());
        let mut source = FastaInputSource::from_file(file, config());
        assert!(!source.stream_is_empty());
        Store::read(opening, &mut source, config(), range, 1, None).expect("the first reading")
    }

    /// The records as `search_cli` reads them without the switch.
    fn old(path: &Path) -> QueryRecords {
        let mut source = FastaInputSource::from_file(File::open(path).unwrap(), config());
        assert!(!source.stream_is_empty());
        read_queries(&mut source)
    }

    /// xorshift64*.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }

        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// A FASTA text of `seed`: LF, CR LF or mixed line ends; records of 0 to about 20,000 letters
    /// (one of 1.1 MiB in some texts: the reader's residue buffer is then taken, not copied), in
    /// lines of 60 or 80 letters or one line (longer than the line reader's 8,191-byte window);
    /// lower case runs, `N`s, IUPAC letters, digits (the reader's warning), records without a title
    /// or without letters, blank and comment lines; letters before the first defline in some; and
    /// the end with or without a line end, with blank lines, or with a record that the reader
    /// refuses.
    fn text(seed: u64) -> Vec<u8> {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        let ends = seed % 3;
        let mut text = Vec::new();
        let eol = |text: &mut Vec<u8>, rng: &mut Rng| match ends {
            0 => text.push(b'\n'),
            1 => text.extend_from_slice(b"\r\n"),
            _ => {
                if rng.below(2) == 0 {
                    text.push(b'\n')
                } else {
                    text.extend_from_slice(b"\r\n")
                }
            }
        };
        if seed.is_multiple_of(7) {
            text.extend_from_slice(b"ACGTACGTTT");
            eol(&mut text, &mut rng);
        }
        let records = 1 + rng.below(40);
        for index in 0..records {
            match rng.below(10) {
                0 => {
                    eol(&mut text, &mut rng);
                }
                1 => {
                    text.extend_from_slice(b";a comment line");
                    eol(&mut text, &mut rng);
                }
                _ => {}
            }
            text.push(b'>');
            if rng.below(8) != 0 {
                text.extend_from_slice(format!("rec{index} some words {seed}").as_bytes());
            }
            eol(&mut text, &mut rng);
            let length = if seed.is_multiple_of(5) && index == 2 {
                (1 << 20) + 100_000
            } else {
                match rng.below(20) {
                    0 => 0,
                    1 => 9_000 + rng.below(11_000) as usize,
                    _ => rng.below(3_000) as usize,
                }
            };
            let width = match rng.below(3) {
                0 => 60,
                1 => 80,
                _ => length.max(1),
            };
            let mut letters: Vec<u8> = (0..length)
                .map(|_| b"ACGTACGTACGTNRY"[rng.below(15) as usize])
                .collect();
            if length > 200 && rng.below(3) == 0 {
                let at = rng.below(length as u64 - 150) as usize;
                letters[at..at + 120].make_ascii_lowercase();
            }
            for (line, chunk) in letters.chunks(width).enumerate() {
                text.extend_from_slice(chunk);
                if line == 1 && rng.below(6) == 0 {
                    text.extend_from_slice(b"12");
                }
                eol(&mut text, &mut rng);
            }
        }
        match seed % 4 {
            0 => {
                // No line end after the last line.
                while matches!(text.last(), Some(b'\n' | b'\r')) {
                    text.pop();
                }
            }
            1 => {
                eol(&mut text, &mut rng);
                eol(&mut text, &mut rng);
            }
            2 => text.extend_from_slice(b">bad\n@@@@\n"),
            _ => {}
        }
        text
    }

    /// `store` with the first reading in small parts on `threads` threads where the file can be cut
    /// (`x_query_scan_par::small_lengths`), so that the second reading reads those parts.
    #[cfg(all(unix, feature = "parallel"))]
    fn store_in_parts(
        path: &Path,
        range: Option<&SequenceRange>,
        threads: usize,
    ) -> (Store, QueryEnd) {
        let file = File::open(path).expect("open the test file");
        let opening = Opening::with_mode(true, path, &file, false, true, NO_RULE);
        assert!(opening.acts());
        let mut source = FastaInputSource::from_file(file, config());
        assert!(!source.stream_is_empty());
        Store::read_with(
            opening,
            &mut source,
            config(),
            range,
            None,
            |file, start| {
                super::super::x_query_scan_par::small_lengths(file, start, config(), threads)
            },
        )
        .expect("the first reading")
    }

    /// The batches of the second reading in `sizes` records (cycled): every record equals the
    /// first's.
    fn read_in_batches(
        path: &Path,
        old: &QueryRecords,
        range: Option<&SequenceRange>,
        sizes: &[usize],
    ) {
        let (store, end) = store(path, range);
        let pool = crate::utils::threading::SearchPool::x_serial();
        check_batches(&store, &end, path, old, range, sizes, &pool);
    }

    /// The batches of `store`'s second reading in `sizes` records (cycled), taken with `pool`:
    /// every record equals the first's.
    fn check_batches(
        store: &Store,
        end: &QueryEnd,
        path: &Path,
        old: &QueryRecords,
        range: Option<&SequenceRange>,
        sizes: &[usize],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) {
        assert_eq!(end_text(end), end_text(&old.end), "{path:?}");
        let cut = range.map(|range| cut_queries(&old.records, range));
        let searched_old: &[FastaRecord] = cut.as_ref().map_or(&old.records, |cut| &cut.records);
        let queries = Queries::Streamed(store);
        let input = store.query_input();
        let count = input.input_lengths.len();
        assert_eq!(count, old.records.len());
        assert_eq!(queries.input_count(), count);
        let lengths: Vec<usize> = searched_old.iter().map(|r| r.seq().len()).collect();
        assert_eq!(queries.searched_lengths(), lengths);
        let mut start = 0;
        let mut turn = 0;
        while start < count {
            let end = (start + sizes[turn % sizes.len()]).min(count);
            turn += 1;
            let searched = input.searched_in(start..end);
            let batch = queries
                .batch_in(start..end, searched.clone(), pool)
                .unwrap();
            assert_eq!(batch.input(), &old.records[start..end], "{path:?}");
            assert_eq!(
                batch.searched(),
                &searched_old[searched.clone()],
                "{path:?}"
            );
            drop(batch);
            queries.done(end);
            start = end;
        }
    }

    /// Texts that exercise the reader's states between records (line ends, buffer windows, long
    /// records, warnings, an error at the end): the second reading reads the records of the
    /// first, in batches of every size, with and without `-query_loc`.
    #[test]
    fn the_second_reading_reads_the_records_of_the_first() {
        let range = crate::blastinput::seq_range::parse_sequence_range(
            "50-1200",
            RangeRole::Query,
            "BLASTN",
        )
        .unwrap();
        let (mut warned, mut failed, mut long) = (0, 0, 0);
        for seed in 1..=40u64 {
            let file = TempFile::new("records", &text(seed));
            let old = old(&file.0);
            warned += usize::from(old.records.iter().any(|r| !r.warnings.is_empty()));
            failed += usize::from(matches!(old.end, QueryEnd::Error { .. }));
            long += usize::from(old.records.iter().any(|r| r.seq().len() > 1 << 20));
            for sizes in [&[1usize][..], &[3, 1, 7], &[usize::MAX]] {
                read_in_batches(&file.0, &old, None, sizes);
                read_in_batches(&file.0, &old, Some(&range), sizes);
            }
        }
        assert!(
            warned >= 10 && failed >= 5 && long >= 5,
            "{warned} {failed} {long}"
        );
    }

    /// The second reading in parts (the first reading in parts: a file without `\r`, whose reading
    /// ends at its end) reads the records of the first, in batches of every size, on a pool of 2 and
    /// 4 threads, with and without `-query_loc`: their local IDs, titles, letters and messages (the
    /// line numbers of the file); a file with a `\r` is read by one reader.
    #[cfg(all(unix, feature = "parallel"))]
    #[test]
    fn the_second_reading_in_parts_reads_the_records_of_the_first() {
        let range = crate::blastinput::seq_range::parse_sequence_range(
            "50-1200",
            RangeRole::Query,
            "BLASTN",
        )
        .unwrap();
        let (mut in_parts, mut warned) = (0, 0);
        for seed in 1..=60u64 {
            let text = text(seed);
            let file = TempFile::new("parts", &text);
            let old = old(&file.0);
            for threads in [2, 4] {
                for ranged in [None, Some(&range)] {
                    for sizes in [&[1usize][..], &[3, 1, 7], &[usize::MAX]] {
                        let (store, end) = store_in_parts(&file.0, ranged, threads);
                        let parts = store.acting().scan_par;
                        assert_eq!(
                            parts,
                            matches!(store.acting().lock().second, Second::Parts { .. })
                        );
                        if text.contains(&b'\r') {
                            assert!(!parts, "seed {seed}: parts with a `\\r`");
                        }
                        in_parts += usize::from(parts);
                        warned += usize::from(
                            parts && old.records.iter().skip(1).any(|r| !r.warnings.is_empty()),
                        );
                        crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                            check_batches(&store, &end, &file.0, &old, ranged, sizes, pool);
                            Ok(())
                        })
                        .unwrap();
                    }
                }
            }
        }
        assert!(in_parts >= 100 && warned >= 20, "{in_parts} {warned}");
    }

    /// Records of about 1,100 bytes each (a title and ten lines of 100 letters), `count` of them,
    /// the letters of record `index` from `letter`.
    #[cfg(all(unix, feature = "parallel"))]
    fn even_records(count: usize, letter: impl Fn(usize) -> u8) -> Vec<u8> {
        let mut bytes = Vec::new();
        for index in 0..count {
            bytes.extend_from_slice(format!(">r{index} a title\n").as_bytes());
            for _ in 0..10 {
                bytes.extend(std::iter::repeat_n(letter(index), 100));
                bytes.push(b'\n');
            }
        }
        bytes
    }

    /// The batches of one record of the second reading in parts of the file `path` changed to
    /// `changed` after the first reading: the records read up to the first error (in order) and
    /// the error, where the reading failed. A second reading of the failed batch is the same error.
    #[cfg(all(unix, feature = "parallel"))]
    fn read_changed(path: &Path, changed: &[u8]) -> (Vec<FastaRecord>, Option<(usize, String)>) {
        let (store, _) = store_in_parts(path, None, 2);
        assert!(store.acting().scan_par, "the first reading in parts");
        std::fs::write(path, changed).unwrap();
        let queries = Queries::Streamed(&store);
        let count = queries.input_count();
        crate::utils::threading::with_search_pool(2, "blastn", |pool| {
            let mut records = Vec::new();
            for index in 0..count {
                match queries.batch_in(index..index + 1, index..index + 1, pool) {
                    Ok(batch) => {
                        records.extend_from_slice(batch.input());
                        drop(batch);
                        queries.done(index + 1);
                    }
                    Err(error) => {
                        let again = queries
                            .batch_in(index..index + 1, index..index + 1, pool)
                            .err()
                            .expect("the error again");
                        assert_eq!(again.to_string(), error.to_string());
                        return Ok((records, Some((index, error.to_string()))));
                    }
                }
            }
            Ok((records, None))
        })
        .unwrap()
    }

    /// The second reading in parts of a file changed after the first reading: a record of another
    /// length is LOSAT's error at the batch that needs it (the parts are read whole, ahead of the
    /// loop's batch: the error waits for that batch, and the records before it are read); a change
    /// that keeps the bytes of every part and its lines is read as the new content; a change that
    /// moves a part is an error; letters added to the last record are an error, a record added after
    /// it is not read.
    #[cfg(all(unix, feature = "parallel"))]
    #[test]
    fn a_file_changed_between_the_readings_in_parts() {
        let bytes = even_records(40, |_| b'A');
        let file = TempFile::new("changed_parts", &bytes);
        let read = |path: &Path| old(path).records;
        let old_records = read(&file.0);
        // Record 30 one letter shorter.
        let mut shorter = bytes.clone();
        let at = shorter
            .windows(5)
            .position(|window| window == b">r30 ")
            .unwrap();
        let line = at + b">r30 a title\n".len();
        shorter.remove(line);
        let (records, error) = read_changed(&file.0, &shorter);
        let (index, error) = error.expect("an error");
        assert_eq!(index, 30);
        assert!(
            error.contains("changed while LOSAT read it (record 31: 999 letters where the first reading read 1000)"),
            "{error}"
        );
        assert_eq!(records, old_records[..30]);
        // Record 25's letters changed, every byte in its place: the new records.
        std::fs::write(&file.0, &bytes).unwrap();
        let same = even_records(40, |index| if index == 25 { b'C' } else { b'A' });
        let (records, error) = read_changed(&file.0, &same);
        assert_eq!(error, None);
        let file_same = TempFile::new("changed_same", &same);
        assert_eq!(records, read(&file_same.0));
        assert_ne!(records, old_records);
        // A longer title in record 5 moves every later byte: its part, which still ends at the old
        // place, reads records of the old lengths (the new title), and the next part no longer
        // starts at a record start: an error there.
        std::fs::write(&file.0, &bytes).unwrap();
        let mut moved = bytes.clone();
        let at = moved
            .windows(5)
            .position(|window| window == b">r5 a")
            .unwrap();
        moved.insert(at + 4, b'x');
        let (records, error) = read_changed(&file.0, &moved);
        let (index, error) = error.expect("an error");
        assert!(
            index > 5 && error.contains("no longer starts at a record start"),
            "{index} {error}"
        );
        assert_eq!(records[..5], old_records[..5]);
        assert_eq!(records[5].title, b"r5 xa title");
        assert_eq!(records[6..], old_records[6..index]);
        // Letters added to the last record: an error at the last batch.
        std::fs::write(&file.0, &bytes).unwrap();
        let mut longer = bytes.clone();
        longer.extend_from_slice(b"ACGT\n");
        let (records, error) = read_changed(&file.0, &longer);
        let (index, error) = error.expect("an error");
        assert_eq!(index, 39);
        assert!(error.contains("record 40: 1004 letters"), "{error}");
        assert_eq!(records, old_records[..39]);
        // A record added after the last: not read.
        std::fs::write(&file.0, &bytes).unwrap();
        let mut added = bytes.clone();
        added.extend_from_slice(b">r40\nACGT\n");
        let (records, error) = read_changed(&file.0, &added);
        assert_eq!(error, None);
        assert_eq!(records, old_records);
    }

    /// The warnings of every record read again (a search that reads its batches and searches
    /// none) are those of the first reading's records.
    #[test]
    fn the_input_warnings_are_those_of_the_records() {
        for seed in [3u64, 6, 12, 18] {
            let file = TempFile::new("warnings", &text(seed));
            let old = old(&file.0);
            let expected: Vec<u8> = old
                .records
                .iter()
                .flat_map(|record| record.warnings.iter().copied())
                .collect();
            let (store, _) = store(&file.0, None);
            let mut written = Vec::new();
            Queries::Streamed(&store)
                .write_input_warnings(&mut written)
                .unwrap();
            assert_eq!(written, expected, "seed {seed}");
            // The same where the second reading reads the parts of the first (the texts that end
            // with a record the reader refuses are read by one reader).
            #[cfg(all(unix, feature = "parallel"))]
            {
                let (store, _) = store_in_parts(&file.0, None, 4);
                assert_eq!(store.acting().scan_par, seed % 4 != 2, "seed {seed}");
                let mut written = Vec::new();
                Queries::Streamed(&store)
                    .write_input_warnings(&mut written)
                    .unwrap();
                assert_eq!(written, expected, "seed {seed} in parts");
            }
        }
    }

    /// A record whose length changed between the two readings ends the search with LOSAT's error.
    #[test]
    fn a_file_changed_between_the_readings_is_an_error() {
        let file = TempFile::new("changed", b">a\nACGTACGT\n>b\nACGT\n");
        let (store, _) = store(&file.0, None);
        std::fs::write(&file.0, b">a\nACGTACG\n>b\nACGT\n").unwrap();
        let error = Queries::Streamed(&store)
            .batch(0..2, 0..2)
            .err()
            .expect("an error")
            .to_string();
        assert!(error.contains("changed while LOSAT read it"), "{error}");
    }

    /// A read error (a record whose length changed) fails every later reading: the loop's next
    /// call is LOSAT's error again, not a batch with records missing (the records before the
    /// failed one are lost, and the reader is past it).
    #[test]
    fn a_read_error_fails_every_later_reading() {
        let mut bytes = Vec::new();
        for index in 0..6 {
            bytes.extend_from_slice(format!(">r{index}\nACGT\n").as_bytes());
        }
        let file = TempFile::new("read_error", &bytes);
        let (store, _) = store(&file.0, None);
        let mut changed = Vec::new();
        for index in 0..6 {
            let letters = if index == 3 { "ACG" } else { "ACGT" };
            changed.extend_from_slice(format!(">r{index}\n{letters}\n").as_bytes());
        }
        std::fs::write(&file.0, changed).unwrap();
        let queries = Queries::Streamed(&store);
        let first = queries.batch(0..1, 0..1).unwrap();
        assert_eq!(first.input().len(), 1);
        drop(first);
        queries.done(1);
        let failed = queries
            .batch(1..4, 1..4)
            .err()
            .expect("the reading's error")
            .to_string();
        assert!(failed.contains("record 4: 3 letters"), "{failed}");
        let error = queries
            .batch(1..4, 1..4)
            .err()
            .expect("the error again")
            .to_string();
        assert_eq!(error, failed);
    }

    /// The loop's batch is released when the loop is done with it: no chunk of the window and no
    /// `current` holds it, and its records are freed with the loop's last handle.
    #[test]
    fn done_releases_the_loops_batch() {
        let mut bytes = Vec::new();
        for index in 0..6 {
            bytes.extend_from_slice(format!(">r{index}\nACGTACGT\n").as_bytes());
        }
        let file = TempFile::new("release", &bytes);
        let (store, _) = store(&file.0, None);
        let queries = Queries::Streamed(&store);
        let batch = queries.batch(0..3, 0..3).unwrap();
        let held = match &batch {
            BatchRecords::Read { input, .. } => Arc::clone(input),
            BatchRecords::Borrowed { .. } | BatchRecords::Taken(_) => panic!("records read again"),
        };
        // The batch's handle, `current` and this one.
        assert_eq!(Arc::strong_count(&held), 3);
        drop(batch);
        queries.done(3);
        {
            let state = store.acting().lock();
            assert!(state.current.is_none());
            assert!(state.chunks.iter().all(|(start, _)| *start >= 3));
        }
        assert_eq!(Arc::strong_count(&held), 1);
    }

    /// Where the switch does not act: off, the reports after the loop, standard input, a directory
    /// (not a regular file).
    #[test]
    fn inputs_that_cannot_be_read_again() {
        let file = TempFile::new("off", b">a\nACGT\n");
        let handle = File::open(&file.0).unwrap();
        assert_eq!(
            Opening::with_mode(false, &file.0, &handle, false, true, NO_RULE).off,
            Some(Off::NotSet)
        );
        assert_eq!(
            Opening::with_mode(true, &file.0, &handle, false, false, NO_RULE).off,
            Some(Off::Reports)
        );
        assert_eq!(
            Opening::with_mode(true, Path::new("-"), &handle, false, true, NO_RULE).off,
            Some(Off::StandardInput)
        );
        assert_eq!(
            Opening::with_mode(true, &file.0, &handle, true, true, NO_RULE).off,
            Some(Off::StandardInput)
        );
        assert!(Opening::with_mode(true, &file.0, &handle, false, true, NO_RULE).acts());
        // The rule turns the switch off; it reads the file's length. The file is 8 bytes: off for a
        // batch of 5 letters, on for 4.
        for (letters, off) in [(5, Some(Off::Batch)), (4, None)] {
            let rule = Rule {
                batch_letters: letters,
            };
            assert_eq!(
                Opening::with_mode(true, &file.0, &handle, false, true, rule).off,
                off
            );
        }
        #[cfg(unix)]
        {
            let dir = File::open(std::env::temp_dir()).unwrap();
            assert_eq!(
                Opening::with_mode(true, &std::env::temp_dir(), &dir, false, true, NO_RULE).off,
                Some(Off::NotAFile)
            );
        }
    }

    /// The first reading's lengths and end are those of `read_queries`.
    #[test]
    fn the_first_reading_is_that_of_read_queries() {
        for seed in 1..=40u64 {
            let file = TempFile::new("scan", &text(seed));
            let old = old(&file.0);
            let mut source = FastaInputSource::from_file(File::open(&file.0).unwrap(), config());
            assert!(!source.stream_is_empty());
            let (lengths, end, _) = scan(&mut source, None);
            let old_lengths: Vec<usize> = old.records.iter().map(|r| r.seq().len()).collect();
            assert_eq!(lengths, old_lengths, "seed {seed}");
            assert_eq!(end_text(&end), end_text(&old.end), "seed {seed}");
        }
    }

    /// A rule that never turns the switch off.
    const NO_RULE: Rule = Rule { batch_letters: 0 };

    /// The rule: off for a file shorter than twice the largest batch; the files of v4-D's ws windows
    /// (1,775 cDNAs 5,398,904 bytes, 3,492 cDNAs 10,779,147, 5,000 cDNAs 15,440,261, 40 MB
    /// 40,005,015, 31,168 cDNAs 96,191,719).
    #[test]
    fn the_rule_follows_the_batch() {
        const BLASTN: u64 = 999_000;
        const DC: u64 = 4_999_000;
        let rule = |batch_letters| Rule { batch_letters };
        assert_eq!(NO_RULE.off(0), None);
        // The batch: blastn off below 1,998,000 bytes, megablast and dc-megablast below 9,998,000,
        // BATCH_SIZE 4,000 below 8,000.
        assert_eq!(rule(BLASTN).off(1_997_999), Some(Off::Batch));
        assert_eq!(rule(BLASTN).off(1_998_000), None);
        for bytes in [5_398_904, 10_779_147, 15_440_261, 40_005_015, 96_191_719] {
            assert_eq!(rule(BLASTN).off(bytes), None);
        }
        assert_eq!(rule(DC).off(5_398_904), Some(Off::Batch));
        assert_eq!(rule(DC).off(9_997_999), Some(Off::Batch));
        assert_eq!(rule(DC).off(9_998_000), None);
        assert_eq!(rule(DC).off(96_191_719), None);
        assert_eq!(rule(4_000).off(7_999), Some(Off::Batch));
        assert_eq!(rule(4_000).off(8_000), None);
        // The ends saturate.
        assert_eq!(rule(u64::MAX).off(u64::MAX - 1), Some(Off::Batch));
        assert_eq!(rule(u64::MAX).off(u64::MAX), None);
        assert_eq!(rule(0).off(0), None);
    }
}
