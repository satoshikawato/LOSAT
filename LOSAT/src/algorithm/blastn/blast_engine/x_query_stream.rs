//! LOSAT_X_QSTREAM (on by default, strict set; shadow LOSAT_X_QSTREAMSHADOW): blastn reads the
//! records of each query batch from the query file when the batch loop of `run_in_pool` takes the
//! batch (or a LOSAT_X_QBATCHPAR slot searches it ahead of the loop), as NCBI reads its query
//! batches one at a time, instead of holding the records of the whole query file for the whole
//! search. `LOSAT_X_QSTREAM=0` or `LOSAT_X_DEFAULT=0` reads them all before the search, as before.
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
//! mixer, LOSAT_X_QBATCHPAR's plan and rule, LOSAT_X_RPTSTREAM's checks, `-query_loc`'s skipped
//! records) and how the reading ends (an `eEOF`, a reader error, LOSAT's rejection of a Seq-id
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
//!   loop takes its batch's records when it reaches the batch (`Queries::batch`), and a slot copies
//!   those of the batch it searches ahead (`Queries::slot_records`). The records read but not yet
//!   taken by the loop wait in a window; the loop's batch is held until the batch is done (a slot
//!   may still copy from it), then dropped.
//!
//! Why the records are the same: the reader is a function of the bytes it reads and of its
//! configuration (the same `ReaderConfig`; the local-ID counter starts at 1 in a new reader and
//! moves one per record read, so the second reading numbers the records as the first). The second
//! reading starts at the same byte, with the same calls per record (`End()`, then
//! `GetNextSequence`), so its line reader refills its buffers at the same bytes (the in-avail
//! window of a regular file depends only on the file's size and the position) and reads the same
//! lines. The first reading's `IsIStreamEmpty` leaves no state but the position (it seeks back and
//! drops its buffered bytes), which the second source starts at. It reads only the records the
//! first reading read without an error, so it never meets the error or the `eEOF`. Each record's
//! length is compared with the first reading's; a difference (the file changed between the two
//! readings) ends the search with LOSAT's error. The batches are cut from the first reading's
//! lengths by the same code as before (`next_ranged_batch`), so they have the same boundaries, the
//! same records in the same order, the same cut letters (`FastaRecord::cut` at the placement of
//! `query_input_of`), the same lower case letters (the `-lcase_masking` masks), the same local ids
//! and titles (every outfmt's ids and titles are made from them), and the same reader warnings per
//! record (written per batch, as before). DUST and the rest of the search read only the batch's
//! records.
//!
//! Where the switch does not act (`Off`; the records are read before the search as without it):
//! where the reports are written after the loop (LOSAT_X_RPTSTREAM off: every record is then
//! read by a writer after the search, so the records would all be held anyway and the second
//! reading would only cost time), standard input (`-query -`, also with `-subject -`), a query
//! that is not a regular file (a named pipe, a character device such as `/dev/stdin`, a socket:
//! `metadata().is_file()`, which follows symbolic links), and a file whose position or second
//! handle cannot be had (`Io`, for example `try_clone` on WASI). LOSAT reads no compressed query
//! file. `run_local` (the Web adapter and ABI v1) gets its records from its caller.
//!
//! What is kept for every query (`query_input_of`): the length of each input record, whether it is
//! skipped, the input position and placement of each searched query and its searched length
//! (about 33 bytes per record without `-query_loc`). The records of every batch are kept as well
//! (cloned when the loop takes them) where a writer reads them after their batch: LOSAT_X_RPTSTREAM
//! on but not acting for the search (one of its refusals, `x_rpt_stream::check`), or its shadow
//! (which writes the off path's reports to compare); there the memory is that of the records read
//! before the search. With LOSAT_X_RPTSTREAM acting without its shadow they are not kept.
//!
//! The shadow (`LOSAT_X_QSTREAMSHADOW`, implies the switch) also reads the whole file the old way
//! (a third source at the same position: `IsIStreamEmpty`, then `read_queries`) and compares: the
//! first reading's lengths and end (its kind, the error's text, kind and line, and the warnings)
//! and the query input made from the lengths (with `cut_queries` or `QueryInput::whole`); then
//! every batch the loop takes (its input records and its searched records) and every slot's copy
//! with the old records (`FastaRecord` as a whole: local id, title, letters with their case, and
//! the reader's warnings). A difference prints `[X_SHADOW] QSTREAM batches= records= DIFFERENT
//! (...)` and panics; at the end of the search `[X_SHADOW] QSTREAM batches= records= slot_copies=
//! off= identical` (batches: the loop's batches compared; records: the records compared, the
//! loop's input records and every slot's copies; off: why the switch did not act, `none` when it
//! did). With `LOSAT_X_STATS` and the switch, `[X_STATS] QSTREAM records= read= batches=
//! slot_copies= slot_records= stale= scan_s= read_s= peak_window= off=` (records: the input
//! records; read: the records read again; stale: slot requests for records the loop had passed,
//! answered with QBATCHPAR's `Cancelled`; scan_s: the first reading; read_s: the second, summed
//! over the threads that read; peak_window: the most records read ahead of the loop at once).
//! Both lines are printed for searches of the command line only.

use std::collections::VecDeque;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::{Duration, Instant};

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::{
    read_queries, FastaInputSource, FastaRecord, ParseErrorCode, QueryEnd, QueryRecords, ReadError,
    ReaderConfig,
};
use crate::blastinput::seq_range::{
    cut_queries, record_interval, Placements, QueryInput, RecordInterval, RecordPlacement,
    SequenceRange,
};

// No NCBI counterpart: reads the switches once. LOSAT_X_QSTREAM is in the strict set (on unless
// `LOSAT_X_QSTREAM=0` or `LOSAT_X_DEFAULT=0`); the shadow is set to any value and implies the
// switch. They choose when LOSAT reads the same records.
/// 0 = off, 1 = LOSAT_X_QSTREAM, 2 = LOSAT_X_QSTREAMSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_QSTREAMSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_QSTREAM") {
            1
        } else {
            0
        }
    })
}

/// Why the records of a search are read before the search (`Opening::new`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The switch is off (`LOSAT_X_QSTREAM=0`, `LOSAT_X_DEFAULT=0`).
    NotSet,
    /// The reports are written after the loop (LOSAT_X_RPTSTREAM off), which reads every record.
    Reports,
    /// `-query -` (standard input), which cannot be read twice.
    StandardInput,
    /// A query that is not a regular file (a named pipe, a device, a socket).
    NotAFile,
    /// The file's position or a second handle on it could not be had.
    Io,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Reports => "reports",
            Off::StandardInput => "standard_input",
            Off::NotAFile => "not_a_file",
            Off::Io => "io",
        }
    }
}

/// Whether the query file of a search can be read a second time: a second handle on the open
/// file and the position its reading starts at.
pub(super) struct Opening {
    mode: u8,
    off: Option<Off>,
    again: Option<(File, u64)>,
}

impl Opening {
    /// `path` is `-query` and `file` the file `run` opened for it, before its reader is made;
    /// `shared_standard_input` is `-query -` with `-subject -`; `reports_streamed` is whether
    /// LOSAT_X_RPTSTREAM is on (`x_rpt_stream::switch_on`).
    pub(super) fn new(
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
    ) -> Opening {
        Opening::with_mode(mode(), path, file, shared_standard_input, reports_streamed)
    }

    /// `new` with the switches' value `mode` (`mode()`; the unit tests set it).
    fn with_mode(
        mode: u8,
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
    ) -> Opening {
        let off = |off| Opening {
            mode,
            off: Some(off),
            again: None,
        };
        if mode == 0 {
            return off(Off::NotSet);
        }
        if !reports_streamed {
            return off(Off::Reports);
        }
        if path.as_os_str() == "-" || shared_standard_input {
            return off(Off::StandardInput);
        }
        match file.metadata() {
            Ok(metadata) if metadata.is_file() => {}
            Ok(_) => return off(Off::NotAFile),
            Err(_) => return off(Off::Io),
        }
        let mut handle = file;
        let Ok(start) = handle.stream_position() else {
            return off(Off::Io);
        };
        match file.try_clone() {
            Ok(again) => Opening {
                mode,
                off: None,
                again: Some((again, start)),
            },
            Err(_) => off(Off::Io),
        }
    }

    /// Whether the records are read batch by batch (`Store::read`).
    pub(super) fn acts(&self) -> bool {
        self.again.is_some()
    }
}

/// The first reading: the loop of `read_queries` (`GetNextSeqBatch` until `End()`), each record's
/// length kept and the record dropped, and how the reading ended (with the warnings of the record
/// whose reading failed, as `read_queries` keeps them).
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
fn scan<R: Read>(source: &mut FastaInputSource<R>) -> (Vec<usize>, QueryEnd) {
    let mut lengths = Vec::new();
    loop {
        if source.end() {
            return (lengths, QueryEnd::Input);
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => lengths.push(record.sequence.len()),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => return (lengths, QueryEnd::EofError),
            Err(error) => return (lengths, QueryEnd::Error { error, warnings }),
        }
    }
}

/// The records of a search's queries read batch by batch, or why they are not (its drop prints
/// the summary lines).
pub(super) struct Store {
    mode: u8,
    off: Option<Off>,
    acting: Option<Acting>,
}

struct Acting {
    query_input: QueryInput,
    /// The searched length of each searched query.
    lengths: Vec<usize>,
    ranged: bool,
    scan: Duration,
    state: Mutex<State>,
    /// The shadow: the records as `read_queries` reads them, and the searched ones (cut with
    /// `-query_loc`).
    old: Option<(Vec<FastaRecord>, Option<Vec<FastaRecord>>)>,
}

/// The second reading.
struct State {
    source: FastaInputSource<File>,
    /// The input records read again.
    next: usize,
    /// The input records `base..next`, read and not taken by the loop.
    base: usize,
    window: VecDeque<FastaRecord>,
    /// The loop's batch: its first input record and its input records (up to `base`).
    current: Option<(usize, Arc<Vec<FastaRecord>>)>,
    stats: Stats,
}

#[derive(Default)]
struct Stats {
    batches: u64,
    slot_copies: u64,
    slot_records: u64,
    stale: u64,
    read: Duration,
    peak_window: usize,
    shadow_batches: u64,
    shadow_records: u64,
}

impl Store {
    /// The store of a search whose records are read before it (`opening` does not act).
    pub(super) fn off(opening: Opening) -> Store {
        debug_assert!(!opening.acts());
        Store {
            mode: opening.mode,
            off: opening.off,
            acting: None,
        }
    }

    /// The first reading, on `source` (the reader `search_cli` made of the file, after its
    /// `IsIStreamEmpty`), where `read_queries` reads without the switch: each record's length is
    /// kept and the record dropped. Returns how the reading ended. `config` is the source's
    /// configuration and `range` `-query_loc`.
    pub(super) fn read(
        opening: Opening,
        source: &mut FastaInputSource<File>,
        config: ReaderConfig,
        range: Option<&SequenceRange>,
    ) -> Result<(Store, QueryEnd)> {
        let (again, start) = opening
            .again
            .expect("LOSAT_X_QSTREAM reads a file it can read again");
        let started = Instant::now();
        let (input_lengths, end) = scan(source);
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
        let old = if opening.mode == 2 {
            let mut old_source = FastaInputSource::from_file(reopen(&again)?, config);
            let empty = old_source.stream_is_empty();
            let QueryRecords {
                records,
                end: old_end,
            } = read_queries(&mut old_source);
            drop(old_source);
            let old_input = match range {
                Some(range) => cut_queries(&records, range),
                None => crate::blastinput::seq_range::RangedQueries {
                    input: QueryInput::whole(&records),
                    records: Vec::new(),
                },
            };
            let difference = if empty {
                Some("the old reading finds the file empty".to_string())
            } else if end_text(&old_end) != end_text(&end) {
                Some(format!(
                    "the reading ends with {} where the old one ends with {}",
                    end_text(&end),
                    end_text(&old_end)
                ))
            } else {
                same_query_input(&query_input, &old_input.input)
                    .err()
                    .or_else(|| {
                        let old_lengths: Vec<usize> = match range {
                            Some(_) => old_input.records.iter().map(|r| r.seq().len()).collect(),
                            None => records.iter().map(|r| r.seq().len()).collect(),
                        };
                        (old_lengths != lengths)
                            .then(|| "the searched lengths differ from the old ones".to_string())
                    })
            };
            if let Some(difference) = difference {
                eprintln!("[X_SHADOW] QSTREAM batches=0 records=0 DIFFERENT (the first reading: {difference})");
                panic!("LOSAT_X_QSTREAMSHADOW: the first reading differs from the old reading");
            }
            Some((records, range.map(|_| old_input.records)))
        } else {
            None
        };
        let source = FastaInputSource::from_file(reopen(&again)?, config);
        drop(again);
        let store = Store {
            mode: opening.mode,
            off: None,
            acting: Some(Acting {
                ranged: range.is_some(),
                query_input,
                lengths,
                scan,
                state: Mutex::new(State {
                    source,
                    next: 0,
                    base: 0,
                    window: VecDeque::new(),
                    current: None,
                    stats: Stats::default(),
                }),
                old,
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

// No NCBI counterpart: the summary lines of the switch (with LOSAT_X_STATS) and of the shadow;
// they print counters and change no value NCBI computes.
impl Drop for Store {
    fn drop(&mut self) {
        // After a panic (a shadow's difference) the lines would claim what the run did not finish.
        if self.mode == 0 || std::thread::panicking() {
            return;
        }
        let off = self.off.map_or("none", Off::name);
        let (records, scan, stats) = match &self.acting {
            Some(acting) => {
                let mut state = acting.lock();
                (
                    acting.query_input.input_lengths.len(),
                    acting.scan,
                    (state.next, std::mem::take(&mut state.stats)),
                )
            }
            None => (0, Duration::ZERO, (0, Stats::default())),
        };
        let (read, stats) = stats;
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!(
                "[X_STATS] QSTREAM records={records} read={read} batches={} slot_copies={} slot_records={} stale={} scan_s={:.4} read_s={:.4} peak_window={} off={off}",
                stats.batches,
                stats.slot_copies,
                stats.slot_records,
                stats.stale,
                scan.as_secs_f64(),
                stats.read.as_secs_f64(),
                stats.peak_window,
            );
        }
        if self.mode == 2 {
            eprintln!(
                "[X_SHADOW] QSTREAM batches={} records={} slot_copies={} off={off} identical",
                stats.shadow_batches, stats.shadow_records, stats.slot_copies
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

    /// The shadow: `got` (records `first..` of `what`) are the old reading's.
    fn compare(&self, state: &mut State, what: &str, got: &[FastaRecord], old: &[FastaRecord]) {
        let difference = if got.len() != old.len() {
            Some(format!(
                "{} records where the old reading has {}",
                got.len(),
                old.len()
            ))
        } else {
            got.iter().zip(old).position(|(a, b)| a != b).map(|at| {
                let (a, b) = (&got[at], &old[at]);
                let field = if a.local_id != b.local_id {
                    "local id"
                } else if a.title != b.title {
                    "title"
                } else if a.sequence != b.sequence {
                    "letters"
                } else {
                    "warnings"
                };
                format!("record {at} of the batch differs in its {field}")
            })
        };
        if let Some(difference) = difference {
            eprintln!(
                "[X_SHADOW] QSTREAM batches={} records={} DIFFERENT ({what}: {difference})",
                state.stats.shadow_batches, state.stats.shadow_records
            );
            panic!("LOSAT_X_QSTREAMSHADOW: the records read batch by batch differ from those read before the search");
        }
        state.stats.shadow_records += got.len() as u64;
    }
}

impl State {
    /// Reads the input records up to `end` into the window (`input_lengths`: the first reading's
    /// lengths).
    fn read_to(&mut self, end: usize, input_lengths: &[usize]) -> Result<()> {
        if self.next >= end {
            return Ok(());
        }
        let started = Instant::now();
        while self.next < end {
            let record = self.read_one(input_lengths)?;
            self.window.push_back(record);
            self.next += 1;
        }
        self.stats.read += started.elapsed();
        self.stats.peak_window = self.stats.peak_window.max(self.window.len());
        Ok(())
    }

    /// The next input record, read as the first reading read it (`End()`, then
    /// `GetNextSequence`), with the length the first reading found.
    fn read_one(&mut self, input_lengths: &[usize]) -> Result<FastaRecord> {
        let index = self.next;
        let changed = |what: String| {
            anyhow::anyhow!(
                "LOSAT_X_QSTREAM: the query file changed while LOSAT read it (record {}: {what})",
                index + 1
            )
        };
        if self.source.end() {
            return Err(changed("the file ends before it".to_string()));
        }
        match self.source.next_sequence(&mut |_| Ok(())) {
            Ok(record) if record.sequence.len() == input_lengths[index] => Ok(record),
            Ok(record) => Err(changed(format!(
                "{} letters where the first reading read {}",
                record.sequence.len(),
                input_lengths[index]
            ))),
            Err(error) => Err(changed(format!("{error}"))),
        }
    }

    /// Input record `ordinal`, read and not yet dropped (in the loop's batch or the window).
    fn record(&self, ordinal: usize) -> &FastaRecord {
        if ordinal >= self.base {
            &self.window[ordinal - self.base]
        } else {
            let (start, records) = self
                .current
                .as_ref()
                .expect("an input record of the loop's batch");
            &records[ordinal - start]
        }
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
}

impl BatchRecords<'_> {
    pub(super) fn input(&self) -> &[FastaRecord] {
        match self {
            BatchRecords::Borrowed { input, .. } => input,
            BatchRecords::Read { input, .. } => input,
        }
    }

    pub(super) fn searched(&self) -> &[FastaRecord] {
        match self {
            BatchRecords::Borrowed { searched, .. } => searched,
            BatchRecords::Read { input, cut } => cut.as_deref().unwrap_or(input),
        }
    }
}

impl<'a> Queries<'a> {
    /// Whether the records are read batch by batch.
    pub(super) fn streamed(&self) -> bool {
        matches!(self, Queries::Streamed(_))
    }

    /// The number of input records.
    pub(super) fn input_count(&self) -> usize {
        match self {
            Queries::Records { input, .. } => input.len(),
            Queries::Streamed(store) => store.acting().query_input.input_lengths.len(),
        }
    }

    /// The searched length of each searched query.
    pub(super) fn searched_lengths(&self) -> Vec<usize> {
        match self {
            Queries::Records { searched, .. } => {
                searched.iter().map(|record| record.seq().len()).collect()
            }
            Queries::Streamed(store) => store.acting().lengths.clone(),
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
            Queries::Streamed(store) => {
                let acting = store.acting();
                let input_lengths = &acting.query_input.input_lengths;
                let mut state = acting.lock();
                while state.next < input_lengths.len() {
                    debug_assert!(state.window.is_empty() && state.current.is_none());
                    let record = state.read_one(input_lengths)?;
                    if let Some((old, _)) = &acting.old {
                        let index = state.next;
                        acting.compare(
                            &mut state,
                            "an input record",
                            std::slice::from_ref(&record),
                            &old[index..index + 1],
                        );
                    }
                    state.next += 1;
                    state.base = state.next;
                    out.write_all(&record.warnings)?;
                }
            }
        }
        Ok(())
    }

    /// The records of the loop's batch of input records `input` and searched queries `searched`.
    /// The loop takes its batches in order, each starting where the one before ended.
    pub(super) fn batch(
        &self,
        input: Range<usize>,
        searched: Range<usize>,
    ) -> Result<BatchRecords<'a>> {
        match *self {
            Queries::Records {
                input: input_records,
                searched: searched_records,
            } => Ok(BatchRecords::Borrowed {
                input: &input_records[input],
                searched: &searched_records[searched],
            }),
            Queries::Streamed(store) => {
                let acting = store.acting();
                let mut state = acting.lock();
                anyhow::ensure!(
                    input.start == state.base,
                    "LOSAT_X_QSTREAM: the batch of the input records {input:?} after those up to {}",
                    state.base
                );
                state.read_to(input.end, &acting.query_input.input_lengths)?;
                let records: Vec<FastaRecord> = state.window.drain(..input.len()).collect();
                state.base = input.end;
                let records = Arc::new(records);
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
                if let Some((old, old_searched)) = &acting.old {
                    acting.compare(
                        &mut state,
                        "the input records of a batch of the loop",
                        &records,
                        &old[input.clone()],
                    );
                    if let (Some(cut), Some(old_searched)) = (&cut, old_searched) {
                        acting.compare(
                            &mut state,
                            "the searched records of a batch of the loop",
                            cut,
                            &old_searched[searched],
                        );
                    }
                    state.stats.shadow_batches += 1;
                }
                Ok(BatchRecords::Read {
                    input: records,
                    cut,
                })
            }
        }
    }

    /// The loop is done with its batch that ends at input record `input_end` (its reports are
    /// written, or kept): its records are dropped (a slot that wants them now searches a batch
    /// the loop has passed).
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

    /// The searched records `searched` for a LOSAT_X_QBATCHPAR slot: borrowed from the records
    /// read before the search, or a copy of those read again. `None` when the loop has passed
    /// their first input record (the loop never takes a batch it has passed, so the slot's search
    /// of it is cancelled).
    pub(super) fn slot_records(
        &self,
        searched: Range<usize>,
    ) -> Result<Option<std::borrow::Cow<'a, [FastaRecord]>>> {
        let store = match *self {
            Queries::Records {
                searched: records, ..
            } => return Ok(Some(std::borrow::Cow::Borrowed(&records[searched]))),
            Queries::Streamed(store) => store,
        };
        if searched.is_empty() {
            return Ok(Some(std::borrow::Cow::Owned(Vec::new())));
        }
        let acting = store.acting();
        let first = acting.query_input.ordinal(searched.start);
        let last = acting.query_input.ordinal(searched.end - 1) + 1;
        let mut state = acting.lock();
        let held_from = state
            .current
            .as_ref()
            .map_or(state.base, |(start, _)| *start);
        if first < held_from {
            state.stats.stale += 1;
            return Ok(None);
        }
        state.read_to(last, &acting.query_input.input_lengths)?;
        let copy: Vec<FastaRecord> = searched
            .clone()
            .map(|index| acting.searched(state.record(acting.query_input.ordinal(index)), index))
            .collect();
        state.stats.slot_copies += 1;
        state.stats.slot_records += copy.len() as u64;
        if let Some((old, old_searched)) = &acting.old {
            let old_searched = old_searched.as_deref().unwrap_or(old);
            acting.compare(
                &mut state,
                "the records of a batch searched by a slot",
                &copy,
                &old_searched[searched],
            );
        }
        Ok(Some(std::borrow::Cow::Owned(copy)))
    }

    /// The searched records of the batches before searched query `start`, for a writer that
    /// writes their reports: the records read before the search, or `kept` (the records of the
    /// batches taken, which the loop keeps where such a writer runs: `run_in_pool`).
    pub(super) fn before<'s>(&'s self, kept: &'s [FastaRecord], start: usize) -> &'s [FastaRecord] {
        match self {
            Queries::Records { searched, .. } => &searched[..start],
            Queries::Streamed(_) => {
                assert!(
                    start <= kept.len(),
                    "LOSAT_X_QSTREAM: a writer reads the records of earlier batches, which were not kept"
                );
                &kept[..start]
            }
        }
    }

    /// Every searched record, for the reports after the loop (`before` with every query).
    pub(super) fn all<'s>(&'s self, kept: &'s [FastaRecord]) -> &'s [FastaRecord] {
        match self {
            Queries::Records { searched, .. } => searched,
            Queries::Streamed(_) => kept,
        }
    }
}

/// The query input and the searched lengths of a search from its input records' lengths: what
/// `cut_queries` makes of the records with `-query_loc` (`range`), else `QueryInput::whole`,
/// without their letters (the shadow compares them; tests below too).
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
    use crate::blastinput::seq_range::RangeRole;

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
    /// with the switches' value `mode`, and how the first reading ended.
    fn store(mode: u8, path: &Path, range: Option<&SequenceRange>) -> (Store, QueryEnd) {
        let file = File::open(path).expect("open the test file");
        let opening = Opening::with_mode(mode, path, &file, false, true);
        assert!(opening.acts());
        let mut source = FastaInputSource::from_file(file, config());
        assert!(!source.stream_is_empty());
        Store::read(opening, &mut source, config(), range).expect("the first reading")
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
        let mut eol = |text: &mut Vec<u8>, rng: &mut Rng| match ends {
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
        if seed % 7 == 0 {
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
            let length = if seed % 5 == 0 && index == 2 {
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

    /// The batches of the second reading in `sizes` records (cycled), with the slots' copies of
    /// the next records and a request the loop has passed: every record equals the first's.
    fn read_in_batches(
        mode: u8,
        path: &Path,
        old: &QueryRecords,
        range: Option<&SequenceRange>,
        sizes: &[usize],
    ) {
        let (store, end) = store(mode, path, range);
        assert_eq!(end_text(&end), end_text(&old.end), "{path:?}");
        let cut = range.map(|range| cut_queries(&old.records, range));
        let searched_old: &[FastaRecord] = cut.as_ref().map_or(&old.records, |cut| &cut.records);
        let queries = Queries::Streamed(&store);
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
            // A slot copies records ahead of the loop (read into the window first).
            let ahead = input.searched_in(end..(end + 3).min(count));
            let copy = queries.slot_records(ahead.clone()).unwrap().unwrap();
            assert_eq!(&*copy, &searched_old[ahead.clone()], "{path:?}");
            let batch = queries.batch(start..end, searched.clone()).unwrap();
            assert_eq!(batch.input(), &old.records[start..end], "{path:?}");
            assert_eq!(
                batch.searched(),
                &searched_old[searched.clone()],
                "{path:?}"
            );
            // A slot of this batch (or a mispredicted one starting here) copies from it.
            let copy = queries.slot_records(searched.clone()).unwrap().unwrap();
            assert_eq!(&*copy, &searched_old[searched.clone()], "{path:?}");
            if searched.start > 0 {
                // A batch before this one: the loop never takes it.
                assert!(queries
                    .slot_records(searched.start - 1..searched.start)
                    .unwrap()
                    .is_none());
            }
            drop(batch);
            queries.done(end);
            start = end;
        }
    }

    /// Texts that exercise the reader's states between records (line ends, buffer windows, long
    /// records, warnings, an error at the end): the second reading reads the records of the
    /// first, in batches of every size, with and without `-query_loc`, in the switch's mode and in
    /// the shadow's (which compares inside too).
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
                for mode in [1, 2] {
                    read_in_batches(mode, &file.0, &old, None, sizes);
                    read_in_batches(mode, &file.0, &old, Some(&range), sizes);
                }
            }
        }
        assert!(
            warned >= 10 && failed >= 5 && long >= 5,
            "{warned} {failed} {long}"
        );
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
            for mode in [1, 2] {
                let (store, _) = store(mode, &file.0, None);
                let mut written = Vec::new();
                Queries::Streamed(&store)
                    .write_input_warnings(&mut written)
                    .unwrap();
                assert_eq!(written, expected, "seed {seed}");
            }
        }
    }

    /// A file that changes between the two readings ends the search with LOSAT's error.
    #[test]
    fn a_file_changed_between_the_readings_is_an_error() {
        let file = TempFile::new("changed", b">a\nACGTACGT\n>b\nACGT\n");
        let (store, _) = store(1, &file.0, None);
        std::fs::write(&file.0, b">a\nACGTACG\n>b\nACGT\n").unwrap();
        let error = Queries::Streamed(&store)
            .batch(0..2, 0..2)
            .err()
            .expect("an error")
            .to_string();
        assert!(error.contains("changed while LOSAT read it"), "{error}");
    }

    /// Where the switch does not act: off, the reports after the loop, standard input, a directory
    /// (not a regular file).
    #[test]
    fn inputs_that_cannot_be_read_again() {
        let file = TempFile::new("off", b">a\nACGT\n");
        let handle = File::open(&file.0).unwrap();
        assert_eq!(
            Opening::with_mode(0, &file.0, &handle, false, true).off,
            Some(Off::NotSet)
        );
        assert_eq!(
            Opening::with_mode(1, &file.0, &handle, false, false).off,
            Some(Off::Reports)
        );
        assert_eq!(
            Opening::with_mode(1, Path::new("-"), &handle, false, true).off,
            Some(Off::StandardInput)
        );
        assert_eq!(
            Opening::with_mode(1, &file.0, &handle, true, true).off,
            Some(Off::StandardInput)
        );
        assert!(Opening::with_mode(1, &file.0, &handle, false, true).acts());
        assert!(Opening::with_mode(2, &file.0, &handle, false, true).acts());
        #[cfg(unix)]
        {
            let dir = File::open(std::env::temp_dir()).unwrap();
            assert_eq!(
                Opening::with_mode(1, &std::env::temp_dir(), &dir, false, true).off,
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
            let (lengths, end) = scan(&mut source);
            let old_lengths: Vec<usize> = old.records.iter().map(|r| r.seq().len()).collect();
            assert_eq!(lengths, old_lengths, "seed {seed}");
            assert_eq!(end_text(&end), end_text(&old.end), "seed {seed}");
        }
    }
}
