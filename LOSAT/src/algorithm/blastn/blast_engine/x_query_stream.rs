//! EXPERIMENT (LOSAT_X_QSTREAM, in the strict set since 2026-10-18: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set; shadow LOSAT_X_QSTREAMSHADOW, off unless set): blastn reads the records of each query batch from the query file when the batch loop of
//! `run_in_pool` takes the batch (or a LOSAT_X_QBATCHPAR slot searches it ahead of the loop), as
//! NCBI reads its query batches one at a time, instead of holding the records of the whole query
//! file for the whole search. It saves memory (the records of the batches not in flight) and costs
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
//!   loop takes its batch's records when it reaches the batch (`Queries::batch`), and a slot gets
//!   those of the batch it searches ahead (`Queries::slot_records`). The records read but not yet
//!   taken by the loop wait in a window of chunks, each a run of consecutive records read at once
//!   (`read_to`): a slot's records are read into a chunk of their own (the records before them
//!   into another), which the slot shares with the window and the loop takes as it is when it
//!   reaches that batch, so a record is held once; a request that is not exactly one chunk (a
//!   mispredicted batch, `-query_loc`'s cut records) gets a copy. The loop's batch is held until
//!   the batch is done (a slot may still read from it), then dropped.
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
//! length is compared with the first reading's: a record whose length changed, a file that ends
//! early or a reader error (the file changed between the two readings) ends the search with
//! LOSAT's error, and every later reading of the search returns it (`State::failed`, also after a
//! slot's reading failed for a batch the loop never takes); a change that keeps every record's
//! length is read as the new content (as NCBI, which reads per batch, would read it). The batches
//! are cut from the first reading's
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
//! The rule (`Rule::off`; not for the shadow): the switch does not act on a regular file where the
//! second reading saves little memory and costs time (v4-D sessions 4 and 5, ws windows w1, w4):
//! - `-task megablast` (`Megablast`): its searches are fast for their query bytes, so the serial
//!   first reading is several % of the time (31,168 cDNAs against 3 Mb at 4 and 16 threads: 1.052
//!   and 1.037 of the time; 1,775 cDNAs at 8 and 16 threads: 1.039 and 1.063), and its max RSS
//!   stays below NCBI's without the switch (31,168 cDNAs at 4 threads: 241 MiB, NCBI 358.5);
//! - a file shorter than twice the largest query batch (`Batch`, `Rule::batch_letters`: `BATCH_SIZE`,
//!   or the `CBatchSizeMixer` maximum, 999,000 letters for blastn and 4,999,000 for dc-megablast):
//!   the batch in flight holds most of the file anyway. The file's bytes are compared with the
//!   batch's letters (titles and line ends make the file a little longer than its letters);
//! - where the slots of LOSAT_X_QBATCHPAR may run (LOSAT_X_QBATCHPAR on, the conditions of
//!   `run.rs`'s `x_qbatch_par_off`: `-subject`, two threads or more, fewer subject units than
//!   threads; and a budget: LOSAT_X_QBATCHPARMEM when set, else the default one of
//!   `x_qbatch_par::ncbi_splits_by_queries`, or a slot count set), a file shorter than
//!   twice what their default budget holds (`Slots`), `SLOT_WINDOW_BYTES_PER_THREAD` bytes per thread
//!   beyond the first: the slots read most of the file ahead of the loop anyway (blastn against 3 Mb
//!   at 16 threads: a 40 MB query file -12 MiB at 1.012 of the time, 15 MB -1.6 MiB at 1.014).
//!
//! Elsewhere the second reading saved 3-77 MiB at 0.994-1.006 of the time (blastn and
//! dc-megablast).
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
//! (...)` and panics; at the end of the search `[X_SHADOW] QSTREAM batches= records= slots=
//! off= identical` (`DIFFERENT` instead if a difference was found, also one whose panic on a slot's
//! thread was dropped with a batch the loop never took; batches: the loop's batches compared;
//! records: the records compared, the
//! loop's input records and every slot's records; slots: the slots' requests; off: why the switch
//! did not act, `none` when it did). With `LOSAT_X_STATS` and the switch, `[X_STATS] QSTREAM
//! records= read= batches= slots= slot_copies= slot_records= kept= stale= scan_s= read_s=
//! peak_window= off=` (records: the input records; read: the records read again; slots: the
//! slots' requests, slot_copies of them answered with a copy of slot_records records, the others
//! with a shared chunk; kept: the searched records the loop kept for a writer after their batch,
//! 0 where LOSAT_X_RPTSTREAM writes each batch's reports; stale: slot requests for records the loop had passed, answered with QBATCHPAR's
//! `Cancelled`; scan_s: the first reading; read_s: the second, summed over the threads that read;
//! peak_window: the most records read ahead of the loop at once).
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

// No NCBI counterpart: reads the switches once (LOSAT_X_QSTREAM is in the strict set; any value
// turns the shadow on, and the shadow implies the switch). They choose when LOSAT reads the same
// records.
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
    /// Neither switch is set.
    NotSet,
    /// The reports are written after the loop (LOSAT_X_RPTSTREAM off), which reads every record.
    Reports,
    /// `-query -` (standard input), which cannot be read twice.
    StandardInput,
    /// A query that is not a regular file (a named pipe, a device, a socket).
    NotAFile,
    /// The file's position or a second handle on it could not be had.
    Io,
    /// `-task megablast` (`Rule::off`).
    Megablast,
    /// The file is shorter than twice the largest query batch (`Rule::off`).
    Batch,
    /// The slots of LOSAT_X_QBATCHPAR would read at least half of the file ahead anyway
    /// (`Rule::off`).
    Slots,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Reports => "reports",
            Off::StandardInput => "standard_input",
            Off::NotAFile => "not_a_file",
            Off::Io => "io",
            Off::Megablast => "megablast",
            Off::Batch => "batch",
            Off::Slots => "slots",
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
    /// LOSAT_X_RPTSTREAM is on (`x_rpt_stream::switch_on`); `rule` is what the rule reads of the
    /// search (`Rule::off`).
    pub(super) fn new(
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
        rule: Rule,
    ) -> Opening {
        Opening::with_mode(
            mode(),
            path,
            file,
            shared_standard_input,
            reports_streamed,
            rule,
        )
    }

    /// `new` with the switches' value `mode` (`mode()`; the unit tests set it).
    fn with_mode(
        mode: u8,
        path: &Path,
        file: &File,
        shared_standard_input: bool,
        reports_streamed: bool,
        rule: Rule,
    ) -> Opening {
        let off = |off| {
            if mode != 0 {
                crate::utils::xreach::x_reach!("LOSAT_X_QSTREAM", false);
            }
            Opening {
                mode,
                off: Some(off),
                again: None,
            }
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
        let bytes = match file.metadata() {
            Ok(metadata) if metadata.is_file() => metadata.len(),
            Ok(_) => return off(Off::NotAFile),
            Err(_) => return off(Off::Io),
        };
        // The rule; the shadow reads the file again wherever it can (a test mode).
        if mode == 1 {
            if let Some(why) = rule.off(bytes) {
                return off(why);
            }
        }
        let mut handle = file;
        let Ok(start) = handle.stream_position() else {
            return off(Off::Io);
        };
        match file.try_clone() {
            Ok(again) => {
                crate::utils::xreach::x_reach!("LOSAT_X_QSTREAM", true);
                Opening {
                    mode,
                    off: None,
                    again: Some((again, start)),
                }
            }
            Err(_) => off(Off::Io),
        }
    }

    /// Whether the records are read batch by batch (`Store::read`).
    pub(super) fn acts(&self) -> bool {
        self.again.is_some()
    }
}

/// No NCBI counterpart (when LOSAT reads the records only): the query bytes that the default
/// budget of LOSAT_X_QBATCHPAR's slots holds in flight per thread beyond the first, about 2 Mb
/// of query in both profiles (`x_qbatch_par::allowance_mib` over the per-letter memory of a batch
/// in flight: 96 MiB at about 50 bytes a blastn letter, 28 MiB at about 15 bytes a megablast or
/// dc-megablast letter).
const SLOT_WINDOW_BYTES_PER_THREAD: u64 = 2_000_000;

/// No NCBI counterpart (when LOSAT reads the records only): what the rule of LOSAT_X_QSTREAM
/// reads of a search (`Rule::off`).
#[derive(Clone, Copy, Debug)]
pub(super) struct Rule {
    /// `-task megablast`.
    pub(super) megablast: bool,
    /// The letters of the largest query batch of the search: `BATCH_SIZE` when it is set, else
    /// the `CBatchSizeMixer` maximum of the task (`run.rs`, `x_query_batch_letters`).
    pub(super) batch_letters: u64,
    /// `-num_threads`.
    pub(super) threads: usize,
    /// Whether the slots of LOSAT_X_QBATCHPAR may run: the switch on, the conditions of
    /// `x_qbatch_par_off` (`run.rs`), and a budget (LOSAT_X_QBATCHPARMEM when it is set, else
    /// `x_qbatch_par::ncbi_splits_by_queries`) or a slot count (LOSAT_X_QBATCHPARSLOTS).
    pub(super) slots_may_run: bool,
}

impl Rule {
    /// Why the switch does not act on a regular query file of `file_bytes` bytes (the module's
    /// rule), or `None` where it acts: `-task megablast`; a file shorter than twice the largest
    /// query batch; where the slots may run, a file shorter than twice the bytes their default
    /// budget holds (`SLOT_WINDOW_BYTES_PER_THREAD` per thread beyond the first).
    fn off(self, file_bytes: u64) -> Option<Off> {
        if self.megablast {
            return Some(Off::Megablast);
        }
        if file_bytes < self.batch_letters.saturating_mul(2) {
            return Some(Off::Batch);
        }
        let window =
            SLOT_WINDOW_BYTES_PER_THREAD.saturating_mul(self.threads.saturating_sub(1) as u64);
        if self.slots_may_run && file_bytes < window.saturating_mul(2) {
            return Some(Off::Slots);
        }
        None
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
fn scan<R: Read>(
    source: &mut FastaInputSource<R>,
    archive: bool,
) -> Result<(Vec<usize>, QueryEnd)> {
    let mut lengths = Vec::new();
    loop {
        if source.end() {
            return Ok((lengths, QueryEnd::Input));
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => {
                if archive {
                    crate::report::archive::check_titles(
                        std::slice::from_ref(&record),
                        "query",
                        "BLASTN",
                    )?;
                }
                lengths.push(record.sequence.len());
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => return Ok((lengths, QueryEnd::EofError)),
            Err(error) => return Ok((lengths, QueryEnd::Error { error, warnings })),
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
    /// The input records `base..next`, read and not taken by the loop, in chunks of consecutive
    /// records (first input record, records), which the slots that search them share.
    base: usize,
    chunks: VecDeque<(usize, Arc<Vec<FastaRecord>>)>,
    /// The loop's batch: its first input record and its input records (up to `base`).
    current: Option<(usize, Arc<Vec<FastaRecord>>)>,
    /// The error of a failed reading: every later reading of the search returns it (the records
    /// of the failed reading are lost, and the reader is past the record that failed).
    failed: Option<String>,
    stats: Stats,
}

#[derive(Default)]
struct Stats {
    batches: u64,
    slots: u64,
    slot_copies: u64,
    slot_records: u64,
    stale: u64,
    /// The searched records the loop kept for a writer after their batch (`Queries::note_kept`).
    kept: u64,
    read: Duration,
    peak_window: usize,
    shadow_batches: u64,
    shadow_records: u64,
    /// The shadow found a difference (its panic, on a slot's thread, may be dropped with a batch
    /// the loop never takes).
    differed: bool,
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
        // The BLAST archive (outfmt 11) refuses a query title it cannot write before the
        // first byte: this reading sees every record.
        let (input_lengths, end) = scan(source, config.keep_defline)?;
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
                    chunks: VecDeque::new(),
                    current: None,
                    failed: None,
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
                "[X_STATS] QSTREAM records={records} read={read} batches={} slots={} slot_copies={} slot_records={} kept={} stale={} scan_s={:.4} read_s={:.4} peak_window={} off={off}",
                stats.batches,
                stats.slots,
                stats.slot_copies,
                stats.slot_records,
                stats.kept,
                stats.stale,
                scan.as_secs_f64(),
                stats.read.as_secs_f64(),
                stats.peak_window,
            );
        }
        if self.mode == 2 {
            eprintln!(
                "[X_SHADOW] QSTREAM batches={} records={} slots={} off={off} {}",
                stats.shadow_batches,
                stats.shadow_records,
                stats.slots,
                if stats.differed {
                    "DIFFERENT"
                } else {
                    "identical"
                }
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
            state.stats.differed = true;
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
    /// Reads the input records up to `end` into one new chunk (`input_lengths`: the first
    /// reading's lengths).
    fn read_to(&mut self, end: usize, input_lengths: &[usize]) -> Result<()> {
        if let Some(failed) = &self.failed {
            return Err(anyhow::anyhow!("{failed}"));
        }
        if self.next >= end {
            return Ok(());
        }
        let started = Instant::now();
        let first = self.next;
        let mut records = Vec::with_capacity(end - first);
        while self.next < end {
            records.push(self.read_one(input_lengths)?);
            self.next += 1;
        }
        self.chunks.push_back((first, Arc::new(records)));
        self.stats.read += started.elapsed();
        self.stats.peak_window = self.stats.peak_window.max(self.next - self.base);
        Ok(())
    }

    /// A chunk boundary at input record `at` (`base <= at <= next`): the chunk across it is cut in
    /// two, its records moved (copied while a slot still shares the chunk).
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
    ) -> Result<Arc<Vec<FastaRecord>>> {
        self.read_to(input.end, input_lengths)?;
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
        let what = if self.source.end() {
            "the file ends before it".to_string()
        } else {
            match self.source.next_sequence(&mut |_| Ok(())) {
                Ok(record) if record.sequence.len() == input_lengths[index] => return Ok(record),
                Ok(record) => format!(
                    "{} letters where the first reading read {}",
                    record.sequence.len(),
                    input_lengths[index]
                ),
                Err(error) => format!("{error}"),
            }
        };
        let message = format!(
            "LOSAT_X_QSTREAM: the query file changed while LOSAT read it (record {}: {what})",
            index + 1
        );
        self.failed = Some(message.clone());
        Err(anyhow::anyhow!(message))
    }

    /// Input record `ordinal`, read and not yet dropped (in the loop's batch or a chunk).
    fn record(&self, ordinal: usize) -> &FastaRecord {
        if ordinal >= self.base {
            let at = self.chunks.partition_point(|(start, _)| *start <= ordinal) - 1;
            let (start, records) = &self.chunks[at];
            &records[ordinal - start]
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

/// The records a LOSAT_X_QBATCHPAR slot searches: those read before the search, a chunk of the
/// records read again (shared with the window and the loop), or a copy.
pub(super) enum SlotRecords<'a> {
    Borrowed(&'a [FastaRecord]),
    Shared(Arc<Vec<FastaRecord>>),
    Owned(Vec<FastaRecord>),
}

impl std::ops::Deref for SlotRecords<'_> {
    type Target = [FastaRecord];

    fn deref(&self) -> &[FastaRecord] {
        match self {
            SlotRecords::Borrowed(records) => records,
            SlotRecords::Shared(records) => records,
            SlotRecords::Owned(records) => records,
        }
    }
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
                    debug_assert!(state.chunks.is_empty() && state.current.is_none());
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
                let records = state.take(input.clone(), &acting.query_input.input_lengths)?;
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

    /// The loop keeps `count` searched records of its batch for a writer after the batch (the
    /// `kept=` count of the `[X_STATS]` line).
    pub(super) fn note_kept(&self, count: usize) {
        if let Queries::Streamed(store) = self {
            store.acting().lock().stats.kept += count as u64;
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
    /// read before the search; read again, the chunk that holds exactly them (shared, without
    /// `-query_loc`: the records gap before them are read into a chunk of their own, so the slots
    /// of consecutive batches each get one), or a copy. `None` when the loop has passed their first
    /// input record (the loop never takes a batch it has passed, so the slot's search of it is
    /// cancelled).
    pub(super) fn slot_records(&self, searched: Range<usize>) -> Result<Option<SlotRecords<'a>>> {
        let store = match *self {
            Queries::Records {
                searched: records, ..
            } => return Ok(Some(SlotRecords::Borrowed(&records[searched]))),
            Queries::Streamed(store) => store,
        };
        if searched.is_empty() {
            return Ok(Some(SlotRecords::Owned(Vec::new())));
        }
        let acting = store.acting();
        let input_lengths = &acting.query_input.input_lengths;
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
        state.read_to(first, input_lengths)?;
        state.read_to(last, input_lengths)?;
        state.stats.slots += 1;
        let shared = if acting.ranged || first < state.base {
            None
        } else {
            state.split(first);
            state.split(last);
            state
                .chunks
                .iter()
                .find(|(start, records)| *start == first && start + records.len() == last)
                .map(|(_, records)| Arc::clone(records))
        };
        let records = match shared {
            Some(records) => SlotRecords::Shared(records),
            None => {
                let copy: Vec<FastaRecord> = searched
                    .clone()
                    .map(|index| {
                        acting.searched(state.record(acting.query_input.ordinal(index)), index)
                    })
                    .collect();
                state.stats.slot_copies += 1;
                state.stats.slot_records += copy.len() as u64;
                SlotRecords::Owned(copy)
            }
        };
        if let Some((old, old_searched)) = &acting.old {
            let old_searched = old_searched.as_deref().unwrap_or(old);
            acting.compare(
                &mut state,
                "the records of a batch searched by a slot",
                &records,
                &old_searched[searched],
            );
        }
        Ok(Some(records))
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
        let opening = Opening::with_mode(mode, path, &file, false, true, NO_RULE);
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
            // A slot gets records ahead of the loop (read into the window first).
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

    /// The slots of consecutive batches, asked out of order, each get one chunk of their own,
    /// shared and not copied; the loop takes the very chunk a slot searched; a request across
    /// chunks gets a copy; every record equals the first reading's.
    #[test]
    fn slots_share_the_chunks_the_loop_takes() {
        let mut bytes = Vec::new();
        for index in 0..12 {
            bytes.extend_from_slice(format!(">r{index} title\nACGTTGCA{index}\nGGCC\n").as_bytes());
        }
        let file = TempFile::new("chunks", &bytes);
        let old = old(&file.0);
        let (store, _) = store(1, &file.0, None);
        let queries = Queries::Streamed(&store);
        let shared = |records: &SlotRecords<'_>| match records {
            SlotRecords::Shared(chunk) => Arc::clone(chunk),
            _ => panic!("a shared chunk"),
        };
        let late = queries.slot_records(6..9).unwrap().unwrap();
        let middle = queries.slot_records(3..6).unwrap().unwrap();
        assert_eq!(&*late, &old.records[6..9]);
        assert_eq!(&*middle, &old.records[3..6]);
        let (late, middle) = (shared(&late), shared(&middle));
        let first = queries.batch(0..3, 0..3).unwrap();
        assert_eq!(first.input(), &old.records[0..3]);
        drop(first);
        queries.done(3);
        let second = queries.batch(3..6, 3..6).unwrap();
        match &second {
            BatchRecords::Read { input, .. } => assert!(Arc::ptr_eq(input, &middle)),
            BatchRecords::Borrowed { .. } => panic!("records read again"),
        }
        drop(second);
        queries.done(6);
        // Across the chunk a slot still shares (6..9) and the next records: a copy.
        let across = queries.slot_records(7..10).unwrap().unwrap();
        assert!(matches!(across, SlotRecords::Owned(_)));
        assert_eq!(&*across, &old.records[7..10]);
        drop(late);
        let third = queries.batch(6..12, 6..12).unwrap();
        assert_eq!(third.input(), &old.records[6..12]);
        let state = store.acting().lock();
        assert_eq!((state.stats.slots, state.stats.slot_copies), (3, 1));
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

    /// A record whose length changed between the two readings ends the search with LOSAT's error.
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

    /// A read error in a slot (a record whose length changed) fails every later reading: the
    /// loop's next batch is LOSAT's error, not a batch with records missing (the slot's records
    /// before the failed one are lost, and the reader is past it).
    #[test]
    fn a_read_error_in_a_slot_fails_the_loops_next_batch() {
        let mut bytes = Vec::new();
        for index in 0..6 {
            bytes.extend_from_slice(format!(">r{index}\nACGT\n").as_bytes());
        }
        let file = TempFile::new("slot_error", &bytes);
        let (store, _) = store(1, &file.0, None);
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
        let slot = queries
            .slot_records(2..5)
            .err()
            .expect("the slot's error")
            .to_string();
        assert!(slot.contains("record 4: 3 letters"), "{slot}");
        let error = queries
            .batch(1..4, 1..4)
            .err()
            .expect("the loop's error")
            .to_string();
        assert_eq!(error, slot);
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
        let (store, _) = store(1, &file.0, None);
        let queries = Queries::Streamed(&store);
        let batch = queries.batch(0..3, 0..3).unwrap();
        let held = match &batch {
            BatchRecords::Read { input, .. } => Arc::clone(input),
            BatchRecords::Borrowed { .. } => panic!("records read again"),
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
            Opening::with_mode(0, &file.0, &handle, false, true, NO_RULE).off,
            Some(Off::NotSet)
        );
        assert_eq!(
            Opening::with_mode(1, &file.0, &handle, false, false, NO_RULE).off,
            Some(Off::Reports)
        );
        assert_eq!(
            Opening::with_mode(1, Path::new("-"), &handle, false, true, NO_RULE).off,
            Some(Off::StandardInput)
        );
        assert_eq!(
            Opening::with_mode(1, &file.0, &handle, true, true, NO_RULE).off,
            Some(Off::StandardInput)
        );
        assert!(Opening::with_mode(1, &file.0, &handle, false, true, NO_RULE).acts());
        assert!(Opening::with_mode(2, &file.0, &handle, false, true, NO_RULE).acts());
        // The rule turns the switch off, not the shadow; it reads the file's length.
        let megablast = Rule {
            megablast: true,
            ..NO_RULE
        };
        assert_eq!(
            Opening::with_mode(1, &file.0, &handle, false, true, megablast).off,
            Some(Off::Megablast)
        );
        assert!(Opening::with_mode(2, &file.0, &handle, false, true, megablast).acts());
        // The file is 8 bytes: off for a batch of 5 letters, on for 4.
        for (letters, off) in [(5, Some(Off::Batch)), (4, None)] {
            let rule = Rule {
                batch_letters: letters,
                ..NO_RULE
            };
            assert_eq!(
                Opening::with_mode(1, &file.0, &handle, false, true, rule).off,
                off
            );
        }
        #[cfg(unix)]
        {
            let dir = File::open(std::env::temp_dir()).unwrap();
            assert_eq!(
                Opening::with_mode(1, &std::env::temp_dir(), &dir, false, true, NO_RULE).off,
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
            let (lengths, end) = scan(&mut source, false).unwrap();
            let old_lengths: Vec<usize> = old.records.iter().map(|r| r.seq().len()).collect();
            assert_eq!(lengths, old_lengths, "seed {seed}");
            assert_eq!(end_text(&end), end_text(&old.end), "seed {seed}");
        }
    }

    /// A rule that never turns the switch off.
    const NO_RULE: Rule = Rule {
        megablast: false,
        batch_letters: 0,
        threads: 1,
        slots_may_run: false,
    };

    /// The rule: off for `-task megablast`, for a file shorter than twice the largest batch, and,
    /// where the slots may run, for a file shorter than twice their window (2,000,000 bytes per
    /// thread beyond the first); the files of v4-D's ws windows (1,775 cDNAs 5,398,904 bytes,
    /// 3,492 cDNAs 10,779,147, 5,000 cDNAs 15,440,261, 40 MB 40,005,015, 31,168 cDNAs 96,191,719).
    #[test]
    fn the_rule_follows_the_task_the_batch_and_the_slots_window() {
        const BLASTN: u64 = 999_000;
        const DC: u64 = 4_999_000;
        let rule = |megablast, batch_letters, threads, slots_may_run| Rule {
            megablast,
            batch_letters,
            threads,
            slots_may_run,
        };
        assert_eq!(NO_RULE.off(0), None);
        // megablast: off at any size, with or without slots.
        for (bytes, threads, slots) in [
            (5_398_904, 16, false),
            (96_191_719, 16, true),
            (1, 1, false),
        ] {
            assert_eq!(
                rule(true, DC, threads, slots).off(bytes),
                Some(Off::Megablast)
            );
        }
        // The batch: blastn off below 1,998,000 bytes, dc-megablast below 9,998,000 (named before
        // the slots), BATCH_SIZE 4,000 below 8,000.
        assert_eq!(
            rule(false, BLASTN, 1, false).off(1_997_999),
            Some(Off::Batch)
        );
        assert_eq!(rule(false, BLASTN, 1, false).off(1_998_000), None);
        assert_eq!(rule(false, BLASTN, 16, false).off(5_398_904), None);
        assert_eq!(rule(false, DC, 16, false).off(5_398_904), Some(Off::Batch));
        assert_eq!(rule(false, DC, 16, true).off(5_398_904), Some(Off::Batch));
        assert_eq!(rule(false, DC, 4, false).off(9_997_999), Some(Off::Batch));
        assert_eq!(rule(false, DC, 4, false).off(9_998_000), None);
        assert_eq!(rule(false, DC, 16, false).off(96_191_719), None);
        assert_eq!(rule(false, 4_000, 4, false).off(7_999), Some(Off::Batch));
        assert_eq!(rule(false, 4_000, 4, false).off(8_000), None);
        // The slots (blastn's batch): 16 threads, a window of 30 MB, off below 60 MB.
        for bytes in [15_440_261, 40_005_015, 59_999_999] {
            assert_eq!(rule(false, BLASTN, 16, true).off(bytes), Some(Off::Slots));
            assert_eq!(rule(false, BLASTN, 16, false).off(bytes), None);
        }
        assert_eq!(rule(false, BLASTN, 16, true).off(60_000_000), None);
        assert_eq!(rule(false, DC, 16, true).off(96_191_719), None);
        // 8 threads: off below 28 MB; 4 threads: below 12 MB; 2 threads: below 4 MB; one thread
        // has no window.
        assert_eq!(rule(false, BLASTN, 8, true).off(40_005_015), None);
        assert_eq!(
            rule(false, BLASTN, 8, true).off(27_999_999),
            Some(Off::Slots)
        );
        assert_eq!(rule(false, BLASTN, 4, true).off(15_440_261), None);
        assert_eq!(
            rule(false, BLASTN, 4, true).off(10_779_147),
            Some(Off::Slots)
        );
        assert_eq!(rule(false, BLASTN, 2, true).off(10_779_147), None);
        assert_eq!(rule(false, BLASTN, 1, true).off(1_998_000), None);
        // The ends saturate.
        assert_eq!(
            rule(false, u64::MAX, 1, false).off(u64::MAX - 1),
            Some(Off::Batch)
        );
        assert_eq!(rule(false, u64::MAX, 1, false).off(u64::MAX), None);
        assert_eq!(rule(false, 0, usize::MAX, true).off(0), Some(Off::Slots));
        assert_eq!(rule(false, 0, usize::MAX, true).off(u64::MAX), None);
    }
}
