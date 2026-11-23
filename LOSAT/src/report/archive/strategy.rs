//! The search strategy file (`-export_search_strategy`): the `Blast4-request` of each query
//! batch, in ASN.1 text.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:635-651
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:366-369
//! The applications call `SaveSearchStrategy` for every query batch, before the batch's
//! search (blastp_app.cpp:262, blastn_app.cpp:282, blastx_app.cpp:264, tblastn_app.cpp:280,
//! tblastx_app.cpp:181): the file holds one `Blast4-request` per batch. `CExportStrategy` also
//! builds the request of the archive (outfmt 11, build_archive.cpp:156-200), after the search.
//! The strategy's request therefore differs from the archive's of the same batch in what the
//! search changes (NCBI BLAST+ 2.17.0, 418 runs of the RG-F1 probe):
//!
//! - no `ident` (the applications pass no client id, blast_app_util.cpp:335-337);
//! - `EffectiveSearchSpace` is the option as it stands before the batch: `-searchsp` as typed in
//!   the first request, else absent there, and the value of the batch before in the later ones
//!   (the options handle is shared by the batches);
//! - blastn's `LCaseMask` parameters hold the lower-case runs of `-lcase_masking` only; the
//!   archive's also hold the DUST intervals that the search added to the query's masks.
//!
//! `ARCHIVE_FORMAT` does not change the file: it is always ASN.1 text.
//!
//! Which requests the file holds when a run stops follows from where they are written: the
//! programs that take their batches one at a time (BLASTN, TBLASTX, TBLASTN, BLASTX) write a
//! batch's request where NCBI does, once the batch is read and is not empty and before any error
//! of its set-up, search or reports is acted on (TBLASTN, TBLASTX and BLASTX may have searched the
//! batch ahead; only the order of the file and of what the run writes matters), so a batch that
//! fails after that point has its request and the batches after it have none; BLASTP, which
//! searches the batches of a call together, releases them where the call ends, up to the batch
//! that ends it. One cursor of the run (`State::next_query`), keyed by the run's index of each
//! batch's first query, keeps every request written once and in order (a batch searched again
//! writes none), and the batch's searched space moves the next request's `EffectiveSearchSpace`
//! option (`note_searched`), with what the batch's set-up sets in the options handle (RPS's
//! database length and `.aux` costs, `OptionsHandle::pass`). A broken order is LOSAT's internal
//! error, not a file that differs. NCBI writes a batch's reports before the next batch's request,
//! while LOSAT may search later batches first (BLASTP's calls, the XML documents that keep every
//! report to the run's end): the run marks where it writes each batch's start (the request's mark
//! among the batch's posts, `diag_log::Posts::request`), and a run that NCBI's error ends there
//! (an exception of the XML stream in a report) has the file cut back to the requests of the
//! batches reached (`end_run`).
//!
//! The applications whose search is another program's write their own requests through the same
//! cursor (`options::RequestKind`): psiblast with a FASTA query one per query (BLASTP's batches
//! of one query, or `psi_round`'s rounds, before the query's first round), rpsblast BLASTP's,
//! rpstblastn BLASTX's. A PSSM query (psiblast's and tblastn's `-in_pssm`) is one query whose
//! request holds the PSSM as the application read it (`requests::RunRequests::with_pssm_query`).
//! NCBI makes that request only with `-db`: with `-subject` its `SaveSearchStrategy` stops the run
//! with "Remote PSI-BL2SEQ is not supported", which the programs raise where NCBI does, after the
//! prolog (`blastp::blast_engine::psiblast_stop_at_pssm_query`, `tblastn::tblastn_app`).

use std::io::Write;
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::Result;

use super::options::OptionInputs;
use super::writer::{write_strategy_request, RequestInputs};
use super::{ArchiveRequest, HandleCarry, OptionsHandle, SearchedContext0};
use crate::report::diag_log::{Posts, ERROR};
use crate::report::run_files::{FileReceiver, RunFileName, RunFiles, StrategyTarget as Target};

/// A writer that counts the bytes written through it to `inner` (the strategy file's ends of
/// requests, `State::unreached`).
struct Counted<'a, W: Write + ?Sized> {
    inner: &'a mut W,
    written: u64,
}

impl<W: Write + ?Sized> Write for Counted<'_, W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(bytes)?;
        self.written += n as u64;
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// The stream of a strategy file: a writer whose written bytes can be cut back (`truncate`) where
/// it is a regular file (`can_cut`).
pub trait StrategySink: Write + Send {
    /// Cuts the file back to its first `len` bytes, after the bytes buffered are written.
    fn truncate(&mut self, len: u64) -> std::io::Result<()>;

    /// Whether `truncate` can cut it: a regular file (not a FIFO, a pipe, a terminal or another
    /// device, whose bytes are gone once written).
    fn can_cut(&self) -> bool;

    /// The run has ended (`end_run`): the file is complete. LOSAT Web's file goes to the run's
    /// receiver here (`ReceivedFile`); a file of the disk is the file system's already.
    fn hand_over(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The sink of a strategy file that LOSAT refuses before any request (`-`, NCBI's standard
/// output): it never gets a byte (`StrategyExport::write_with`).
struct NoSink;

impl Write for NoSink {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl StrategySink for NoSink {
    fn truncate(&mut self, _len: u64) -> std::io::Result<()> {
        Ok(())
    }

    fn can_cut(&self) -> bool {
        true
    }
}

impl StrategySink for std::io::BufWriter<std::fs::File> {
    fn truncate(&mut self, len: u64) -> std::io::Result<()> {
        self.flush()?;
        self.get_ref().set_len(len)
    }

    fn can_cut(&self) -> bool {
        self.get_ref()
            .metadata()
            .is_ok_and(|metadata| metadata.file_type().is_file())
    }
}

/// The strategy file of a caller without files (LOSAT Web, `RunFiles::Receiver`): the file in
/// memory under the name that the command line writes it to, a regular file that can be cut, which
/// the receiver gets whole where the run ends (`end_run`, `hand_over`: its error is the run's
/// internal error), or, for a run that stops before its export is made, when the made file is
/// dropped: the receiver gets the file that the command line leaves then, empty
/// (`create_strategy_file`). (No NCBI counterpart: NCBI writes a file.)
pub struct ReceivedFile {
    name: Vec<u8>,
    bytes: Vec<u8>,
    receiver: Arc<dyn FileReceiver>,
    /// The receiver has got the file (`hand_over`).
    handed: bool,
}

impl std::fmt::Debug for ReceivedFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReceivedFile")
            .field("name", &String::from_utf8_lossy(&self.name))
            .field("len", &self.bytes.len())
            .finish_non_exhaustive()
    }
}

impl Write for ReceivedFile {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        debug_assert!(
            !self.handed,
            "a write of the strategy file after its run ended"
        );
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl StrategySink for ReceivedFile {
    fn truncate(&mut self, len: u64) -> std::io::Result<()> {
        self.bytes
            .truncate(usize::try_from(len).unwrap_or(usize::MAX));
        Ok(())
    }

    fn can_cut(&self) -> bool {
        true
    }

    fn hand_over(&mut self) -> std::io::Result<()> {
        if std::mem::replace(&mut self.handed, true) {
            return Ok(());
        }
        self.receiver.receive(&self.name, &self.bytes)
    }
}

impl Drop for ReceivedFile {
    fn drop(&mut self) {
        // The file of a run that stopped before its export was made (`create_strategy_file_in`):
        // the receiver keeps it (LOSAT Web's run files, emitted after the run). A failure has no
        // run to fail any more; in a debug build it is a loss that a test sees.
        let received = self.hand_over();
        debug_assert!(
            received.is_ok() || std::thread::panicking(),
            "the strategy file {} did not reach the run's files: {received:?}",
            String::from_utf8_lossy(&self.name)
        );
    }
}

struct State {
    sink: Box<dyn StrategySink>,
    /// The requests made for a sink that cannot be cut (`StrategySink::can_cut`), by the run's
    /// first query of each one's batch, in order, until the run reaches their batches' starts
    /// (`release_through`) or ends (`end_run`): such a sink gets only the requests NCBI writes.
    /// `None` for a regular file, which gets each request where it is made and is cut back.
    held: Option<std::collections::VecDeque<(usize, Vec<u8>)>>,
    /// The bytes written to `sink`.
    written: u64,
    /// The end in the file of the last request whose batch's start the run reached
    /// (`release_through`).
    reached: u64,
    /// What the searches set in the options handle after the batches noted so far
    /// (`note_searched`).
    carry: HandleCarry,
    /// The run's searched queries whose batches' requests are written: the next request starts
    /// here.
    next_query: usize,
    /// The run's searched queries whose batches' searched spaces moved `carry`.
    noted_through: usize,
    /// No request starts at or after this query (`stop_after`; `usize::MAX` without a stop).
    stop: usize,
    /// The request of each batch written (`write_with`), by the run's first query of the batch,
    /// with NCBI's posts of the bytes that ASN.1 text replaced in it, until the run takes them
    /// (`take_posts`) to write them where the batch starts: a request mark (`Posts::request`) and
    /// the posts.
    pending: Vec<(usize, Posts)>,
    /// The end in the file of each request written that the run has not reached yet
    /// (`release_through`), by the run's first query of its batch, in order.
    unreached: std::collections::VecDeque<(usize, u64)>,
}

struct Shared {
    request: ArchiveRequest,
    /// The run's options handle, shared with the run's archive (`OptionsHandle`).
    options: OptionsHandle,
    /// LOSAT's refusal of the file (`create_strategy_file`, `check_report_stream`).
    refusal: Option<String>,
    /// The report streams' bytes, each stream's held apart from the run's start until NCBI's
    /// first request (`ReportHolds`): only for a file that LOSAT refuses.
    hold: Option<ReportHolds>,
    state: Mutex<State>,
}

/// What a run with `-export_search_strategy` holds: the request of its options and the file.
#[derive(Clone)]
pub struct StrategyExport(Arc<Shared>);

impl std::fmt::Debug for StrategyExport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StrategyExport").finish_non_exhaustive()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3560-3569
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
/// A strategy file made where NCBI's `SetOptions` makes it (`create_strategy_file`,
/// `open_strategy_file_in`): the file and LOSAT's refusal of it, if any, which the run raises
/// where NCBI writes its first request (`check_report_stream`).
#[derive(Debug)]
pub struct StrategyFile {
    target: StrategyTarget,
    refusal: Option<String>,
}

/// Where the requests of a `StrategyFile` go.
#[derive(Debug)]
enum StrategyTarget {
    /// Nowhere: a file that LOSAT refuses before any request (`-`, NCBI's standard output; the
    /// standard output of LOSAT Web's run, or a file there that is not one of the run's), or the
    /// null device of LOSAT Web's run (`run_files::StrategyTarget::Nowhere`).
    None,
    /// The file on the file system (the command line).
    Disk(std::fs::File),
    /// The file in memory that the caller's receiver gets (LOSAT Web).
    Received(ReceivedFile),
}

/// Creates (and empties) the strategy file where NCBI's `SetOptions` does: after the argument
/// parser accepted the command line and before any option handler runs, so before the
/// `-outfmt` text is read, the database list is resolved, and `-subject`, `-query` and `-out`
/// are opened. A file that cannot be created is NCBI's argument error (exit status 1).
///
/// `out` is the `-out` argument. LOSAT refuses a strategy file that is the standard output
/// (`-`) or the same file as the standard output or `-out`: NCBI writes each request to the
/// output's stream before the batch's report (`-`), or to a second stream on the same file,
/// whose bytes then replace or are replaced by the report's as the two buffers are flushed.
/// NCBI writes its first request only after its formatter is made, its prolog written and the
/// first batch read (the first `SaveSearchStrategy`, blastp_app.cpp:262), so every run that NCBI
/// ends before that point (an empty query, a database error, an option or format error, a first
/// batch it cannot read) writes no request, and its bytes are NCBI's. The refusal is therefore
/// only recorded here and raised where NCBI makes its first request (`check_report_stream`);
/// until then the run goes on as NCBI's, its report stream held back (`hold_report`).
pub fn create_strategy_file(path: &Path, out: Option<&Path>, label: &str) -> Result<StrategyFile> {
    if path.as_os_str() == "-" {
        return Ok(standard_output_refusal(label));
    }
    let file = std::fs::File::create(path)
        .map_err(|_| crate::cli::inaccessible("export_search_strategy", path))?;
    let refusal = is_same_as_output(&file, path, out).then(|| same_as_output(label));
    Ok(StrategyFile {
        target: StrategyTarget::Disk(file),
        refusal,
    })
}

/// `-` (`create_strategy_file`): no file, and LOSAT's refusal.
fn standard_output_refusal(label: &str) -> StrategyFile {
    StrategyFile {
        target: StrategyTarget::None,
        refusal: Some(format!(
            "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's {label}; give a file name"
        )),
    }
}

/// LOSAT's refusal of a strategy file that is the output's (`create_strategy_file`).
fn same_as_output(label: &str) -> String {
    format!(
        "an -export_search_strategy file that is also the output (-out or the standard output) is not supported by LOSAT's {label}: NCBI writes both through separate streams and the file's bytes follow the order of their flushes"
    )
}

/// LOSAT Web's refusal of the strategy file `path`, which the command line's run writes where the
/// run has no file (`what`, `create_strategy_file_in`).
fn not_in_the_browser(path: &Path, what: &str) -> String {
    format!(
        "-export_search_strategy {} ({what}) is not supported in the browser; give a file name of the run's own",
        path.display()
    )
}

/// The strategy file `path` of a run whose files go to `files` (`RunFiles`): on `Disk` the
/// command line's (`create_strategy_file`); for a `Receiver` (LOSAT Web, which has no file system
/// and no `-out`), the file that the command line on Linux makes of the name
/// (`RunFiles::strategy_target`, the model of `-logfile`'s names), in the run's file in memory
/// (`ReceivedFile`) under the name as typed, with the command line's error and refusals: NCBI's
/// error at once for a name that the command line cannot open; `-` and the standard output
/// (`/dev/stdout`) refused as the command line refuses them, at NCBI's first request; no file for
/// the null device (`/dev/null`, `/dev/stdin`). LOSAT Web refuses, at the same point, a name whose
/// command line run writes where the run has no file: the standard error, a device whose writes
/// fail or wait for a reader, a file of the system (procfs, a message queue), or the file of
/// another option of the run (`others`: `-query`, `-subject`, `-logfile`, `-window_masker_db`,
/// the files that the search reads from the run's folder, `web_adapter::open_strategy_file`),
/// whose bytes the command line's strategy file replaces or mixes with its own.
pub fn create_strategy_file_in(
    files: &RunFiles,
    path: &Path,
    out: Option<&Path>,
    label: &str,
    others: &[RunFileName<'_>],
) -> Result<StrategyFile> {
    let RunFiles::Receiver(receiver) = files else {
        return create_strategy_file(path, out, label);
    };
    if path.as_os_str() == "-" {
        return Ok(standard_output_refusal(label));
    }
    let browser = |what: &str| StrategyFile {
        target: StrategyTarget::None,
        refusal: Some(not_in_the_browser(path, what)),
    };
    Ok(match RunFiles::strategy_target(path.as_os_str(), others) {
        Target::File => StrategyFile {
            target: StrategyTarget::Received(ReceivedFile {
                name: path.as_os_str().as_encoded_bytes().to_vec(),
                bytes: Vec::new(),
                receiver: receiver.clone(),
                handed: false,
            }),
            refusal: None,
        },
        Target::Nowhere => StrategyFile {
            target: StrategyTarget::None,
            refusal: None,
        },
        Target::Inaccessible => {
            return Err(crate::cli::inaccessible("export_search_strategy", path))
        }
        Target::StandardOutput => StrategyFile {
            target: StrategyTarget::None,
            refusal: Some(same_as_output(label)),
        },
        Target::System => {
            browser("a device or a file of the system, where the run has no file")
        }
        Target::OtherFile(option) => browser(&format!(
            "the file of -{option} too, whose bytes the command line's strategy file replaces or mixes with its own"
        )),
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:250-262 (the same order in
// blastn_app.cpp:255-282, blastx_app.cpp:251-264, tblastn_app.cpp:267-280, tblastx_app.cpp:170-181,
// psiblast_app.cpp:503-517 and :340, rpsblast_app.cpp:188-199, rpstblastn_app.cpp:185-196)
/// LOSAT's refusal of the run's strategy file (`create_strategy_file`: the standard output, or
/// the file of the output), raised where NCBI writes its first request: after its formatter and
/// its prolog, once the first batch is read and before it is searched (`SaveSearchStrategy`).
/// Every stop of NCBI's before that point (an empty query, a database or option error, a query
/// that cannot be read, a bad `BATCH_SIZE`, psiblast's `-in_pssm` with `-subject`) comes first
/// with NCBI's own bytes: the run's report stream holds what the run writes meanwhile
/// (`hold_report`), and `end_run` writes it before the error. Here the held bytes are dropped
/// and the run refuses, so the stream has no byte of the report before the refusal. The
/// programs call it where they make the first batch's request (`StrategyExport::write_with`
/// calls it too, so no request is ever written to such a file) and, where they write that
/// request after searching, at NCBI's point (BLASTP's `search`).
pub fn check_report_stream(strategy: Option<&StrategyExport>) -> Result<()> {
    let Some(strategy) = strategy else {
        return Ok(());
    };
    match strategy.0.refusal.as_deref() {
        Some(refusal) => {
            if let Some(hold) = &strategy.0.hold {
                hold.drop_held();
            }
            Err(anyhow::anyhow!("{refusal}"))
        }
        None => Ok(()),
    }
}

/// What a report stream gets of the bytes held back until NCBI's first request
/// (`check_report_stream`).
enum HoldState {
    /// The run's bytes so far.
    Holding(Vec<u8>),
    /// The run ended without the request (`end_run`): the bytes wait to be written, once, before
    /// the next write or flush of the stream.
    Released(Vec<u8>),
    /// The request point refused the run: the stream takes no byte.
    Dropped,
}

/// The held state of one report stream of a run whose strategy file LOSAT refuses.
#[derive(Clone)]
struct ReportHold(Arc<Mutex<HoldState>>);

/// What a write of the held stream does (`ReportHold::step`).
enum HoldStep {
    /// The bytes are kept (or, dropped, thrown away).
    Kept,
    /// Write these bytes first, then the write itself.
    Through(Vec<u8>),
}

impl ReportHold {
    fn with_state<R>(&self, f: impl FnOnce(&mut HoldState) -> R) -> R {
        f(&mut self
            .0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner))
    }

    /// A write of `bytes` (a flush writes none).
    fn step(&self, bytes: &[u8]) -> HoldStep {
        self.with_state(|state| match state {
            HoldState::Holding(held) => {
                held.extend_from_slice(bytes);
                HoldStep::Kept
            }
            HoldState::Released(held) => HoldStep::Through(std::mem::take(held)),
            HoldState::Dropped => HoldStep::Kept,
        })
    }

    /// The run's end or its request point (`ReportHolds::end`): the held bytes are released to
    /// the stream (`keep`) or dropped.
    fn end(&self, keep: bool) {
        self.with_state(|state| {
            if let HoldState::Holding(held) = state {
                *state = if keep {
                    HoldState::Released(std::mem::take(held))
                } else {
                    HoldState::Dropped
                };
            }
        });
    }
}

/// Where a run's holds are (`ReportHolds`): before NCBI's first request and the run's end, then
/// released or dropped.
#[derive(Clone, Copy, PartialEq, Eq)]
enum HoldPhase {
    Holding,
    Released,
    Dropped,
}

/// The held report streams of a run whose strategy file LOSAT refuses: one hold for each stream
/// (`hold_report`), so that each stream gets back only its own bytes, as the command line's one
/// report stream gets its own (LOSAT Web writes several formats, each to its own stream). The
/// request point drops them all, the run's end releases or drops them all; a stream held after
/// that gets the run's phase at once.
#[derive(Default)]
struct ReportHolds(Mutex<(Option<HoldPhase>, Vec<ReportHold>)>);

impl ReportHolds {
    fn with_holds<R>(&self, f: impl FnOnce(&mut (Option<HoldPhase>, Vec<ReportHold>)) -> R) -> R {
        f(&mut self
            .0
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner))
    }

    /// The hold of one more report stream.
    fn hold(&self) -> ReportHold {
        self.with_holds(|(phase, holds)| {
            let hold = ReportHold(Arc::new(Mutex::new(
                match phase.unwrap_or(HoldPhase::Holding) {
                    HoldPhase::Holding => HoldState::Holding(Vec::new()),
                    HoldPhase::Released => HoldState::Released(Vec::new()),
                    HoldPhase::Dropped => HoldState::Dropped,
                },
            )));
            holds.push(hold.clone());
            hold
        })
    }

    fn drop_held(&self) {
        self.end(false);
    }

    /// The run ended: the held bytes are each stream's, as NCBI's error left them, unless `keep`
    /// is false (LOSAT's refusal, which has no byte before it). Only the first end counts.
    fn end(&self, keep: bool) {
        self.with_holds(|(phase, holds)| {
            if phase.is_some() {
                return;
            }
            *phase = Some(if keep {
                HoldPhase::Released
            } else {
                HoldPhase::Dropped
            });
            for hold in holds.iter() {
                hold.end(keep);
            }
        });
    }
}

/// The destination of a report stream (`cli::report_destination`) behind its hold.
struct HeldReport<W: Write> {
    inner: W,
    hold: ReportHold,
}

impl<W: Write> Write for HeldReport<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        match self.hold.step(bytes) {
            HoldStep::Kept => Ok(bytes.len()),
            HoldStep::Through(held) => {
                self.inner.write_all(&held)?;
                self.inner.write(bytes)
            }
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self.hold.step(&[]) {
            HoldStep::Kept => Ok(()),
            HoldStep::Through(held) => {
                self.inner.write_all(&held)?;
                self.inner.flush()
            }
        }
    }
}

impl<W: Write> Drop for HeldReport<W> {
    fn drop(&mut self) {
        let _ = Write::flush(self);
    }
}

/// The destination `inner` of the command line's report stream (the `-out` file or standard
/// output) behind the hold of a run whose strategy file LOSAT refuses (`strategy`): from the
/// run's start until NCBI's first request (`check_report_stream`) the stream's bytes are held,
/// not written. NCBI writes nothing to the stream before its formatter, and its first request
/// only after its prolog, so what is held is the prolog (and the reports of the stops that come
/// first); the run's end writes them before its error (`end_run`), the request point drops them
/// and refuses. Each call holds its own stream (`ReportHolds`): a run with more than one report
/// stream (LOSAT Web's formats) gets each stream's bytes back to that stream. The destination as
/// it is, without a hold or a cost.
pub(crate) fn hold_report<'a>(
    strategy: Option<&StrategyExport>,
    inner: Box<dyn Write + 'a>,
) -> Box<dyn Write + 'a> {
    match strategy.and_then(|strategy| strategy.0.hold.as_ref()) {
        Some(holds) => Box::new(HeldReport {
            inner,
            hold: holds.hold(),
        }),
        None => inner,
    }
}

/// `hold_report` of a destination that can move between threads.
pub(crate) fn hold_report_send<'a>(
    strategy: Option<&StrategyExport>,
    inner: Box<dyn Write + Send + 'a>,
) -> Box<dyn Write + Send + 'a> {
    match strategy.and_then(|strategy| strategy.0.hold.as_ref()) {
        Some(holds) => Box::new(HeldReport {
            inner,
            hold: holds.hold(),
        }),
        None => inner,
    }
}

/// Whether the created strategy `file` is the file of `-out` (which exists then: NCBI creates
/// it later) or of the standard output.
#[cfg(unix)]
fn is_same_as_output(file: &std::fs::File, _path: &Path, out: Option<&Path>) -> bool {
    use std::os::fd::AsFd;
    use std::os::unix::fs::{FileTypeExt, MetadataExt};
    let Ok(own) = file.metadata() else {
        return false;
    };
    // The null device: what NCBI writes there is lost whatever the order of the writes.
    if own.file_type().is_char_device()
        && std::fs::metadata("/dev/null").is_ok_and(|null| null.rdev() == own.rdev())
    {
        return false;
    }
    let same = |other: std::io::Result<std::fs::Metadata>| {
        other.is_ok_and(|other| other.dev() == own.dev() && other.ino() == own.ino())
    };
    let stdout_is_report = out.is_none_or(|out| out.as_os_str() == "-");
    match out.filter(|out| out.as_os_str() != "-") {
        Some(out) if same(std::fs::metadata(out)) => return true,
        _ => {}
    }
    stdout_is_report
        && std::io::stdout()
            .as_fd()
            .try_clone_to_owned()
            .is_ok_and(|fd| same(std::fs::File::from(fd).metadata()))
}

/// Where there are no inode numbers (Windows, WASI): the strategy file and `-out` (which must
/// exist, as above) have the same path form (`utils::same_file`). The standard output is not
/// asked about here.
///
/// No NCBI counterpart: NCBI never asks whether its strategy file is its `-out` file (see
/// `create_strategy_file`).
#[cfg(not(unix))]
fn is_same_as_output(_file: &std::fs::File, path: &Path, out: Option<&Path>) -> bool {
    use crate::utils::same_file::path_identity;
    let (Some(out), Some(own)) = (
        out.filter(|out| out.as_os_str() != "-"),
        path_identity(path),
    ) else {
        return false;
    };
    path_identity(out) == Some(own)
}

/// The arguments that make the `HitlistSize` option of a run's requests besides its format
/// (`StrategyExport::for_run`).
#[derive(Clone, Copy, Debug, Default)]
pub struct HitlistArgs {
    pub max_target_seqs: Option<usize>,
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
}

/// What a program's arguments give its search strategy's requests (`StrategyArgs`): the recorded
/// options, the counts of its `HitlistSize` and the database of `-db` as NCBI keeps it.
pub(crate) struct StrategyInputs<'a> {
    pub(crate) inputs: OptionInputs,
    pub(crate) counts: HitlistArgs,
    pub(crate) database: Option<&'a str>,
}

/// The arguments of a program that writes `-export_search_strategy` (`StrategyExport::for_args`):
/// each program's in its archive module, with its archive's recorded options.
pub(crate) trait StrategyArgs {
    fn strategy_inputs(&self) -> StrategyInputs<'_>;
}

/// The `HitlistSize` option of a run of format `outfmt` (`report_counts`): the archive's is that
/// of outfmt 11, while outfmt 0 to 4 take the larger of `-num_descriptions` and `-num_alignments`
/// (blast_args.cpp:2910-2936, probe `nd` of RG-F1).
pub(crate) fn hitlist_size(outfmt: i32, counts: HitlistArgs) -> i32 {
    crate::blastinput::report_counts::report_counts(
        outfmt,
        counts.max_target_seqs,
        counts.num_descriptions,
        counts.num_alignments,
    )
    .hitlist_size as i32
}

/// The strategy file of a run's `-export_search_strategy` argument (`path`), created where NCBI
/// creates it (`create_strategy_file`), or `None` without the argument. A file name that is not
/// UTF-8 is refused first (`check_utf8_file_name`). `out` is the `-out` argument; `label` names
/// the program in a refusal.
///
/// With a strategy file NCBI skips its early reading of `-outfmt` (blast_args.cpp:3624-3627,
/// cited at `create_strategy_file`'s callers): a value its parser refuses (`99`, `abc`) fails
/// only in the formatting handler, the last of the handlers, so `-subject`, `-subject_loc`,
/// `-query`, `-out` and `-query_loc` come first, in the order of the other arguments. The
/// callers read the value with `app::early_format_choice`, which leaves it to their own
/// formatting handler, and make `-out` where they always do (line RG F/d36 review 2c B1).
pub fn open_strategy_file(
    path: Option<&Path>,
    out: Option<&Path>,
    label: &str,
) -> Result<Option<StrategyFile>> {
    open_strategy_file_in(&RunFiles::Disk, path, out, label, &[])
}

/// `open_strategy_file` of a run whose files go to `files` (`create_strategy_file_in`), whose
/// other files are `others`: LOSAT Web's adapter opens the file here, where the command line does.
pub fn open_strategy_file_in(
    files: &RunFiles,
    path: Option<&Path>,
    out: Option<&Path>,
    label: &str,
    others: &[RunFileName<'_>],
) -> Result<Option<StrategyFile>> {
    let Some(path) = path else {
        return Ok(None);
    };
    crate::blastinput::input_files::check_utf8_file_name(path, "export_search_strategy", label)?;
    create_strategy_file_in(files, path, out, label, others).map(Some)
}

/// The options handle of a run with the search strategy `strategy` (its archive's too,
/// `StrategyExport::options`), or a new one without it.
pub fn options_of(strategy: Option<&StrategyExport>) -> OptionsHandle {
    strategy
        .map(|strategy| strategy.options().clone())
        .unwrap_or_default()
}

impl StrategyExport {
    /// The export of a run of `inputs.kind` into `sink` (`create_strategy_file`). `database`
    /// is the `-db` argument as NCBI keeps it (`None` for `-subject`); `options` is the run's
    /// options handle, which its archive shares. The request is the archive's
    /// (`ArchiveRequest::for_inputs`).
    pub fn new(
        sink: Box<dyn StrategySink>,
        inputs: &OptionInputs,
        database: Option<&str>,
        options: OptionsHandle,
    ) -> Self {
        Self::with_refusal(sink, inputs, database, options, None)
    }

    fn with_refusal(
        sink: Box<dyn StrategySink>,
        inputs: &OptionInputs,
        database: Option<&str>,
        options: OptionsHandle,
        refusal: Option<String>,
    ) -> Self {
        Self(Arc::new(Shared {
            request: ArchiveRequest::for_inputs(inputs, database),
            options,
            hold: refusal.is_some().then(ReportHolds::default),
            refusal,
            state: Mutex::new(State {
                held: (!sink.can_cut()).then(std::collections::VecDeque::new),
                sink,
                written: 0,
                reached: 0,
                carry: HandleCarry::default(),
                next_query: 0,
                noted_through: 0,
                stop: usize::MAX,
                pending: Vec::new(),
                unreached: std::collections::VecDeque::new(),
            }),
        }))
    }

    /// The export of a run into its strategy file (`open_strategy_file`), buffered: the program's
    /// recorded options `inputs` with the `HitlistSize` of the run's format `outfmt` and `counts`
    /// (the archive's is outfmt 11's), and the run's options handle `options` (made once per run
    /// and given to its archive too, `ArchiveHandle::for_program`). `outfmt` is none for a
    /// `-outfmt` value that was not read (`app::EarlyFormat::Deferred`): the run stops at its
    /// formatting handler before any request, so the options keep their own hitlist size.
    pub fn for_run(
        file: StrategyFile,
        mut inputs: OptionInputs,
        outfmt: Option<i32>,
        counts: HitlistArgs,
        database: Option<&str>,
        options: &OptionsHandle,
    ) -> Self {
        if let Some(outfmt) = outfmt {
            inputs.hitlist_size = hitlist_size(outfmt, counts);
        }
        // A refused file gets no request (`check_report_stream`, `write_with`); `-` has none.
        let sink: Box<dyn StrategySink> = match file.target {
            StrategyTarget::Disk(file) => Box::new(std::io::BufWriter::new(file)),
            StrategyTarget::Received(file) => Box::new(file),
            StrategyTarget::None => Box::new(NoSink),
        };
        Self::with_refusal(sink, &inputs, database, options.clone(), file.refusal)
    }

    /// The export of the run of the arguments `args` into `file` (`for_run`), with the format
    /// `outfmt` and the run's options handle `options`: the one maker of the command line's runs
    /// and of LOSAT Web's (`web_adapter::give_strategy`).
    pub(crate) fn for_args(
        args: &dyn StrategyArgs,
        file: StrategyFile,
        outfmt: Option<i32>,
        options: &OptionsHandle,
    ) -> Self {
        let StrategyInputs {
            inputs,
            counts,
            database,
        } = args.strategy_inputs();
        Self::for_run(file, inputs, outfmt, counts, database, options)
    }

    /// The run's options handle.
    pub fn options(&self) -> &OptionsHandle {
        &self.0.options
    }

    /// Whether the run's report streams are held until NCBI's first request (`hold_report`): the
    /// file has LOSAT's refusal.
    pub(crate) fn holds_report(&self) -> bool {
        self.0.hold.is_some()
    }

    fn with_state<R>(&self, f: impl FnOnce(&mut State) -> R) -> R {
        let mut state = self
            .0
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        f(&mut state)
    }

    /// Whether the request of the batch of the run's queries `queries` is still to be written:
    /// not one already written (a batch that a run searches again, or a later part of a batch the
    /// search split) and none at or after the run's stop (`stop_after`).
    pub fn takes(&self, queries: Range<usize>) -> bool {
        self.with_state(|state| state.takes(&queries))
    }

    /// The run's searched queries whose batches' requests are written.
    pub fn written_through(&self) -> usize {
        self.with_state(|state| state.next_query)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-288
    /// Makes the request of the batch of the run's searched queries `queries`, where NCBI writes
    /// it: after the batch is read and made a query factory, before its search. `make` builds the
    /// request (from the batch's records, `requests::RunRequests`) and hands it to its argument,
    /// which writes it with what the batches before set in the options handle (`note_searched`),
    /// all under one lock, and the file is complete after it (the command line's exit does not run
    /// destructors). NCBI writes a batch's request, searches the batch and writes its reports
    /// before the next batch's request, so a run that its reports end (an exception of the XML
    /// stream, a failed write) has the requests up to that batch; LOSAT may search several batches
    /// before their reports (BLASTP's calls, the XML documents that hold every report until the
    /// run's end). So the run marks where it reaches each batch's start, with NCBI's posts of the
    /// bytes that ASN.1 text replaced in the batch's request (`replaced_posts`): the request mark
    /// that `take_posts` hands out with them (`Posts::write`, `release_through`); `end_run` cuts the
    /// file back to the requests reached when NCBI's error ends the run. (Holding the requests'
    /// bytes until then instead would keep every request of a BLASTN XML run in memory: 2.6 MB on a
    /// run of 82 MB, review 1 of K4h gap 2; so only a sink that cannot be cut holds them,
    /// `State::held`.)
    /// Returns false, calling nothing, when the request is not to be written (`takes`). The
    /// requests are made in the run's order, each starting where the one before ended, after
    /// the batch before was noted; anything else is LOSAT's internal error.
    pub fn write_with(
        &self,
        queries: Range<usize>,
        make: impl FnOnce(&mut dyn FnMut(&RequestInputs<'_>) -> Result<()>) -> Result<()>,
    ) -> Result<bool> {
        check_report_stream(Some(self))?;
        self.with_state(|state| {
            if !state.takes(&queries) {
                return Ok(false);
            }
            anyhow::ensure!(
                queries.start == state.next_query && state.noted_through == state.next_query,
                "LOSAT internal error: the search strategy's request of the queries {queries:?} does not follow the requests written ({}) and noted ({})",
                state.next_query,
                state.noted_through
            );
            let shared = &self.0;
            let mut posts = Posts::request(self.clone(), queries.start);
            let mut made_one = false;
            let made = make(&mut |request| {
                let sets = shared.options.request_sets(&state.carry);
                // A big integer of -2^63 (RPS-BLAST's `DbLength` that a batch before set,
                // `-dbsize`), which NCBI writes as `-0`, is written as its digits
                // (PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits).
                // A sink that cannot be cut gets the request where the run reaches its batch.
                let mut held = Vec::new();
                let mut counted = Counted {
                    inner: match state.held {
                        Some(_) => &mut held as &mut dyn Write,
                        None => &mut *state.sink as &mut dyn Write,
                    },
                    written: 0,
                };
                let written = write_strategy_request(
                    &mut counted,
                    &shared.request,
                    request,
                    &sets,
                    &shared.options.strategy_database_params(),
                    |replaced| {
                        posts.append(&replaced_posts(shared.request.application, replaced));
                        Ok(())
                    },
                );
                let flushed = match state.held.as_mut() {
                    Some(requests) => {
                        requests.push_back((queries.start, held));
                        Ok(())
                    }
                    None => {
                        state.written += counted.written;
                        state.sink.flush()
                    }
                };
                state.unreached.push_back((queries.start, state.written));
                state.next_query = queries.end;
                made_one = true;
                written?;
                flushed?;
                Ok(())
            });
            if made_one {
                state.pending.push((queries.start, posts));
            }
            made?;
            Ok(true)
        })
    }

    /// The run reaches the start of the batch whose first query (the run's) is `first`
    /// (`Posts::write` of its request mark), after the reports of the batches before it: the
    /// requests of the batches through it are NCBI's whatever ends the run (`end_run`); a sink
    /// that cannot be cut gets them here.
    pub(crate) fn release_through(&self, first: usize) -> std::io::Result<()> {
        self.with_state(|state| {
            while let Some(&(at, end)) = state.unreached.front() {
                if at > first {
                    break;
                }
                state.reached = end;
                state.unreached.pop_front();
            }
            state.write_held(|at| at <= first)
        })
    }

    /// Whether `other` is this export (`Posts`'s request marks compare by it).
    pub(crate) fn same(&self, other: &StrategyExport) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// `write_with` of a request already made.
    #[cfg(test)]
    pub fn write_request(
        &self,
        queries: Range<usize>,
        request: &RequestInputs<'_>,
    ) -> Result<bool> {
        self.write_with(queries, |write| write(request))
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-282
    /// Takes the posts of the requests written for the batches whose first query (the run's) is
    /// below `before`, in the run's order, each with that query: the request's mark (which marks the
    /// request reached where the run writes it, `release_through`) and NCBI's posts of its replaced
    /// bytes.
    /// NCBI writes the request and posts them after it read the batch (the reader's messages of
    /// its records) and before the batch's search, so after the reports of the batch before; the
    /// run writes them there, with the posts before the first query's report (`QueryWarnings`),
    /// or before the error that stops the batch.
    pub fn take_posts(&self, before: usize) -> Vec<(usize, Posts)> {
        self.with_state(|state| {
            let taken = state.pending.partition_point(|&(first, _)| first < before);
            state.pending.drain(..taken).collect()
        })
    }

    /// Every post that `take_posts` would take, as one sequence.
    pub fn take_all_posts(&self) -> Posts {
        let mut all = Posts::new();
        for (_, posts) in self.take_posts(usize::MAX) {
            all.append(&posts);
        }
        all
    }

    /// Writes every post waiting (`take_all_posts`) to the run's diagnostics `diagnostics`, each
    /// whole post a message of the run's archive `archive`, and their requests reached
    /// (`Posts::write`): where the run writes the start of the batch whose request was written
    /// last, after its reader's messages and before its error.
    pub fn write_posts(
        &self,
        diagnostics: &mut dyn Write,
        archive: Option<&super::ArchiveHandle>,
    ) -> std::io::Result<()> {
        self.take_all_posts().write(diagnostics, archive)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:177-181
    /// Moves what the requests after the batch of the run's searched queries `queries` hold of the
    /// options handle past that batch, once the batch's set-up has written it into the handle
    /// (above; `OptionsHandle::pass`: the set-up's parameters, then the search space): `searched`
    /// is the batch's context 0 when NCBI searched the batch (the archive's
    /// `BatchArchive::searched_space`, from the same values), `None` for a batch NCBI does not
    /// search. A batch noted already (a search of the batches again, psiblast's later rounds of a
    /// query) is not passed again; a batch noted out of the run's order is LOSAT's internal error.
    pub fn note_searched(
        &self,
        queries: Range<usize>,
        searched: Option<SearchedContext0>,
    ) -> Result<()> {
        self.with_state(|state| {
            if queries.end <= state.noted_through {
                return Ok(());
            }
            anyhow::ensure!(
                queries.start == state.noted_through,
                "LOSAT internal error: the search strategy's batch of the queries {queries:?} is noted after the batches through {}",
                state.noted_through
            );
            self.0.options.pass(&mut state.carry, searched);
            state.noted_through = queries.end;
            Ok(())
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
    /// NCBI stops the run when it formats the batch whose queries end at `queries_end`, after
    /// that batch's request (the BLAST XML report data of a matrix name it does not know,
    /// `blast_xml::matrix::XmlMatrixFailure`, `QueryBatches::xml_matrix_request_end`): no request
    /// of a later batch is written.
    pub fn stop_after(&self, queries_end: usize) {
        self.with_state(|state| state.stop = state.stop.min(queries_end));
    }
}

// No NCBI counterpart: NCBI posts while it writes the request, so no post can wait there.
/// The end of a run of a program with `-export_search_strategy` (`strategy`; `None` without
/// it). The requests written (`StrategyExport::write_with`) whose batches' starts the run did not
/// reach: a run that ends with one of NCBI's errors (`cli::native_error_of`) ended before those starts,
/// where NCBI's run ends before it writes them (its reports' exception, or a later batch's error,
/// whose batch's start the run writes before the error), so the file is cut back to the requests
/// reached; a run that ends with `Ok` keeps them all, as NCBI writes every request (and so does a
/// refusal of LOSAT's, whose file is not NCBI's). A strategy file that cannot be cut (a FIFO, a
/// pipe, a terminal: `StrategySink::can_cut`) has got only the requests reached
/// (`StrategyExport::release_through`): the others are dropped at NCBI's error, written at any
/// other end. A regular file that fails to be cut is LOSAT's error, not a file with requests that
/// NCBI does not write. A run that ends as NCBI's (`Ok` or NCBI's error)
/// has written every post of its requests at its batch's start (`StrategyExport::take_posts`), so
/// a post still waiting is LOSAT's internal error instead of output that silently lacks it.
/// LOSAT's refusals (any other error) keep their own text. Each program calls this once, where
/// its run's result is known (`blastn`, `blastp` for psiblast and rpsblast too, `blastx` for
/// rpstblastn too, `tblastn`, `tblastx`; LOSAT Web's records searches, `ReportOutputs::run_records`).
/// The posts still waiting go then: each holds the export (its request's mark), so that the export
/// would outlive the run, and a sink that ends with the export (LOSAT Web's file of the run,
/// `ReceivedFile`) would never reach its receiver.
pub fn end_run(strategy: Option<&StrategyExport>, result: Result<()>) -> Result<()> {
    let Some(strategy) = strategy else {
        return result;
    };
    let ended = end_export(strategy, result);
    let waiting = strategy.with_state(|state| std::mem::take(&mut state.pending));
    drop(waiting);
    // The file is complete: LOSAT Web's goes to the run's files now (`ReceivedFile`), and a file
    // that does not get there is LOSAT's internal error instead of a run without it. A refusal of
    // LOSAT's keeps its own text.
    let handed = strategy.with_state(|state| state.sink.hand_over());
    match handed {
        Err(error) if !matches!(&ended, Err(refusal) if crate::cli::native_error_of(refusal).is_none()) => {
            Err(anyhow::anyhow!(
                "LOSAT internal error: the -export_search_strategy file did not reach the run's files ({error})"
            ))
        }
        _ => ended,
    }
}

/// `end_run` of a run with the export `strategy`, before its waiting posts go.
fn end_export(strategy: &StrategyExport, result: Result<()>) -> Result<()> {
    // NCBI's error (`cli::native_error_of`, which the command line's exit reads too): a
    // `NativeError`, or an exception of NCBI's that a report's writer carries (the `-html`
    // links' registry, which NCBI's prolog reads before its first request).
    let native = matches!(&result, Err(error) if crate::cli::native_error_of(error).is_some());
    // The report stream's held bytes (`hold_report`): where the run ends before NCBI's first
    // request they are NCBI's, written before its error (or its success); LOSAT's refusal has
    // none before it.
    if let Some(hold) = &strategy.0.hold {
        hold.end(result.is_ok() || native);
    }
    let ended = strategy.with_state(|state| -> std::io::Result<()> {
        state.unreached.clear();
        if !native {
            return state.write_held(|_| true);
        }
        if let Some(held) = state.held.as_mut() {
            held.clear();
        }
        if state.written > state.reached {
            state.sink.truncate(state.reached)?;
            state.written = state.reached;
        }
        Ok(())
    });
    if let Err(error) = ended {
        // A refusal of LOSAT's keeps its own text.
        if result.is_err() && !native {
            return result;
        }
        return Err(anyhow::anyhow!(
            "LOSAT could not end the -export_search_strategy file with the requests NCBI BLAST+ writes ({error})"
        ));
    }
    if result.is_err() && !native {
        return result;
    }
    let waiting = strategy.with_state(|state| {
        state
            .pending
            .iter()
            .find(|(_, posts)| posts.has_text())
            .map(|&(first, _)| first)
    });
    match waiting {
        Some(first) => Err(anyhow::anyhow!(
            "LOSAT internal error: the run ended ({}) with NCBI's posts of the search strategy request of the batch from query {first} not written",
            match &result {
                Ok(()) => "exit status 0".to_string(),
                Err(error) => error.to_string().trim_end().to_string(),
            }
        )),
        None => result,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
/// NCBI's posts of the bytes that ASN.1 text replaced in a request of a run of the application
/// `application` (`encoder::ReplacedChar`, with its prefix, `ArchiveRequest::application`): one
/// error post each, which the
/// application's handler also saves as a message of the batch's archive (outfmt 11) when it
/// writes one (`Posts::write`; the archive's own posts are not saved, `super::post_replaced`).
/// NCBI posts them while it writes the request, after it read the batch (whose reader messages
/// come first) and before the batch's search: blastn_app.cpp:277-282.
fn replaced_posts(
    application: &str,
    replaced: &[crate::report::serial::encoder::ReplacedChar],
) -> Posts {
    let mut posts = Posts::new();
    for byte in replaced {
        posts.post(ERROR, byte.post(Some(application)));
    }
    posts
}

impl State {
    // No NCBI counterpart: NCBI writes each request to its file before the batch's search
    // (blastp_app.cpp:259-262, cited at `StrategyExport::write_with`); LOSAT, which may search
    // later batches before an earlier batch's report ends the run, holds them for a sink it
    // cannot cut back.
    /// Writes to a sink that cannot be cut its held requests whose batches' first queries pass
    /// `through` (`State::held`), in order, then flushes it.
    fn write_held(&mut self, through: impl Fn(usize) -> bool) -> std::io::Result<()> {
        let Some(held) = self.held.as_mut() else {
            return Ok(());
        };
        let mut wrote = false;
        while held.front().is_some_and(|(at, _)| through(*at)) {
            if let Some((_, bytes)) = held.pop_front() {
                self.sink.write_all(&bytes)?;
                wrote = true;
            }
        }
        if wrote {
            self.sink.flush()?;
        }
        Ok(())
    }

    /// `StrategyExport::takes`.
    fn takes(&self, queries: &Range<usize>) -> bool {
        queries.end > self.next_query && queries.start < self.stop
    }
}

#[cfg(test)]
#[path = "strategy_tests.rs"]
mod tests;
