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
use super::writer::{encode_strategy_request, EncodedStrategyRequest, RequestInputs};
use super::{ArchiveHandle, ArchiveRequest, HandleCarry, OptionsHandle, SearchedContext0};
use crate::report::diag_log::{Posts, ERROR};
use crate::report::ostream::{emit_strategy_posts, StrategyStream};
use crate::report::run_files::{FileReceiver, RunFileName, RunFiles, StrategyTarget as Target};
use crate::report::serial::encoder::BufferPiece;

// No NCBI counterpart: Rust abstraction over `ostream::StrategyStream` (CLI disk/device sinks),
// `ReceivedFile` (LOSAT Web in-memory sink), and `NoSink` (refused/null sinks).
/// The stream of a strategy file: writes each encoded `Blast4-request` (with NCBI's object-stream
/// `BufferPiece` writes and `803.7` posts) as NCBI's `CNcbiOfstream` writes it
/// (`ostream::StrategyStream`), or into memory for LOSAT Web (`ReceivedFile`).
pub trait StrategySink: Send {
    /// Writes one complete `Blast4-request`, interleaving `posts` (`(written_before_post, line)`)
    /// with NCBI's buffer flushes and posting any write failure to `diagnostics` (and `archive`
    /// for outfmt 11).
    fn write_request(
        &mut self,
        bytes: &[u8],
        pieces: &[BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()>;

    /// Writes the partial bytes flushed before a `NullElement` throw (`encoder::NullElement`).
    fn write_partial(
        &mut self,
        bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()>;

    /// The run has ended (`end_run`): the file is complete. LOSAT Web's file goes to the run's
    /// receiver here (`ReceivedFile`); a file of the disk is the file system's already.
    fn hand_over(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// No NCBI counterpart: sink for refused targets (`-`) and LOSAT Web's `/dev/null` that discards
// request bytes while still emitting any `803.7` diagnostics.
/// The sink of a strategy target with no output file (`-`, refused before any request; or LOSAT
/// Web's `/dev/null`, which drops bytes while still posting any `803.7` diagnostics).
struct NoSink;

impl StrategySink for NoSink {
    fn write_request(
        &mut self,
        _bytes: &[u8],
        _pieces: &[BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        Ok(())
    }

    fn write_partial(
        &mut self,
        _bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        Ok(())
    }
}

impl StrategySink for StrategyStream {
    fn write_request(
        &mut self,
        bytes: &[u8],
        pieces: &[BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        self.write_object(bytes, pieces, posts, diagnostics, archive)
    }

    fn write_partial(
        &mut self,
        bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        StrategyStream::write_partial(self, bytes, posts, diagnostics, archive)
    }
}

/// The strategy file of a caller without files (LOSAT Web, `RunFiles::Receiver`): the file in
/// memory under the name that the command line writes it to, which the receiver gets whole where
/// the run ends (`end_run`, `hand_over`: its error is the run's internal error), or, for a run
/// that stops before its export is made, when the made file is dropped: the receiver gets the
/// file that the command line leaves then, empty (`create_strategy_file`). (No NCBI counterpart:
/// NCBI writes a file.)
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

impl StrategySink for ReceivedFile {
    fn write_request(
        &mut self,
        bytes: &[u8],
        _pieces: &[BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        debug_assert!(
            !self.handed,
            "a write of the strategy file after its run ended"
        );
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }

    fn write_partial(
        &mut self,
        bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        debug_assert!(
            !self.handed,
            "a write of the strategy file after its run ended"
        );
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        self.bytes.extend_from_slice(bytes);
        Ok(())
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

// No NCBI counterpart: tracks a batch's pending `Posts` and whether its request carried `803.7`
// diagnostics so `end_export` can detect untaken posts without storing a flag on `RequestMark`.
struct PendingBatch {
    first: usize,
    has_posts: bool,
    posts: Posts,
}

struct State {
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
    /// The request mark of each batch encoded (`write_with`), by the run's first query of the
    /// batch, until the run takes it (`take_posts`) to release the request where the batch starts.
    pending: Vec<PendingBatch>,
}

struct Shared {
    request: ArchiveRequest,
    /// The run's options handle, shared with the run's archive (`OptionsHandle`).
    options: OptionsHandle,
    /// LOSAT's refusal of the file (`create_strategy_file`, `create_strategy_file_in`,
    /// `check_report_stream`).
    refusal: Option<String>,
    /// The run's report stream (`report::ostream::ReportStream::of_cli`), whose writes of NCBI's
    /// before its first request `check_report_stream` makes: only for a file that LOSAT refuses.
    report: Mutex<Option<crate::report::ostream::RequestPoint>>,
    state: Mutex<State>,
    /// The strategy sink and its encoded requests held until the run reaches each batch's start
    /// (`release_through`). Its own lock, taken after `state`'s where both are, so that the panic
    /// hook can wait for it (`hook_lock`).
    held: Mutex<HeldSink>,
}

// No NCBI counterpart: holds an encoded `Blast4-request` (or partial bytes before a `NullElement`
// throw) until its batch's start is reached in report order; NCBI writes directly in its batch loop.
/// One encoded request waiting in `HeldSink`.
enum HeldRequest {
    Complete(EncodedStrategyRequest),
    Partial(Vec<u8>),
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-288
/// The strategy sink with the requests encoded for it, by the run's first query of each batch,
/// in order, until the run reaches each batch's start (`StrategyExport::release_through`) or ends
/// (`end_run`, or a panic: `install_panic_hook`). Releasing at each batch's start guarantees
/// every sink (shared `-out` / `/dev/stdout` files, failing devices such as `/dev/full` or
/// `RLIMIT_FSIZE`, FIFOs, pipes, terminals, and LOSAT Web's `ReceivedFile`) receives each request
/// at NCBI's exact `SaveSearchStrategy` point after prior batches' reports.
struct HeldSink {
    sink: Box<dyn StrategySink>,
    requests: std::collections::VecDeque<(usize, HeldRequest)>,
}

impl HeldSink {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-288
    /// Writes the held requests whose batches' first queries pass `through`, in order.
    fn write_through(
        &mut self,
        through: impl Fn(usize) -> bool,
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        while self.requests.front().is_some_and(|(at, _)| through(*at)) {
            if let Some((_, request)) = self.requests.pop_front() {
                match request {
                    HeldRequest::Complete(encoded) => self.sink.write_request(
                        &encoded.bytes,
                        &encoded.pieces,
                        &encoded.posts,
                        diagnostics,
                        archive,
                    )?,
                    HeldRequest::Partial(bytes) => {
                        self.sink.write_partial(&bytes, &[], diagnostics, archive)?;
                    }
                }
            }
        }
        Ok(())
    }
}

impl Shared {
    // No NCBI counterpart: acquires the panic-hook-aware mutex around `HeldSink::write_through`.
    /// `HeldSink::write_through` under the panic-hook-aware lock.
    fn write_held(
        &self,
        through: impl Fn(usize) -> bool,
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        hook_lock(&self.held).write_through(through, diagnostics, archive)
    }
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:777-798
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
/// Creates (and empties) the strategy file where NCBI's `SetOptions` does: after the argument
/// parser accepted the command line and before any option handler runs, so before the
/// `-outfmt` text is read, the database list is resolved, and `-subject`, `-query` and `-out`
/// are opened. A file that cannot be created is NCBI's argument error (exit status 1).
///
/// Any non-`-` file or device path (including `/dev/stdout`, the same file as `-out`, symlinks,
/// and failing devices such as `/dev/full`) is opened here on its own descriptor and written
/// through `ostream::StrategyStream` at each batch's start. LOSAT still refuses `-` (standard
/// output interleaved into NCBI's report stream, reserved for its own follow-up): the refusal is
/// recorded here and raised where NCBI makes its first request (`check_report_stream`), so every
/// stop before that point ends as NCBI's.
pub fn create_strategy_file(path: &Path, label: &str) -> Result<StrategyFile> {
    if path.as_os_str() == "-" {
        return Ok(standard_output_refusal(label));
    }
    let file = std::fs::File::create(path)
        .map_err(|_| crate::cli::inaccessible("export_search_strategy", path))?;
    Ok(StrategyFile {
        target: StrategyTarget::Disk(file),
        refusal: None,
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

/// LOSAT Web's refusal of a strategy file that is the output's (`create_strategy_file_in`,
/// `Target::StandardOutput`).
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
/// (`/dev/stdout`) refused at NCBI's first request; no file for
/// the null device (`/dev/null`, `/dev/stdin`). LOSAT Web refuses, at the same point, a name whose
/// command line run writes where the run has no file: the standard error, a device whose writes
/// fail or wait for a reader, a file of the system (procfs, a message queue), or the file of
/// another option of the run (`others`: `-query`, `-subject`, `-logfile`, `-window_masker_db`,
/// the files that the search reads from the run's folder, `web_adapter::open_strategy_file`),
/// whose bytes the command line's strategy file replaces or mixes with its own.
pub fn create_strategy_file_in(
    files: &RunFiles,
    path: &Path,
    label: &str,
    others: &[RunFileName<'_>],
) -> Result<StrategyFile> {
    let RunFiles::Receiver(receiver) = files else {
        return create_strategy_file(path, label);
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
/// LOSAT's refusal of the run's strategy file (`create_strategy_file`: `-` on the command line;
/// `create_strategy_file_in`: `-`, standard output, or non-run files in LOSAT Web), raised where
/// NCBI writes its first request: after its formatter and its prolog, once the first batch is read
/// and before it is searched (`SaveSearchStrategy`).
/// Every stop of NCBI's before that point (an empty query, a database or option error, a query
/// that cannot be read, a bad `BATCH_SIZE`, psiblast's `-in_pssm` with `-subject`, a device
/// that fails at the prolog's flush) comes first with NCBI's own bytes, because the run's report
/// stream is NCBI's until here. Here the stream makes NCBI's writes so far while the formatter
/// lives (`RequestPoint::make_writes`): a write that fails ends the run as NCBI's failed write
/// ends it (`BLAST failed to write output`, exit code 6, on `/dev/full`), before any request.
/// Otherwise the run refuses; what the stream holds then is written before the refusal, as at
/// LOSAT's other refusals (NCBI's prolog, which NCBI writes before its request too). The
/// programs call it where they make the first batch's request (`StrategyExport::write_with`
/// calls it too, so no request is ever written to such a file) and, where they write that
/// request after searching, at NCBI's point (BLASTP's `search`).
pub fn check_report_stream(strategy: Option<&StrategyExport>) -> Result<()> {
    let Some(strategy) = strategy else {
        return Ok(());
    };
    match strategy.0.refusal.as_deref() {
        Some(refusal) => {
            let report = strategy
                .0
                .report
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .take();
            if let Some(report) = report {
                report.make_writes()?;
            }
            Err(anyhow::anyhow!("{refusal}"))
        }
        None => Ok(()),
    }
}

impl StrategyExport {
    /// Whether LOSAT refuses the file (`check_report_stream`). NCBI writes its prolog before its
    /// first request, so a program that writes its prolog with its first reports writes it before
    /// that refusal: a device that fails at the prolog's flush then ends the run as NCBI's.
    pub fn refused(&self) -> bool {
        self.0.refusal.is_some()
    }

    /// The run's report stream (`report::ostream::ReportStream::of_cli`), whose writes before
    /// NCBI's first request `check_report_stream` makes where LOSAT refuses the file.
    pub(crate) fn attach_report(&self, report: crate::report::ostream::RequestPoint) {
        if self.0.refusal.is_some() {
            *self
                .0
                .report
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(report);
        }
    }
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
/// UTF-8 is refused first (`check_utf8_file_name`). `label` names the program in a refusal.
///
/// With a strategy file NCBI skips its early reading of `-outfmt` (blast_args.cpp:3624-3627,
/// cited at `create_strategy_file`'s callers): a value its parser refuses (`99`, `abc`) fails
/// only in the formatting handler, the last of the handlers, so `-subject`, `-subject_loc`,
/// `-query`, `-out` and `-query_loc` come first, in the order of the other arguments. The
/// callers read the value with `app::early_format_choice`, which leaves it to their own
/// formatting handler, and make `-out` where they always do (line RG F/d36 review 2c B1).
pub fn open_strategy_file(path: Option<&Path>, label: &str) -> Result<Option<StrategyFile>> {
    open_strategy_file_in(&RunFiles::Disk, path, label, &[])
}

/// `open_strategy_file` of a run whose files go to `files` (`create_strategy_file_in`), whose
/// other files are `others`: LOSAT Web's adapter opens the file here, where the command line does.
pub fn open_strategy_file_in(
    files: &RunFiles,
    path: Option<&Path>,
    label: &str,
    others: &[RunFileName<'_>],
) -> Result<Option<StrategyFile>> {
    let Some(path) = path else {
        return Ok(None);
    };
    crate::blastinput::input_files::check_utf8_file_name(path, "export_search_strategy", label)?;
    create_strategy_file_in(files, path, label, others).map(Some)
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
        let requests = std::collections::VecDeque::new();
        let export = Self(Arc::new(Shared {
            request: ArchiveRequest::for_inputs(inputs, database),
            options,
            report: Mutex::new(None),
            refusal,
            state: Mutex::new(State {
                carry: HandleCarry::default(),
                next_query: 0,
                noted_through: 0,
                stop: usize::MAX,
                pending: Vec::new(),
            }),
            held: Mutex::new(HeldSink { sink, requests }),
        }));
        let mut open = hook_lock(&HOLDING);
        open.retain(|shared| shared.strong_count() > 0);
        open.push(Arc::downgrade(&export.0));
        export
    }

    /// The export of a run into its strategy file (`open_strategy_file`), buffered through
    /// `ostream::StrategyStream`: the program's recorded options `inputs` with the `HitlistSize`
    /// of the run's format `outfmt` and `counts` (the archive's is outfmt 11's), and the run's
    /// options handle `options` (made once per run and given to its archive too,
    /// `ArchiveHandle::for_program`). `outfmt` is none for a `-outfmt` value that was not read
    /// (`app::EarlyFormat::Deferred`): the run stops at its formatting handler before any request,
    /// so the options keep their own hitlist size.
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
            StrategyTarget::Disk(file) => Box::new(StrategyStream::new(
                Box::new(file),
                inputs.kind.application(),
            )),
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
    /// run's end). Every sink therefore holds encoded requests in `Shared::held` and releases each
    /// batch at its start mark (`take_posts` / `Posts::write` -> `release_through`), interleaving
    /// `803.7` posts and write-failure diagnostics with NCBI's buffer flushes; `end_run` drops
    /// unreached requests when NCBI's error ends the run.
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
            let mut posts = Posts::new();
            let mut has_posts = false;
            let mut made_one = false;
            let made = make(&mut |request| {
                let sets = shared.options.request_sets(&state.carry);
                // A big integer of -2^63 (RPS-BLAST's `DbLength` that a batch before set,
                // `-dbsize`), which NCBI writes as `-0`, is written as its digits
                // (PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits).
                let encoded = encode_strategy_request(
                    &shared.request,
                    request,
                    &sets,
                    &shared.options.strategy_database_params(),
                )?;
                state.next_query = queries.end;
                made_one = true;
                has_posts |= !encoded.posts.is_empty();
                // NCBI's stream threw while it wrote this request (an element stated as `NULL` in a
                // PSSM query, `encoder::NullElement`): NCBI's error ends the run there, after the
                // bytes its stream flushed, so those partial bytes go to the sink now.
                if let Some(thrown) = encoded.thrown {
                    for (_, text) in encoded.posts {
                        posts.post(ERROR, text);
                    }
                    hook_lock(&shared.held)
                        .requests
                        .push_back((queries.start, HeldRequest::Partial(encoded.bytes)));
                    shared.write_held(|at| at <= queries.start, &mut std::io::sink(), None)?;
                    return Err(thrown.into());
                }
                posts.append(&Posts::request(self.clone(), queries.start));
                hook_lock(&shared.held)
                    .requests
                    .push_back((queries.start, HeldRequest::Complete(encoded)));
                Ok(())
            });
            if made_one {
                state.pending.push(PendingBatch {
                    first: queries.start,
                    has_posts: has_posts || posts.has_text(),
                    posts,
                });
            }
            made?;
            Ok(true)
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-288
    /// The run reaches the start of the batch whose first query (the run's) is `first`
    /// (`Posts::write` of its request mark), after the reports of the batches before it: releases
    /// the encoded requests through `first` to the strategy stream, interleaving any `803.7`
    /// posts and strategy write-failure diagnostics on `diagnostics` (and `archive` for outfmt
    /// 11).
    pub(crate) fn release_through(
        &self,
        first: usize,
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        self.0.write_held(|at| at <= first, diagnostics, archive)
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
            let taken = state.pending.partition_point(|batch| batch.first < before);
            state
                .pending
                .drain(..taken)
                .map(|batch| (batch.first, batch.posts))
                .collect()
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

// No NCBI counterpart: NCBI's process does not panic; a panic is LOSAT's own failure.
/// The exports whose sinks cannot be cut (`Shared::held`: their requests wait in memory until the
/// run reaches their batches), for the panic hook (`install_panic_hook`).
static HOLDING: Mutex<Vec<std::sync::Weak<Shared>>> = Mutex::new(Vec::new());

thread_local! {
    /// The locks that the panic hook takes (`hook_lock`) that this thread holds.
    static HOOK_LOCKS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// A lock that the panic hook takes too (`HOLDING`, `Shared::held`), counted for its thread: a
/// thread that panics while it holds one does not wait for it in the hook, which would never end
/// (`write_held_at_panic`). Each is held only for steps that make no request: `HOLDING` to add
/// an export, `Shared::held` to add a request made or to write requests to the sink, which takes
/// as long as the sink takes them (a FIFO whose reader stops reading holds it).
struct HookLock<'a, T>(std::sync::MutexGuard<'a, T>);

/// `HookLock` of `mutex`, which waits for it.
fn hook_lock<T>(mutex: &Mutex<T>) -> HookLock<'_, T> {
    let guard = mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    HOOK_LOCKS.set(HOOK_LOCKS.get() + 1);
    HookLock(guard)
}

impl<T> std::ops::Deref for HookLock<'_, T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.0
    }
}

impl<T> std::ops::DerefMut for HookLock<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

impl<T> Drop for HookLock<'_, T> {
    fn drop(&mut self) {
        HOOK_LOCKS.set(HOOK_LOCKS.get() - 1);
    }
}

// No NCBI counterpart (see `HOLDING`).
/// Makes a panic of the command line end each strategy file as LOSAT's refusals end it
/// (`end_run`): the requests held for a sink that cannot be cut (a FIFO, a pipe) are written, so
/// it gets every request made, as a regular file has them (each request is written and flushed to
/// such a file when it is made). The release profile ends the process at a panic (`panic =
/// "abort"`), so no `end_run` and no destructor runs after it, only the panic hook. A build that
/// unwinds may catch a panic and go on (the first phase of `blastx`'s batches, which raises it
/// again where the unit is written): there the hook could write requests of a run that then ends
/// with NCBI's error, so it is installed only where a panic ends the process. The hook writes the
/// panic's message first (the hook before it, Rust's own by default), then the requests; it waits
/// for a thread that is adding or writing requests (`hook_lock`), so a FIFO whose reader stops
/// reading stops it as it stops the run's own writes. The command line's entry
/// (`cli::main_entry`) calls it once, before the program runs.
pub fn install_panic_hook() {
    if !cfg!(panic = "abort") {
        return;
    }
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        previous(info);
        write_held_at_panic();
    }));
}

/// The panic hook's part (`install_panic_hook`): the exports of `HOLDING` (`write_held_of`),
/// unless the panicking thread holds one of the hook's locks (a panic inside such a step: the
/// requests of every export are left as they are).
fn write_held_at_panic() {
    if HOOK_LOCKS.get() > 0 {
        return;
    }
    let open = hook_lock(&HOLDING);
    write_held_of(&open);
}

/// Writes every held request of each live export of `open`, after the steps of other threads on
/// its sink (`hook_lock`); a request still being made (`StrategyExport::write_with`) is not held
/// yet.
fn write_held_of(open: &[std::sync::Weak<Shared>]) {
    for shared in open.iter().filter_map(std::sync::Weak::upgrade) {
        let _ = shared.write_held(|_| true, &mut std::io::sink(), None);
    }
}

// No NCBI counterpart: NCBI posts while it writes the request, so no post can wait there.
/// The end of a run of a program with `-export_search_strategy` (`strategy`; `None` without
/// it). The requests encoded (`StrategyExport::write_with`) whose batches' starts the run did not
/// reach: a run that ends with one of NCBI's errors (`cli::NativeError::caught`, including a
/// `StrategyWriteFault`) ended before those starts, so those unreached requests are dropped; a
/// run that ends with `Ok` (or a refusal of LOSAT's) flushes any remaining held requests. A run
/// that ends as NCBI's (`Ok` or NCBI's error) has written every `803.7` post of its requests at
/// its batch's start (`StrategyExport::take_posts`), so a post still waiting is LOSAT's internal
/// error instead of output that silently lacks it. LOSAT's refusals (any other error) keep their
/// own text. Each program calls this once, where its run's result is known (`blastn`, `blastp`
/// for psiblast and rpsblast too, `blastx` for rpstblastn too, `tblastn`, `tblastx`; LOSAT Web's
/// records searches, `ReportOutputs::run_records`). The posts still waiting go then: each holds
/// the export (its request's mark), so that the export would outlive the run, and a sink that
/// ends with the export (LOSAT Web's file of the run, `ReceivedFile`) would never reach its
/// receiver.
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
    let handed = hook_lock(&strategy.0.held).sink.hand_over();
    let refused = matches!(&ended, Err(error) if crate::cli::NativeError::caught(error).is_none());
    match handed {
        Err(error) if !refused => {
            Err(anyhow::anyhow!(
                "LOSAT internal error: the -export_search_strategy file did not reach the run's files ({error})"
            ))
        }
        _ => ended,
    }
}

/// `end_run` of a run with the export `strategy`, before its waiting posts go.
fn end_export(strategy: &StrategyExport, result: Result<()>) -> Result<()> {
    // NCBI's error (`cli::NativeError::caught`, which the command line's exit reads too): a
    // `NativeError`, a strategy write fault, or an exception of NCBI's that a report's writer
    // carries (the `-html` links' registry, which NCBI's prolog reads before its first request).
    let native = matches!(&result, Err(error) if crate::cli::NativeError::caught(error).is_some());
    let ended = if native {
        hook_lock(&strategy.0.held).requests.clear();
        Ok(())
    } else {
        strategy.0.write_held(|_| true, &mut std::io::sink(), None)
    };
    if let Err(error) = ended {
        // A refusal of LOSAT's keeps its own text.
        if result.is_err() && !native {
            return result;
        }
        if let Some(native) =
            crate::cli::NativeError::of_writer_exception::<dyn std::error::Error>(&error)
        {
            return Err(native.into());
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
            .find(|batch| batch.has_posts)
            .map(|batch| batch.first)
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

impl State {
    /// `StrategyExport::takes`.
    fn takes(&self, queries: &Range<usize>) -> bool {
        queries.end > self.next_query && queries.start < self.stop
    }
}

#[cfg(test)]
#[path = "strategy_tests.rs"]
mod tests;
