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

struct State {
    /// The sink of a regular file (`StrategySink::can_cut`), which gets each request where it is
    /// made and is cut back (`end_run`); `None` for a sink that cannot be cut, whose requests wait
    /// in `Shared::held`.
    sink: Option<Box<dyn StrategySink>>,
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
    /// The run's report stream (`report::ostream::ReportStream::of_cli`), whose writes of NCBI's
    /// before its first request `check_report_stream` makes: only for a file that LOSAT refuses.
    report: Mutex<Option<crate::report::ostream::RequestPoint>>,
    state: Mutex<State>,
    /// A sink that cannot be cut, with its requests held (`None` for a regular file, whose sink is
    /// `State::sink`). Its own lock, taken after `state`'s where both are, so that the panic hook
    /// can wait for it (`hook_lock`).
    held: Option<Mutex<HeldSink>>,
}

// No NCBI counterpart: NCBI writes each request to its file before the batch's search
// (blastp_app.cpp:259-262, cited at `StrategyExport::write_with`); LOSAT, which may search later
// batches before an earlier batch's report ends the run, holds them for a sink it cannot cut back.
/// A sink that cannot be cut (`StrategySink::can_cut`: a FIFO, a pipe, a terminal) with the
/// requests made for it, by the run's first query of each one's batch, in order, until the run
/// reaches their batches' starts (`StrategyExport::release_through`) or ends (`end_run`, or a
/// panic: `install_panic_hook`): such a sink gets only the requests NCBI writes. Its lock is held
/// only to add a request made (outside the lock) or to write and drop requests, never while a
/// request is made.
struct HeldSink {
    sink: Box<dyn StrategySink>,
    requests: std::collections::VecDeque<(usize, Vec<u8>)>,
}

impl HeldSink {
    /// Writes the held requests whose batches' first queries pass `through`, in order, then
    /// flushes the sink.
    fn write_through(&mut self, through: impl Fn(usize) -> bool) -> std::io::Result<()> {
        let mut wrote = false;
        while self.requests.front().is_some_and(|(at, _)| through(*at)) {
            if let Some((_, bytes)) = self.requests.pop_front() {
                self.sink.write_all(&bytes)?;
                wrote = true;
            }
        }
        if wrote {
            self.sink.flush()?;
        }
        Ok(())
    }
}

impl Shared {
    /// `HeldSink::write_through` of a sink that cannot be cut (nothing for a regular file).
    fn write_held(&self, through: impl Fn(usize) -> bool) -> std::io::Result<()> {
        match &self.held {
            Some(held) => hook_lock(held).write_through(through),
            None => Ok(()),
        }
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
/// A strategy file made where NCBI's `SetOptions` makes it (`create_strategy_file`): the file
/// (`None` for `-`, NCBI's standard output) and LOSAT's refusal of it, if any, which the run
/// raises where NCBI writes its first request (`check_report_stream`).
#[derive(Debug)]
pub struct StrategyFile {
    file: Option<std::fs::File>,
    refusal: Option<String>,
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
/// until then the run goes on as NCBI's.
pub fn create_strategy_file(path: &Path, out: Option<&Path>, label: &str) -> Result<StrategyFile> {
    if path.as_os_str() == "-" {
        return Ok(StrategyFile {
            file: None,
            refusal: Some(format!(
                "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's {label}; give a file name"
            )),
        });
    }
    let file = std::fs::File::create(path)
        .map_err(|_| crate::cli::inaccessible("export_search_strategy", path))?;
    let refusal = is_same_as_output(&file, path, out).then(|| {
        format!(
            "an -export_search_strategy file that is also the output (-out or the standard output) is not supported by LOSAT's {label}: NCBI writes both through separate streams and the file's bytes follow the order of their flushes"
        )
    });
    Ok(StrategyFile {
        file: Some(file),
        refusal,
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:250-262 (the same order in
// blastn_app.cpp:255-282, blastx_app.cpp:251-264, tblastn_app.cpp:267-280, tblastx_app.cpp:170-181,
// psiblast_app.cpp:503-517 and :340, rpsblast_app.cpp:188-199, rpstblastn_app.cpp:185-196)
/// LOSAT's refusal of the run's strategy file (`create_strategy_file`: the standard output, or
/// the file of the output), raised where NCBI writes its first request: after its formatter and
/// its prolog, once the first batch is read and before it is searched (`SaveSearchStrategy`).
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
    let Some(path) = path else {
        return Ok(None);
    };
    crate::blastinput::input_files::check_utf8_file_name(path, "export_search_strategy", label)?;
    create_strategy_file(path, out, label).map(Some)
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
        let can_cut = sink.can_cut();
        let (sink, held) = if can_cut {
            (Some(sink), None)
        } else {
            let requests = std::collections::VecDeque::new();
            (None, Some(Mutex::new(HeldSink { sink, requests })))
        };
        let export = Self(Arc::new(Shared {
            request: ArchiveRequest::for_inputs(inputs, database),
            options,
            report: Mutex::new(None),
            refusal,
            state: Mutex::new(State {
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
            held,
        }));
        if !can_cut {
            let mut open = hook_lock(&HOLDING);
            open.retain(|shared| shared.strong_count() > 0);
            open.push(Arc::downgrade(&export.0));
        }
        export
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
        let sink: Box<dyn StrategySink> = match file.file {
            Some(file) => Box::new(std::io::BufWriter::new(file)),
            None => Box::new(NoSink),
        };
        Self::with_refusal(sink, &inputs, database, options.clone(), file.refusal)
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
    /// run's end). So the run marks where it reaches each batch's start, with NCBI's posts of the
    /// bytes that ASN.1 text replaced in the batch's request (`replaced_posts`): the request mark
    /// that `take_posts` hands out with them (`Posts::write`, `release_through`); `end_run` cuts the
    /// file back to the requests reached when NCBI's error ends the run. (Holding the requests'
    /// bytes until then instead would keep every request of a BLASTN XML run in memory: 2.6 MB on a
    /// run of 82 MB, review 1 of K4h gap 2; so only a sink that cannot be cut holds them,
    /// `Shared::held`.)
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
                    inner: match state.sink.as_mut() {
                        Some(sink) => &mut **sink as &mut dyn Write,
                        None => &mut held as &mut dyn Write,
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
                let wrote = counted.written;
                let flushed = match state.sink.as_mut() {
                    Some(sink) => {
                        state.written += wrote;
                        sink.flush()
                    }
                    // A sink that cannot be cut is `Shared::held` (`with_refusal` makes one of
                    // the two from the sink).
                    None => match &shared.held {
                        Some(waiting) => {
                            hook_lock(waiting).requests.push_back((queries.start, held));
                            Ok(())
                        }
                        None => Err(std::io::Error::other(
                            "LOSAT internal error: the search strategy export has no sink",
                        )),
                    },
                };
                state.unreached.push_back((queries.start, state.written));
                state.next_query = queries.end;
                made_one = true;
                // NCBI's stream threw while it wrote this request (an element stated as `NULL` in a
                // PSSM query, `encoder::NullElement`): NCBI's error ends the run there, after the
                // bytes its stream flushed, so the request is reached.
                if matches!(&written, Err(error) if error.is::<crate::cli::NativeError>()) {
                    flushed?;
                    state.reach_through(queries.start);
                    shared.write_held(|at| at <= queries.start)?;
                    return written;
                }
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
            state.reach_through(first);
            self.0.write_held(|at| at <= first)
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
        let _ = shared.write_held(|_| true);
    }
}

// No NCBI counterpart: NCBI posts while it writes the request, so no post can wait there.
/// The end of a run of a program with `-export_search_strategy` (`strategy`; `None` without
/// it). The requests written (`StrategyExport::write_with`) whose batches' starts the run did not
/// reach: a run that ends with one of NCBI's errors (`cli::NativeError`) ended before those starts,
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
/// rpstblastn too, `tblastn`, `tblastx`).
pub fn end_run(strategy: Option<&StrategyExport>, result: Result<()>) -> Result<()> {
    let Some(strategy) = strategy else {
        return result;
    };
    let native = matches!(&result, Err(error) if error.is::<crate::cli::NativeError>());
    let ended = strategy.with_state(|state| -> std::io::Result<()> {
        state.unreached.clear();
        if !native {
            return strategy.0.write_held(|_| true);
        }
        if let Some(held) = &strategy.0.held {
            hook_lock(held).requests.clear();
        }
        if let Some(sink) = state.sink.as_mut() {
            if state.written > state.reached {
                sink.truncate(state.reached)?;
                state.written = state.reached;
            }
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
    /// The run reaches the batches whose first queries are at most `first`
    /// (`StrategyExport::release_through`).
    fn reach_through(&mut self, first: usize) {
        while let Some(&(at, end)) = self.unreached.front() {
            if at > first {
                break;
            }
            self.reached = end;
            self.unreached.pop_front();
        }
    }

    /// `StrategyExport::takes`.
    fn takes(&self, queries: &Range<usize>) -> bool {
        queries.end > self.next_query && queries.start < self.stop
    }
}

#[cfg(test)]
#[path = "strategy_tests.rs"]
mod tests;
