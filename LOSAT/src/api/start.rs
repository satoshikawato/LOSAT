//! The NCBI application start of a library call (line RG item E1w-2).
//!
//! NCBI's application starts before it reads its arguments: it opens the file of `-logfile`,
//! loads its configuration (`<program>.ini` or the file of `-conffile`, then `.ncbirc`) and,
//! after the arguments, ends a `-dryrun` (`ncbi_environment::ApplicationStart`, cited there).
//! The command line's entry runs one start per process (`cli::main_entry`); the library runs one
//! per search call, at its beginning (`CallStart::new`), so that a call refuses what the command
//! line refuses there (a usage report setting that NCBI aborts on, a registry entry that LOSAT
//! does not read), whatever its inputs and options.

use std::ffi::OsStr;
use std::io::Write;

use super::error::Error;
use super::input::Inputs;
use super::search::Search;
use crate::blastinput::ncbi_environment::{
    log_file, run_messages_follow_log, ApplicationSettings, ApplicationStart, DiagTarget, TextPlace,
};
use crate::cli::Commands;
use crate::report::run_files::RunFiles;

/// A start whose `-logfile` file is open (`open`), before its configuration load (`load`).
pub(super) struct Opened {
    pub(super) start: ApplicationStart,
    pub(super) target: DiagTarget,
}

/// The application start of `program` with the arguments `words` (after the program's name), of
/// this process (its environment and its registry files), up to the configuration load: the
/// usage report's messages go to `diagnostics` (NCBI writes them to `cerr` at once), the file of
/// `-logfile` is opened (appended to) and gets them (`ApplicationStart::write_log_setup`). `Err`
/// with LOSAT's refusal of the start.
pub(super) fn open<S: AsRef<OsStr> + Clone>(
    program: &str,
    words: &[S],
    diagnostics: &mut dyn Write,
) -> Result<Opened, Error> {
    open_with(
        ApplicationStart::of_process(program, words),
        words,
        &RunFiles::Disk,
        diagnostics,
    )
}

/// `open` with the start `start` (made from the arguments `words`), whose file of `-logfile` goes
/// to `files` (`RunFiles::log_target`).
fn open_with<S: AsRef<OsStr> + Clone>(
    start: ApplicationStart,
    words: &[S],
    files: &RunFiles,
    diagnostics: &mut dyn Write,
) -> Result<Opened, Error> {
    match start.usage_report() {
        Ok(text) => {
            let _ = diagnostics.write_all(text.as_bytes());
        }
        Err(reason) => return Err(Error::Search(reason.into())),
    }
    let mut target = files.log_target(log_file(words));
    start.write_log_setup(&mut target, diagnostics);
    Ok(Opened { start, target })
}

impl Opened {
    /// The configuration load (`ApplicationStart::write_config_load`): its messages go where the
    /// diagnostics go, the file of `-logfile` or `diagnostics`. A load that ends the program
    /// (exit status 2 or 255 on the command line) is an error whose message is what the command
    /// line writes for it without a `-logfile` (the file has its part, as there); `Err` with
    /// LOSAT's refusal too.
    pub(super) fn load(mut self, diagnostics: &mut dyn Write) -> Result<Started, Error> {
        self.write_load(diagnostics)?;
        Ok(Started {
            settings: self.start.settings(),
            target: self.target,
        })
    }

    /// `load`, which leaves the start and the target here.
    fn write_load(&mut self, diagnostics: &mut dyn Write) -> Result<(), Error> {
        let mut stderr = Vec::new();
        match self.start.write_config_load(&mut self.target, &mut stderr) {
            Ok(None) => {
                let _ = diagnostics.write_all(&stderr);
                Ok(())
            }
            Ok(Some(_status)) => {
                let mut text = Vec::new();
                let _ = self
                    .start
                    .write_config_load(&mut DiagTarget::Stderr, &mut text);
                Err(Error::written(text))
            }
            Err(reason) => Err(Error::Search(reason.into())),
        }
    }
}

/// What a library call's application start left.
pub(super) struct Started {
    /// The settings, or LOSAT's refusal of them (`ApplicationStart::settings`).
    pub(super) settings: Result<ApplicationSettings, String>,
    /// Where the diagnostics go: the file of `-logfile`, `cerr` (the call's diagnostics) or
    /// nowhere.
    pub(super) target: DiagTarget,
}

/// The whole start (`open`, then `Opened::load`).
pub(super) fn begin<S: AsRef<OsStr> + Clone>(
    program: &str,
    words: &[S],
    diagnostics: &mut dyn Write,
) -> Result<Started, Error> {
    open(program, words, diagnostics)?.load(diagnostics)
}

/// The application start of one search call of the library, made at the call's beginning. The
/// searches of a call share it (two of `run`'s groups, the batches of `run_batches`).
pub(crate) struct CallStart {
    /// The start's settings (`ApplicationStart::settings`), which LOSAT does not refuse.
    settings: ApplicationSettings,
    /// The `-logfile` target of a program whose run's messages go there
    /// (`run_messages_follow_log`: tblastx), where it is not `cerr`: the file of `-logfile`, or
    /// nowhere. With `cerr` the call's diagnostics take the run's messages, as every other
    /// program's.
    log: Option<DiagTarget>,
    /// Whether the argv has `-dryrun`.
    dry_run: bool,
    /// The call's argv as the command line's parser read it with the call's kind of inputs
    /// (`Search::command`), once per call as the command line parses once per run: the call's
    /// steps take it from here (`CallStart::command`) instead of parsing the argv again.
    /// `AssertUnwindSafe` keeps `Batches`, which holds the start, unwind safe as it was before
    /// the start kept the parse: the only fields of `Commands` that are not are LOSAT Web's
    /// hooks (`BlastnArgs::window_masker_counts`, `BlastxArgs::rpstblastn`), which the parser
    /// skips, so a parse holds `None` there, and the start only reads its parse or clones it.
    command: std::panic::AssertUnwindSafe<Commands>,
    /// The diagnostic stream of NCBI's main thread of the call's run (`ncbi_exception::MainDiag`),
    /// which the call's code reads where it runs (`CallStart::enter`).
    diag: std::sync::Arc<crate::ncbi_exception::MainDiag>,
}

impl CallStart {
    /// The calling thread runs this call until the guard is dropped: the call's posts read its
    /// own diagnostic stream (`ncbi_exception::main_diag_failbit`), not another call's.
    pub(super) fn enter(&self) -> crate::ncbi_exception::MainDiagGuard {
        crate::ncbi_exception::MainDiag::enter(Some(self.diag.clone()))
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (cited at
    // `cli::main_entry`: the BLAST programs keep the default `DryRun`, which returns 0 after the
    // arguments have passed their definitions)
    /// The start of a call of `search` with `inputs`, in the command line's order
    /// (`cli::main_entry`): the start, with its messages on `diagnostics`; the parse of the argv
    /// with the kind of the call's subject (`Search::command`: `-subject_loc` with a database fails
    /// here, dry run or not; it reads no record); then the settings that LOSAT refuses
    /// (`ApplicationStart::settings`), the call's error, as the command line refuses them after its
    /// parse (`start.settings()?`).
    pub(super) fn new(
        search: &Search,
        inputs: Inputs<'_>,
        diagnostics: &mut dyn Write,
    ) -> Result<CallStart, Error> {
        // The call is a run of its own: NCBI's application of it starts with a diagnostic stream
        // that no other run in the process fails or reads.
        let diag: std::sync::Arc<crate::ncbi_exception::MainDiag> = Default::default();
        let _run = crate::ncbi_exception::MainDiag::enter(Some(diag.clone()));
        let program = search.program().name();
        let started = begin(program, search.argv(), diagnostics)?;
        let command = search.command(inputs, &started.target)?;
        let settings = started
            .settings
            .map_err(|reason| Error::Search(reason.into()))?;
        Ok(CallStart {
            settings,
            log: run_messages_follow_log(program)
                .then_some(started.target)
                .filter(|target| !matches!(target, DiagTarget::Stderr)),
            dry_run: search.dry_run(),
            command: std::panic::AssertUnwindSafe(command),
            diag,
        })
    }

    // No NCBI counterpart: NCBI's application parses its arguments once per run (cited at
    // `cli::main_entry`); a library call's steps share the one parse of its start instead of
    // parsing the argv again (the parse is deterministic: the same argv and inputs give the
    // same arguments).
    /// The call's parsed argv (`Search::command` with the call's inputs, read at the call's
    /// start), for a step that reads it.
    pub(super) fn command_ref(&self) -> &Commands {
        &self.command.0
    }

    /// The call's parsed argv, for a step that takes its own (each group of `run`, each batch).
    pub(super) fn command(&self) -> Commands {
        self.command.0.clone()
    }

    /// Whether the call is a dry run, which ends here: the start, the parse with the call's kind
    /// of subject and the settings have passed (`CallStart::new`). Neither the records nor
    /// `-outfmt` are read (NCBI reads them in `Run`).
    pub(super) fn dry_run(&self) -> bool {
        self.dry_run
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:345-347
    /// The application settings of the call's start. The records of a library call are read before
    /// its start, with NCBI's data loaders (`read_fasta`, `read_inputs`); settings that turn them
    /// off (`[BLAST] DATA_LOADERS`) read the same records: the loaders only make NCBI's reader take
    /// a record's first line as a Seq-id, which LOSAT's reader refuses
    /// (`fasta_reader::FastaReader::read_one_seq`), so a record that it read is read the same
    /// without them.
    pub(crate) fn settings(&self) -> ApplicationSettings {
        self.settings.clone()
    }

    /// The `-logfile` target of the call's run messages (tblastx's), which the caller holds for
    /// the run (`RunLog`).
    pub(super) fn take_log(&mut self) -> RunLog {
        RunLog(self.log.take())
    }
}

/// The `-logfile` target of a run whose messages go there (`CallStart::take_log`: tblastx's), or
/// none.
pub(super) struct RunLog(Option<DiagTarget>);

impl RunLog {
    /// The writer of one search's messages. With `to_log` they go where the run's messages go
    /// (the `-logfile` target; `cerr` is `diagnostics`); without, to `diagnostics` alone: the
    /// messages of a later search of the call, which repeat the first's and which the command line
    /// writes once (`run`, `run_collect`), or messages gathered to be cut and then written with
    /// `RunLog::route` (`run_batches`). A failure goes to the target either way
    /// (`RunSink::failed`).
    pub(super) fn sink<'a>(
        &'a mut self,
        diagnostics: &'a mut (dyn Write + Send),
        to_log: bool,
    ) -> RunSink<'a> {
        RunSink {
            log: self,
            diagnostics,
            to_log,
        }
    }

    /// Where the messages of a sink with `to_log` end up (`RunSink::write`): the target's place
    /// (`DiagTarget::place`), or the caller's diagnostics without a target.
    pub(super) fn place(&self) -> TextPlace {
        self.0
            .as_ref()
            .map_or(TextPlace::StandardError, DiagTarget::place)
    }

    /// Writes `messages`, gathered from a search, where the run's messages go, and gives back what
    /// goes to the call's diagnostics: all of it without a target or for `cerr`, nothing for the
    /// file of `-logfile` or no output.
    pub(super) fn route(&mut self, messages: Vec<u8>) -> Vec<u8> {
        match &mut self.0 {
            Some(target) => {
                let mut rest = Vec::new();
                target.write_text(&messages, &mut rest);
                rest
            }
            None => messages,
        }
    }

    /// The error of a failed run (`Error::search`), posted (`RunLog::post`).
    pub(super) fn failed(&mut self, error: &anyhow::Error) -> Error {
        self.post(error);
        Error::search(error)
    }

    /// Writes the text of `error`, a `NativeError` of the run, to the `-logfile` file, where the
    /// command line writes it (its standard error is the file there); the caller's error keeps
    /// it. LOSAT's own refusals stay off the file, as on the command line (`cli::main_entry`).
    pub(super) fn post(&mut self, error: &anyhow::Error) {
        if let Some(native) = error.downcast_ref::<crate::cli::NativeError>() {
            self.post_text(&native.message);
        }
    }

    /// Writes `text`, the message of a `NativeError` of the run, to the `-logfile` file (on disk or
    /// in memory), when the run's messages go there (`RunLog::post`).
    fn post_text(&mut self, text: &[u8]) {
        if let Some(target) = &mut self.0 {
            target.write_to_file(text);
        }
    }
}

/// The messages of one search of a call (`RunLog::sink`).
pub(super) struct RunSink<'a> {
    log: &'a mut RunLog,
    diagnostics: &'a mut (dyn Write + Send),
    to_log: bool,
}

impl RunSink<'_> {
    /// `RunLog::failed`.
    pub(super) fn failed(&mut self, error: &anyhow::Error) -> Error {
        self.log.failed(error)
    }

    /// Where the messages written to this sink end up (`RunSink::write`): where the run's messages
    /// go (`RunLog::place`) with `to_log`, else the caller's diagnostics.
    fn place(&self) -> TextPlace {
        match self.to_log {
            true => self.log.place(),
            false => TextPlace::StandardError,
        }
    }
}

impl Write for RunSink<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match (&mut self.log.0, self.to_log) {
            (Some(target), true) => {
                target.write_text(buf, &mut *self.diagnostics);
                Ok(buf.len())
            }
            _ => self.diagnostics.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.diagnostics.flush()
    }
}

/// Whose error a caller's run writes (`CallerSink::error`), which decides where the command line
/// writes it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorSource {
    /// NCBI's message (a `NativeError`'s), which NCBI writes where the run's messages go
    /// (`CallerSink::add`).
    Ncbi,
    /// LOSAT's own refusal (`Error: ` and the error, as `main` writes it), which NCBI has no
    /// message for.
    Losat,
}

/// The writer of the messages that a caller's run writes for `cerr` (`CallerStart::sink`): they go
/// where the run's messages go (`RunLog::sink`: the file of `-logfile` in TBLASTX, else the caller's
/// diagnostics), as they come. It also writes what the caller's run adds to the command line's
/// messages for its arguments (`CallerSink::add`), among them the error of an output that the run
/// could not write (`CallerSink::error`), and says where each went.
pub struct CallerSink<'a>(RunSink<'a>);

impl CallerSink<'_> {
    /// Writes `text`, what the caller's run adds to the command line's messages for its arguments
    /// (LOSAT Web writes outfmt 0, 6 and 7 besides the argv's `-outfmt`, and gives what the command
    /// line writes for those formats), and gives back where that is: the caller's diagnostics, the
    /// standard error, also where the run's messages go to the file of `-logfile`, which keeps what
    /// the command line's file holds for the arguments (line WQA's sweep, D3 review r2,
    /// 2026-11-19: the file had the other formats' lines); nowhere where the run's messages go
    /// nowhere (`-logfile ""` and `-logfile /dev/null` in TBLASTX, `DiagTarget::Discard`), as the
    /// command line's run of any format writes them nowhere there (D3 review r3, 2026-11-19).
    pub fn add(&mut self, text: &[u8]) -> TextPlace {
        match self.0.place() {
            TextPlace::Nowhere => TextPlace::Nowhere,
            TextPlace::StandardError | TextPlace::LogFile => {
                let _ = self.0.diagnostics.write_all(text);
                TextPlace::StandardError
            }
        }
    }

    /// Writes `text`, the error of an output that the run could not write (one that the command
    /// line's run of the arguments does not write: that run's own error ends the caller's run), and
    /// gives back where that is: NCBI's message (`ErrorSource::Ncbi`) where the run's other added
    /// messages go (`CallerSink::add`), LOSAT's own refusal (`ErrorSource::Losat`) to the caller's
    /// diagnostics, the standard error, which the command line writes it to whatever `-logfile` is
    /// (`cli::main_entry`). The place is the writer's own, so it stays where the text is (line
    /// WQA's sweep, D7, 2026-11-20).
    pub fn error(&mut self, text: &[u8], source: ErrorSource) -> TextPlace {
        match source {
            ErrorSource::Ncbi => self.add(text),
            ErrorSource::Losat => {
                let _ = self.0.diagnostics.write_all(text);
                TextPlace::StandardError
            }
        }
    }
}

impl Write for CallerSink<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// The application start of a run that its caller makes (LOSAT Web's adapter, which runs the start
/// of every run, as the command line does: `web_adapter::start_run`), with the files of the run
/// going to the caller's `RunFiles`: what `begin_with` left.
pub struct CallerStart {
    /// The error that ends the run at the configuration load (NCBI's report when the program's
    /// `.ini` or the file of `-conffile` ends it, LOSAT's refusal of the load), as `begin`'s.
    ended: Option<String>,
    /// The settings, or LOSAT's refusal of them, which the caller reports after it parsed the
    /// arguments (`ApplicationStart::settings`).
    settings: Result<ApplicationSettings, String>,
    /// Whether the arguments have `-dryrun` (`preparse`).
    dry_run: bool,
    /// The version text that the run prints in place of the rest (`version`).
    version: Option<String>,
    /// The `-logfile` target when the run's messages go there (`run_messages_follow_log`).
    log: RunLog,
    /// The `-logfile` target otherwise.
    target: Option<DiagTarget>,
    /// The file of `-logfile` (`log_file`) and where it goes, for `finish`.
    log_name: Option<std::ffi::OsString>,
    files: RunFiles,
}

/// The start `start`, made from the arguments `words` (after the program's name) by the caller,
/// up to its configuration load (`open`, `Opened::load`): the usage report's messages and the
/// configuration's that go to `cerr` are written to `diagnostics`, the file of `-logfile` goes to
/// `files`. `Err` when the start ends before its `-logfile` target is open (LOSAT's refusal of the
/// usage report); a start that ends at the load is a
/// `CallerStart` whose `ended` is the error, and whose file of `-logfile` has its messages. A
/// version word (`-version` or a full version) is printed where the command line's entry prints
/// it, after the `-logfile` target is open and before the load, which the start then does not
/// make (`CallerStart::version`; line WQA's D14, 2026-11-20).
pub fn begin_with<S: AsRef<OsStr> + Clone>(
    start: ApplicationStart,
    words: &[S],
    files: &RunFiles,
    diagnostics: &mut dyn Write,
) -> Result<CallerStart, String> {
    use crate::blastinput::ncbi_environment::{preparse, Preparsed};
    // A run of its own (LOSAT Web's), whose diagnostic stream never fails: LOSAT Web refuses
    // `-in_pssm` (`ncbi_exception::main_diag_failbit`), and no call it runs reads another's.
    let follows = run_messages_follow_log(start.program());
    let program = start.program();
    let version = crate::cli::ncbi_version_request(
        [OsStr::new("losat"), OsStr::new(program)]
            .into_iter()
            .chain(words.iter().map(AsRef::as_ref)),
    );
    // LOSAT Web's errors are text (its JSON): the start's messages are made from text
    // (`ApplicationStart`'s messages), so `Display` replaces no byte of them: the one message that
    // can hold bytes that are not UTF-8, a `-conffile` name's, is made from a name that its host
    // gives as text.
    let mut opened =
        open_with(start, words, files, diagnostics).map_err(|error| error.to_string())?;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001
    // NCBI's application reads and prints the versions before it loads the configuration.
    let ended = match version {
        Some(_) => None,
        None => opened
            .write_load(diagnostics)
            .err()
            .map(|error| error.to_string()),
    };
    let (log, target) = match follows {
        true => (RunLog(Some(opened.target)), None),
        false => (RunLog(None), Some(opened.target)),
    };
    Ok(CallerStart {
        ended,
        settings: opened.start.settings(),
        dry_run: matches!(preparse(words), Preparsed::Words { dryrun: true, .. }),
        version,
        log,
        target,
        log_name: log_file(words).map(OsStr::to_os_string),
        files: files.clone(),
    })
}

impl CallerStart {
    /// The error that ended the start, which ends the run, before its arguments are parsed.
    pub fn ended(&self) -> Option<&str> {
        self.ended.as_deref()
    }

    /// The settings of the start, or LOSAT's refusal of them (the error of the run after its
    /// arguments are parsed).
    pub fn settings(&self) -> Result<&ApplicationSettings, &str> {
        self.settings.as_ref().map_err(String::as_str)
    }

    /// Whether the run is a dry run, which ends after its arguments and settings.
    pub fn dry_run(&self) -> bool {
        self.dry_run
    }

    /// The text that the run prints to the command line's standard output in place of everything
    /// after the start's opening, where the arguments have a version word (`-version`,
    /// `-version-full`, `-version-full-xml`, `-version-full-json`): NCBI's application prints it
    /// before the configuration load (which this start then did not make) and the arguments, exit
    /// status 0 (`cli::ncbi_version_request`, which the command line's entry prints).
    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }

    /// The writer of the messages that the run writes for `cerr`: they go where the run's
    /// messages go (`RunLog::sink`: the file of `-logfile` in TBLASTX, else `diagnostics`), as
    /// they come (`CallerSink`).
    pub fn sink<'a>(&'a mut self, diagnostics: &'a mut (dyn Write + Send)) -> CallerSink<'a> {
        CallerSink(self.log.sink(diagnostics, true))
    }

    /// Where the messages written to `sink` end up (`RunLog::place`): the file of `-logfile`, or
    /// nowhere, where the run's messages go to its target (TBLASTX's, `run_messages_follow_log`),
    /// else the caller's diagnostics. LOSAT Web tells its host where the run's messages go (line
    /// WQA's sweep, D7, 2026-11-20).
    pub fn messages_place(&self) -> TextPlace {
        self.log.place()
    }

    /// Writes `text`, the message of a `NativeError` that ended the run, to the `-logfile` file
    /// when the run's messages go there, where the command line writes it (`RunLog::post`).
    pub fn post(&mut self, text: &[u8]) {
        self.log.post_text(text);
    }

    /// The run's `-logfile` target, where the posts of an argument error go
    /// (`web_adapter::try_parse_cli`), as the command line's go to its own: `cerr` without one.
    pub(crate) fn target(&self) -> &DiagTarget {
        match (&self.log.0, &self.target) {
            (Some(target), _) | (None, Some(target)) => target,
            (None, None) => &DiagTarget::Stderr,
        }
    }

    /// Hands the file of `-logfile` to the run's files (`RunFiles::deliver_log`) at the end of the
    /// run.
    pub fn finish(self) -> std::io::Result<()> {
        match self.log.0.or(self.target) {
            Some(target) => self.files.deliver_log(self.log_name.as_deref(), target),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
#[path = "start_tests.rs"]
mod start_tests;
