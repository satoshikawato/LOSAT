//! The NCBI application start of a library call (line RG item E1w-2).
//!
//! NCBI's application starts before it reads its arguments: it opens the file of `-logfile`,
//! loads its configuration (`<program>.ini` or the file of `-conffile`, then `.ncbirc`) and,
//! after the arguments, ends a `-dryrun` (`ncbi_environment::ApplicationStart`, quoted there).
//! The command line's entry runs one start per process (`cli::main_entry`); the library runs one
//! per call that needs it: at the call's beginning when its argv has `-logfile`, `-conffile` or
//! `-dryrun` (`Search::starts`), else when the call first reads NCBI's settings, where NCBI's
//! library looks them up itself (a database search, blastn's `-window_masker_taxid`, the organism
//! report of outfmt 18, the file of a `-matrix` that NCBI's tables do not have). A search of
//! records without those reads no NCBI setting (the contract of `api`).

use std::ffi::OsStr;
use std::io::Write;

use super::error::Error;
use super::input::Inputs;
use super::search::Search;
use crate::blastinput::ncbi_environment::{
    log_file, run_messages_follow_log, ApplicationSettings, ApplicationStart, DiagTarget,
};
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
/// to `files` (`RunFiles::log_target`: LOSAT Web's names that it refuses are refusals here).
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
        Err(reason) => return Err(Error::Search(reason.to_string())),
    }
    let mut target = files.log_target(log_file(words)).map_err(Error::Search)?;
    start.write_log_setup(&mut target);
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
                Err(Error::Search(
                    String::from_utf8_lossy(&text)
                        .trim_end_matches('\n')
                        .to_string(),
                ))
            }
            Err(reason) => Err(Error::Search(reason)),
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

/// The application start of one search call of the library: made at the call's beginning when its
/// argv has a word that the start reads (`Search::starts`), else the first time the call reads
/// NCBI's settings (`CallStart::settings`). The searches of a call share it (two of `run`'s
/// groups, the batches of `run_batches`).
pub(crate) struct CallStart {
    program: &'static str,
    /// The start's settings, once it has run.
    settings: Option<Result<ApplicationSettings, String>>,
    /// The `-logfile` target of a program whose run's messages go there
    /// (`run_messages_follow_log`: tblastx).
    log: Option<DiagTarget>,
    /// Whether the argv has `-dryrun`.
    dry_run: bool,
}

impl CallStart {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (quoted at
    // `cli::main_entry`: the BLAST programs keep the default `DryRun`, which returns 0 after the
    // arguments have passed their definitions)
    /// The start of a call of `search` with `inputs`, in the command line's order
    /// (`cli::main_entry`): the start, run now with its messages on `diagnostics` when the argv
    /// has `-logfile`, `-conffile` or `-dryrun`; the parse of the argv with the kind of the call's
    /// subject (`Search::command`: `-subject_loc` with a database fails here, dry run or not; it
    /// reads no record); then the settings that LOSAT refuses (`ApplicationStart::settings`), the
    /// call's error when the start has run, as the command line refuses them after its parse
    /// (`start.settings()?`).
    pub(super) fn new(
        search: &Search,
        inputs: Inputs<'_>,
        diagnostics: &mut dyn Write,
    ) -> Result<CallStart, Error> {
        let program = search.program().name();
        let mut call = CallStart {
            program,
            settings: None,
            log: None,
            dry_run: search.dry_run(),
        };
        let started = match search.starts() {
            true => Some(begin(program, search.argv(), diagnostics)?),
            false => None,
        };
        search.command(inputs)?;
        if let Some(started) = started {
            if let Err(reason) = &started.settings {
                return Err(Error::Search(reason.clone()));
            }
            call.settings = Some(started.settings);
            if run_messages_follow_log(program) {
                call.log = Some(started.target);
            }
        }
        Ok(call)
    }

    /// Whether the call is a dry run, which ends here: the start, the parse with the call's kind
    /// of subject and the settings have passed (`CallStart::new`). Neither the records nor
    /// `-outfmt` are read (NCBI reads them in `Run`).
    pub(super) fn dry_run(&self) -> bool {
        self.dry_run
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:345-347
    // ```c++
    //     if (m_Config.GetDataLoaderConfig().UseDataLoaders()) {
    //         m_InputReader.reset
    //             (new CBlastInputReader(m_Config.GetDataLoaderConfig(),
    // ```
    /// The application settings of the call, after running the start now (its messages on
    /// `diagnostics`) when it has not run. The records of a library call are read before its
    /// start, with NCBI's data loaders (`read_fasta`, `read_inputs`); settings that turn them off
    /// (`[BLAST] DATA_LOADERS`) read the same records: the loaders only make NCBI's reader take a
    /// record's first line as a Seq-id, which LOSAT's reader refuses
    /// (`fasta_reader::FastaReader::read_one_seq`), so a record that it read is read the same without
    /// them.
    pub(crate) fn settings(
        &mut self,
        diagnostics: &mut dyn Write,
    ) -> anyhow::Result<ApplicationSettings> {
        if self.settings.is_none() {
            let started = begin(self.program, &[] as &[&str], diagnostics)
                .map_err(|error| anyhow::Error::msg(error.message().to_string()))?;
            self.settings = Some(started.settings);
        }
        self.settings
            .clone()
            .expect("the start has run")
            .map_err(anyhow::Error::msg)
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
            self.post_text(native.message.as_bytes());
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

/// The application start of a run that its caller makes (LOSAT Web's adapter, which runs the start
/// of every run, as the command line does: `web_adapter::start_run`), with the files of the run
/// going to `files`: what `begin_with` left.
pub struct CallerStart {
    /// The error that ends the run at the configuration load (NCBI's report when the program's
    /// `.ini` or the file of `-conffile` ends it, LOSAT's refusal of the load), as `begin`'s.
    ended: Option<String>,
    /// The settings, or LOSAT's refusal of them, which the caller reports after it parsed the
    /// arguments (`ApplicationStart::settings`).
    settings: Result<ApplicationSettings, String>,
    /// The `-logfile` target when the run's messages go there (`run_messages_follow_log`).
    log: RunLog,
    /// The `-logfile` target otherwise.
    target: Option<DiagTarget>,
}

/// The start `start`, made from the arguments `words` (after the program's name) by the caller,
/// up to its configuration load (`open`, `Opened::load`): the usage report's messages and the
/// configuration's that go to `cerr` are written to `diagnostics`, the file of `-logfile` goes to
/// `files`. `Err` when the start ends before its `-logfile` target is open (LOSAT's refusal of the
/// usage report, a `-logfile` name that `files` refuses); a start that ends at the load is a
/// `CallerStart` whose `ended` is the error, and whose file of `-logfile` has its messages.
pub fn begin_with<S: AsRef<OsStr> + Clone>(
    start: ApplicationStart,
    words: &[S],
    files: &RunFiles,
    diagnostics: &mut dyn Write,
) -> Result<CallerStart, String> {
    let follows = run_messages_follow_log(start.program());
    let mut opened =
        open_with(start, words, files, diagnostics).map_err(|error| error.message().to_string())?;
    let ended = opened
        .write_load(diagnostics)
        .err()
        .map(|error| error.message().to_string());
    let (log, target) = match follows {
        true => (RunLog(Some(opened.target)), None),
        false => (RunLog(None), Some(opened.target)),
    };
    Ok(CallerStart {
        ended,
        settings: opened.start.settings(),
        log,
        target,
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

    /// Writes `messages`, which the run wrote for `cerr`, where the run's messages go
    /// (`RunLog::route`), and gives back what goes to `cerr`.
    pub fn route(&mut self, messages: Vec<u8>) -> Vec<u8> {
        self.log.route(messages)
    }

    /// Writes `text`, the message of a `NativeError` that ended the run, to the `-logfile` file
    /// when the run's messages go there, where the command line writes it (`RunLog::post`).
    pub fn post(&mut self, text: &[u8]) {
        self.log.post_text(text);
    }

    /// Hands the file of `-logfile` to `files` (`RunFiles::deliver_log`) at the end of the run;
    /// `words` are the start's.
    pub fn finish<S: AsRef<OsStr>>(self, words: &[S], files: &RunFiles) -> std::io::Result<()> {
        match self.log.0.or(self.target) {
            Some(target) => files.deliver_log(log_file(words), target),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
#[path = "start_tests.rs"]
mod start_tests;
