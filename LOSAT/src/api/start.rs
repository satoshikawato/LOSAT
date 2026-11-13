//! The NCBI application start of a library call (line RG item E1w-2).
//!
//! NCBI's application starts before it reads its arguments: it opens the file of `-logfile`,
//! loads its configuration (`<program>.ini` or the file of `-conffile`, then `.ncbirc`) and,
//! after the arguments, ends a `-dryrun` (`ncbi_environment::ApplicationStart`, quoted there).
//! The command line's entry runs one start per process (`cli::main_entry`); the library runs one
//! per call that needs it: at the call's beginning when its argv has `-logfile`, `-conffile` or
//! `-dryrun` (`Search::starts`), else when the call first reads NCBI's settings (a database
//! search, blastn's `-window_masker_taxid`). A search of records without those reads no NCBI
//! setting (the contract of `api`).

use std::ffi::OsStr;
use std::io::Write;

use super::error::Error;
use super::search::Search;
use crate::blastinput::ncbi_environment::{
    diag_target, log_file, run_messages_follow_log, ApplicationSettings, ApplicationStart,
    DiagStream, DiagTarget,
};

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
    let start = ApplicationStart::of_process(program, words);
    match start.usage_report() {
        Ok(text) => {
            let _ = diagnostics.write_all(text.as_bytes());
        }
        Err(reason) => return Err(Error::Search(reason.to_string())),
    }
    let mut target = diag_target(log_file(words));
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
        let mut stderr = Vec::new();
        match self.start.write_config_load(&mut self.target, &mut stderr) {
            Ok(None) => {
                let _ = diagnostics.write_all(&stderr);
                Ok(Started {
                    settings: self.start.settings(),
                    target: self.target,
                })
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
    /// The start of a call of `search`: run now, with its messages on `diagnostics`, when the
    /// argv has `-logfile`, `-conffile` or `-dryrun`.
    pub(super) fn new(search: &Search, diagnostics: &mut dyn Write) -> Result<CallStart, Error> {
        let program = search.program().name();
        let mut call = CallStart {
            program,
            settings: None,
            log: None,
            dry_run: search.dry_run(),
        };
        if search.starts() {
            let started = begin(program, search.argv(), diagnostics)?;
            call.settings = Some(started.settings);
            if run_messages_follow_log(program) {
                call.log = Some(started.target);
            }
        }
        Ok(call)
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (quoted at
    // `cli::main_entry`: the BLAST programs keep the default `DryRun`, which returns 0 after the
    // arguments have passed their definitions)
    /// Whether the call is a dry run, which ends here: the start has run and the arguments have
    /// passed the parser (`validate`), and LOSAT refuses the settings it refuses in a run, as the
    /// command line's dry run does. Neither the records nor `-outfmt` are read (NCBI reads them in
    /// `Run`).
    pub(super) fn dry_run(&self) -> Result<bool, Error> {
        if !self.dry_run {
            return Ok(false);
        }
        match &self.settings {
            Some(Err(reason)) => Err(Error::Search(reason.clone())),
            Some(Ok(_)) | None => Ok(true),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:345-347
    // ```c++
    //     if (m_Config.GetDataLoaderConfig().UseDataLoaders()) {
    //         m_InputReader.reset
    //             (new CBlastInputReader(m_Config.GetDataLoaderConfig(),
    // ```
    /// The application settings of the call, after running the start now (its messages on
    /// `diagnostics`) when it has not run. The records of a library call are read before its
    /// start, with NCBI's data loaders (`read_fasta`, `read_inputs`): settings that turn them off
    /// (`[BLAST] DATA_LOADERS`) would read a record's first line otherwise, and are refused.
    pub(crate) fn settings(
        &mut self,
        diagnostics: &mut dyn Write,
    ) -> anyhow::Result<ApplicationSettings> {
        if self.settings.is_none() {
            let started = begin(self.program, &[] as &[&str], diagnostics)
                .map_err(|error| anyhow::Error::msg(error.message().to_string()))?;
            self.settings = Some(started.settings);
        }
        let settings = self
            .settings
            .clone()
            .expect("the start has run")
            .map_err(anyhow::Error::msg)?;
        if !settings.data_loaders {
            anyhow::bail!(
                "the registry turns NCBI BLAST+'s data loaders off ([BLAST] DATA_LOADERS), with which its reader takes a record's first line otherwise; LOSAT's library reads the records before the call's application start, with the data loaders, and this is not supported by LOSAT"
            );
        }
        Ok(settings)
    }

    /// The `-logfile` target of the call's run messages (tblastx's), which the caller holds for
    /// the run (`RunLog`).
    pub(super) fn take_log(&mut self) -> RunLog {
        RunLog(self.log.take())
    }
}

/// The `-logfile` target of a run whose messages go there (`CallStart::take_log`), or none.
pub(super) struct RunLog(Option<DiagTarget>);

impl RunLog {
    /// No target: the run's messages go to its diagnostics.
    pub(super) fn none() -> RunLog {
        RunLog(None)
    }

    /// The writer of the run's messages: `diagnostics`, or the `-logfile` target (`cerr` is
    /// `diagnostics`).
    pub(super) fn diagnostics<'a>(
        &'a mut self,
        diagnostics: &'a mut (dyn Write + Send),
    ) -> Box<dyn Write + Send + 'a> {
        match &mut self.0 {
            Some(target) => Box::new(DiagStream {
                target,
                stderr: diagnostics,
            }),
            None => Box::new(diagnostics),
        }
    }

    /// The error of a failed run (`Error::search`). The text of its `NativeError` goes to the
    /// `-logfile` file too, where the command line writes it (its standard error is the file
    /// there); the error keeps it.
    pub(super) fn failed(&mut self, error: &anyhow::Error) -> Error {
        if let (Some(DiagTarget::File(file)), Some(native)) =
            (&mut self.0, error.downcast_ref::<crate::cli::NativeError>())
        {
            let _ = file.write_all(native.message.as_bytes());
        }
        Error::search(error)
    }
}
