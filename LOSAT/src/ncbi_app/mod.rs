//! NCBI's application (`CNcbiApplication::AppMain`, c++/src/corelib/ncbiapp.cpp) for the programs
//! that LOSAT reads with NCBI's argument reader rather than clap's: `losat blastdbcmd`, `losat
//! dustmasker`, `losat segmasker`, `losat convert2blastmask`, `losat makembindex` and `losat
//! windowmasker`. Each program is a `Tool` (its name, version, frozen `-help` text; its argument
//! XML is NCBI's) and a body, its `Run`; `run_with` does what NCBI's application does around it,
//! in NCBI's order: the diagnostics of `-logfile`, the words taken before the reader (`-version`,
//! `-conffile`, `-dryrun`), the configuration load, the arguments (`args`, NCBI's
//! `CArgDescriptions`), the usage and the report of an argument error, and LOSAT's refusals of
//! the NCBI C++ Toolkit options that it does not serve.

pub(crate) mod args;
pub(crate) mod files;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::{
    diag_target, log_file, preparse, ApplicationSettings, ApplicationStart, DiagTarget, Preparsed,
};
use crate::report::serial::encoder::ReplacedChar;

use args::{ArgException, Args, ArgsError};

/// A program run by NCBI's application and argument reader. Its arguments are those of its NCBI
/// `-xmlhelp` XML (`args::arguments` of `program`).
pub(crate) struct Tool {
    /// The program's name (`GetProgramBasename`): the start of its `-version` line, its usage,
    /// its `<program>.ini` and its settings (`ApplicationStart::settings`).
    pub(crate) program: &'static str,
    /// The program's name in LOSAT's refusals (`BLASTDBCMD`).
    pub(crate) refusal_name: &'static str,
    /// NCBI's `-help` text (`CArgDescriptions::PrintUsage` of the descriptions with the details),
    /// frozen. Its part before the details (`REQUIRED ARGUMENTS` or `OPTIONAL ARGUMENTS`) is
    /// the usage that `-h` and the argument errors print (`short_usage`); its synopsis lists the
    /// arguments in the order of the program's `Init` (`args::usage_position`).
    pub(crate) help: &'static str,
    /// The prefix of the program's diagnostics (`SetDiagPostPrefix`), if it sets one.
    pub(crate) post_prefix: Option<&'static str>,
    /// The components of its full versions (`CVersionAPI::AddComponentVersion`; windowmasker's).
    pub(crate) components: &'static [crate::ncbi_version::Component],
    /// NCBI's `-help-full` text where it is not `help` (makeblastdb's, which lists its hidden
    /// `-randomize`), frozen; `None`: `-help-full` writes `help`, as NCBI's programs without a
    /// hidden argument do.
    pub(crate) help_full: Option<&'static str>,
    /// The description that `-h` writes where it is not the one of `help` (rpsbproc's short
    /// description, `SetUsageContext` of its `Init`), frozen: the lines after `DESCRIPTION`, each
    /// with its indentation and line feed (`short_usage`).
    pub(crate) short_description: Option<&'static str>,
    /// The exit status of a run whose write of a help text (on the standard output), a usage or an
    /// argument error's report (on the standard error) fails: rpsbproc's 4, as its failed write of
    /// a report (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 9). `None`: a failed write is dropped,
    /// as NCBI's programs drop a failed write to `cout` and `cerr`.
    pub(crate) write_failure: Option<i32>,
    /// What runs in place of `Run` with `-dryrun` (`DryRun`, after `Init`): the work of the
    /// program's `Init` that LOSAT does in its `Run` (makeprofiledb's parameters and files), or
    /// nothing (`None`: exit status 0).
    pub(crate) dry_run: Option<Body>,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3829-3834
/// The text of `-h`: the usage and the description of `tool`'s `-help` text (`Tool::help`; its
/// short description in place of the description where the tool has one,
/// `Tool::short_description`), then NCBI's pointer to `-help`. The details start with the
/// mandatory arguments (`REQUIRED ARGUMENTS`, convert2blastmask's) or the optional ones. The
/// description states LOSAT's version (`help_text`).
pub(crate) fn short_usage(tool: &Tool) -> String {
    let help = tool.help;
    let end = ["REQUIRED ARGUMENTS", "OPTIONAL ARGUMENTS"]
        .iter()
        .filter_map(|heading| help.find(heading))
        .min()
        .unwrap_or(help.len());
    const DESCRIPTION: &str = "\nDESCRIPTION\n";
    let usage = match tool.short_description {
        Some(description) => {
            let start = help[..end]
                .find(DESCRIPTION)
                .map_or(end, |at| at + DESCRIPTION.len());
            format!("{}{description}\n", &help[..start])
        }
        None => help[..end].to_string(),
    };
    format!(
        "{}Use '-help' to print detailed descriptions of command line arguments\n",
        help_text(&usage)
    )
}

/// LOSAT's text of NCBI's frozen `-help` text `help`: its description states LOSAT's version
/// (`BLAST database client, version LOSAT 0.2.0`; `wording::losat_text`).
pub(crate) fn help_text(help: &str) -> std::borrow::Cow<'_, str> {
    crate::wording::losat_text(help)
}

/// Where a program's diagnostics go: the target of `-logfile` (`stderr` for `cerr`), and the
/// standard error stream, on which the usage of an argument error is printed.
pub(crate) struct Diag<'a> {
    target: DiagTarget,
    stderr: &'a mut dyn Write,
    tool: &'static Tool,
    /// A copy of the texts written since `start_record` (`take_record`).
    record: Option<Vec<u8>>,
    /// Whether a write to the standard error failed (`Tool::write_failure`).
    failed: bool,
}

impl<'a> Diag<'a> {
    /// The diagnostics of `tool` that go to `target` (`stderr` for `cerr`), with the standard error
    /// stream `stderr`.
    pub(crate) fn new(target: DiagTarget, stderr: &'a mut dyn Write, tool: &'static Tool) -> Self {
        Diag {
            target,
            stderr,
            tool,
            record: None,
            failed: false,
        }
    }

    /// The target of the diagnostics, for what the program writes after its arguments.
    pub(crate) fn into_target(self) -> DiagTarget {
        self.target
    }
}

/// A standard error stream that notes a failed write (`Diag::failed`).
struct Noted<'w> {
    stream: &'w mut dyn Write,
    failed: &'w mut bool,
}

impl Write for Noted<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.stream.write(buf);
        *self.failed |= written.is_err();
        written
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let flushed = self.stream.flush();
        *self.failed |= flushed.is_err();
        flushed
    }
}

impl Diag<'_> {
    /// The standard error stream, which notes a failed write.
    fn noted(&mut self) -> Noted<'_> {
        Noted {
            stream: &mut *self.stderr,
            failed: &mut self.failed,
        }
    }

    /// Posts one diagnostic line (`LOG_POST`, `ERR_POST`).
    pub(crate) fn post(&mut self, line: &str) {
        self.post_bytes(line.as_bytes());
    }

    /// `post` of a line that is not UTF-8.
    pub(crate) fn post_bytes(&mut self, line: &[u8]) {
        let mut text = line.to_vec();
        text.push(b'\n');
        self.write_text(&text);
    }

    /// Writes `text`, diagnostic lines that end with their line feeds (an NCBI report made by
    /// `crate::ncbi_exception`: `err_post_x`, `NcbiException::report_default_ex`).
    pub(crate) fn write_text(&mut self, text: &[u8]) {
        if let Some(record) = &mut self.record {
            record.extend_from_slice(text);
        }
        let mut stderr = Noted {
            stream: &mut *self.stderr,
            failed: &mut self.failed,
        };
        self.target.write_text(text, &mut stderr);
    }

    /// Starts keeping a copy of the texts that `write_text` writes (the posts of a read that
    /// NCBI repeats where LOSAT reads once: windowmasker's `-mk_counts`).
    pub(crate) fn start_record(&mut self) {
        self.record = Some(Vec::new());
    }

    /// The texts written since `start_record`, which stops keeping them.
    pub(crate) fn take_record(&mut self) -> Vec<u8> {
        self.record.take().unwrap_or_default()
    }

    /// Writes `text` to the standard error stream itself (an insertion into `cerr`), not to the
    /// target of `-logfile`: segmasker's report of an exception its `Run` catches.
    pub(crate) fn cerr(&mut self, text: &[u8]) {
        let mut stderr = self.noted();
        let _ = stderr.write_all(text);
        let _ = stderr.flush();
    }

    /// LOSAT's refusal: a line on the standard error stream itself, not on the target of
    /// `-logfile` (NCBI has no text for it, and `-logfile /dev/null` must not hide it: as the
    /// searches' refusals, line RG DECISIONS 32), exit status 1.
    pub(crate) fn refuse(&mut self, reason: &str) -> i32 {
        self.cerr(format!("Error: {reason}\n").as_bytes());
        1
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1175-1177
    /// An `ERR_POST` of an `Error` with `message` inside `Run`: the severity, the prefix of the
    /// application (`Tool::post_prefix`), the message.
    pub(crate) fn err_post(&mut self, message: &str) {
        self.err_post_bytes(message.as_bytes());
    }

    /// An `ERR_POST` of a `Warning` with `message` inside `Run`: as `err_post`, of severity
    /// `Warning`.
    pub(crate) fn warning_post(&mut self, message: &str) {
        self.write_text(&crate::ncbi_exception::post_line(
            crate::ncbi_exception::Severity::Warning,
            None,
            self.tool.post_prefix,
            message.as_bytes(),
        ));
    }

    /// `err_post` of a message that is not UTF-8.
    pub(crate) fn err_post_bytes(&mut self, message: &[u8]) {
        self.write_text(&crate::ncbi_exception::error_post_line(
            None,
            self.tool.post_prefix,
            message,
        ));
    }

    /// The posts of the bytes that NCBI's object stream replaced while it wrote an object
    /// (`ReplacedChar`: an `ERR_POST_X` of subcode 7 with a message, of error code 803 for each
    /// byte, as `Error: (803.7) <prefix> <message>`).
    pub(crate) fn post_replaced(&mut self, replaced: &[ReplacedChar]) {
        for byte in replaced {
            debug_assert!(!byte.trace.is_empty(), "a byte replaced through `Traced`");
            self.write_text(&byte.post(self.tool.post_prefix));
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:704-710
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:274-285
    /// A `CArgException` that the program's `Run` throws (`files::open_input_file`,
    /// `Args::output_file`): the application adds `Application's execution failed` to it and
    /// reports it as an argument error, after the usage; exit status 1.
    pub(crate) fn run_arg_exception(&mut self, error: ArgException) -> i32 {
        self.report_arg_exception(&error, Some("Application's execution failed"))
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:675-690
    /// A `CArgException` that the program's `Init` throws (makeprofiledb's `AsInputFile`,
    /// `AsOutputFile`): reported as the reader's argument errors are, without a backlog; exit
    /// status 1.
    pub(crate) fn init_arg_exception(&mut self, error: ArgException) -> i32 {
        self.report_arg_exception(&error, None)
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7877-7888
    // With `eDPF_Severity` alone each exception of the chain is posted as `Error: <text>`, oldest
    // first, then the title's line, `Error: ` and the empty title, ` (CArgException::<code>) `
    // and the text of the newest exception that has its main text: the one the reader threw (a
    // backlog has none).
    /// NCBI's report of the argument error `error` (the backlog `backlog` added, when `Run`
    /// threw it): the usage, the line of `=` and the posts; exit status 1. Its callers are the
    /// reader's entry, `read`, and `run_arg_exception`.
    fn report_arg_exception(&mut self, error: &ArgException, backlog: Option<&str>) -> i32 {
        if error.code != "eNoValue" {
            // One write of the usage and the line of `=`: a caller's writer that fails keeps the
            // whole text (`api::tools`'s record of the diagnostics), as `write!`'s pieces would not.
            let mut text = short_usage(self.tool);
            // The diagnostics of a `-logfile` go to a stream that is not `cerr`: NCBI then says
            // where they are.
            if !matches!(self.target, DiagTarget::Stderr) {
                text.push_str(
                    "Error in command-line arguments. See error logs for more details.\n",
                );
            }
            text.push_str(&"=".repeat(72));
            text.push_str("\n\n");
            let mut stderr = self.noted();
            let _ = stderr.write_all(text.as_bytes());
            let _ = stderr.flush();
        }
        let messages = error.messages.iter().map(Vec::as_slice);
        for message in messages.chain(backlog.map(str::as_bytes)) {
            self.post_bytes(&[b"Error: ".as_slice(), message].concat());
        }
        // NCBI's `Error:  (CArgException::<code>) <message>`, in LOSAT's words (`crate::wording`:
        // without the class and code, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12).
        self.post_bytes(&[b"Error:  ".as_slice(), error.message()].concat());
        1
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-697
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112 (cited at
// `Diag::report_arg_exception`)
/// NCBI's reader of `tool`'s arguments `words` (`args::parse`, the words that the application
/// leaves to it: `ncbi_environment::preparse`) with what NCBI's application writes for each
/// outcome; the arguments, or the exit status of the command line that the reader ends. The one
/// producer of these bytes for every caller of the reader.
///
/// - `-h`: the short usage (`short_usage`) on `out`, exit status 0.
/// - `-help` and `-help-full`: the program's `-help` text (`help_text`; `-help-full` the text of
///   `Tool::help_full` where the program has one) on `out`, exit status 0.
/// - `-xmlhelp`: NCBI's XML of the program's arguments on `out`, made by `xml_help` (the
///   command line's; LOSAT's XML leaves out the arguments that LOSAT refuses,
///   `ncbi_app::args::refusal`, `cli::xml_help_text`), exit status 0; a caller without it (the
///   library API) refuses the word (`ReadEnd::Refused`: `cli::xml_help_refusal`, line RG DECISIONS
///   E2x-1).
/// - An argument error (`ArgException`): on `diag`'s standard error, the short usage, then, when
///   the diagnostics go to the file of `-logfile`, the line `Error in command-line arguments. See
///   error logs for more details.`, then a line of 72 `=` and an empty line; then NCBI's posts of
///   the exception on `diag`'s target (the file of `-logfile`, else the standard error): `Error:
///   <message>` for each message of the chain, oldest first, and `Error:  <message>` of the
///   exception's own (`Diag::report_arg_exception`); exit status 1.
///
/// A failed write is dropped, or, for a tool with a status of its own (`Tool::write_failure`:
/// rpsbproc), ends the run with that status: a failed write of the text on `out` is reported on
/// the standard error first (`<program> failed to write output: <error>`).
pub(crate) fn read(
    tool: &'static Tool,
    words: &[OsString],
    out: &mut dyn Write,
    diag: &mut Diag,
    xml_help: Option<fn(&str) -> Option<String>>,
) -> Result<Args, ReadEnd> {
    let text = match args::parse(tool, words) {
        Ok(args) => return Ok(args),
        Err(ArgsError::Usage) => short_usage(tool),
        Err(ArgsError::Help) => help_text(tool.help).into_owned(),
        Err(ArgsError::HelpFull) => help_text(tool.help_full.unwrap_or(tool.help)).into_owned(),
        Err(ArgsError::XmlHelp) => match xml_help.and_then(|xml| xml(tool.program)) {
            Some(xml) => xml,
            None => {
                return Err(ReadEnd::Refused(crate::cli::xml_help_refusal(
                    tool.refusal_name,
                )))
            }
        },
        // `eNoValue` (which prints no usage) is not raised by the reader.
        Err(ArgsError::Invalid(error)) => {
            let status = diag.report_arg_exception(&error, None);
            return Err(ReadEnd::Exit(match (diag.failed, tool.write_failure) {
                (true, Some(failure)) => failure,
                _ => status,
            }));
        }
    };
    if let Err(error) = out.write_all(text.as_bytes()).and_then(|()| out.flush()) {
        if let Some(failure) = tool.write_failure {
            diag.cerr(format!("{} failed to write output: {error}\n", tool.program).as_bytes());
            return Err(ReadEnd::Exit(failure));
        }
    }
    Err(ReadEnd::Exit(0))
}

/// How a read of the arguments (`read`) ends where it gives none.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum ReadEnd {
    /// The run ends with this exit status: `read` has written what NCBI writes (the help on the
    /// standard output, an argument error's report).
    Exit(i32),
    /// LOSAT's refusal of a word that the caller does not run (`-xmlhelp` in the library API),
    /// not written: the caller writes it as its refusals (`Diag::refuse` on the command line).
    Refused(String),
}

/// Runs `tool` with `words` (the command line after `losat <program>`), after the command line's
/// application start `start` of `tool`, and its body `body` on the process's standard streams;
/// the exit status.
pub(crate) fn run(
    tool: &'static Tool,
    words: &[OsString],
    start: &ApplicationStart,
    body: Body,
) -> i32 {
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    debug_assert_eq!(
        start.program(),
        tool.program,
        "the start of another program"
    );
    let status = run_app(
        tool,
        words,
        (&mut std::io::stdin(), &mut out, &mut err),
        body,
        // The command line's `-xmlhelp`: NCBI's XML without the arguments that LOSAT refuses
        // (of the programs of `ncbi_app`, makeprofiledb's `-binary` at its run, `ncbi_xmlhelp_tests`).
        (start, Some(crate::cli::xml_help_text)),
    );
    let _ = out.flush();
    let _ = err.flush();
    status
}

/// A program's `Run`: it gets the arguments, the application's settings, the standard input
/// (`files::Cin`, what an `eInputFile` argument `-` reads: `files::open_input_file`), the
/// standard output and the diagnostics, and returns the exit status.
pub(crate) type Body = fn(
    args: &Args,
    settings: &ApplicationSettings,
    input: &mut files::Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32;

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit` loads the
// configuration before the arguments are read; cited at `ApplicationStart`)
/// The run of `tool` with `words`, its standard input `input`, its output on `out` and its
/// errors on `err`; the exit status. The library API's run, which runs NCBI's application start
/// of this process as the command line's entry does (`ApplicationStart::of_process`): the usage
/// report's messages go to `err` at once.
///
/// The order is NCBI's: the diagnostic file of `-logfile` is opened first (append; stderr when
/// it cannot be opened), then `-version` and the full versions are read (the first one is
/// printed, in the library API too, whatever the bytes of the other words, `ncbi_version::print`),
/// then the configuration is loaded (the file of `-conffile` in place of `<program>.ini`; its
/// messages go where the diagnostics go; a `<program>.ini` that NCBI reports as malformed ends the
/// run with exit status 2: g_reg probes `o2`, `inibad_null`), then the arguments (`read`, which
/// writes what NCBI writes for a help word and an argument error); `-xmlhelp`, whose XML only the
/// command line writes, is refused where the reader takes it. `-dryrun` ends the run with exit
/// status 0 after the arguments and the settings. A value that is not UTF-8 is refused where the
/// body would start (`Args::text_values_refusal`). Then `body` runs.
pub(crate) fn run_with(
    tool: &'static Tool,
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
    body: Body,
) -> i32 {
    let start = ApplicationStart::of_process(tool.program, words);
    match start.usage_report() {
        Ok(text) => {
            let _ = err.write_all(text.as_bytes());
        }
        Err(reason) => {
            let _ = writeln!(err, "Error: {reason}");
            return 1;
        }
    }
    run_app(tool, words, (input, out, err), body, (&start, None))
}

/// The run of `run` and `run_with` after the application start `start` (whose usage report's
/// messages the caller has written), with the writer of NCBI's XML of `-xmlhelp` that the command
/// line's `run` has (`xml_help`); the library API's run has none and refuses `-xmlhelp`
/// (`cli::xml_help_refusal`, line RG DECISIONS E2x-1), so the library does not carry the XML.
fn run_app(
    tool: &'static Tool,
    words: &[OsString],
    (input, out, err): (&mut dyn Read, &mut dyn Write, &mut dyn Write),
    body: Body,
    (start, xml_help): (&ApplicationStart, Option<fn(&str) -> Option<String>>),
) -> i32 {
    // NCBI's application start (`ncbi_environment::ApplicationStart`): the file of `-logfile`
    // gets the usage report's messages when NCBI opens it, before `-version`.
    let mut target = diag_target(log_file(words));
    start.write_log_setup(&mut target, err);
    let (words, dryrun) = match preparse(words) {
        // `-version` and the full versions (`ncbi_version::print`, line RG item E3), before the
        // arguments are read, whatever their bytes.
        Preparsed::Version(request) => {
            let _ = out.write_all(
                crate::ncbi_version::print_with_components(request, tool.program, tool.components)
                    .as_bytes(),
            );
            return 0;
        }
        Preparsed::Words { words, dryrun, .. } => (words, dryrun),
    };
    let mut diag = Diag::new(target, err, tool);
    // The configuration (the file of `-conffile` in place of `<program>.ini`,
    // `ApplicationStart::of_process`).
    match start.write_config_load(&mut diag.target, diag.stderr) {
        Ok(None) => {}
        Ok(Some(status)) => return status,
        Err(reason) => return diag.refuse(&reason),
    }
    // NCBI's XML, which the command line's run writes (`cli::xml_help_text`); the library API
    // writes none.
    let args = match read(tool, &words, out, &mut diag, xml_help) {
        Ok(args) => args,
        Err(ReadEnd::Exit(status)) => return status,
        Err(ReadEnd::Refused(reason)) => return diag.refuse(&reason),
    };
    let settings = match start.settings() {
        Ok(settings) => settings,
        Err(reason) => return diag.refuse(&reason),
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (cited at
    // `cli::main_entry`: blastdbcmd and dustmasker keep the default `DryRun`, which returns 0
    // after the arguments have passed their definitions; LOSAT refuses the settings it refuses
    // in a run first). NCBI BLAST+ 2.17.0 dustmasker: `-in missing.fa -dryrun` and
    // `-infmt blastdb -in missing -dryrun` give exit status 0 and no output.
    if dryrun {
        return match tool.dry_run {
            Some(init) => init(
                &args,
                &settings,
                &mut files::Cin::new(input),
                out,
                &mut diag,
            ),
            None => 0,
        };
    }
    // The programs read their values as text (`Args::value`): one that is not UTF-8 is refused
    // where NCBI's `Run` starts, after the reader's errors and `-dryrun`, which runs no body.
    if let Some(refusal) = args.text_values_refusal(tool) {
        return diag.refuse(&refusal);
    }
    body(
        &args,
        &settings,
        &mut files::Cin::new(input),
        out,
        &mut diag,
    )
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
