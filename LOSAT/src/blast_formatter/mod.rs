//! `losat blast_formatter`: NCBI BLAST+ 2.17.0 blast_formatter (`CBlastFormatterApp`,
//! c++/src/app/blast/blast_formatter.cpp), which writes the reports of a search again from its
//! BLAST archive (`-outfmt 11`, `-archive`) without searching again. Its command line is NCBI's
//! (`ncbi_app`); `Run` opens the archive and hands it to `format_archive`, which reads its objects
//! and writes each one's report. NCBI's `-rid`, which fetches the results of a search from NCBI's
//! servers, is refused (`refusals`), as the search programs refuse `-remote`.

mod archive;
mod formatting;
mod pairwise;
pub(crate) mod refusals;
mod search_program;
mod tabular;

use std::ffi::{OsStr, OsString};
use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{self, Cin};
use crate::ncbi_app::{Diag, Tool};
use crate::report::run_files::RunFiles;

use refusals::Refusal;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLAST_FORMATTER";

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:112-120
/// The text of `blast_formatter -help` (NCBI 2.17.0's, frozen; line RG item F3, probe P2.201):
/// `-logfile`, `-conffile`, `-version-full`, `-xmlhelp` and `-dryrun` are hidden. Its part before
/// `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument errors print (probe P2.202).
const HELP: &str = include_str!("help.txt");

/// NCBI's DESCRIPTION line of `-help` (the line after `DESCRIPTION` in `HELP`), which the top
/// level's `-help` lists in LOSAT's words (`cli::entry_program_description`).
pub(crate) fn description() -> &'static str {
    HELP.split_once("DESCRIPTION\n")
        .and_then(|(_, rest)| rest.lines().next())
        .map_or("", str::trim)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:57-69
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:122-144
/// blast_formatter as NCBI's application runs it (`ncbi_app`): BLAST's version, BLAST's usage
/// report (`ncbi_environment`), the arguments and dependencies of its `Init` (the table made of
/// its `-xmlhelp` XML, `ncbi_xmlhelp/blast_formatter.xml`; `CDebugArgs` adds nothing to a release
/// build) and no prefix on its diagnostics (its `Run` sets none: `Warning: ...`, line RG item F3
/// probe P2.002).
static TOOL: Tool = Tool {
    program: "blast_formatter",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: &[],
    help_full: None,
    short_description: None,
    write_failure: None,
    dry_run: None,
};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:588-592
/// Runs `blast_formatter` with `words` (the command line after `losat blast_formatter`), after
/// its application start `start`, on the process's standard streams; the exit status. A negative
/// status is NCBI's `std::terminate` (`cli::NativeError::ios_failure_abort`, its report already
/// written): the process aborts, as NCBI's does (SIGABRT).
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    let status = crate::ncbi_app::run(&TOOL, words, start, run_body);
    if status < 0 {
        std::process::abort();
    }
    status
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around blast_formatter's `Run`:
/// `ncbi_app::run_with`). The library API's `run_blast_formatter`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, run_body)
}

/// `run_with` of LOSAT Web's run (`ncbi_app::run_in`): the files that the run writes (`-out FILE`,
/// `-logfile FILE`) go to `files`, and `-xmlhelp` prints NCBI's XML, as on the command line.
pub(crate) fn run_in(
    words: &[OsString],
    files: &RunFiles,
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_in(&TOOL, words, files, (input, out, err), run_body)
}

/// The checks of a blast_formatter run before it opens `-archive` (LOSAT Web's validation of the
/// words, `web_adapter::check_blast_formatter`): `run_in` with `check_body`, which stops where
/// `Run` would open `-archive`. The standard input is empty; the exit status is 0 for a valid argv
/// (including help/version words), with errors on `err` and the `-logfile` target in `files`.
pub(crate) fn check_in(
    words: &[OsString],
    run_files: &RunFiles,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_in(
        &TOOL,
        words,
        run_files,
        (&mut std::io::empty(), out, err),
        check_body,
    )
}

// No NCBI counterpart: LOSAT Web checks `-logfile` against the input names the run opens for
// reading, so a `-logfile` naming `-archive` (even when missing) is refused instead of creating
// the file before `-archive` reads it. `-out` is opened only after `-archive` is read, so
// `for_out` is empty (existing files in the run directory are still refused by `OpenFiles`).
/// The input file names that a `blast_formatter` run with `words` opens for reading:
/// `(for_logfile, for_out)`.
pub(crate) fn read_files<S: AsRef<OsStr>>(words: &[S]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut for_logfile = Vec::new();
    crate::ncbi_app::args::scan_argv(TOOL.program, words, |name, value| {
        if name == "archive" {
            if let Some(path) = value.filter(|v| !v.is_empty() && *v != OsStr::new("-")) {
                for_logfile.push(path.as_encoded_bytes().to_vec());
            }
        }
    });
    (for_logfile, Vec::new())
}

/// The command line of a run that formats an archive (`format_archive`), as `PrintFormattedOutput`
/// reads it for each archive object (blast_formatter.cpp:329-343): the formatting arguments
/// (`formatting::extract_formatting` reads `-outfmt` and the counts there, after the object), the
/// `-out` file (opened when the first object is formatted) and `-parse_deflines`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Options {
    /// `-outfmt` as given (raw bytes, preserving non-UTF-8 bytes as NCBI does), or its default `0`.
    pub(crate) outfmt: Vec<u8>,
    pub(crate) show_gis: bool,
    pub(crate) num_descriptions: Option<usize>,
    pub(crate) num_alignments: Option<usize>,
    pub(crate) line_length: Option<usize>,
    pub(crate) html: bool,
    /// `-sorthits` and `-sorthsps` (0 to 4).
    pub(crate) sort_hits: Option<u8>,
    pub(crate) sort_hsps: Option<u8>,
    pub(crate) max_target_seqs: Option<usize>,
    /// `-out` (its bytes, as NCBI opens the file): `-` is the standard output
    /// (`files::open_output_file`).
    pub(crate) out: OsString,
    pub(crate) parse_deflines: bool,
}

impl Options {
    /// The options of the arguments `args`, which passed their constraints: the counts are 0 or
    /// more (`-line_length` and `-max_target_seqs` 1 or more) and the sort options 0 to 4.
    fn of(args: &Args) -> Options {
        let count = |name: &str| {
            args.integer(name)
                .and_then(|number| usize::try_from(number).ok())
        };
        let sort = |name: &str| {
            args.integer(name)
                .and_then(|number| u8::try_from(number).ok())
        };
        Options {
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2801-2851
            outfmt: args
                .bytes("outfmt")
                .map_or_else(|| b"0".to_vec(), std::borrow::Cow::into_owned),
            show_gis: args.flag("show_gis"),
            num_descriptions: count("num_descriptions"),
            num_alignments: count("num_alignments"),
            line_length: count("line_length"),
            html: args.flag("html"),
            sort_hits: sort("sorthits"),
            sort_hsps: sort("sorthsps"),
            max_target_seqs: count("max_target_seqs"),
            out: args.os_value("out").unwrap_or(OsStr::new("-")).to_owned(),
            parse_deflines: args.flag("parse_deflines"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:499-545
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:879-884
/// blast_formatter's `Run`: the archive of `-archive` (`-` is the standard input `input`), read
/// and formatted by `format_archive`; without `-archive`, `-rid`'s value (`rid`). The errors are
/// `CATCH_ALL`'s (`catch_all`).
fn run_body(
    args: &Args,
    settings: &ApplicationSettings,
    input: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    if let Some(archive) = args.os_value("archive") {
        let mut archive = match files::open_input_file("archive", archive, input) {
            Ok(stream) => stream,
            Err(error) => {
                return catch_all(
                    diag,
                    crate::cli::command_line_argument_error(error.message()),
                )
            }
        };
        return archive::format_archive(
            &Options::of(args),
            settings,
            &mut archive,
            out,
            args.files(),
            diag,
        );
    }
    rid(args, diag)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:499-545
/// The checks of `run_body` before opening `-archive` (LOSAT Web's validation, `check_in`):
/// when `-archive` is given, only checks whether `-out` is a name LOSAT Web refuses to overwrite;
/// without `-archive`, `run_body`'s own `-rid` steps (`rid`), which open no files.
fn check_body(
    args: &Args,
    _settings: &ApplicationSettings,
    _input: &mut Cin<'_>,
    _out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    if args.os_value("archive").is_some() {
        let out = args.os_value("out").unwrap_or(OsStr::new("-"));
        if out != OsStr::new("-") {
            if let Some(reason) = args.files().refused_output(out) {
                return diag.refuse(&crate::seqmasks_io::run::unsupported(&reason, PROGRAM));
            }
        }
        return 0;
    }
    rid(args, diag)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:499-545
// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:879-884
/// `Run` without `-archive` (`run_body`, and `check_body` of LOSAT Web's validation, which opens no
/// file here either): `-rid`'s value. None is NCBI's error of the argument without a value (NCBI
/// BLAST+ 2.17.0 `blast_formatter` alone, probe P2.205), the empty one NCBI's error of
/// `CRemoteBlast` (P2.221), any other LOSAT's refusal where NCBI asks its servers for the search
/// (`Refusal::Rid`). The errors are `CATCH_ALL`'s (`catch_all`).
fn rid(args: &Args, diag: &mut Diag) -> i32 {
    let rid = match args.as_bytes("rid") {
        Ok(rid) => rid,
        Err(error) => {
            return catch_all(
                diag,
                crate::cli::command_line_argument_error(error.message()),
            )
        }
    };
    if rid.is_empty() {
        return catch_all(
            diag,
            crate::blastinput::app::engine_error("Empty RID string specified"),
        );
    }
    Refusal::Rid.refuse(diag)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:166-180
/// `CATCH_ALL` of `Run` for `error`: the text of a `NativeError` (made by the shared producers:
/// `cli::command_line_argument_error`, `blastinput::app::options_error`, `engine_error`) posted
/// where the diagnostics go (`-logfile`), and its exit status; any other error is LOSAT's refusal.
fn catch_all(diag: &mut Diag, error: anyhow::Error) -> i32 {
    match error.downcast_ref::<NativeError>() {
        Some(native) => {
            diag.write_text(&native.message);
            native.exit
        }
        None => diag.refuse(&format!("{error:#}")),
    }
}

#[cfg(test)]
#[path = "run_tests.rs"]
mod run_tests;
