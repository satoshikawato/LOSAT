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
/// and formatted by `format_archive`; without `-archive`, `-rid`'s value: none is NCBI's error of
/// the argument without a value (NCBI BLAST+ 2.17.0 `blast_formatter` alone, probe P2.205), the
/// empty one NCBI's error of `CRemoteBlast` (P2.221), any other LOSAT's refusal where NCBI asks its
/// servers for the search (`Refusal::Rid`). The errors are `CATCH_ALL`'s (`catch_all`).
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
