//! `losat blast_formatter`: NCBI BLAST+ 2.17.0 blast_formatter (`CBlastFormatterApp`,
//! c++/src/app/blast/blast_formatter.cpp), which writes the reports of a search again from its
//! BLAST archive (`-outfmt 11`, `-archive`) without searching again. Its command line is NCBI's
//! (`ncbi_app`); `Run` opens the archive and hands it to `format_archive`, which reads its objects
//! and writes each one's report. NCBI's `-rid`, which fetches the results of a search from NCBI's
//! servers, is refused (`refusals`), as the search programs refuse `-remote`.

mod archive;
mod formatting;
pub(crate) mod refusals;
mod search_program;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::ncbi_app::args::{flag, key, Args, Constraint, Dep, Desc, Kind};
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:122-141
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2656-2737
// `CDebugArgs::SetArgumentDescriptions` (blast_formatter.cpp:143-144) adds nothing to a release
// build (its `-help` lists no debug option).
/// The arguments of `CBlastFormatterApp::Init` in the order it describes them: `-sorthits` and
/// `-sorthsps` take 0 to 4 (`CAlignFormatUtil::eEvalue` to `eQueryCoverage`, `eHspEvalue` to
/// `eSubjectStart`).
const DESCS: &[Desc] = &[
    key("rid", Kind::Text, None, None),
    key("archive", Kind::InputFile, None, None),
    key("outfmt", Kind::Text, Some("0"), None),
    flag("show_gis"),
    key(
        "num_descriptions",
        Kind::Integer,
        None,
        Some(Constraint::AtLeast(0)),
    ),
    key(
        "num_alignments",
        Kind::Integer,
        None,
        Some(Constraint::AtLeast(0)),
    ),
    key(
        "line_length",
        Kind::Integer,
        None,
        Some(Constraint::AtLeast(1)),
    ),
    flag("html"),
    key(
        "sorthits",
        Kind::Integer,
        None,
        Some(Constraint::Between(0, 4)),
    ),
    key(
        "sorthsps",
        Kind::Integer,
        None,
        Some(Constraint::Between(0, 4)),
    ),
    key(
        "max_target_seqs",
        Kind::Integer,
        None,
        Some(Constraint::AtLeast(1)),
    ),
    key("out", Kind::OutputFile, Some("-"), None),
    flag("parse_deflines"),
];

/// The `SetDependency` calls of `Init` (quoted at `DESCS`), in their order.
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    vec![
        ("rid", Dep::Excludes, "archive"),
        ("max_target_seqs", Dep::Excludes, "num_descriptions"),
        ("max_target_seqs", Dep::Excludes, "num_alignments"),
    ]
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:57-69
/// blast_formatter as NCBI's application runs it (`ncbi_app`): BLAST's version, BLAST's usage
/// report (`ncbi_environment`), its dependencies and no prefix on its diagnostics (its `Run` sets
/// none: `Warning: ...`, line RG item F3 probe P2.002).
static TOOL: Tool = Tool {
    program: "blast_formatter",
    refusal_name: PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies,
    post_prefix: None,
    components: &[],
};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:588-592
/// Runs `blast_formatter` with `words` (the command line after `losat blast_formatter`), after
/// its application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, run_body)
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
    /// `-outfmt` as given, or its default `0`.
    pub(crate) outfmt: String,
    pub(crate) show_gis: bool,
    pub(crate) num_descriptions: Option<usize>,
    pub(crate) num_alignments: Option<usize>,
    pub(crate) line_length: Option<usize>,
    pub(crate) html: bool,
    /// `-sorthits` and `-sorthsps` (0 to 4).
    pub(crate) sort_hits: Option<u8>,
    pub(crate) sort_hsps: Option<u8>,
    pub(crate) max_target_seqs: Option<usize>,
    /// `-out`: `-` is the standard output (`files::open_output_file`).
    pub(crate) out: String,
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
            outfmt: args.value("outfmt").unwrap_or("0").to_string(),
            show_gis: args.flag("show_gis"),
            num_descriptions: count("num_descriptions"),
            num_alignments: count("num_alignments"),
            line_length: count("line_length"),
            html: args.flag("html"),
            sort_hits: sort("sorthits"),
            sort_hsps: sort("sorthsps"),
            max_target_seqs: count("max_target_seqs"),
            out: args.value("out").unwrap_or("-").to_string(),
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
    if let Some(archive) = args.value("archive") {
        let mut archive = match files::open_input_file("archive", archive, input) {
            Ok(stream) => stream,
            Err(error) => {
                return catch_all(
                    diag,
                    crate::cli::command_line_argument_error(error.message()),
                )
            }
        };
        return archive::format_archive(&Options::of(args), settings, &mut archive, out, diag);
    }
    let Some(rid) = args.value("rid") else {
        return catch_all(
            diag,
            crate::cli::command_line_argument_error(&crate::ncbi_app::args::no_value_message(
                "rid",
            )),
        );
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
