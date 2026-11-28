//! `losat dustmasker`: NCBI BLAST+ 2.17.0 dustmasker 1.0.0 (`CDustMaskApplication`), the masker
//! of low-complexity regions of nucleotide sequences by the symmetric DUST algorithm. Its
//! command line is NCBI's (`ncbi_app`, NCBI's application and argument reader): `main` hands it
//! the words after `dustmasker`. The work of its `Run` is in `app`.

mod app;
mod fasta;
mod masker;
mod work;

use std::ffi::{OsStr, OsString};
use std::io::{Read, Write};

use crate::ncbi_app::Tool;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "DUSTMASKER";

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:76-85
/// The text of `dustmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions of `CDustMaskApplication::Init`, dust_mask_app.cpp:86-117, after `HideStdArgs`:
/// `-logfile`, `-conffile`, `-version` and `-dryrun` are not shown). Its part before
/// `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.hpp:51-55
/// dustmasker as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's states its
/// own (`ncbi_version::print`), the arguments of its XML (`ncbi_xmlhelp/dustmasker.xml`) and no
/// prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "dustmasker",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: &[],
    help_full: None,
    short_description: None,
    write_failure: None,
    dry_run: None,
};

// NCBI reference (598d8ae6): c++/src/app/dustmasker/main.cpp:39-42
/// Runs `dustmasker` with `words` (the command line after `losat dustmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, app::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around dustmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_dustmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, app::run)
}

/// `run_with` of LOSAT Web's run (`ncbi_app::run_in`; line NP WW4c, 2026-11-28): the files that
/// the run writes (`-out FILE`, `-logfile FILE`) go to `files`, and `-xmlhelp` prints NCBI's XML,
/// as on the command line.
pub(crate) fn run_in(
    words: &[OsString],
    files: &crate::report::run_files::RunFiles,
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_in(&TOOL, words, files, (input, out, err), app::run)
}

/// The checks of a dustmasker run before it opens a file (LOSAT Web's validation of the words,
/// `web_adapter::check_dustmasker`; line NP WW4c): `run_in` with the body `check`. The standard input
/// is empty; the exit status, with the standard error on `err` and the file of `-logfile` to
/// `files`.
pub(crate) fn check_in(
    words: &[OsString],
    files: &crate::report::run_files::RunFiles,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_in(
        &TOOL,
        words,
        files,
        (&mut std::io::empty(), out, err),
        check,
    )
}

/// `Run` up to where it opens its first file (`seqmasks_io::run::check_output`).
fn check(
    args: &crate::ncbi_app::args::Args,
    _settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    _stdin: &mut crate::ncbi_app::files::Cin<'_>,
    _out: &mut dyn Write,
    diag: &mut crate::ncbi_app::Diag,
) -> i32 {
    crate::seqmasks_io::run::check_output(args, diag, PROGRAM)
}

// No NCBI counterpart: LOSAT Web reads `-logfile` and `-out` against the input names the run
// opens for reading, so a write name equal to a missing input is refused instead of diverging
// from the command line, which makes the file and then reads it (line WQA D19; line NP WW4c).
/// The input file names that a `dustmasker` run with `words` opens for reading: `(for_logfile,
/// for_out)` (`log_inputs::masker_read_parts`: `-in`, or the files that reading the database of
/// `-infmt blastdb` checks or opens, `log_inputs::tool_read_files`).
pub(crate) fn read_files<S: AsRef<OsStr> + Clone>(words: &[S]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    crate::report::log_inputs::tool_read_files(
        TOOL.program,
        crate::seqdb::SeqType::Nucleotide,
        words,
        crate::report::log_inputs::masker_read_parts(TOOL.program, words),
    )
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
