//! `losat makeprofiledb`: NCBI BLAST+ 2.17.0 makeprofiledb (`CMakeProfileDBApp`,
//! c++/src/app/blastdb/makeprofiledb.cpp) with `-dbtype rps` and `-dbtype delta`, which writes the
//! database of profiles that rpsblast, rpstblastn and deltablast search: the PSSMs (`.rps`), their
//! frequency ratios (`.freq`), the table of their words (`.loo`), their statistics (`.aux`), the
//! protein BLAST database of their queries (`crate::writedb`), and for `-dbtype delta` the
//! weighted residue frequencies (`.wcounts`) and effective observations (`.obsr`). Its command
//! line is NCBI's (`ncbi_app`); the work of its `Run` is in `work`, the files of a volume in
//! `rps_files` and `delta`, the checks and the statistics of a profile in `profile`.

mod delta;
mod profile;
mod rps_files;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::Tool;

/// makeprofiledb's name in LOSAT's refusals (`<what> is not supported by LOSAT's
/// MAKEPROFILEDB`, `crate::seqmasks_io::run::unsupported`).
const PROGRAM: &str = "MAKEPROFILEDB";

/// The text of `makeprofiledb -help` (NCBI 2.17.0's, frozen). Its part before `REQUIRED
/// ARGUMENTS` is the usage that `-h` and the argument errors print.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:379-484
const HELP: &str = include_str!("help.txt");

/// makeprofiledb as NCBI's application runs it (`ncbi_app`): no prefix on its diagnostics; its
/// `Init` reads its parameters (`work::dry_run` with `-dryrun`).
static TOOL: Tool = Tool {
    program: "makeprofiledb",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: &[],
    help_full: None,
    short_description: None,
    dry_run: Some(work::dry_run),
};

/// Runs `makeprofiledb` with `words` (the command line after `losat makeprofiledb`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around makeprofiledb's `Run`:
/// `ncbi_app::run_with`).
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}
