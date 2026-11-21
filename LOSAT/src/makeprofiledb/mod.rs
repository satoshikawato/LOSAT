//! `losat makeprofiledb`: NCBI BLAST+ 2.17.0 makeprofiledb (`CMakeProfileDBApp`,
//! c++/src/app/blastdb/makeprofiledb.cpp) with `-dbtype rps`, which writes the database of
//! profiles that rpsblast and rpstblastn search: the PSSMs (`.rps`), their frequency ratios
//! (`.freq`), the table of their words (`.loo`), their statistics (`.aux`) and the protein BLAST
//! database of their queries (`crate::writedb`). Its command line is NCBI's (`ncbi_app`); the
//! work of its `Run` is in `work`, the files of a volume in `rps_files`, the checks and the
//! statistics of a profile in `profile`.

mod profile;
mod rps_files;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, mandatory_key, Constraint, Dep, Desc, Kind};
use crate::ncbi_app::Tool;

/// The text of `makeprofiledb -help` (NCBI 2.17.0's, frozen). Its part before `REQUIRED
/// ARGUMENTS` is the usage that `-h` and the argument errors print.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:379-484
const HELP: &str = include_str!("help.txt");

/// The arguments of `CMakeProfileDBApp::Init` in the order it describes them, each with the
/// constraint its `SetConstraint` gives; the values of `CArgAllow_Strings` in the order of its
/// `std::set`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:379-484
const DESCS: &[Desc] = &[
    mandatory_key("in", Kind::InputFile),
    flag("binary"),
    key("title", Kind::Text, None, None),
    key("threshold", Kind::Double, Some("9.82"), None),
    flag("force"),
    key("out", Kind::Text, None, None),
    key(
        "blastdb_version",
        Kind::Integer,
        Some("5"),
        Some(Constraint::Integers(4, 5)),
    ),
    key("max_smp_vol", Kind::Integer, Some("2500"), None),
    key(
        "dbtype",
        Kind::Text,
        Some("rps"),
        Some(Constraint::Strings(&["cobalt", "delta", "rps"])),
    ),
    key("index", Kind::Boolean, Some("true"), None),
    key("gapopen", Kind::Integer, None, None),
    key("gapextend", Kind::Integer, None, None),
    key("scale", Kind::Double, Some("100"), None),
    key(
        "matrix",
        Kind::Text,
        Some("BLOSUM62"),
        Some(Constraint::Strings(&[
            "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM250", "PAM30", "PAM70",
        ])),
    ),
    key("obsr_threshold", Kind::Double, Some("6"), None),
    key("exclude_invalid", Kind::Boolean, Some("true"), None),
    key("taxid", Kind::Integer, None, Some(Constraint::AtLeast(0))),
    key("taxid_map", Kind::InputFile, None, None),
];

/// The `SetDependency` calls of `Init`.
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:476
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    vec![("taxid", Dep::Excludes, "taxid_map")]
}

/// makeprofiledb as NCBI's application runs it (`ncbi_app`): no prefix on its diagnostics; its
/// `Init` reads its parameters (`work::dry_run` with `-dryrun`).
static TOOL: Tool = Tool {
    program: "makeprofiledb",
    refusal_name: "MAKEPROFILEDB",
    help: HELP,
    descs: DESCS,
    dependencies,
    post_prefix: None,
    components: &[],
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
