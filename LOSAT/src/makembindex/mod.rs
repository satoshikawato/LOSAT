//! `losat makembindex`: NCBI BLAST+ 2.17.0 makembindex (`CMkIndexApplication`,
//! c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp), which writes the megablast index of a
//! BLAST database (`.NN.idx` volumes and the `.shd` super header of each database volume) or,
//! with `-old_style_index true`, of a FASTA file or a database. Its command line is NCBI's
//! (`ncbi_app`); the index is `crate::dbindex`'s. The work of its `Run` is in `work`.

mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Dep, Desc, Kind};
use crate::ncbi_app::Tool;

/// The text of `makembindex -help` (NCBI 2.17.0's, frozen). Its part before `OPTIONAL
/// ARGUMENTS` is the usage that `-h` and the argument errors print.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:61-69
const HELP: &str = include_str!("help.txt");

/// The `int` limit of `CArgAllow_Integers(1, kMax_Int)`.
const INT_MAX: i64 = i32::MAX as i64;

/// The arguments of `CMkIndexApplication::Init` in the order it describes them, each with the
/// constraint its `SetConstraint` gives; the values that `-verbosity` and `-iformat` allow are
/// listed in the order of their `std::set`.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:65-147
const DESCS: &[Desc] = &[
    key("input", Kind::Text, None, None),
    key("output", Kind::Text, None, None),
    key(
        "verbosity",
        Kind::Text,
        Some("normal"),
        Some(Constraint::Strings(&["normal", "quiet", "verbose"])),
    ),
    key(
        "iformat",
        Kind::Text,
        None,
        Some(Constraint::Strings(&["blastdb", "fasta"])),
    ),
    key("legacy", Kind::Boolean, Some("true"), None),
    key("idmap", Kind::Boolean, Some("false"), None),
    key("db_mask", Kind::Text, None, None),
    flag("show_filters"),
    key(
        "nmer",
        Kind::Integer,
        None,
        Some(Constraint::Integers(8, 15)),
    ),
    key(
        "ws_hint",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key(
        "volsize",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("stat", Kind::Text, None, None),
    key(
        "stride",
        Kind::Integer,
        None,
        Some(Constraint::Integers(1, INT_MAX)),
    ),
    key("old_style_index", Kind::Boolean, Some("false"), None),
];

/// The `SetDependency` calls of `Init`.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:143-146
fn dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    vec![
        ("show_filters", Dep::Excludes, "output"),
        ("db_mask", Dep::Requires, "input"),
    ]
}

/// makembindex as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's prints
/// the toolkit's release of the BLAST package (`2.17.0`; makembindex sets no version of its own),
/// no usage report and no prefix on its diagnostics.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/main.cpp:51-54
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.hpp:41-56
static TOOL: Tool = Tool {
    program: "makembindex",
    refusal_name: crate::dbindex::PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies,
    post_prefix: None,
    components: &[],
    dry_run: None,
};

/// Runs `makembindex` with `words` (the command line after `losat makembindex`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around makembindex's `Run`:
/// `ncbi_app::run_with`). The library API's `run_makembindex`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}

#[cfg(test)]
mod tests;
