//! `losat dustmasker`: NCBI BLAST+ 2.17.0 dustmasker 1.0.0 (`CDustMaskApplication`), the masker
//! of low-complexity regions of nucleotide sequences by the symmetric DUST algorithm. Its
//! command line is NCBI's (`ncbi_app`, NCBI's application and argument reader), not clap's:
//! `main` hands it the words after `dustmasker`. The work of its `Run` is in `app`.

mod app;
mod fasta;
mod masker;
mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "DUSTMASKER";

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:76-85
/// The text of `dustmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `DESCS` after `HideStdArgs`: `-logfile`, `-conffile`, `-version` and `-dryrun`
/// are not shown). Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument
/// errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:86-117
// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_cmdline_args.cpp:40-51
/// The arguments of `CDustMaskApplication::Init` in the order it describes them; NCBI's reader
/// adds the help flags before them and `-logfile` and `-conffile` after them (`ncbi_app::args`).
/// The values that `-outfmt` allows are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key("window", Kind::Integer, Some("64"), None),
    key("level", Kind::Integer, Some("20"), None),
    key("linker", Kind::Integer, Some("1"), None),
    key("infmt", Kind::Text, Some("fasta"), None),
    key(
        "outfmt",
        Kind::Text,
        Some("interval"),
        Some(Constraint::Strings(&[
            "acclist",
            "fasta",
            "interval",
            "maskinfo_asn1_bin",
            "maskinfo_asn1_text",
            "maskinfo_xml",
            "seqloc_asn1_bin",
            "seqloc_asn1_text",
            "seqloc_xml",
        ])),
    ),
    flag("parse_seqids"),
    flag("hard_masking"),
];

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.hpp:51-55
/// dustmasker as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's states its
/// own (`ncbi_version::print`), no dependencies between its arguments and no prefix on its
/// diagnostics.
static TOOL: Tool = Tool {
    program: "dustmasker",
    refusal_name: PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
    components: &[],
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

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
