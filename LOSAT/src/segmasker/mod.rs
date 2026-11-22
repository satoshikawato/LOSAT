//! `losat segmasker`: NCBI BLAST+ 2.17.0 segmasker 1.0.0 (`SegMaskerApplication`,
//! c++/src/app/segmasker/segmasker.cpp), the masker of low-complexity regions of protein
//! sequences by the SEG algorithm. Its command line is NCBI's (`ncbi_app`), as dustmasker's;
//! its input and output are the masking programs' (`seqmasks_io`), and its SEG is BLAST's
//! (`utils::seg::SegMasker`, without merging overlapping segments). The work of its `Run` is in
//! `work`.

mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "SEGMASKER";

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:96-104
/// The text of `segmasker -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`, `-version`
/// and `-dryrun` are hidden. Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the
/// argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:110-150
// (`kInputFormats` and the eight `kOutputFormats`: mask_cmdline_args.cpp:40-51, cited at
// dustmasker's descriptions; segmasker allows no `acclist`. The defaults `12`, `2.2`, `2.5` are
// those of segmasker 2.17.0's `-help`.)
/// The arguments of `SegMaskerApplication::Init` in the order it describes them; the values that
/// `-infmt` and `-outfmt` allow are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key(
        "infmt",
        Kind::Text,
        Some("fasta"),
        Some(Constraint::Strings(&["blastdb", "fasta"])),
    ),
    flag("parse_seqids"),
    key(
        "outfmt",
        Kind::Text,
        Some("interval"),
        Some(Constraint::Strings(&[
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
    key("window", Kind::Integer, Some("12"), None),
    key("locut", Kind::Double, Some("2.2"), None),
    key("hicut", Kind::Double, Some("2.5"), None),
];

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:69-74
/// segmasker as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's states its
/// own (`ncbi_version::print`), no dependencies between its arguments and no prefix on its
/// diagnostics.
static TOOL: Tool = Tool {
    program: "segmasker",
    refusal_name: PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
    components: &[],
};

// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:265-272
/// Runs `segmasker` with `words` (the command line after `losat segmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around segmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_segmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
