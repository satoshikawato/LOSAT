//! `losat convert2blastmask`: NCBI BLAST+ 2.17.0 convert2blastmask
//! (`CConvert2BlastMaskApplication`, c++/src/app/blastdb/convert2blastmask.cpp), which turns the
//! lower case letters of FASTA input, and its lines `from - to`, into the masking information
//! that makeblastdb's `-mask_data` reads. Its command line is NCBI's (`ncbi_app`); its reader is
//! NCBI's `CFastaReader` with the flags of `CMaskFromFasta`
//! (`writedb::fasta_source::reader::ReaderFlags::mask_from_fasta`) and its writer the masking
//! programs' `CMaskWriterBlastDbMaskInfo` (`seqmasks_io::writer`). The work of its `Run` is in
//! `work`.

mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, mandatory_key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "CONVERT2BLASTMASK";

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:141-145
/// The text of `convert2blastmask -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`,
/// `-version-full`, `-xmlhelp` and `-dryrun` are hidden. Its part before `REQUIRED ARGUMENTS`
/// is the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:154-186
/// The arguments of `CConvert2BlastMaskApplication::Init` in the order it describes them; the
/// values that `-outfmt` allows are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key(
        "outfmt",
        Kind::Text,
        Some("maskinfo_asn1_text"),
        Some(Constraint::Strings(&[
            "interval",
            "maskinfo_asn1_bin",
            "maskinfo_asn1_text",
            "maskinfo_xml",
        ])),
    ),
    flag("parse_seqids"),
    mandatory_key("masking_algorithm", Kind::Text),
    mandatory_key("masking_options", Kind::Text),
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:111-120
/// convert2blastmask as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's
/// states BLAST's (`ncbi_version::print`), BLAST's usage report (`ncbi_environment`), no
/// dependencies between its arguments and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "convert2blastmask",
    refusal_name: PROGRAM,
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
    components: &[],
    dry_run: None,
};

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:259-263
/// Runs `convert2blastmask` with `words` (the command line after `losat convert2blastmask`),
/// after its application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around convert2blastmask's `Run`:
/// `ncbi_app::run_with`). The library API's `run_convert2blastmask`.
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
