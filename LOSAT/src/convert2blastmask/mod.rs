//! `losat convert2blastmask`: NCBI BLAST+ 2.17.0 convert2blastmask
//! (`CConvert2BlastMaskApplication`, c++/src/app/blastdb/convert2blastmask.cpp), which turns the
//! lower case letters of FASTA input, and its lines `from - to`, into the masking information
//! that makeblastdb's `-mask_data` reads. Its command line is NCBI's (`ncbi_app`); its reader is
//! NCBI's `CFastaReader` with the flags of `CMaskFromFasta`
//! (`writedb::fasta_source::reader::ReaderFlags::mask_from_fasta`) and its writer the masking
//! programs' `CMaskWriterBlastDbMaskInfo` (`seqmasks_io::writer`). The work of its `Run` is in
//! `work`.

mod work;

use std::ffi::{OsStr, OsString};
use std::io::{Read, Write};

use crate::ncbi_app::Tool;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "CONVERT2BLASTMASK";

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:141-145
/// The text of `convert2blastmask -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`,
/// `-version-full`, `-xmlhelp` and `-dryrun` are hidden. Its part before `REQUIRED ARGUMENTS`
/// is the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:111-120
/// convert2blastmask as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's
/// states BLAST's (`ncbi_version::print`), BLAST's usage report (`ncbi_environment`), the
/// arguments of its XML (`ncbi_xmlhelp/convert2blastmask.xml`) and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "convert2blastmask",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: &[],
    help_full: None,
    short_description: None,
    write_failure: None,
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
    crate::ncbi_app::run_in(&TOOL, words, files, (input, out, err), work::run)
}

/// The checks of a convert2blastmask run up to its opening of `-out` (LOSAT Web's validation of the
/// words, `web_adapter::check_convert2blastmask`; line NP WW4c, WW4f): `run_in` with the body
/// `work::check`. The standard input is empty; the exit status, with the standard error on `err`
/// and the file of `-logfile` to `files`.
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
        work::check,
    )
}

// No NCBI counterpart: LOSAT Web reads `-logfile` against the input names the run opens for
// reading (line WQA D19; line NP WW4c). `-in` is opened before `-out`, so a missing input ends the
// run before `-out` opens and `for_out` is empty.
/// The input file names that a `convert2blastmask` run with `words` opens for reading:
/// `(for_logfile, for_out)`: `-in`, unless it is the standard input.
pub(crate) fn read_files<S: AsRef<OsStr>>(words: &[S]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut for_logfile = Vec::new();
    crate::ncbi_app::args::scan_argv(TOOL.program, words, |name, value| {
        if name == "in" {
            if let Some(path) = value.filter(|v| !v.is_empty() && *v != OsStr::new("-")) {
                for_logfile.push(path.as_encoded_bytes().to_vec());
            }
        }
    });
    (for_logfile, Vec::new())
}

#[cfg(test)]
mod tests;
