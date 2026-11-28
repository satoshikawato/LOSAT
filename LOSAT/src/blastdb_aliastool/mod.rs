//! `losat blastdb_aliastool`: NCBI BLAST+ 2.17.0 blastdb_aliastool (`CBlastDBAliasApp`,
//! c++/src/app/blastdb/blastdb_aliastool.cpp), every mode of it: alias files of a database and a
//! list (`-db` with `-gilist`, `-seqidlist`, `-taxidlist` or `-oid_masks`), of several databases
//! (`-dblist`, `-dblist_file`), of the volumes of one (`-num_volumes`) and of VDB names
//! (`-vdblist`, `-vdblist_file`); the binary GI list of a text one (`-gi_file_in`,
//! `-gi_file_out`); the binary seqid list of a text one (`-seqid_file_in`, `-seqid_file_out`,
//! `-seqid_title`, `-seqid_db`, `-seqid_dbtype`) and what one says of itself
//! (`-seqid_file_info`). Its command line is NCBI's (`ncbi_app`); the files are written by
//! LOSAT's writer (`writedb::alias_file`, `writedb::binary_list`, `writedb::seqidlist_writer`) in
//! the formats its reader reads (`seqdb`). The work of its `Run` is in `work`.

mod work;

use std::ffi::OsString;

use crate::ncbi_app::Tool;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLASTDB_ALIASTOOL";

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:146-340
/// The text of `blastdb_aliastool -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`,
/// `-version-full`, `-xmlhelp` and `-dryrun` are hidden. Its part before `OPTIONAL ARGUMENTS` is
/// the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:56-68
/// blastdb_aliastool as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's
/// states BLAST's (`CBlastVersion`), BLAST's usage report, the arguments of its XML
/// (`ncbi_xmlhelp/blastdb_aliastool.xml`, whose `-gi_file_out` is opened as it is read: `fPreOpen`)
/// and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "blastdb_aliastool",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: &[],
    help_full: None,
    short_description: None,
    dry_run: None,
    closed_log_ends_run: false,
};

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp:704-708
/// Runs `blastdb_aliastool` with `words` (the command line after `losat blastdb_aliastool`),
/// after its application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around blastdb_aliastool's `Run`:
/// `ncbi_app::run_with`), which the tests of NCBI's readers run in the process.
#[cfg(test)]
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn std::io::Read,
    out: &mut dyn std::io::Write,
    err: &mut dyn std::io::Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}
