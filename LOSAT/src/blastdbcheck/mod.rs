//! `losat blastdbcheck`: NCBI BLAST+ 2.17.0 blastdbcheck (`CBlastDbCheckApplication`,
//! c++/src/app/blastdb/blastdbcheck.cpp), the check of the databases that `-db` names or that a
//! directory holds (`-dir`, `-recursive`): each volume's metadata, the sequences that `-full`,
//! `-stride`, `-random` or `-ends` choose (their definition lines, sequence, Bioseq, taxid with
//! `-must_have_taxids` and the lookup of each Seq-id unless `-no_isam`), and each alias file's
//! lines, reported at the `-verbosity` asked. Its command line is NCBI's (`ncbi_app`); the
//! databases are read by LOSAT's reader (`seqdb`). The work of its `Run` is in `work`, the tests
//! of a volume in `volume`, the test of an alias file in `alias`.

mod alias;
mod log;
mod path;
mod volume;
mod work;

use std::ffi::OsString;

use crate::ncbi_app::Tool;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLASTDBCHECK";

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:142-252
/// The text of `blastdbcheck -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`,
/// `-version-full`, `-xmlhelp` and `-dryrun` are hidden (the program's own `-logfile` too, which
/// has the standard argument's name). Its part before `OPTIONAL ARGUMENTS` is the usage that `-h`
/// and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1380-1386
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcheck.cpp:1546-1559
/// blastdbcheck as NCBI's application runs it (`ncbi_app`): BLAST's usage report, the arguments
/// of its XML (`ncbi_xmlhelp/blastdbcheck.xml`) and the prefix `blastdbcheck` on the diagnostics
/// of its `Run` (`SetDiagPostPrefix`).
static TOOL: Tool = Tool {
    program: "blastdbcheck",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: Some("blastdbcheck"),
    components: &[],
    help_full: None,
    short_description: None,
    dry_run: None,
    closed_log_ends_run: false,
};

/// Runs `blastdbcheck` with `words` (the command line after `losat blastdbcheck`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around blastdbcheck's `Run`:
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

/// `run_with` whose `Run` draws its random sample with the seed `seed` (the clock's seconds that
/// NCBI's sample test takes): the tests of the sampled rows.
#[cfg(test)]
pub(crate) fn run_seeded(
    words: &[OsString],
    seed: u32,
    out: &mut dyn std::io::Write,
    err: &mut dyn std::io::Write,
) -> i32 {
    work::SEED.with(|cell| cell.set(Some(seed)));
    let status =
        crate::ncbi_app::run_with(&TOOL, words, &mut std::io::empty(), out, err, work::run);
    work::SEED.with(|cell| cell.set(None));
    status
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
