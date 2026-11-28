//! `losat windowmasker`: NCBI BLAST+ 2.17.0 windowmasker 1.0.0 (`CWinMaskApplication`), the
//! masker of over-represented and low-complexity regions of nucleotide sequences by the counts
//! of their units (WindowMasker, `utils::winmask`) and, with `-dust`, the symmetric DUST
//! algorithm. Its command line is NCBI's (`ncbi_app`): `main` hands it the words
//! after `windowmasker`. Its three modes are the making of a counts file (`-mk_counts`, with
//! `-checkdup` its check for duplicated sequences), the conversion of a counts file to another
//! format (`-convert`) and the masking with a counts file (`-ustat`). The work of its `Run` is in
//! `app` (the configuration, `CWinMaskConfig`), `counts` and `checkdup` (the counts) and `work`
//! (the masking).

mod app;
mod checkdup;
pub(crate) mod counts;
mod ids;
mod sdust;
mod work;

use std::ffi::{OsStr, OsString};
use std::io::{Read, Write};

use crate::ncbi_app::Tool;
use crate::ncbi_version::Component;

/// The name in LOSAT's refusals.
const PROGRAM: &str = "WINDOWMASKER";

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:76,99-110
/// The text of `windowmasker -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions of `CWinMaskConfig::AddWinMaskArgs` after `HideStdArgs`), frozen. Its part before
/// `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:79-84
// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:50-254
/// windowmasker as NCBI's application runs it (`ncbi_app`): LOSAT's version where NCBI's states its
/// own (`ncbi_version::print`), the arguments of its XML (`ncbi_xmlhelp/windowmasker.xml`: those
/// of `CWinMaskConfig::AddWinMaskArgs` for `eAny`, the input determined, and its `SetDependency`
/// calls) and no prefix on its diagnostics.
static TOOL: Tool = Tool {
    program: "windowmasker",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: None,
    components: COMPONENTS,
    help_full: None,
    short_description: None,
    dry_run: None,
    closed_pipe_ends_run: true,
};

/// The components that `CWinMaskApplication`'s constructor adds: the algorithm, the counts
/// algorithm and the four counts formats, which LOSAT's masker and counts readers implement.
///
/// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:84-95
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:57-60,
/// seq_masker_ostat.cpp:41-47 (`STAT_ALGO_*`), seq_masker_ostat_opt_bin.cpp:39-43 (and the other
/// formats' `STAT_FMT_*`)
const COMPONENTS: &[Component] = &[
    component("window-masker-algorithm", "", 0),
    component("windowmasker-statistics-algorithm-version", "", 0),
    component("windowmasker-statistics-format-version", "ascii ", 1),
    component("windowmasker-statistics-format-version", "binary ", 1),
    component("windowmasker-statistics-format-version", "oascii ", 1),
    component("windowmasker-statistics-format-version", "obinary ", 1),
];

const fn component(name: &'static str, prefix: &'static str, minor: u32) -> Component {
    Component {
        name,
        prefix,
        major: 1,
        minor,
        patch: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/app/winmasker/main.cpp
/// Runs `windowmasker` with `words` (the command line after `losat windowmasker`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub(crate) fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, app::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around windowmasker's `Run`:
/// `ncbi_app::run_with`). The library API's `run_windowmasker`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, app::run)
}

/// `run_with` of LOSAT Web's run (`ncbi_app::run_in`): the files that the run writes (`-out FILE`,
/// `-logfile FILE`) go to `files`, and `-xmlhelp` prints NCBI's XML, as on the command line.
pub(crate) fn run_in(
    words: &[OsString],
    files: &crate::report::run_files::RunFiles,
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_in(&TOOL, words, files, (input, out, err), app::run)
}

/// The checks of a windowmasker run before it reads its input or opens a file (LOSAT Web's
/// validation of the words, `web_adapter::check_windowmasker`; line NP WW3a, D-WW3-7, 2026-11-22,
/// and D-WW3-10): `run_in` with the body `app::check`, which ends where `Run` would read or open
/// its first file. The standard input is
/// empty and the standard output dropped (a help or a version is printed in place of a run, which
/// is valid); the exit status, with the standard error on `err` and the file of `-logfile` to
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
        app::check,
    )
}

// No NCBI counterpart: LOSAT Web reads `-logfile` and `-out` against the input names the run
// opens for reading, so a write name equal to a missing input is refused instead of diverging
// from the command line, which creates the file and then reads it back (line WQA D19 round R2,
// D19b).
/// The input file names that a `windowmasker` run with `words` opens for reading:
/// `(for_logfile, for_out)`.
/// - `for_logfile`: every file that the run opens after the application start opens `-logfile`:
///   `-in` (split via `counts::file_list` for `-mk_counts`, including the `-fa_list true` list file
///   and every entry read from it; for `-infmt blastdb`, the database lookup files from
///   `seqdb::alias::database_lookup_files`), `-ustat` (in masking mode), `-ids` and `-exclude_ids`
///   (when not `-convert`).
/// - `for_out`: the subset that the command line opens and reads after opening `-out` in a way
///   that changes the run when the file is missing: every input/database entry of `-mk_counts`
///   (comma-separated or read from `-fa_list true`, FASTA or `-infmt blastdb`), and `-ustat` of
///   masking.
pub(crate) fn read_files<S: AsRef<OsStr> + Clone>(words: &[S]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let (log_paths, out_paths, databases, databases_for_out) =
        crate::report::log_inputs::windowmasker_read_parts(words);
    let encode = |paths: Vec<std::path::PathBuf>| -> Vec<Vec<u8>> {
        paths
            .into_iter()
            .map(|p| p.into_os_string().into_encoded_bytes())
            .collect()
    };
    let mut for_log = encode(log_paths);
    let mut for_out = encode(out_paths);
    if !databases.is_empty() {
        let blastdb =
            crate::blastinput::ncbi_environment::ApplicationStart::of_process(PROGRAM, words)
                .settings()
                .ok()
                .and_then(|s| s.blastdb);
        for db in databases {
            for path in crate::seqdb::alias::database_lookup_files(
                &crate::utils::os_bytes::os_str_bytes(&db),
                crate::seqdb::SeqType::Nucleotide,
                blastdb.as_deref(),
            ) {
                if !path.is_empty() {
                    if databases_for_out {
                        for_out.push(path.clone());
                    }
                    for_log.push(path);
                }
            }
        }
    }
    (for_log, for_out)
}

#[cfg(test)]
mod tests;
