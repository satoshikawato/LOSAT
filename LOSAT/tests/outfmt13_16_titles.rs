//! The refusal of a title that the XML2 and JSON documents (outfmt 13 to 16) cannot write as NCBI
//! does (`report::blast_xml::names`: a byte outside US-ASCII or a control byte) looks only at the
//! reports that NCBI BLAST+ 2.17.0 writes, on whole runs of the binary of each search program with
//! `-subject` and of the library API (line NP, series SM).
//!
//! NCBI writes a `-subject` search's report of a query when the next query's results arrive, or
//! in the epilog (`m_IsBl2Seq`, blast_format.cpp:2362-2404, 2206-2219). A run that a later batch
//! stops (a reader error, a batch that fails its set-up) never reaches the epilog: the report of
//! the last query of the batches before is never made, so NCBI never fails on its titles, and
//! LOSAT writes NCBI's output instead of its refusal. The `xml2_title*` rows of
//! `tests/fasta_input_fixtures.py` hold NCBI's bytes of BLASTN, TBLASTX, TBLASTN and BLASTP; this
//! test adds BLASTX, the side files of outfmt 13 and 14, a subject's title, a batch that fails
//! its set-up and the library API, which that harness does not run.
//!
//! The query files (`tests/fixtures/fasta_input/`): `batch_title_late` (record 5 has `café` in its
//! title, record 6 a data line that NCBI's reader refuses) and `batch_title_first` (`café` on
//! record 1, the same error). A nucleotide record has 180 letters and a protein record 60:
//! `BATCH_SIZE=100` (50) makes a batch of each record, 200 (100) of two.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Without the caller's experiment switches (`LOSAT_X_*`) and NCBI's variables that change a
/// search (`tests/search_environment.rs`); `LOSAT_BLASTX_WORKER_LOG` stays, which a build with
/// every feature needs once it starts a pool.
fn clear_environment(command: &mut Command) {
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT"
                    | "BATCH_SIZE"
                    | "BL2SEQ_LEGACY"
                    | "BLASTDB"
                    | "CHUNK_SIZE"
                    | "OLD_FSC"
                    | "OVERLAP_CHUNK_SIZE"
                    | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/fasta_input")
        .join(name)
}

/// The five programs: the kind of their queries and of their subjects.
const PROGRAMS: [(&str, &str, &str); 5] = [
    ("blastn", "nuc", "nuc"),
    ("blastp", "prot", "prot"),
    ("blastx", "nuc", "prot"),
    ("tblastn", "prot", "nuc"),
    ("tblastx", "nuc", "nuc"),
];

fn extension(kind: &str) -> &'static str {
    if kind == "prot" {
        "faa"
    } else {
        "fa"
    }
}

fn file(name: &str, kind: &str) -> PathBuf {
    fixture(&format!("{name}.{kind}.{}", extension(kind)))
}

/// The batch size of one record a batch.
fn each(kind: &str) -> &'static str {
    if kind == "prot" {
        "50"
    } else {
        "100"
    }
}

/// A temporary directory of this test, removed at its end.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-outfmt13-16-titles-{tag}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create the temporary directory");
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

struct Run {
    status: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
    /// The side files of outfmt 13 and 14 (`o_<n>.json`, `o_<n>.xml` of `-out o.json`, `o.xml`),
    /// by name.
    side_files: Vec<String>,
}

/// Runs `command` (a search with `-outfmt outfmt`; 13 and 14 write to `dir/o.<ext>`).
fn run(mut command: Command, outfmt: &str, dir: &Path) -> Run {
    command.args(["-outfmt", outfmt]);
    if matches!(outfmt, "13" | "14") {
        let ext = if outfmt == "13" { "json" } else { "xml" };
        for entry in std::fs::read_dir(dir).expect("list the directory") {
            let path = entry.expect("an entry").path();
            if path.file_name().is_some_and(|name| {
                let name = name.to_string_lossy();
                name.starts_with("o.") || name.starts_with("o_")
            }) {
                std::fs::remove_file(path).expect("remove an old output");
            }
        }
        command.arg("-out").arg(dir.join(format!("o.{ext}")));
    }
    let output = command.output().expect("the search runs");
    let mut side_files: Vec<String> = std::fs::read_dir(dir)
        .expect("list the directory")
        .map(|entry| {
            entry
                .expect("an entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .filter(|name| name.starts_with("o_"))
        .collect();
    side_files.sort();
    Run {
        status: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        side_files,
    }
}

fn search(
    program: &str,
    query: &Path,
    subject: &Path,
    env: &[(&str, &str)],
    threads: &str,
) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    clear_environment(&mut command);
    command
        .arg(program)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-num_threads", threads]);
    for (name, value) in env {
        command.env(name, value);
    }
    command
}

const REFUSAL: &str = "not supported by LOSAT";

/// The reports that a run wrote: the documents of outfmt 15 and 16 on the output, the side files
/// of 13 and 14.
fn reports(run: &Run, outfmt: &str) -> usize {
    let stdout = String::from_utf8_lossy(&run.stdout);
    match outfmt {
        "16" => stdout.matches("\n<BlastOutput2>").count(),
        "15" => stdout.matches("\"report\": {").count(),
        _ => run.side_files.len(),
    }
}

/// A title in a report that NCBI never writes: LOSAT ends as NCBI does (`status`, the error text
/// `error`), after the `written` reports before.
fn not_refused(what: &str, run: &Run, outfmt: &str, status: i32, error: &str, written: usize) {
    assert!(
        !run.stderr.contains(REFUSAL),
        "{what}: refused: {}",
        run.stderr
    );
    assert_eq!(
        run.status,
        Some(status),
        "{what}: status, stderr {}",
        run.stderr
    );
    assert!(run.stderr.contains(error), "{what}: stderr {}", run.stderr);
    assert_eq!(reports(run, outfmt), written, "{what}: reports");
}

/// A title that NCBI writes: LOSAT refuses it before the first byte.
fn refused(what: &str, run: &Run) {
    assert_eq!(run.status, Some(1), "{what}: status, stderr {}", run.stderr);
    assert!(
        run.stderr.contains(REFUSAL),
        "{what}: stderr {}",
        run.stderr
    );
    assert!(
        run.stdout.is_empty() && run.side_files.is_empty(),
        "{what}: wrote {} bytes and {} side files",
        run.stdout.len(),
        run.side_files.len()
    );
}

const READER_ERROR: &str = "BLAST query error: CFastaReader: Near line";

#[test]
fn the_title_of_the_report_that_ncbi_holds_back_and_never_writes_is_not_refused() {
    let dir = TempDir::new("held");
    for (program, query_kind, subject_kind) in PROGRAMS {
        let each = each(query_kind);
        let subject = fixture(&format!(
            "base.s.{subject_kind}.{}",
            extension(subject_kind)
        ));
        let late = file("batch_title_late", query_kind);
        let first = file("batch_title_first", query_kind);
        for outfmt in ["13", "14", "15", "16"] {
            for threads in ["1", "4"] {
                let what = |case: &str| format!("{program} {case} outfmt {outfmt} n{threads}");
                // Record 5 is the last query of the batches before record 6's: its report is held
                // back and never written. Records 1 to 4 are.
                let result = run(
                    search(program, &late, &subject, &[("BATCH_SIZE", each)], threads),
                    outfmt,
                    &dir.0,
                );
                not_refused(&what("late each"), &result, outfmt, 1, READER_ERROR, 4);
                // Record 1's report is written when record 2's results arrive.
                let result = run(
                    search(program, &first, &subject, &[("BATCH_SIZE", each)], threads),
                    outfmt,
                    &dir.0,
                );
                refused(&what("first each"), &result);
                // A run that reads to the end of the file writes every report in the epilog: the
                // second record's title (`batch_title.s` read as the query) is refused.
                let result = run(
                    search(
                        program,
                        &file("batch_title.s", query_kind),
                        &subject,
                        &[("BATCH_SIZE", each)],
                        threads,
                    ),
                    outfmt,
                    &dir.0,
                );
                refused(&what("last record, no error"), &result);
            }
        }
    }
}
