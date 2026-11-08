//! The refusal of a title that the BLAST archive (outfmt 11) cannot write as NCBI does
//! (`report::archive::check_titles`, a byte outside `0x20..=0x7E` in text) looks only at the
//! records that NCBI BLAST+ 2.17.0 writes, on a whole run of the binary of each search program
//! with `-subject` (line NP, series SM, DECISIONS 24).
//!
//! NCBI writes the archive of each query batch that it completes (`WriteArchive`: the batch's
//! queries and every `-subject` record); the batch that stops the run (a reader error, a last
//! batch without letters, `Empty CBlastQueryVector` with `BATCH_SIZE=0`) writes neither, and the
//! error archive holds no sequence. Where NCBI writes the title (with `#` and its `Bad char`
//! messages), LOSAT refuses it before the first byte; where NCBI does not, LOSAT writes NCBI's
//! output (the `arch_title*` rows of `tests/fasta_input_fixtures.py` and `db.archtitle.*` of
//! `tests/db_search_manifest.tsv` hold NCBI's bytes; this test adds BLASTX's `-subject`, which
//! the first harness does not run, LOSAT's experimental reading and report paths of BLASTN, and
//! the library API's `-subject` searches with a batch size, which no harness runs).
//!
//! The query files (`tests/fixtures/fasta_input/`): `batch_title_late` (record 5 has `café` in its
//! title, record 6 a data line that NCBI's reader refuses), `batch_title_lastempty` (records 5 to 7
//! without letters, record 7 with `café`), `batch_late_bad` and `batch_clean`; the subject file
//! `batch_title.s` has `café` on its second record. A nucleotide record has 180 letters and a
//! protein record 60: `BATCH_SIZE=100` (50) makes a batch of each record, 200 (100) of two.

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

fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    clear_environment(&mut command);
    command
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

fn file(name: &str, kind: &str) -> PathBuf {
    fixture(&format!(
        "{name}.{kind}.{}",
        if kind == "prot" { "faa" } else { "fa" }
    ))
}

/// The batch sizes of one record a batch and of two (`each`, `two`).
fn sizes(kind: &str) -> (&'static str, &'static str) {
    if kind == "prot" {
        ("50", "100")
    } else {
        ("100", "200")
    }
}

struct Run {
    status: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
}

fn run(program: &str, query: &Path, subject: &Path, env: &[(&str, &str)], threads: &str) -> Run {
    let mut command = losat();
    command
        .arg(program)
        .args(["-query"])
        .arg(query)
        .arg("-subject")
        .arg(subject);
    command.args(["-outfmt", "11", "-num_threads", threads]);
    for (name, value) in env {
        command.env(name, value);
    }
    let output = command.output().expect("losat runs");
    Run {
        status: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

const REFUSAL: &str = "has a character the archive cannot write as NCBI does";

/// A title that NCBI never writes: LOSAT ends as NCBI does (`status`, the error text `error`),
/// after the archives of the batches before (`written`: whether there are some).
fn not_refused(what: &str, run: &Run, status: i32, error: &str, written: bool) {
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
    let archives = String::from_utf8_lossy(&run.stdout)
        .matches("Blast4-archive ::=")
        .count();
    // The archives of the batches before and the error archive.
    assert_eq!(archives > 1, written, "{what}: {archives} archives");
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
        run.stdout.is_empty(),
        "{what}: wrote {} bytes",
        run.stdout.len()
    );
}

const READER_ERROR: &str = "BLAST query error: CFastaReader: Near line";

#[test]
fn a_title_in_a_batch_that_ncbi_does_not_write_is_not_refused() {
    for (program, query_kind, subject_kind) in PROGRAMS {
        let (each, two) = sizes(query_kind);
        let subject = fixture(&format!(
            "base.s.{subject_kind}.{}",
            if subject_kind == "prot" { "faa" } else { "fa" }
        ));
        let bad_subject = file("batch_title.s", subject_kind);
        for threads in ["1", "4"] {
            let what = |case: &str| format!("{program} {case} threads {threads}");
            // Record 5 shares its batch with record 6's reader error.
            let late = file("batch_title_late", query_kind);
            let result = run(program, &late, &subject, &[("BATCH_SIZE", two)], threads);
            not_refused(&what("late two"), &result, 1, READER_ERROR, true);
            // The default batch holds every record up to the error.
            let result = run(program, &late, &subject, &[], threads);
            not_refused(&what("late default"), &result, 1, READER_ERROR, false);
            // The last batch has no letters: NCBI fails it.
            let last = file("batch_title_lastempty", query_kind);
            let result = run(program, &last, &subject, &[("BATCH_SIZE", each)], threads);
            not_refused(
                &what("lastempty each"),
                &result,
                3,
                "BLAST engine error: Warning: Sequence contains no data",
                true,
            );
            // No batch is completed: the subjects are never written.
            let bad = file("batch_late_bad", query_kind);
            let result = run(program, &bad, &bad_subject, &[], threads);
            not_refused(
                &what("subject, reader error"),
                &result,
                1,
                READER_ERROR,
                false,
            );
            // A batch size of 0 reads no query (BLASTN's is the mixer's, which completes the batch).
            let clean = file("batch_clean", query_kind);
            let result = run(
                program,
                &clean,
                &bad_subject,
                &[("BATCH_SIZE", "0")],
                threads,
            );
            if program == "blastn" {
                refused(&what("subject, BATCH_SIZE=0"), &result);
            } else {
                not_refused(
                    &what("subject, BATCH_SIZE=0"),
                    &result,
                    3,
                    "BLAST engine error: Empty CBlastQueryVector",
                    false,
                );
            }
        }
    }
}

#[test]
fn a_title_in_a_batch_that_ncbi_writes_is_refused_before_the_first_byte() {
    for (program, query_kind, subject_kind) in PROGRAMS {
        let (each, _) = sizes(query_kind);
        let subject = fixture(&format!(
            "base.s.{subject_kind}.{}",
            if subject_kind == "prot" { "faa" } else { "fa" }
        ));
        for threads in ["1", "4"] {
            let what = |case: &str| format!("{program} {case} threads {threads}");
            // Record 5's batch is completed after four others: BLASTN refuses it before the first
            // batch's archive, whether its reading keeps the records or reads them again batch by
            // batch, and whether its reports are written batch by batch or after the loop.
            let late = file("batch_title_late", query_kind);
            let paths: &[&[(&str, &str)]] = if program == "blastn" {
                &[
                    &[],
                    &[("LOSAT_X_QSTREAM", "0")],
                    &[("LOSAT_X_RPTSTREAM", "0")],
                ]
            } else {
                &[&[]]
            };
            for path in paths {
                let mut env = vec![("BATCH_SIZE", each)];
                env.extend_from_slice(path);
                let result = run(program, &late, &subject, &env, threads);
                refused(&what(&format!("late each {path:?}")), &result);
            }
            // The subjects are written with the first batch.
            let clean = file("batch_clean", query_kind);
            let result = run(
                program,
                &clean,
                &file("batch_title.s", subject_kind),
                &[],
                threads,
            );
            refused(&what("subject"), &result);
        }
    }
}

/// An example of the library API, built (with the profile of the binary) when a single test
/// target was asked for.
fn example_binary(name: &str) -> PathBuf {
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_losat"))
        .parent()
        .unwrap()
        .to_path_buf();
    let example = bin_dir
        .join("examples")
        .join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    if !example.exists() {
        let profile: &[&str] = if bin_dir.ends_with("release") {
            &["--release"]
        } else {
            &["--profile", "test"]
        };
        let status = Command::new(env!("CARGO"))
            .args(["build", "--example", name])
            .args(profile)
            .arg("--target-dir")
            .arg(bin_dir.parent().unwrap())
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .status()
            .expect("build the example");
        assert!(status.success(), "cargo build --example {name}");
    }
    example
}

/// The library API (`examples/library_command_line`, the command line's arguments through
/// `LOSAT::api::run`; a failure ends with status 1) reads the whole query before it searches, so
/// its batches meet no reader error; they stop where the command line's do otherwise, and its
/// refusal looks at the same records. (Its `-subject` entries lost their own checks of every
/// record.)
#[test]
fn the_library_refuses_only_the_titles_that_ncbi_writes() {
    let example = example_binary("library_command_line");
    let library = |program: &str, query: &Path, subject: &Path, env: &[(&str, &str)]| -> Run {
        let mut command = Command::new(&example);
        clear_environment(&mut command);
        command
            .arg(program)
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(subject);
        command.args(["-outfmt", "11", "-num_threads", "1"]);
        for (name, value) in env {
            command.env(name, value);
        }
        let output = command.output().expect("the example runs");
        Run {
            status: output.status.code(),
            stdout: output.stdout,
            stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        }
    };
    for (program, query_kind, subject_kind) in PROGRAMS {
        let (each, _) = sizes(query_kind);
        let what = |case: &str| format!("library {program} {case}");
        let subject = fixture(&format!(
            "base.s.{subject_kind}.{}",
            if subject_kind == "prot" { "faa" } else { "fa" }
        ));
        let bad_subject = file("batch_title.s", subject_kind);
        let clean = file("batch_clean", query_kind);
        let last = file("batch_title_lastempty", query_kind);
        // The last batch has no letters: NCBI fails it.
        let result = library(program, &last, &subject, &[("BATCH_SIZE", each)]);
        not_refused(
            &what("lastempty each"),
            &result,
            1,
            "BLAST engine error: Warning: Sequence contains no data",
            true,
        );
        // One batch with letters: its archive holds the title.
        refused(
            &what("lastempty default"),
            &library(program, &last, &subject, &[]),
        );
        // A batch size of 0 reads no query (BLASTN's is the mixer's).
        let result = library(program, &clean, &bad_subject, &[("BATCH_SIZE", "0")]);
        if program == "blastn" {
            refused(&what("subject, BATCH_SIZE=0"), &result);
        } else {
            not_refused(
                &what("subject, BATCH_SIZE=0"),
                &result,
                1,
                "BLAST engine error: Empty CBlastQueryVector",
                false,
            );
        }
        refused(
            &what("subject"),
            &library(program, &clean, &bad_subject, &[]),
        );
    }
}
