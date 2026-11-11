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

/// The reports that a run wrote: the Iterations of outfmt 5, the documents of outfmt 15 and 16 on
/// the output, the side files of 13 and 14.
fn reports(run: &Run, outfmt: &str) -> usize {
    let stdout = String::from_utf8_lossy(&run.stdout);
    match outfmt {
        "5" => stdout.matches("<Iteration>").count(),
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
        let queries = Records::of(query_kind);
        let last = write(
            &dir.0,
            &format!("{program}.last"),
            &format!(">qa plain\n{}\n>qb caf\u{e9}\n{}\n", queries.a, queries.b),
        );
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
                // Without an error: record 2's report (`batch_title.s` read as the query) is
                // written when record 3's results arrive, and the last record's in the epilog.
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
                refused(&what("a later record, no error"), &result);
                let result = run(
                    search(program, &last, &subject, &[("BATCH_SIZE", each)], threads),
                    outfmt,
                    &dir.0,
                );
                refused(&what("the last record, no error"), &result);
            }
        }
    }
}

/// `count` proteins of `length` residues (without X or a stop) from a fixed linear congruential
/// generator, so that each hits only itself.
fn random_proteins(count: usize, length: usize) -> Vec<String> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state: u64 = 20_261_111;
    (0..count)
        .map(|_| {
            (0..length)
                .map(|_| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    RESIDUES[(state >> 33) as usize % RESIDUES.len()] as char
                })
                .collect()
        })
        .collect()
}

/// One codon of each residue of `protein`.
fn back_translate(protein: &str) -> String {
    protein
        .chars()
        .map(|residue| match residue {
            'A' => "GCT",
            'C' => "TGT",
            'D' => "GAT",
            'E' => "GAA",
            'F' => "TTT",
            'G' => "GGT",
            'H' => "CAT",
            'I' => "ATT",
            'K' => "AAA",
            'L' => "CTG",
            'M' => "ATG",
            'N' => "AAT",
            'P' => "CCG",
            'Q' => "CAG",
            'R' => "CGT",
            'S' => "TCT",
            'T' => "ACC",
            'V' => "GTT",
            'W' => "TGG",
            'Y' => "TAT",
            _ => unreachable!("a residue of `random_proteins`"),
        })
        .collect()
}

/// The records of the tests below, as `kind` (`prot` or `nuc`) writes them: `a` and `b` (120
/// residues each, or their 360 nucleotides) hit only themselves; `$` is a data line that NCBI's
/// reader refuses (`READER_ERROR`).
struct Records {
    a: String,
    b: String,
}

impl Records {
    fn of(kind: &str) -> Self {
        let proteins = random_proteins(2, 120);
        let letters = |protein: &str| {
            if kind == "prot" {
                protein.to_string()
            } else {
                back_translate(protein)
            }
        };
        Records {
            a: letters(&proteins[0]),
            b: letters(&proteins[1]),
        }
    }
}

const STOP: &str = ">e9 stop\n$$$$$$$$$$\nACGT\n";

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).expect("write a FASTA file");
    path
}

#[test]
fn a_subject_title_that_only_the_held_report_shows_is_not_refused() {
    let dir = TempDir::new("subject");
    for (program, query_kind, subject_kind) in PROGRAMS {
        let each = each(query_kind);
        let (queries, subjects) = (Records::of(query_kind), Records::of(subject_kind));
        // Subject s2's title cannot be written; only query b hits s2.
        let subject = write(
            &dir.0,
            &format!("{program}.s"),
            &format!(">s1 plain\n{}\n>s2 caf\u{e9}\n{}\n", subjects.a, subjects.b),
        );
        let held = write(
            &dir.0,
            &format!("{program}.held"),
            &format!(">qa plain\n{}\n>qb plain\n{}\n{STOP}", queries.a, queries.b),
        );
        let written = write(
            &dir.0,
            &format!("{program}.written"),
            &format!(">qb plain\n{}\n>qa plain\n{}\n{STOP}", queries.b, queries.a),
        );
        for outfmt in ["14", "15", "16"] {
            for threads in ["1", "4"] {
                let what = |case: &str| format!("{program} {case} outfmt {outfmt} n{threads}");
                // Query b's report is held back and never written: s2 is not shown. (The
                // e-value cuts query a's random hits on s2.)
                let mut command =
                    search(program, &held, &subject, &[("BATCH_SIZE", each)], threads);
                command.args(["-evalue", "1e-10"]);
                let result = run(command, outfmt, &dir.0);
                not_refused(&what("held"), &result, outfmt, 1, READER_ERROR, 1);
                // Query b's report is written when query a's results arrive.
                let mut command = search(
                    program,
                    &written,
                    &subject,
                    &[("BATCH_SIZE", each)],
                    threads,
                );
                command.args(["-evalue", "1e-10"]);
                let result = run(command, outfmt, &dir.0);
                refused(&what("written"), &result);
            }
        }
    }
}

const NO_DATA: &str = "BLAST engine error: Warning: Sequence contains no data";

#[test]
fn the_title_before_a_last_batch_without_letters_is_not_refused() {
    let dir = TempDir::new("lastempty");
    for (program, query_kind, subject_kind) in PROGRAMS {
        let each = each(query_kind);
        let subject = fixture(&format!(
            "base.s.{subject_kind}.{}",
            extension(subject_kind)
        ));
        let queries = Records::of(query_kind);
        // The last batch has no letters: NCBI fails it, and the report of qb before it is held
        // back and never written.
        let query = write(
            &dir.0,
            &format!("{program}.q"),
            &format!(
                ">qa plain\n{}\n>qb caf\u{e9}\n{}\n>qc no letters\n",
                queries.a, queries.b
            ),
        );
        for outfmt in ["14", "15", "16"] {
            for threads in ["1", "4"] {
                let what = format!("{program} outfmt {outfmt} n{threads}");
                let result = run(
                    search(program, &query, &subject, &[("BATCH_SIZE", each)], threads),
                    outfmt,
                    &dir.0,
                );
                not_refused(&what, &result, outfmt, 3, NO_DATA, 1);
            }
        }
        // A title in the batch that fails (`batch_title_lastempty`: records 5 to 7 without
        // letters, record 7 with `café`) is never written by any format: outfmt 5 writes records
        // 1 to 4 and the error, outfmt 16 records 1 to 3 (record 4's report is held back).
        for outfmt in ["5", "16"] {
            let result = run(
                search(
                    program,
                    &file("batch_title_lastempty", query_kind),
                    &subject,
                    &[("BATCH_SIZE", each)],
                    "1",
                ),
                outfmt,
                &dir.0,
            );
            let what = format!("{program} title in the failing batch, outfmt {outfmt}");
            let written = if outfmt == "5" { 4 } else { 3 };
            not_refused(&what, &result, outfmt, 3, NO_DATA, written);
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
/// its batches meet no reader error; a last batch without letters stops it as it stops the
/// command line. Its refusals look at the same reports and come before the first byte (BLASTX's
/// library wrote the start of the document, or the reports of the batches before, first).
#[test]
fn the_library_refuses_only_the_titles_that_ncbi_writes_before_the_first_byte() {
    let example = example_binary("library_command_line");
    let dir = TempDir::new("library");
    for (program, query_kind, subject_kind) in PROGRAMS {
        let each = each(query_kind);
        let base = fixture(&format!(
            "base.s.{subject_kind}.{}",
            extension(subject_kind)
        ));
        let (queries, subjects) = (Records::of(query_kind), Records::of(subject_kind));
        let lastempty = write(
            &dir.0,
            &format!("{program}.lastempty"),
            &format!(
                ">qa plain\n{}\n>qb caf\u{e9}\n{}\n>qc no letters\n",
                queries.a, queries.b
            ),
        );
        let bad_subject = write(
            &dir.0,
            &format!("{program}.s"),
            &format!(">s1 plain\n{}\n>s2 caf\u{e9}\n{}\n", subjects.a, subjects.b),
        );
        let clean = write(
            &dir.0,
            &format!("{program}.clean"),
            &format!(">qa plain\n{}\n>qb plain\n{}\n", queries.a, queries.b),
        );
        let library = |query: &Path, subject: &Path, outfmt: &str| {
            let mut command = Command::new(&example);
            clear_environment(&mut command);
            command
                .arg(program)
                .arg("-query")
                .arg(query)
                .arg("-subject")
                .arg(subject)
                .args(["-num_threads", "1"])
                .env("BATCH_SIZE", each);
            run(command, outfmt, &dir.0)
        };
        for outfmt in ["15", "16"] {
            let what = |case: &str| format!("library {program} {case} outfmt {outfmt}");
            not_refused(
                &what("lastempty"),
                &library(&lastempty, &base, outfmt),
                outfmt,
                1,
                NO_DATA,
                1,
            );
            // The second batch's title, written in the epilog.
            refused(
                &what("second batch"),
                &library(&file("batch_title.s", query_kind), &base, outfmt),
            );
            // A subject's title, shown by query b's report.
            refused(&what("subject"), &library(&clean, &bad_subject, outfmt));
        }
    }
}

/// The library API takes several formats in one call. Outfmt 5 writes the report that a
/// `-subject` XML2 or JSON document holds back, so with outfmt 5 among the formats that report's
/// title is refused before either format writes a byte; with outfmt 16 alone it is not.
/// (BLASTP's default batches of 10,000 letters make a batch of each 10,000-residue query; the
/// last batch has no letters, which NCBI fails.)
#[test]
fn a_library_call_with_outfmt_5_refuses_the_held_title_before_any_format_writes() {
    use LOSAT::api::{read_inputs, run, validate, Inputs, Program};
    let proteins = random_proteins(3, 10_000);
    let query = format!(
        ">qa plain\n{}\n>qb caf\u{e9}\n{}\n>qc no letters\n",
        proteins[0], proteins[1]
    );
    let subject = format!(">s1 plain\n{}\n", &proteins[0][..300]);
    let search = validate(Program::Blastp, &[]).expect("blastp");
    let (queries, subjects) =
        read_inputs(&search, query.as_bytes(), subject.as_bytes()).expect("the records");
    let call = |formats: &[&str]| {
        let mut sinks: Vec<Vec<u8>> = vec![Vec::new(); formats.len()];
        let mut diagnostics = Vec::new();
        let result = {
            let mut reports: Vec<(&str, &mut (dyn std::io::Write + Send))> = formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(&format, sink)| (format, sink as &mut (dyn std::io::Write + Send)))
                .collect();
            run(
                &search,
                Inputs::new(&queries, &subjects),
                &mut reports,
                &mut diagnostics,
            )
        };
        (result, sinks)
    };
    let (result, sinks) = call(&["16"]);
    let message = result.expect_err("NCBI's error").message().to_string();
    assert!(message.contains("Sequence contains no data"), "{message}");
    assert_eq!(
        String::from_utf8_lossy(&sinks[0])
            .matches("\n<BlastOutput2>")
            .count(),
        1,
        "16 alone: qa's report"
    );
    for formats in [["16", "5"], ["5", "16"]] {
        let (result, sinks) = call(&formats);
        let message = result.expect_err("the refusal").message().to_string();
        assert!(message.contains(REFUSAL), "{formats:?}: {message}");
        assert!(
            sinks.iter().all(Vec::is_empty),
            "{formats:?}: wrote {:?} bytes",
            sinks.iter().map(Vec::len).collect::<Vec<_>>()
        );
    }
}
