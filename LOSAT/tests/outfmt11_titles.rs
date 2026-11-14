//! Consistency checks of outfmt 11 (the BLAST archive) with a title that NCBI's serial stream
//! replaces bytes of (`café`: `#` for each byte in the archive and a `Bad char` report on
//! standard error per byte, as NCBI BLAST+ 2.17.0 writes it; line NP, series SM, DECISIONS 24
//! and 25, whose refusal of such a title is gone, line RG item K4b).
//!
//! NCBI's bytes of the plain command-line runs are frozen elsewhere: the `arch_title*` rows of
//! `tests/fasta_input_fixtures.py` (`-subject`, BLASTN, BLASTP, TBLASTN, TBLASTX), `db.archtitle.*`
//! of `tests/db_search_manifest/` (`-db`, all five programs), `tests/fixtures/strategy/cases.tsv`
//! (the search strategy) and `tests/archive_replaced_bytes.rs`. This file keeps what those do not
//! run, each as the same bytes (standard output, standard error and exit status) as the plain
//! command-line run of the same arguments:
//!
//! - BLASTX `-subject` (no frozen row has it): the archive of each batch that NCBI completes holds
//!   the title as `caf##` with its two `Bad char` reports, the batch that stops the run (a reader
//!   error, a last batch without letters, `Empty CBlastQueryVector` with `BATCH_SIZE=0`) writes
//!   neither, with `-num_threads` 1 and 4 alike;
//! - BLASTN's experimental reading and report paths (`LOSAT_X_QSTREAM=0`, `LOSAT_X_RPTSTREAM=0`)
//!   and BLASTP's block search (`LOSAT_X_BPQBLOCK`) give the plain run's bytes;
//! - the library API's `-subject` searches with a batch size (`examples/library_command_line`,
//!   which reads the whole query before it searches and ends a failure with status 1) give the
//!   command line's bytes where the command line meets no reader error;
//! - `-db` without its header file, whose batch with hits NCBI fails (DECISIONS 25; the message
//!   holds the absolute path, so no frozen row has it): a title in the batch that fails, or after
//!   it, is never written and changes no byte of the run, and a title in the first batch is written.
//!
//! The query files (`tests/fixtures/fasta_input/`): `batch_title_late` (record 5 has `café` in its
//! title, record 6 a data line that NCBI's reader refuses), `batch_title_first` (`café` on record
//! 1), `batch_title_lastempty` (records 5 to 7 without letters, record 7 with `café`),
//! `batch_late_bad` and `batch_clean`; the subject file `batch_title.s` has `café` on its second
//! record. A nucleotide record has 180 letters and a protein record 60: `BATCH_SIZE=100` (50)
//! makes a batch of each record, 200 (100) of two.

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::library_examples::example_binary;

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

/// A title with `café` (`utf8`: the bytes 0xC3 0xA9) or with DEL (`del`: 0x7F) in a text archive
/// is written as NCBI writes it (line RG item K4b; the other outfmt 11 tests call this): the
/// archive holds `#` for each byte and standard error one `Error: (803.7)` / `Bad char` report
/// for each byte, with the string. NCBI's bytes are frozen in `tests/archive_replaced/`.
pub fn assert_replaced(kind: &str, output: &std::process::Output, label: &str) {
    let (written, bytes): (&str, &[&str]) = if kind == "utf8" {
        ("caf##", &["0xC3", "0xA9"])
    } else {
        ("ti#tle", &["0x7F"])
    };
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stdout.contains(written),
        "{label}: no `{written}` in the archive"
    );
    assert_eq!(
        stderr.matches("Bad char [").count(),
        bytes.len(),
        "{label}: reports: {stderr}"
    );
    for byte in bytes {
        assert!(
            stderr.contains(&format!("Bad char [{byte}] in string at line ")),
            "{label}: {stderr}"
        );
    }
    assert!(
        stderr.contains("Error: (803.7) ["),
        "{label}: no report: {stderr}"
    );
    assert!(
        !stderr.contains("cannot write as NCBI does") && !stderr.contains("visible ASCII"),
        "{label}: refused: {stderr}"
    );
}

/// A `-db` name with a byte that the stream replaces names a database that does not exist: the
/// error archive and the error line of any missing database (status 2), not a refusal of the name.
pub fn assert_missing_database(output: &std::process::Output, label: &str) {
    assert_eq!(output.status.code(), Some(2), "{label}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with("BLAST Database error: No alias or index file found for"),
        "{label}: {stderr}"
    );
    assert!(
        output.stdout.starts_with(b"Blast4-archive ::= {"),
        "{label}: standard output"
    );
    assert!(!stderr.contains("-db name"), "{label}: {stderr}");
}

/// The standard output, standard error and exit status of two runs are the same.
fn same(what: &str, left: &Run, right: &Run) {
    assert_eq!(left.status, right.status, "{what}: status");
    assert!(
        left.stdout == right.stdout,
        "{what}: standard output ({} and {} bytes)",
        left.stdout.len(),
        right.stdout.len()
    );
    assert_eq!(left.stderr, right.stderr, "{what}: standard error");
}

/// How many times the archive of `run` holds the title `café` (`caf##`) and how many `Bad char`
/// reports the run made (two for each written title, one for each byte of `é`).
fn titles(run: &Run) -> (usize, usize) {
    (
        String::from_utf8_lossy(&run.stdout)
            .matches("caf##")
            .count(),
        run.stderr.matches("Bad char [").count(),
    )
}

/// The batch sizes and files of a run: `(case, query file, subject file, env)`.
type Case<'a> = (&'a str, PathBuf, PathBuf, Vec<(&'a str, &'a str)>);

/// The runs of one program that put the title of `batch_title_*` or of `batch_title.s` in a batch
/// that NCBI completes (`written`) or does not, with NCBI's exit status (the frozen rows of
/// `tests/fasta_input_fixtures/` have them): `(case, status, titles written)`.
const BLASTX_EXPECTED: [(&str, i32, usize); 11] = [
    ("late two", 1, 0),
    ("late each", 1, 1),
    ("late default", 1, 0),
    ("first two", 1, 1),
    ("first default", 1, 0),
    ("first zero", 3, 0),
    ("lastempty each", 3, 0),
    ("lastempty default", 0, 1),
    ("subject after error", 1, 0),
    ("subject clean", 0, 1),
    ("subject zero", 3, 0),
];

fn cases<'a>(
    query_kind: &str,
    subject_kind: &str,
    subject: &Path,
    each: &'a str,
    two: &'a str,
) -> Vec<Case<'a>> {
    let bad_subject = file("batch_title.s", subject_kind);
    let subject = subject.to_path_buf();
    let q = |name: &str| file(name, query_kind);
    vec![
        (
            "late two",
            q("batch_title_late"),
            subject.clone(),
            vec![("BATCH_SIZE", two)],
        ),
        (
            "late each",
            q("batch_title_late"),
            subject.clone(),
            vec![("BATCH_SIZE", each)],
        ),
        (
            "late default",
            q("batch_title_late"),
            subject.clone(),
            vec![],
        ),
        (
            "first two",
            q("batch_title_first"),
            subject.clone(),
            vec![("BATCH_SIZE", two)],
        ),
        (
            "first default",
            q("batch_title_first"),
            subject.clone(),
            vec![],
        ),
        (
            "first zero",
            q("batch_title_first"),
            subject.clone(),
            vec![("BATCH_SIZE", "0")],
        ),
        (
            "lastempty each",
            q("batch_title_lastempty"),
            subject.clone(),
            vec![("BATCH_SIZE", each)],
        ),
        (
            "lastempty default",
            q("batch_title_lastempty"),
            subject,
            vec![],
        ),
        (
            "subject after error",
            q("batch_late_bad"),
            bad_subject.clone(),
            vec![],
        ),
        (
            "subject clean",
            q("batch_clean"),
            bad_subject.clone(),
            vec![],
        ),
        (
            "subject zero",
            q("batch_clean"),
            bad_subject,
            vec![("BATCH_SIZE", "0")],
        ),
    ]
}

fn base_subject(subject_kind: &str) -> PathBuf {
    fixture(&format!(
        "base.s.{subject_kind}.{}",
        if subject_kind == "prot" { "faa" } else { "fa" }
    ))
}

#[test]
fn blastx_subject_writes_the_titles_of_the_batches_that_ncbi_completes() {
    let (_, query_kind, subject_kind) = PROGRAMS[2];
    let (each, two) = sizes(query_kind);
    let subject = base_subject(subject_kind);
    let cases = cases(query_kind, subject_kind, &subject, each, two);
    assert_eq!(cases.len(), BLASTX_EXPECTED.len());
    for ((name, query, subject, env), (expected, status, written)) in
        cases.iter().zip(BLASTX_EXPECTED)
    {
        assert_eq!(*name, expected);
        let what = |threads: &str| format!("blastx {name} threads {threads}");
        let one = run("blastx", query, subject, env, "1");
        assert_eq!(one.status, Some(status), "{}: {}", what("1"), one.stderr);
        assert_eq!(titles(&one), (written, 2 * written), "{}", what("1"));
        assert!(
            !one.stderr.contains("cannot write as NCBI does"),
            "{}: {}",
            what("1"),
            one.stderr
        );
        same(&what("4"), &one, &run("blastx", query, subject, env, "4"));
    }
}

/// The experimental reading and report paths of BLASTN write the plain run's bytes (the title of
/// record 5 is written after four other batches, whether the reading keeps the records or reads
/// them again batch by batch, and whether the reports are written batch by batch or after the
/// loop), and so does BLASTP's block search (`LOSAT_X_BPQBLOCK`, here a block of each batch,
/// which writes each block's archive before it searches the next).
#[test]
fn experimental_paths_write_the_plain_runs_bytes() {
    let (program, query_kind, subject_kind) = PROGRAMS[0];
    let (each, two) = sizes(query_kind);
    let subject = base_subject(subject_kind);
    for threads in ["1", "4"] {
        for (name, query, subject, env) in cases(query_kind, subject_kind, &subject, each, two) {
            let plain = run(program, &query, &subject, &env, threads);
            for switch in ["LOSAT_X_QSTREAM", "LOSAT_X_RPTSTREAM"] {
                let mut with = env.clone();
                with.push((switch, "0"));
                same(
                    &format!("{program} {name} {switch}=0 threads {threads}"),
                    &plain,
                    &run(program, &query, &subject, &with, threads),
                );
            }
        }
    }
    let (program, query_kind, subject_kind) = PROGRAMS[1];
    let (each, _) = sizes(query_kind);
    let subject = base_subject(subject_kind);
    let query = file("batch_title.s", query_kind);
    for threads in ["1", "4"] {
        let plain = run(program, &query, &subject, &[("BATCH_SIZE", each)], threads);
        let blocks = run(
            program,
            &query,
            &subject,
            &[("BATCH_SIZE", each), ("LOSAT_X_BPQBLOCK", "60")],
            threads,
        );
        assert_eq!(titles(&plain).0, 1, "blastp blocks: the title is written");
        same(&format!("blastp blocks threads {threads}"), &plain, &blocks);
    }
}

/// The library API (`examples/library_command_line`, the command line's arguments through
/// `losat::api::run`; a failure ends with status 1) reads the whole query before it searches, so
/// its batches meet no reader error; where the command line meets none either, the two give the
/// same bytes (and a failing command line a failing library run).
#[test]
fn the_library_writes_the_command_lines_bytes() {
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
    // The cases whose command line meets no reader error.
    const NO_READER_ERROR: [&str; 5] = [
        "lastempty each",
        "lastempty default",
        "subject clean",
        "subject zero",
        "second batch",
    ];
    for (program, query_kind, subject_kind) in PROGRAMS {
        let (each, two) = sizes(query_kind);
        let subject = base_subject(subject_kind);
        let mut cases = cases(query_kind, subject_kind, &subject, each, two);
        // The second batch's title (`batch_title.s` read as the query, one record a batch).
        cases.push((
            "second batch",
            file("batch_title.s", query_kind),
            subject.clone(),
            vec![("BATCH_SIZE", each)],
        ));
        for (name, query, subject, env) in cases {
            if !NO_READER_ERROR.contains(&name) {
                continue;
            }
            let what = format!("library {program} {name}");
            let line = run(program, &query, &subject, &env, "1");
            let api = library(program, &query, &subject, &env);
            assert_eq!(
                api.status,
                Some(if line.status == Some(0) { 0 } else { 1 }),
                "{what}: status ({:?} on the command line), stderr {}",
                line.status,
                api.stderr
            );
            assert!(
                api.stdout == line.stdout,
                "{what}: standard output ({} and {} bytes)",
                api.stdout.len(),
                line.stdout.len()
            );
            assert_eq!(api.stderr, line.stderr, "{what}: standard error");
        }
    }
}

/// A temporary directory of this test, removed at its end.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-outfmt11-titles-{tag}-{}",
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

/// `count` proteins of `length` residues (without X or a stop) from a fixed linear congruential
/// generator, so that each hits only itself.
fn random_proteins(count: usize, length: usize) -> Vec<String> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state: u64 = 20_261_108;
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

/// Writes the FASTA file `path` of the records `(id, title, sequence)`.
fn write_fasta(path: &Path, records: &[(&str, String, String)]) {
    let text: String = records
        .iter()
        .map(|(id, title, sequence)| format!(">{id} {title}\n{sequence}\n"))
        .collect();
    std::fs::write(path, text).expect("write a FASTA file");
}

/// The version 4 database `name` of the sequences `sequences` (`s0`, `s1`, ...) made by
/// `LOSAT makeblastdb` in `dir`, without its header file (`.phr` or `.nhr`).
fn database_without_headers(dir: &Path, name: &str, dbtype: &str, sequences: &[String]) -> PathBuf {
    let fasta = dir.join(format!("{name}.fa"));
    let ids = ["s0", "s1", "s2"];
    let records: Vec<(&str, String, String)> = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (ids[index], format!("subject {index}"), sequence.clone()))
        .collect();
    write_fasta(&fasta, &records);
    let db = dir.join(name);
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(&fasta)
        .args(["-dbtype", dbtype, "-blastdb_version", "4", "-out"])
        .arg(&db)
        .output()
        .expect("run LOSAT makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let headers = dir.join(format!("{name}.{}hr", &dbtype[..1]));
    std::fs::remove_file(&headers).expect("remove the header file");
    db
}

/// `-db` without its header file (line NP, series SM, DECISIONS 25): NCBI reads a definition
/// line first where it makes the Seq-aligns of a batch with hits, and fails that batch (`BLAST
/// Database error: Error: File (<path>) not found.`, status 2; the path is absolute, so no row of
/// `tests/db_search_manifest/` freezes it) after the archives of the batches before. The query
/// files hold two records of X (protein) or N (nucleotide) only, which have no valid context (NCBI
/// does not search their batch, reads no definition line for it and writes its archive), then two
/// records that hit the database, the first with a hyphen (the reader's warning, read with its
/// batch); `BATCH_SIZE` 60 (180) makes a batch of the first two and one of each record after. A
/// title that the archive cannot write in the batch that fails, or in the batch after it, is never
/// written: the run's bytes are those of the run without that title, also where the run holds its
/// output for the later title (`archive::ArchiveHold`: the first batch's archive comes before the
/// warning in a stream that both go to). One in the first batch is written, with its `Bad char`
/// reports. The command line (`-num_threads` 1 and 4, apart and merged) and the library API
/// (`examples/library_command_line`, whose formats and diagnostics are separate streams of the
/// caller).
#[test]
fn a_title_in_the_batch_that_a_missing_header_file_fails_is_not_written() {
    let dir = TempDir::new("nohdr");
    let proteins = random_proteins(2, 150);
    let nucleotides: Vec<String> = proteins
        .iter()
        .map(|protein| back_translate(protein))
        .collect();
    let protein_db = database_without_headers(&dir.0, "sp", "prot", &proteins);
    let nucleotide_db = database_without_headers(&dir.0, "sn", "nucl", &nucleotides);
    // The query files of each kind: `café` on no record, on u1 (the first batch) or on q3 (the
    // batch that fails).
    let query = |kind: &str, bad: &str| -> PathBuf {
        let (letter, length, hits) = if kind == "prot" {
            ('X', 30, &proteins)
        } else {
            ('N', 90, &nucleotides)
        };
        let title = |id: &str| {
            if id == bad {
                "query café".to_string()
            } else {
                "query".to_string()
            }
        };
        let unsearched: String = std::iter::repeat_n(letter, length).collect();
        // A hyphen in the middle of q3: the reader warns and ignores it.
        let (head, tail) = hits[0].split_at(hits[0].len() / 2);
        let hyphened = format!("{head}-{tail}");
        let path = dir.0.join(format!("q.{kind}.{bad}.fa"));
        write_fasta(
            &path,
            &[
                ("u1", title("u1"), unsearched.clone()),
                ("u2", title("u2"), unsearched),
                ("q3", title("q3"), hyphened),
                ("q4", title("q4"), hits[1].clone()),
            ],
        );
        path
    };
    let example = example_binary("library_command_line");
    for (program, query_kind, subject_kind) in PROGRAMS {
        let db = if subject_kind == "prot" {
            &protein_db
        } else {
            &nucleotide_db
        };
        let batch_size = if query_kind == "prot" { "60" } else { "180" };
        for arm in ["1", "4", "library"] {
            let what = |case: &str| format!("{program} {case} ({arm})");
            let command = |bad: &str| -> Command {
                let mut command = if arm == "library" {
                    let mut command = Command::new(&example);
                    clear_environment(&mut command);
                    command
                } else {
                    losat()
                };
                command
                    .arg(program)
                    .arg("-query")
                    .arg(query(query_kind, bad))
                    .arg("-db")
                    .arg(db)
                    .args(["-outfmt", "11", "-num_threads"])
                    .arg(if arm == "library" { "1" } else { arm })
                    .env("BATCH_SIZE", batch_size);
                command
            };
            let search = |bad: &str| -> Run {
                let output = command(bad).output().expect("losat runs");
                Run {
                    status: output.status.code(),
                    stdout: output.stdout,
                    stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
                }
            };
            // Standard output and standard error into one file, as `2>&1` writes them.
            let merged = |bad: &str| -> Vec<u8> {
                let path = dir.0.join(format!("merged.{program}.{arm}.{bad}"));
                let file = std::fs::File::create(&path).expect("create the merged output");
                let status = command(bad)
                    .stdout(file.try_clone().expect("clone the merged output"))
                    .stderr(file)
                    .status()
                    .expect("losat runs");
                assert!(
                    !status.success(),
                    "{program} {bad} ({arm}): the run succeeded"
                );
                std::fs::read(&path).expect("read the merged output")
            };
            let clean = search("none");
            // The library API ends a failure with status 1.
            let status = if arm == "library" { 1 } else { 2 };
            assert_eq!(
                clean.status,
                Some(status),
                "{}: {}",
                what("no title"),
                clean.stderr
            );
            assert!(
                clean.stderr.contains("BLAST Database error: Error: File ("),
                "{}: stderr {}",
                what("no title"),
                clean.stderr
            );
            assert_eq!(titles(&clean), (0, 0), "{}", what("no title"));
            // The first batch's archive and the error archive.
            assert!(
                String::from_utf8_lossy(&clean.stdout)
                    .matches("Blast4-archive ::=")
                    .count()
                    > 1,
                "{}: no archive of the first batch",
                what("no title")
            );
            assert!(
                clean.stderr.contains("Hyphens are invalid"),
                "{}: no warning of the reader: {}",
                what("no title"),
                clean.stderr
            );
            for bad in ["q3", "q4"] {
                let failed = search(bad);
                assert_eq!(
                    (failed.status, &failed.stdout, &failed.stderr),
                    (clean.status, &clean.stdout, &clean.stderr),
                    "{}: a title in the batch that fails, or after it, changed the run",
                    what(bad)
                );
                if arm != "library" {
                    assert!(
                        merged(bad) == merged("none"),
                        "{}: the merged stream changed",
                        what(bad)
                    );
                }
            }
            let first = search("u1");
            assert_eq!(first.status, clean.status, "{}: status", what("u1"));
            assert_eq!(titles(&first), (1, 2), "{}: {}", what("u1"), first.stderr);
        }
    }
}
