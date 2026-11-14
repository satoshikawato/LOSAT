//! The refusal of a title that the BLAST archive (outfmt 11) cannot write as NCBI does
//! (`report::archive::check_titles`, a byte outside `0x20..=0x7E` in text) looks only at the
//! records that NCBI BLAST+ 2.17.0 writes, on a whole run of the binary of each search program
//! with `-subject` (line NP, series SM, DECISIONS 24), and with `-db` without its header file,
//! whose batch with hits NCBI fails (DECISIONS 25).
//!
//! NCBI writes the archive of each query batch that it completes (`WriteArchive`: the batch's
//! queries and every `-subject` record); the batch that stops the run (a reader error, a last
//! batch without letters, `Empty CBlastQueryVector` with `BATCH_SIZE=0`) writes neither, and the
//! error archive holds no sequence. Where NCBI writes the title (with `#` and its `Bad char`
//! messages), LOSAT refuses it before the first byte; where NCBI does not, LOSAT writes NCBI's
//! output (the `arch_title*` rows of `tests/fasta_input_fixtures.py` and `db.archtitle.*` of
//! `tests/db_search_manifest/` hold NCBI's bytes; this test adds BLASTX's `-subject`, which
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
            // BLASTP's block search (LOSAT_X_BPQBLOCK, here a block of each batch) writes each
            // block's archive before it searches the next: a run with such a title (the second
            // record of `batch_title.s`, read without a reader error) takes the search of all
            // queries, which refuses it before the first byte.
            if program == "blastp" {
                let result = run(
                    program,
                    &file("batch_title.s", query_kind),
                    &subject,
                    &[("BATCH_SIZE", each), ("LOSAT_X_BPQBLOCK", "60")],
                    threads,
                );
                refused(&what("blocks"), &result);
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

/// The library API (`examples/library_command_line`, the command line's arguments through
/// `losat::api::run`; a failure ends with status 1) reads the whole query before it searches, so
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
        // The second batch's title (`batch_title.s` read as the query, one record a batch): the
        // first batch's archive is held until the refusal drops it.
        refused(
            &what("second batch"),
            &library(
                program,
                &file("batch_title.s", query_kind),
                &subject,
                &[("BATCH_SIZE", each)],
            ),
        );
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
/// written: LOSAT gives the run's bytes without that title, also where the run holds its output
/// for the later title (`archive::ArchiveHold`: the first batch's archive comes before the
/// warning in a stream that both go to). One in the first batch is written: LOSAT refuses it
/// before the first byte. The command line (`-num_threads` 1 and 4, apart and merged) and the
/// library API (`examples/library_command_line`, whose formats and diagnostics are separate
/// streams of the caller).
#[test]
fn a_title_in_the_batch_that_a_missing_header_file_fails_is_not_refused() {
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
            not_refused(
                &what("no title"),
                &clean,
                status,
                "BLAST Database error: Error: File (",
                true,
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
            refused(&what("u1"), &search("u1"));
        }
    }
}
