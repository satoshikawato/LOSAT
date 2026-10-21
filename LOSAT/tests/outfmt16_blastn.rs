#![allow(warnings, clippy::all)]

//! BLAST XML2 and JSON (outfmt 16 and 15) of BLASTN: what the frozen-fixture rows cannot show.
//! The `-subject` report that NCBI holds back one query (a failed batch loses it), the documents
//! of a run with no query, and the searches LOSAT refuses before the first byte (a byte >= 0x80 in
//! a title or the `-db` name, DECISIONS 19 Q2; a taxonomy database scientific name, DECISIONS 30)
//! and the sink of a file. The bytes of NCBI BLAST+ 2.17.0 are in `tests/fixtures/outfmt16`
//! (`badbatch_blastn_{sub,db}_{15,16}.out`: `BATCH_SIZE=100 blastn -query qnucl_bad.fna -subject
//! snucl.fna | -db db/snucl -num_threads 1 -outfmt N`; `unsearched_blastn_*`: the same with
//! `qnucl_unsearched.fna`; `emptypipe_blastn_*` and `emptyfile_blastn_*`, run in that directory).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2363-2404, :2206-2219
//! (the `-subject` result is written when the next query's arrives, or by `PrintEpilog`, which a
//! failed batch does not reach).

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/outfmt16")
}

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt16_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn command(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_X_DEFAULT" {
            command.env_remove(name);
        }
    }
    // NCBI's variables change a report (and the search path): none leaks in
    for name in ["BATCH_SIZE", "CHUNK_SIZE", "BLASTDB"] {
        command.env_remove(name);
    }
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt16_worker.log"),
    );
    command.envs(env.iter().copied());
    command
}

fn run_in(dir: &Path, args: &[&str], env: &[(&str, &str)]) -> Output {
    command(env)
        .current_dir(dir)
        .args(args)
        .output()
        .expect("run LOSAT")
}

fn run_with_stdin(args: &[&str], stdin: &[u8]) -> Output {
    use std::io::Write;
    let mut child = command(&[])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    child.stdin.take().unwrap().write_all(stdin).unwrap();
    child.wait_with_output().expect("wait for LOSAT")
}

const SUBJECT: &str = "snucl.fna";
const DB: &str = "db/snucl";

/// The command of a frozen case: `blastn -query QUERY (-subject snucl.fna | -db db/snucl)`, run in
/// `tests/fixtures/outfmt16` (the db name is part of the bytes).
fn args<'a>(query: &'a str, mode: &str, threads: &'a str, format: &'a str) -> Vec<&'a str> {
    let mut args = vec!["blastn", "-query", query];
    args.extend(if mode == "sub" {
        ["-subject", SUBJECT]
    } else {
        ["-db", DB]
    });
    args.extend(["-num_threads", threads, "-outfmt", format]);
    args
}

fn frozen(name: &str) -> String {
    String::from_utf8(fs::read(fixtures().join(name)).unwrap()).unwrap()
}

/// Design §5.2 row 22 for BLASTN (`badbatch_blastn_*`): with `BATCH_SIZE=100` every query is a
/// batch and the third one fails to read. `-db` has written the reports of the two queries
/// before it; `-subject` the first one's only (the second is held back and lost), and neither
/// has an epilog. The document of a failed run is the same in a `-out` file.
#[test]
fn a_failed_batch_keeps_the_db_report_and_loses_the_held_subject_report() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("badbatch_blastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let output = run_in(
                    &dir,
                    &args("qnucl_bad.fna", mode, threads, format),
                    &[("BATCH_SIZE", "100")],
                );
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 13"),
                    "{mode} {format}: {stderr}"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format} threads {threads}"
                );
            }
            // a file sink: the same document (its start is written with the first report, or
            // by the error path when no report was written)
            let scratch = temp_dir("failed_out");
            let out = scratch.join("report.out");
            let mut args = args("qnucl_bad.fna", mode, "1", format);
            args.extend(["-out", out.to_str().unwrap()]);
            let output = run_in(&dir, &args, &[("BATCH_SIZE", "100")]);
            assert_eq!(output.status.code(), Some(1), "{mode} {format} -out");
            assert!(output.stdout.is_empty());
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                expected,
                "{mode} {format} -out"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

/// A first batch that fails to read (the whole file is one batch, `firstbad_blastn_*`): NCBI has
/// written the prolog (`PrintProlog` comes before the first batch is read) and nothing else, in
/// both modes and with or without `-out`.
#[test]
fn a_first_batch_that_fails_to_read_leaves_the_prolog_alone() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("firstbad_blastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let output = run_in(&dir, &args("qnucl_bad.fna", mode, threads, format), &[]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 13"),
                    "{mode} {format}: {stderr}"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format} threads {threads}"
                );
            }
            let scratch = temp_dir("firstbad_out");
            let out = scratch.join("report.out");
            let mut args = args("qnucl_bad.fna", mode, "1", format);
            args.extend(["-out", out.to_str().unwrap()]);
            let output = run_in(&dir, &args, &[]);
            assert_eq!(output.status.code(), Some(1), "{mode} {format} -out");
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                expected,
                "{mode} {format} -out"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

/// A query of N only is a batch NCBI does not search when `BATCH_SIZE` puts it alone: its `Search`
/// has no `hits` member and ` \nNo hits found` after the message (`unsearched_blastn_*`:
/// `BATCH_SIZE=100`, queries `qn1`, `inv1`, `qn2`, `inv2`), while a query of N in a batch with a
/// valid one has the message alone (the frozen row `sa.x16.blastn.invalid`).
#[test]
fn a_query_of_an_unsearched_batch_has_no_hits_member_and_the_no_hits_suffix() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("unsearched_blastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let output = run_in(
                    &dir,
                    &args("qnucl_unsearched.fna", mode, threads, format),
                    &[("BATCH_SIZE", "100")],
                );
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{mode} {format}: {stderr}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format} threads {threads}"
                );
            }
        }
    }
}

/// A run with no query formatted writes the prolog and the epilog (`PrintProlog`, `PrintEpilog`;
/// `emptypipe_blastn_*`); an empty query file is `Query is Empty!` and no bytes
/// (`emptyfile_blastn_*`).
#[test]
fn an_empty_query_pipe_writes_the_prolog_and_the_epilog() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("emptypipe_blastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let (subject, db) = (
                    dir.join(SUBJECT).to_str().unwrap().to_string(),
                    dir.join(DB).to_str().unwrap().to_string(),
                );
                let mut args = vec!["blastn", "-query", "-"];
                args.extend(if mode == "sub" {
                    ["-subject", subject.as_str()]
                } else {
                    ["-db", db.as_str()]
                });
                args.extend(["-num_threads", threads, "-outfmt", format]);
                let output = run_with_stdin(&args, b"");
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{mode} {format}: {stderr}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format}"
                );
            }
            // an empty query file
            let scratch = temp_dir("emptyfile");
            let empty = scratch.join("empty.fna");
            fs::write(&empty, b"").unwrap();
            let output = run_in(&dir, &args(empty.to_str().unwrap(), mode, "1", format), &[]);
            assert_eq!(output.status.code(), Some(0));
            assert!(output.stdout.is_empty(), "Query is Empty!: no bytes");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                frozen(&format!("emptyfile_blastn_{mode}_{format}.err"))
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

fn fasta(path: &Path, records: &[(&[u8], &[u8])]) {
    let mut text = Vec::new();
    for (title, sequence) in records {
        text.extend_from_slice(b">");
        text.extend_from_slice(title);
        text.push(b'\n');
        text.extend_from_slice(sequence);
        text.push(b'\n');
    }
    fs::write(path, text).expect("write FASTA");
}

const DNA: &[u8] = b"AACAGCACATACCGCGAATTGTTAGTCCTAATGTGGATGGAACACTCCGTGCTAAAAGCAAATTTTGCCGTGCTTCGTAAAGTGATGGACGACCCGAAAATACGGGATGGGCCCATGGCTGTAATGAAGGGCTGGGAGAACTATTGCTGGCC";
const OTHER_DNA: &[u8] = b"TCATGCTATCTATGGCACAAATGGACATTATTTTTGGAAGGTCGTTAGTGTCCGAGGGTCTCATAGATCCAGTTTAACGCCTTGGATTGATCATGCTATCTATGGCACAAATGGACATTATTTT";

fn refused(output: &Output, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTN"
        )),
        "{stderr}"
    );
}

/// Q2: a byte >= 0x80 in a title or the `-db` name is refused for 15 and 16 as for 5, before the
/// first byte (NCBI cuts the 16 document at the byte, and writes the 15 one in UTF-8 after
/// guessing the encoding).
#[test]
fn a_title_or_database_name_outside_us_ascii_is_refused_in_15_and_16() {
    let dir = temp_dir("q2");
    let (query, subject) = (dir.join("q.fna"), dir.join("s.fna"));
    for format in ["15", "16"] {
        for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
            for position in 0..2 {
                let mut records: Vec<(&[u8], &[u8])> = vec![(b"q1 ok", DNA), (b"q2 ok", DNA)];
                records[position] = (title, DNA);
                fasta(&query, &records);
                fasta(&subject, &[(b"s1 subject", DNA)]);
                let output = run_in(
                    &dir,
                    &[
                        "blastn",
                        "-query",
                        query.to_str().unwrap(),
                        "-subject",
                        subject.to_str().unwrap(),
                        "-outfmt",
                        format,
                    ],
                    &[],
                );
                refused(&output, "a query title");
            }
        }
        // a subject title with a hit
        fasta(&query, &[(b"q1 ok", DNA)]);
        fasta(
            &subject,
            &[(b"s1 ok", OTHER_DNA), ("s2 caf\u{e9}".as_bytes(), DNA)],
        );
        let output = run_in(
            &dir,
            &[
                "blastn",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        refused(&output, "a subject title");
    }
    // the -db name
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("x16_ntaxid_caf\u{e9}");
    for extension in ["nhr", "nin", "nsq"] {
        fs::copy(
            tests.join(format!("x16_ntaxid.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for format in ["15", "16"] {
        let output = run_in(
            &dir,
            &[
                "blastn",
                "-query",
                tests.join("x16_ntaxid_q.fna").to_str().unwrap(),
                "-db",
                copy.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        refused(&output, "a database name");
    }
    fs::remove_dir_all(dir).ok();
}

/// A `-out` file holds the bytes of standard output for 15 and 16, in a search that succeeds
/// (the prolog comes with the first report, `Xml2Document`) at every thread count; the frozen
/// row `sa.x16.blastn.minus` is the same search on standard output.
#[test]
fn a_file_sink_holds_the_document_of_standard_output() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let query = tests.join("fasta/outfmt0/x16n_q.fna");
    let subject = tests.join("fasta/outfmt0/x16n_s.fna");
    for format in ["15", "16"] {
        let expected = fs::read(
            tests
                .join("fixtures/outfmt0")
                .join(format!("sa.x16.blastn.minus.{format}.out")),
        )
        .unwrap();
        for threads in ["1", "4"] {
            let scratch = temp_dir("file_sink");
            let out = scratch.join("report.out");
            let output = run_in(
                &scratch,
                &[
                    "blastn",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                    "-out",
                    out.to_str().unwrap(),
                ],
                &[],
            );
            assert_eq!(output.status.code(), Some(0), "{format}");
            assert!(output.stdout.is_empty(), "{format}");
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                String::from_utf8_lossy(&expected),
                "{format} threads {threads}"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}
