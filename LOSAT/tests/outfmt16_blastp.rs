#![allow(warnings, clippy::all)]

//! BLAST XML2 and JSON (outfmt 16 and 15) of BLASTP: what the frozen-fixture rows cannot show.
//! The `-subject` report that NCBI holds back one query (a failed batch loses it), the documents
//! of a run with no query, and the searches LOSAT refuses before the first byte (a byte >= 0x80 in
//! a title or the `-db` name, DECISIONS 19 Q2; a `-db` taxid with `taxdb.bti` on the search path,
//! Q4). The bytes of NCBI BLAST+ 2.17.0 are in `tests/fixtures/outfmt16`
//! (`badbatch_{sub,db}_{15,16}.out`: `BATCH_SIZE=100 blastp -query qprot_bad.faa -subject
//! sprot.faa | -db db/sprot -num_threads 1 -outfmt N`, run in that directory).
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

/// Design §5.2 row 22 (`blastp_{sub,db}_badbatch_fmt{15,16}`): with `BATCH_SIZE=100` the third
/// query's batch fails to read. `-db` has written the reports of the two queries before it;
/// `-subject` only the first one's (the second is held back and lost), and neither has an epilog.
#[test]
fn a_failed_batch_keeps_the_db_report_and_loses_the_held_subject_report() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for (mode, subject_args) in [
            ("sub", ["-subject", "sprot.faa"]),
            ("db", ["-db", "db/sprot"]),
        ] {
            let mut args = vec!["blastp", "-query", "qprot_bad.faa"];
            args.extend(subject_args);
            args.extend(["-num_threads", "1", "-outfmt", format]);
            let expected = fs::read(dir.join(format!("badbatch_{mode}_{format}.out"))).unwrap();
            for threads in ["1", "4"] {
                let mut args = args.clone();
                let at = args.iter().position(|a| *a == "-num_threads").unwrap();
                args[at + 1] = threads;
                let output = run_in(&dir, &args, &[("BATCH_SIZE", "100")]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 10"),
                    "{mode} {format}: {stderr}"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&expected),
                    "{mode} {format} threads {threads}"
                );
            }
        }
    }
}

/// A run with no query formatted writes the prolog and the epilog (`PrintProlog`, `PrintEpilog`;
/// the probes `blastp_sub_emptypipe_fmt{15,16}`); an empty query file is `Query is Empty!` and no
/// bytes (the probe `blastp_sub_emptyfile_fmt{15,16}`).
#[test]
fn an_empty_query_pipe_writes_the_prolog_and_the_epilog() {
    let subject = fixtures().join("sprot.faa");
    let expected16 = "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2\n    xmlns=\"http://www.ncbi.nlm.nih.gov\"\n    xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n    xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n></BlastXML2>\n";
    let expected15 = "{\n\"BlastOutput2\": [\n]\n}\n";
    let db = fixtures().join("db/sprot");
    for (format, expected) in [("16", expected16), ("15", expected15)] {
        for subject_args in [
            ["-subject", subject.to_str().unwrap()],
            ["-db", db.to_str().unwrap()],
        ] {
            for threads in ["1", "4"] {
                let mut args = vec!["blastp", "-query", "-"];
                args.extend(subject_args);
                args.extend(["-num_threads", threads, "-outfmt", format]);
                let output = run_with_stdin(&args, b"");
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{args:?}: {stderr}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{args:?}"
                );
            }
        }
        let dir = temp_dir("emptyfile");
        let empty = dir.join("empty.fa");
        fs::write(&empty, b"").unwrap();
        let output = run_in(
            &dir,
            &[
                "blastp",
                "-query",
                empty.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        assert_eq!(output.status.code(), Some(0));
        assert!(output.stdout.is_empty(), "Query is Empty!: no bytes");
        fs::remove_dir_all(dir).ok();
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

const PROTEIN: &[u8] = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
const OTHER_PROTEIN: &[u8] = b"WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM";

fn refused(output: &Output, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTP"
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
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    for format in ["15", "16"] {
        // a query title, first or second; UTF-8 and Latin-1
        for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
            for position in 0..2 {
                let mut records: Vec<(&[u8], &[u8])> =
                    vec![(b"q1 ok", PROTEIN), (b"q2 ok", PROTEIN)];
                records[position] = (title, PROTEIN);
                fasta(&query, &records);
                fasta(&subject, &[(b"s1 subject", PROTEIN)]);
                let output = run_in(
                    &dir,
                    &[
                        "blastp",
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
        fasta(&query, &[(b"q1 ok", PROTEIN)]);
        fasta(
            &subject,
            &[
                (b"s1 ok", OTHER_PROTEIN),
                ("s2 caf\u{e9}".as_bytes(), PROTEIN),
            ],
        );
        let output = run_in(
            &dir,
            &[
                "blastp",
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
    let copy = dir.join("residues_caf\u{e9}");
    for extension in ["phr", "pin", "psq"] {
        fs::copy(
            tests.join(format!("residues.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for format in ["15", "16"] {
        let output = run_in(
            &dir,
            &[
                "blastp",
                "-query",
                tests.join("residues_q.faa").to_str().unwrap(),
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

/// Q3 and Q4: a reported `-db` sequence with a taxid is written (`HitDescr_taxid`, the frozen rows
/// `db.sa.x16_taxid.{15,16}`), but NCBI also writes its scientific name when `taxdb.bti` is on the
/// search path (the current directory, `BLASTDB`), which LOSAT cannot: that search is refused
/// before the first byte. A taxdb with a database whose sequences have no taxid changes nothing.
#[test]
fn a_taxid_database_is_refused_only_with_a_taxdb_on_the_search_path() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let taxid_db = tests.join("x16_taxid");
    let query = tests.join("x16_taxid_q.faa");
    let plain_db = tests.join("residues");
    let plain_query = tests.join("residues_q.faa");
    let dir = temp_dir("taxdb");
    fs::write(dir.join("taxdb.bti"), b"").unwrap();
    fn as_str(args: &[String]) -> Vec<&str> {
        args.iter().map(String::as_str).collect()
    }
    for format in ["15", "16"] {
        let args = |db: &Path, query: &Path| {
            vec![
                "blastp".to_string(),
                "-query".to_string(),
                query.to_str().unwrap().to_string(),
                "-db".to_string(),
                db.to_str().unwrap().to_string(),
                "-outfmt".to_string(),
                format.to_string(),
            ]
        };
        // without a taxdb: the taxid is written
        let work = temp_dir("nottaxdb");
        let output = run_in(&work, &as_str(&args(&taxid_db, &query)), &[]);
        assert_eq!(output.status.code(), Some(0));
        let text = String::from_utf8_lossy(&output.stdout);
        assert!(text.contains("9606"), "{text}");
        assert!(!text.contains("sciname"), "{text}");
        fs::remove_dir_all(work).ok();
        // taxdb.bti in BLASTDB, or in the current directory
        for (cwd, env) in [
            (
                temp_dir("nottaxdb2"),
                vec![("BLASTDB", dir.to_str().unwrap())],
            ),
            (dir.clone(), vec![]),
        ] {
            let output = run_in(&cwd, &as_str(&args(&taxid_db, &query)), &env);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(1), "{format} {env:?}: {stderr}");
            assert!(output.stdout.is_empty(), "nothing is written");
            assert!(
                stderr.contains("with a taxonomy ID searched with a taxonomy database (taxdb.bti)"),
                "{stderr}"
            );
            assert!(
                stderr.contains("LOSAT's BLASTP in outfmt 13 to 16"),
                "{stderr}"
            );
            if cwd != dir {
                fs::remove_dir_all(cwd).ok();
            }
            // a database without taxids is not touched by the taxdb
            let output = run_in(&dir, &as_str(&args(&plain_db, &plain_query)), &[]);
            assert_eq!(
                output.status.code(),
                Some(0),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }
    fs::remove_dir_all(dir).ok();
}
