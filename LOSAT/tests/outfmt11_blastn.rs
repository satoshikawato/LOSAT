//! `blastn -outfmt 11` (the BLAST archive) against the bytes that NCBI BLAST+ 2.17.0 `blastn
//! -num_threads 1` writes, for the cases that need an environment variable or a failing run (the
//! searches without them are the `sa.arch.bn.*` rows of `outfmt0_manifest.tsv` and the
//! `db.sa.arch.bn.*` rows of `db_search_manifest.tsv`): `BATCH_SIZE` (one archive per query
//! batch, `EffectiveSearchSpace` per batch), a reader error in the second batch (the archive of
//! the first batch, then an error archive, rc 1), `ARCHIVE_FORMAT` (`bin...` in any case:
//! ASN.1 binary; any value except `xml` that is not empty: nothing is written, rc 0), `-outfmt
//! "11 qseqid"` (the specification is ignored), a `-subject` file that does not exist (an error
//! archive, rc 1), a `-query` file that does not exist (an error archive, rc 1, with `-subject` and with
//! `-db`), a run that ends with messages no archive carried (an empty pipe as the query with `-max_target_seqs 1`:
//! an error archive without an error, rc 0) and an empty query file (`Query is Empty!`, nothing written, rc 0).
//!
//! Frozen bytes: `tests/fixtures/outfmt11/<case>.out` and `.err`, written by NCBI's `blastn
//! -num_threads 1 <args> -outfmt 11` run from `LOSAT/` with the case's environment
//! (S/of8/scripts/freeze_outfmt11_tests.sh of the line OF8 records). Every case runs at 1 and 4
//! threads, and the batch cases without report streaming (`LOSAT_X_RPTSTREAM=0`, the archives
//! of every batch written after the search from the batch ranges) and with its shadow
//! (`LOSAT_X_RPTSTREAMSHADOW=1`, which compares the streamed bytes with those of the old path).
//! A database that does not exist (rc 2) carries the search path in its message, which depends
//! on where the test runs: only the structure is compared.

use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn frozen(name: &str, ext: &str) -> Vec<u8> {
    std::fs::read(
        engine()
            .join("tests/fixtures/outfmt11")
            .join(format!("{name}.{ext}")),
    )
    .unwrap_or_default()
}

fn input(name: &str) -> String {
    format!("{INPUTS}/{name}")
}

/// One run of `losat blastn` from `LOSAT/` in a clean environment plus `env`.
fn run(env: &[(&str, &str)], args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(engine()).arg("blastn").args(args);
    command.args(["-num_threads", threads]);
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE"
            )
        {
            command.env_remove(name);
        }
    }
    for (name, value) in env {
        command.env(name, value);
    }
    command.output().expect("run losat")
}

/// A case: the frozen NCBI bytes of `name`, the run's environment and arguments, the exit status.
struct Case {
    name: &'static str,
    env: &'static [(&'static str, &'static str)],
    /// Arguments after `blastn`, `-query` and `-subject` first; files by name in `tests/fasta/outfmt11`.
    query: &'static str,
    subject: &'static str,
    outfmt: &'static str,
    rc: i32,
}

const BATCH: (&str, &str) = ("BATCH_SIZE", "500");
const BIN_FORMAT: (&str, &str) = ("ARCHIVE_FORMAT", "bin");

const CASES: &[Case] = &[
    Case {
        name: "batches_sub",
        env: &[BATCH],
        query: "arch_batches.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "batches_sub_bin",
        env: &[BIN_FORMAT, BATCH],
        query: "arch_batches.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "bad_sub",
        env: &[BATCH],
        query: "arch_bad.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 1,
    },
    Case {
        name: "bad_sub_bin",
        env: &[BIN_FORMAT, BATCH],
        query: "arch_bad.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 1,
    },
    Case {
        name: "env_bin_upper",
        env: &[("ARCHIVE_FORMAT", "BIN")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_binary",
        env: &[("ARCHIVE_FORMAT", "binary")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_empty",
        env: &[("ARCHIVE_FORMAT", "")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "spec",
        env: &[],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11 qseqid",
        rc: 0,
    },
    Case {
        name: "nosubject",
        env: &[],
        query: "arch_q.fna",
        subject: "does_not_exist.fa",
        outfmt: "11",
        rc: 1,
    },
    Case {
        name: "noquery_sub",
        env: &[],
        query: "no_such_query.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 1,
    },
    Case {
        name: "empty_query",
        env: &[],
        query: "arch_empty.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_foo",
        env: &[("ARCHIVE_FORMAT", "foo")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_spbin",
        env: &[("ARCHIVE_FORMAT", " bin")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_asn",
        env: &[("ARCHIVE_FORMAT", "asn")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
    Case {
        name: "env_text",
        env: &[("ARCHIVE_FORMAT", "text")],
        query: "arch_q.fna",
        subject: "arch_s.fna",
        outfmt: "11",
        rc: 0,
    },
];

/// The standard error without the lines of the report-streaming shadow.
fn without_shadow(stderr: &[u8]) -> Vec<u8> {
    stderr
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| !line.starts_with(b"[X_SHADOW]"))
        .flatten()
        .copied()
        .collect()
}

fn check(case: &Case, extra_env: &[(&str, &str)], threads: &str) {
    let mut env: Vec<(&str, &str)> = case.env.to_vec();
    env.extend_from_slice(extra_env);
    let query = input(case.query);
    let subject = input(case.subject);
    let output = run(
        &env,
        &[
            "-query",
            &query,
            "-subject",
            &subject,
            "-outfmt",
            case.outfmt,
        ],
        threads,
    );
    let label = format!("{} threads {threads} {extra_env:?}", case.name);
    assert_eq!(output.status.code(), Some(case.rc), "{label}: exit status");
    let want = frozen(case.name, "out");
    assert!(
        output.stdout == want,
        "{label}: stdout differs ({} bytes, NCBI wrote {})",
        output.stdout.len(),
        want.len()
    );
    assert_eq!(
        without_shadow(&output.stderr),
        frozen(case.name, "err"),
        "{label}: stderr"
    );
}

#[test]
fn archives_equal_ncbis_at_one_and_four_threads() {
    for case in CASES {
        for threads in ["1", "4"] {
            check(case, &[], threads);
        }
    }
}

#[test]
fn batch_cases_do_not_depend_on_report_streaming() {
    for case in CASES.iter().filter(|case| case.env.contains(&BATCH)) {
        for extra in [
            &[("LOSAT_X_RPTSTREAM", "0")][..],
            &[("LOSAT_X_RPTSTREAM", "1")][..],
            &[("LOSAT_X_RPTSTREAMSHADOW", "1")][..],
        ] {
            for threads in ["1", "4"] {
                check(case, extra, threads);
            }
        }
    }
}

/// `-out FILE` receives the bytes of the standard output; `ARCHIVE_FORMAT` that NCBI does not
/// write still creates the (empty) file.
#[test]
fn the_out_file_receives_the_archive() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_out_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = input("arch_q.fna");
    let subject = input("arch_s.fna");
    for (format, want) in [
        ("", frozen("env_empty", "out")),
        ("foo", Vec::new()),
        ("bin", frozen("env_bin_upper", "out")),
    ] {
        let file = dir.join(format!("out_{format}.asn"));
        let _ = std::fs::remove_file(&file);
        let output = run(
            &[("ARCHIVE_FORMAT", format)],
            &[
                "-query",
                &query,
                "-subject",
                &subject,
                "-outfmt",
                "11",
                "-out",
                file.to_str().unwrap(),
            ],
            "1",
        );
        assert_eq!(output.status.code(), Some(0), "{format}");
        assert!(
            output.stdout.is_empty(),
            "{format}: nothing on standard output"
        );
        assert_eq!(
            std::fs::read(&file).unwrap(),
            want,
            "-out file with ARCHIVE_FORMAT={format:?}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `ARCHIVE_FORMAT=xml` is NCBI's XML archive, which LOSAT does not write yet: it stops
/// before the first byte and says so.
#[test]
fn xml_archives_are_refused_before_the_first_byte() {
    let query = input("arch_q.fna");
    let subject = input("arch_s.fna");
    for value in ["xml", "XML", "Xml"] {
        let output = run(
            &[("ARCHIVE_FORMAT", value)],
            &["-query", &query, "-subject", &subject, "-outfmt", "11"],
            "1",
        );
        assert_ne!(output.status.code(), Some(0), "{value}");
        assert!(output.stdout.is_empty(), "{value}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("ARCHIVE_FORMAT=xml with -outfmt 11 is not supported by LOSAT"),
            "{value}: {stderr}"
        );
    }
}

/// A title with a byte outside `0x20..=0x7E` (0x80 and above, and 0x7F) is written as `#` by
/// NCBI, with a diagnostic per byte: LOSAT refuses it (query, `-subject` and a `-db` name)
/// before the first byte.
#[test]
fn non_ascii_titles_and_database_names_are_refused() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_title_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let plain = engine().join(input("arch_s.fna"));
    let plain = plain.to_str().unwrap();
    for (bytes, name) in [
        (
            &b">t caf\xc3\xa9\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"[..],
            "utf8",
        ),
        (
            &b">t ti\x7ftle\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"[..],
            "del",
        ),
    ] {
        let titled = dir.join(format!("titled_{name}.fna"));
        std::fs::write(&titled, bytes).unwrap();
        let titled = titled.to_str().unwrap();
        for (args, role) in [
            (vec!["-query", titled, "-subject", plain], "query"),
            (vec!["-query", plain, "-subject", titled], "subject"),
        ] {
            let mut full = args.clone();
            full.extend_from_slice(&["-outfmt", "11"]);
            let output = run(&[], &full, "1");
            assert_ne!(output.status.code(), Some(0), "{role}");
            assert!(output.stdout.is_empty(), "{role}: nothing is written");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains("visible ASCII"), "{name} {role}: {stderr}");
            assert!(stderr.contains("-outfmt 11"), "{name} {role}: {stderr}");
        }
        // The same title is fine in outfmt 8: only the archive refuses it.
        let output = run(
            &[],
            &["-query", titled, "-subject", plain, "-outfmt", "8"],
            "1",
        );
        assert_eq!(output.status.code(), Some(0), "{name}");
    }
    for database in ["caf\u{e9}", "d\u{7f}b"] {
        let output = run(
            &[],
            &["-query", plain, "-db", database, "-outfmt", "11"],
            "1",
        );
        assert_ne!(output.status.code(), Some(0), "{database:?}");
        assert!(output.stdout.is_empty(), "{database:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("-db name"),
            "{database:?}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// A database that does not exist: rc 2 and an error archive whose message is the line on
/// standard error (the message wraps at 78 columns, so the comparison joins the lines).
#[test]
fn a_missing_database_gives_an_error_archive() {
    let query = input("arch_q.fna");
    let output = run(
        &[],
        &[
            "-query",
            &query,
            "-db",
            "tests/fasta/outfmt11/no_such_db",
            "-outfmt",
            "11",
        ],
        "1",
    );
    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).unwrap();
    let line = stderr.trim_end();
    assert!(
        line.starts_with(
            "BLAST Database error: No alias or index file found for nucleotide database ["
        ),
        "{line}"
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let head = "Blast4-archive ::= {\n  request {\n    body get-request-info {\n      request-id \"Error\"\n    }\n  },\n  results {\n  },\n  messages {\n    {\n      code 2,\n      message \"";
    assert!(text.starts_with(head), "{text}");
    let rest = &text[head.len()..];
    let end = rest
        .find("\"\n    }\n  }\n}\n")
        .expect("end of the message");
    assert_eq!(rest[..end].replace('\n', ""), line);
    assert_eq!(&rest[end..], "\"\n    }\n  }\n}\n");
    // Every line of the message but the first starts with a space (the wrap goes before the blank).
    assert!(rest[..end]
        .lines()
        .skip(1)
        .all(|line| line.starts_with(' ')));
}

/// The archive written to a pipe equals the file (`Write` of a pipe in pieces).
#[test]
fn an_archive_read_from_standard_input_equals_the_file_run() {
    use std::io::Write;
    use std::process::Stdio;
    let query = std::fs::read(engine().join(input("arch_q.fna"))).unwrap();
    let subject = input("arch_s.fna");
    let mut child = Command::new(BIN)
        .current_dir(engine())
        .args([
            "blastn",
            "-query",
            "-",
            "-subject",
            &subject,
            "-outfmt",
            "11",
            "-num_threads",
            "1",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&query).unwrap();
    let output = child.wait_with_output().unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert!(output.stdout == frozen("env_empty", "out"));
}

/// A `-query` file that does not exist, with a database: the error archive of NCBI, rc 1.
#[test]
fn a_missing_query_file_gives_an_error_archive_with_a_database() {
    for threads in ["1", "4"] {
        let output = run(
            &[],
            &[
                "-query",
                &input("no_such_query.fna"),
                "-db",
                "tests/db_search/arch_snucl",
                "-outfmt",
                "11",
            ],
            threads,
        );
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(output.stdout, frozen("noquery_db", "out"));
        assert_eq!(output.stderr, frozen("noquery_db", "err"));
    }
}

/// A query that is an empty pipe is not NCBI's `Query is Empty!` (it cannot tell): the search
/// ends with the message of `-max_target_seqs 1` and no archive, so an error archive without an
/// error follows (`PrintErrorArchive` when messages are left, rc 0).
#[test]
fn messages_left_at_the_end_make_an_error_archive() {
    use std::io::Write;
    use std::process::Stdio;
    for threads in ["1", "4"] {
        let mut command = Command::new(BIN);
        command
            .current_dir(engine())
            .args(["blastn", "-query", "-", "-subject"])
            .arg(input("arch_s.fna"))
            .args([
                "-outfmt",
                "11",
                "-max_target_seqs",
                "1",
                "-num_threads",
                threads,
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        for (name, _) in std::env::vars() {
            if name.starts_with("LOSAT_X_")
                || matches!(
                    name.as_str(),
                    "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE"
                )
            {
                command.env_remove(name);
            }
        }
        let mut child = command.spawn().expect("run losat");
        child.stdin.take().expect("stdin").flush().unwrap();
        let output = child.wait_with_output().expect("wait for losat");
        assert_eq!(output.status.code(), Some(0), "threads {threads}");
        assert_eq!(
            output.stdout,
            frozen("empty_pipe_mts1", "out"),
            "threads {threads}"
        );
        assert_eq!(
            output.stderr,
            frozen("empty_pipe_mts1", "err"),
            "threads {threads}"
        );
    }
}
