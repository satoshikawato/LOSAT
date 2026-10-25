//! `blastp -outfmt 11` (the BLAST archive) against the bytes that NCBI BLAST+ 2.17.0 `blastp
//! -num_threads 1` writes, for the cases that need an environment variable, standard input or a
//! failing run (the searches without them are the `sa.arch.bp.*` rows of `outfmt0_manifest.tsv`
//! and the `db.sa.arch.bp.*` rows of `db_search_manifest.tsv`): `BATCH_SIZE` (one archive per
//! query batch, `EffectiveSearchSpace` per batch, an unsearched batch repeating the last value), a
//! reader error in the second batch (the archive of the first batch, then an error archive, rc 1),
//! `ARCHIVE_FORMAT` (`bin...` in any case: ASN.1 binary; any value except `xml` that is not empty:
//! nothing is written, rc 0), `-outfmt "11 qseqid"` (the specification is ignored), the errors
//! NCBI turns into an error archive (a `-subject` or `-query` file that does not exist, with
//! `-subject` and with `-db`; `-query_loc` and `-subject_loc` past the end; an option error inside
//! `Run`: `-ungapped` with composition-based statistics), a query that has no letters at all, a run
//! that ends with messages no archive carried (an empty pipe as the query with `-max_target_seqs
//! 1`: an error archive without an error, rc 0) and an empty query file (`Query is Empty!`,
//! nothing written, rc 0).
//!
//! Frozen bytes: `tests/fixtures/outfmt11/bp_<case>.out` and `.err`, written by NCBI's `blastp
//! -num_threads 1 <args> -outfmt 11` run from `LOSAT/` with the case's environment
//! (S/of8/scripts/freeze_outfmt11_bp_tests.sh of the line OF8 records). Every case runs at 1 and 4
//! threads and, with the query-block search forced (`LOSAT_X_BPQBLOCK=1`, a block of one batch),
//! must give the same bytes: the archive of each batch is the same document in or out of a block,
//! and the run-level state (the last `EffectiveSearchSpace` and the saved messages) is the run's.
//! A database that does not exist (rc 2) carries the search path in its message, which depends
//! on where the test runs: only the structure is compared.

use std::path::PathBuf;
use std::process::{Command, Output};

mod reach_support;

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn frozen(name: &str, ext: &str) -> Vec<u8> {
    std::fs::read(
        engine()
            .join("tests/fixtures/outfmt11")
            .join(format!("bp_{name}.{ext}")),
    )
    .unwrap_or_default()
}

fn input(name: &str) -> String {
    format!("{INPUTS}/{name}")
}

/// One run of `losat blastp` from `LOSAT/` in a clean environment plus `env`.
fn run(env: &[(&str, &str)], args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(engine()).arg("blastp").args(args);
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
    /// Files by name in `tests/fasta/outfmt11`; `subject` may be a name that does not exist.
    query: &'static str,
    subject: &'static str,
    /// Arguments after `-query` and `-subject`.
    extra: &'static [&'static str],
    outfmt: &'static str,
    rc: i32,
}

const BATCH: (&str, &str) = ("BATCH_SIZE", "150");
const BIN_FORMAT: (&str, &str) = ("ARCHIVE_FORMAT", "bin");

const fn case(
    name: &'static str,
    env: &'static [(&'static str, &'static str)],
    query: &'static str,
    extra: &'static [&'static str],
    rc: i32,
) -> Case {
    Case {
        name,
        env,
        query,
        subject: "arch_bp_s.faa",
        extra,
        outfmt: "11",
        rc,
    }
}

const CASES: &[Case] = &[
    case("batches_sub", &[BATCH], "arch_bp_batches.faa", &[], 0),
    case(
        "batches_sub_bin",
        &[BIN_FORMAT, BATCH],
        "arch_bp_batches.faa",
        &[],
        0,
    ),
    case(
        "warn_sub",
        &[("BATCH_SIZE", "60")],
        "arch_bp_warn.faa",
        &["-max_target_seqs", "2"],
        0,
    ),
    case("bad_sub", &[BATCH], "arch_bp_bad.faa", &[], 1),
    case(
        "bad_sub_bin",
        &[BIN_FORMAT, BATCH],
        "arch_bp_bad.faa",
        &[],
        1,
    ),
    case(
        "invbatch_sub",
        &[("BATCH_SIZE", "40")],
        "arch_bp_invbatch.faa",
        &[],
        0,
    ),
    case(
        "env_bin_upper",
        &[("ARCHIVE_FORMAT", "BIN")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    case(
        "env_binary",
        &[("ARCHIVE_FORMAT", "binary")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    case(
        "env_empty",
        &[("ARCHIVE_FORMAT", "")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    Case {
        outfmt: "11 qseqid",
        ..case("spec", &[], "arch_bp_q.faa", &[], 0)
    },
    Case {
        subject: "does_not_exist.fa",
        ..case("nosubject", &[], "arch_bp_q.faa", &[], 1)
    },
    case("noquery_sub", &[], "no_such_query.faa", &[], 1),
    case("empty_query", &[], "arch_bp_empty.faa", &[], 0),
    case("emptyonly", &[], "arch_bp_emptyonly.faa", &[], 3),
    case(
        "qloc_past",
        &[],
        "arch_bp_q.faa",
        &["-query_loc", "500-600"],
        3,
    ),
    case(
        "qloc_rev",
        &[],
        "arch_bp_q.faa",
        &["-query_loc", "50-10"],
        3,
    ),
    case(
        "sloc_past",
        &[],
        "arch_bp_q.faa",
        &["-subject_loc", "500-600"],
        1,
    ),
    case(
        "errcbs",
        &[],
        "arch_bp_q.faa",
        &["-ungapped", "-num_descriptions", "5"],
        1,
    ),
    case(
        "env_foo",
        &[("ARCHIVE_FORMAT", "foo")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    case(
        "env_spbin",
        &[("ARCHIVE_FORMAT", " bin")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    case(
        "env_asn",
        &[("ARCHIVE_FORMAT", "asn")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
    case(
        "env_text",
        &[("ARCHIVE_FORMAT", "text")],
        "arch_bp_q.faa",
        &[],
        0,
    ),
];

/// The cases that the query-block search can take: more than one batch, no reader `stop`.
const MULTI_BATCH: &[&str] = &["batches_sub", "batches_sub_bin", "warn_sub", "invbatch_sub"];

fn check(case: &Case, extra_env: &[(&str, &str)], threads: &str) {
    let mut env: Vec<(&str, &str)> = case.env.to_vec();
    env.extend_from_slice(extra_env);
    let query = input(case.query);
    let subject = input(case.subject);
    let mut args = vec!["-query", &query, "-subject", &subject];
    args.extend_from_slice(case.extra);
    args.extend_from_slice(&["-outfmt", case.outfmt]);
    let output = run(&env, &args, threads);
    let label = format!("{} threads {threads} {extra_env:?}", case.name);
    assert_eq!(output.status.code(), Some(case.rc), "{label}: exit status");
    let want = frozen(case.name, "out");
    assert!(
        output.stdout == want,
        "{label}: stdout differs ({} bytes, NCBI wrote {}): {}",
        output.stdout.len(),
        want.len(),
        crate::outfmt11_blastx::first_difference(&output.stdout, &want)
    );
    assert_eq!(output.stderr, frozen(case.name, "err"), "{label}: stderr");
}

#[test]
fn archives_equal_ncbis_at_one_and_four_threads() {
    for case in CASES {
        for threads in ["1", "4"] {
            check(case, &[], threads);
        }
    }
}

/// The query-block search takes outfmt 11 under the switch of outfmt 8, 9 and 12
/// (`LOSAT_X_BPQBLOCKSA`) and writes the same bytes: with a block of one query letter, every
/// batch is a block of its own (the reach rows of the run show that the block search took the
/// archive: `LOSAT_X_REACH`, `tests/reach_support`); `LOSAT_X_BPQBLOCKSA=0` keeps the search of
/// all queries.
#[test]
fn query_blocks_do_not_change_the_archives() {
    // Only these cases run more than one batch to a searched end: the others have one batch, or
    // stop or fail before a search, so the plan of blocks refuses them.
    for case in CASES.iter().filter(|case| MULTI_BATCH.contains(&case.name)) {
        for threads in ["1", "4"] {
            let reach = reach_support::Reach::new("outfmt11_bp");
            check(case, &[("LOSAT_X_BPQBLOCK", "1"), reach.var()], threads);
            for switch in ["LOSAT_X_BPQBLOCK", "LOSAT_X_BPQBLOCKSA"] {
                let sides = reach.sides(switch);
                assert!(
                    sides.0,
                    "{} threads {threads}: {switch} did not take the block search: {sides:?}",
                    case.name
                );
            }
            check(
                case,
                &[("LOSAT_X_BPQBLOCK", "1"), ("LOSAT_X_BPQBLOCKSA", "0")],
                threads,
            );
        }
    }
}

/// `-out FILE` receives the bytes of the standard output; `ARCHIVE_FORMAT` that NCBI does not
/// write still creates the (empty) file.
#[test]
fn the_out_file_receives_the_archive() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_bp_out_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = input("arch_bp_q.faa");
    let subject = input("arch_bp_s.faa");
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
    let query = input("arch_bp_q.faa");
    let subject = input("arch_bp_s.faa");
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
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_bp_title_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let plain = engine().join(input("arch_bp_s.faa"));
    let plain = plain.to_str().unwrap();
    for (bytes, name) in [
        (
            &b">t caf\xc3\xa9\nMKVLAAGIVALLLAAGCSSKVLAAGIVALLLAAGCSS\n"[..],
            "utf8",
        ),
        (
            &b">t ti\x7ftle\nMKVLAAGIVALLLAAGCSSKVLAAGIVALLLAAGCSS\n"[..],
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
    let query = input("arch_bp_q.faa");
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
            "BLAST Database error: No alias or index file found for protein database ["
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
    let query = std::fs::read(engine().join(input("arch_bp_q.faa"))).unwrap();
    let subject = input("arch_bp_s.faa");
    let mut child = Command::new(BIN)
        .current_dir(engine())
        .args([
            "blastp",
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
                &input("no_such_query.faa"),
                "-db",
                "tests/db_search/arch_sprot",
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
            .args(["blastp", "-query", "-", "-subject"])
            .arg(input("arch_bp_s.faa"))
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

/// NCBI's `PrintErrorArchive` converts the trimmed `-outfmt` text to an integer and writes
/// nothing when that throws (`11 qseqid`); `11` with white space around it converts. NCBI 2.17.0
/// with a missing `-query`: rc 1 and no bytes for `11 qseqid` and ` 11 qseqid`, the error
/// archive of `noquery_sub` for ` 11 `.
#[test]
fn an_outfmt_spec_gets_no_error_archive() {
    for threads in ["1", "4"] {
        for (outfmt, written) in [("11 qseqid", false), (" 11 qseqid", false), (" 11 ", true)] {
            let output = run(
                &[],
                &[
                    "-query",
                    &input("no_such_query.faa"),
                    "-subject",
                    &input("arch_bp_s.faa"),
                    "-outfmt",
                    outfmt,
                ],
                threads,
            );
            assert_eq!(output.status.code(), Some(1), "{outfmt:?}");
            if written {
                assert_eq!(output.stdout, frozen("noquery_sub", "out"), "{outfmt:?}");
            } else {
                assert!(
                    output.stdout.is_empty(),
                    "{outfmt:?}: {} bytes",
                    output.stdout.len()
                );
            }
            assert_eq!(output.stderr, frozen("noquery_sub", "err"), "{outfmt:?}");
        }
    }
}
