//! `tblastx -outfmt 11` (the BLAST archive) against the bytes that NCBI BLAST+ 2.17.0 `tblastx
//! -num_threads 1` writes, for the cases that need an environment variable or a failing run (the
//! searches without them are the `sa.arch.tbx.*` rows of `outfmt0_manifest.tsv` and the
//! `db.sa.arch.tbx.*` rows of `db_search_manifest.tsv`): `BATCH_SIZE` (one archive per query
//! batch, `EffectiveSearchSpace` per batch, an unsearched batch repeating the last value, a
//! batch whose first query has no valid context 0, queries of different lengths in a batch, the
//! batch's own SEG masks), a reader error in the second batch (the archive of the first batch,
//! then an error archive, rc 1), `ARCHIVE_FORMAT` (`bin...` in any case: ASN.1 binary; `xml`:
//! NCBI's XML archive; any other value that is not empty: nothing is written, rc 0), `-outfmt "11 qseqid"` (the specification
//! is ignored), the errors NCBI turns into an error archive (a `-subject` or `-query` file that
//! does not exist, a `-subject_loc` past the end of a subject, a `-query_loc` past the end of
//! every query, a subject file with no letters), a query file with no letters (rc 3), an empty
//! query file (`Query is Empty!`, nothing written, rc 0), and the places where NCBI writes no
//! Bioseq, so that a `-subject` title outside `0x20..=0x7E` is no reason to refuse (the
//! `u8_*` cases).
//!
//! Frozen bytes: `tests/fixtures/outfmt11/tbx_<case>.out` and `.err`, written by NCBI's `tblastx
//! -num_threads 1 <args> -outfmt 11` run from `LOSAT/` with the case's environment
//! (S/of8/scripts/freeze_outfmt11_tbx_tests.sh and ..._tests2.sh of the line OF8 records). Every
//! case runs at 1 and 4 threads, and the multi-batch case also with each strict-set switch whose
//! new side the default run reaches turned off (`LOSAT_X_TBX...=0`). A database that does not
//! exist (rc 2) carries the search path in its message, which depends on where the test runs:
//! only the structure is compared.

use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";
const DB: &str = "tests/db_search/arch_tbx_snucl";

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn frozen(name: &str, ext: &str) -> Vec<u8> {
    std::fs::read(
        engine()
            .join("tests/fixtures/outfmt11")
            .join(format!("tbx_{name}.{ext}")),
    )
    .unwrap_or_default()
}

fn input(name: &str) -> String {
    format!("{INPUTS}/{name}")
}

/// One run of `losat tblastx` from `LOSAT/` in a clean environment plus `env`.
fn run(env: &[(&str, &str)], args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(engine()).arg("tblastx").args(args);
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
    name: String,
    env: &'static [(&'static str, &'static str)],
    /// Files by name in `tests/fasta/outfmt11`; `subject` may be a name that does not exist.
    query: &'static str,
    subject: &'static str,
    /// Arguments after `-query` and the subject (or `-db`).
    extra: &'static [&'static str],
    outfmt: String,
    rc: i32,
    /// `-db tests/db_search/arch_tbx_snucl` in place of `-subject`.
    db: bool,
}

const BATCH: (&str, &str) = ("BATCH_SIZE", "200");
const BIN_FORMAT: (&str, &str) = ("ARCHIVE_FORMAT", "bin");
const XML_FORMAT: (&str, &str) = ("ARCHIVE_FORMAT", "xml");

fn case(
    name: &str,
    env: &'static [(&'static str, &'static str)],
    query: &'static str,
    extra: &'static [&'static str],
    rc: i32,
) -> Case {
    Case {
        name: name.to_string(),
        env,
        query,
        subject: "arch_tbx_s.fna",
        extra,
        outfmt: "11".to_string(),
        rc,
        db: false,
    }
}

fn cases() -> Vec<Case> {
    vec![
        case("default_sub_xml", &[XML_FORMAT], "arch_tbx_q.fna", &[], 0),
        case("batches_sub", &[BATCH], "arch_tbx_batches.fna", &[], 0),
        case(
            "batches_sub_bin",
            &[BIN_FORMAT, BATCH],
            "arch_tbx_batches.fna",
            &[],
            0,
        ),
        case(
            "batch1_sub",
            &[("BATCH_SIZE", "1")],
            "arch_tbx_batches.fna",
            &[],
            0,
        ),
        case(
            "warn_sub",
            &[("BATCH_SIZE", "60")],
            "arch_tbx_warn.fna",
            &["-max_target_seqs", "2"],
            0,
        ),
        case("bad_sub", &[BATCH], "arch_tbx_bad.fna", &[], 1),
        case(
            "bad_sub_bin",
            &[BIN_FORMAT, BATCH],
            "arch_tbx_bad.fna",
            &[],
            1,
        ),
        case(
            "bad_first_sub",
            &[("BATCH_SIZE", "15000")],
            "arch_tbx_bad.fna",
            &[],
            1,
        ),
        case(
            "invbatch_sub",
            &[("BATCH_SIZE", "40")],
            "arch_tbx_invalid.fna",
            &[],
            0,
        ),
        case(
            "invfirst_sub",
            &[("BATCH_SIZE", "1")],
            "arch_tbx_invalid.fna",
            &[],
            0,
        ),
        case(
            "invonly_sub",
            &[("BATCH_SIZE", "40")],
            "arch_tbx_invonly.fna",
            &[],
            0,
        ),
        case(
            "stop1batch_sub",
            &[("BATCH_SIZE", "130")],
            "arch_tbx_stop1.fna",
            &[],
            0,
        ),
        case(
            "lenbatch_sub",
            &[("BATCH_SIZE", "300")],
            "arch_tbx_len.fna",
            &[],
            0,
        ),
        case(
            "segbatch_sub",
            &[("BATCH_SIZE", "100")],
            "arch_tbx_seg.fna",
            &[],
            0,
        ),
        case(
            "qlocbatch_sub",
            &[BATCH],
            "arch_tbx_batches.fna",
            &["-query_loc", "100-120"],
            0,
        ),
        case(
            "zero_sub",
            &[("BATCH_SIZE", "0")],
            "arch_tbx_batches.fna",
            &[],
            3,
        ),
        case(
            "env_bin_upper",
            &[("ARCHIVE_FORMAT", "BIN")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_binary",
            &[("ARCHIVE_FORMAT", "binary")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_empty",
            &[("ARCHIVE_FORMAT", "")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_foo",
            &[("ARCHIVE_FORMAT", "foo")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_spbin",
            &[("ARCHIVE_FORMAT", " bin")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_asn",
            &[("ARCHIVE_FORMAT", "asn")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        case(
            "env_text",
            &[("ARCHIVE_FORMAT", "text")],
            "arch_tbx_q.fna",
            &[],
            0,
        ),
        Case {
            outfmt: "11 qseqid".to_string(),
            ..case("spec", &[], "arch_tbx_q.fna", &[], 0)
        },
        Case {
            subject: "does_not_exist.fa",
            ..case("nosubject", &[], "arch_tbx_q.fna", &[], 1)
        },
        case("noquery_sub", &[], "no_such_query.fna", &[], 1),
        case("empty_query", &[], "arch_tbx_empty.fna", &[], 0),
        case("emptyonly", &[], "arch_tbx_emptyonly.fna", &[], 3),
        Case {
            subject: "arch_tbx_s_emptyonly.fna",
            ..case("emptysubjfile", &[], "arch_tbx_q.fna", &[], 3)
        },
        case(
            "qlocbad",
            &[],
            "arch_tbx_q.fna",
            &["-query_loc", "500-600"],
            3,
        ),
        case(
            "slocbad",
            &[],
            "arch_tbx_q.fna",
            &["-subject_loc", "5000-6000"],
            1,
        ),
        Case {
            db: true,
            ..case("noquery_db", &[], "no_such_query.fna", &[], 1)
        },
        // Q3: a `-subject` title outside `0x20..=0x7E` is no reason to refuse where NCBI writes no
        // Bioseq (an error archive, or nothing): these are NCBI's bytes.
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case("u8_empty_query", &[], "arch_tbx_empty.fna", &[], 0)
        },
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case("u8_noquery", &[], "no_such_query.fna", &[], 1)
        },
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case(
                "u8_bad_first",
                &[("BATCH_SIZE", "15000")],
                "arch_tbx_bad.fna",
                &[],
                1,
            )
        },
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case(
                "u8_zero",
                &[("BATCH_SIZE", "0")],
                "arch_tbx_batches.fna",
                &[],
                3,
            )
        },
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case(
                "u8_qlocbad",
                &[],
                "arch_tbx_q.fna",
                &["-query_loc", "500-600"],
                3,
            )
        },
        Case {
            subject: "arch_tbx_s_utf8.fna",
            ..case("u8_emptyonly", &[], "arch_tbx_emptyonly.fna", &[], 3)
        },
        // A query title outside the range in the batch that a reader error interrupts is not
        // written by NCBI (nor refused).
        case("u8_interrupted", &[BATCH], "arch_tbx_bad_u8.fna", &[], 1),
    ]
}

/// The runs of the strict-set test: (case, switch). `=0` is the reference side of each switch,
/// so a pair is listed only when the default run takes the new side (4 threads, `LOSAT_X_REACH`
/// on `batches_sub`) or the switch acts on every run (class `always`); the `=0` run then takes a
/// path other than the default one. The seed stage's switches (`PAIRPARGATE`, `SEEDCOMPACT`,
/// `TBXHITBF`, `TBXKMERGE`, `TBXSCANBF`), the pool's (`TBXNEXT`, `TBXQFRAMES`, `TBXRPTPAR`,
/// `TBXPIPE`: the next batch starts while the batch before links) and the always-on ones
/// (`SUMPMEMO`) change what the search does, not what the
/// archive reads: the first query's block, search space and masks and the final hits.
/// `DISPTABLE` is not here: it acts only where a format reads displayed counts, which outfmt 11
/// never builds, so its `=0` run is the default run. `lenbatch_sub` x `TBXPIPE` is not here for
/// the same reason (4 threads take the reference side on that case). `TBXEXTBATCH`, `TBXLUT`,
/// `TBXPREFETCH`, `TBXQSLICE`, `TBXUNITMASK` and `SEEDUNMAP` are not needed: they change only
/// the seed scan of the PAIRPAR units (value switches reach their new side with a small value on
/// a 200 kb pair), and the archive reads the final hits and the batch set-up's query statistics,
/// which are made before the seed stage. The review sweep (26 switches, `ERFABSORB`, `TBXFRAMES`
/// and `TBXUNGAP` since graduated, against 26 archive
/// cases, `S/of8/logs/RESULTS_impl_b3.md`) covers every switch.
const STRICT_SET_RUNS: &[(&str, &str)] = &[
    ("batches_sub", "LOSAT_X_PAIRPARGATE"),
    ("batches_sub", "LOSAT_X_SEEDCOMPACT"),
    ("batches_sub", "LOSAT_X_TBXHITBF"),
    ("batches_sub", "LOSAT_X_TBXKMERGE"),
    ("batches_sub", "LOSAT_X_TBXNEXT"),
    ("batches_sub", "LOSAT_X_TBXPIPE"),
    ("batches_sub", "LOSAT_X_TBXQFRAMES"),
    ("batches_sub", "LOSAT_X_TBXRPTPAR"),
    ("batches_sub", "LOSAT_X_TBXSCANBF"),
    ("segbatch_sub", "LOSAT_X_TBXQFRAMES"),
    ("stop1batch_sub", "LOSAT_X_TBXPIPE"),
    ("batches_sub", "LOSAT_X_SUMPMEMO"),
];

fn check(case: &Case, extra_env: &[(&str, &str)], threads: &str) {
    let mut env: Vec<(&str, &str)> = case.env.to_vec();
    env.extend_from_slice(extra_env);
    let query = input(case.query);
    let subject = input(case.subject);
    let mut args = vec!["-query", &query];
    if case.db {
        args.extend_from_slice(&["-db", DB]);
    } else {
        args.extend_from_slice(&["-subject", &subject]);
    }
    args.extend_from_slice(case.extra);
    args.extend_from_slice(&["-outfmt", &case.outfmt]);
    let output = run(&env, &args, threads);
    let label = format!("{} threads {threads} {extra_env:?}", case.name);
    assert_eq!(output.status.code(), Some(case.rc), "{label}: exit status");
    let want = frozen(&case.name, "out");
    assert!(
        output.stdout == want,
        "{label}: stdout differs ({} bytes, NCBI wrote {})",
        output.stdout.len(),
        want.len()
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&frozen(&case.name, "err")),
        "{label}: stderr"
    );
}

#[test]
fn archives_equal_ncbis_at_one_and_four_threads() {
    for case in &cases() {
        for threads in ["1", "4"] {
            check(case, &[], threads);
        }
    }
}

/// The paths of the strict set write the archives their `=0` side writes: the runs of
/// `STRICT_SET_RUNS` are compared with the switch turned off.
#[test]
fn the_strict_set_does_not_change_the_archives() {
    let all = cases();
    for &(name, switch) in STRICT_SET_RUNS {
        let case = all
            .iter()
            .find(|case| case.name == name)
            .unwrap_or_else(|| panic!("no case {name}"));
        check(case, &[(switch, "0")], "4");
    }
}

/// `-out FILE` receives the bytes of the standard output; `ARCHIVE_FORMAT` that NCBI does not
/// write still creates the (empty) file.
#[test]
fn the_out_file_receives_the_archive() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_tbx_out_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = input("arch_tbx_q.fna");
    let subject = input("arch_tbx_s.fna");
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

/// A title with a byte outside `0x20..=0x7E` (0x80 and above, and 0x7F) is written as `#` by
/// NCBI, with a diagnostic per byte: LOSAT refuses it (query, `-subject` and a `-db` name)
/// before the first byte, where NCBI would write that Bioseq: a query of a batch that NCBI
/// completes (a batch before a reader error too), and the `-subject` records whenever a batch
/// is completed. Where NCBI writes no Bioseq the `u8_*` cases above carry NCBI's bytes.
#[test]
fn non_ascii_titles_and_database_names_are_refused() {
    let plain_query = input("arch_tbx_q.fna");
    let plain_subject = input("arch_tbx_s.fna");
    let utf8_query = input("arch_tbx_q_utf8.fna");
    let utf8_subject = input("arch_tbx_s_utf8.fna");
    let del = std::env::temp_dir().join(format!("losat_outfmt11_tbx_title_{}", std::process::id()));
    std::fs::create_dir_all(&del).unwrap();
    // 0x7F is outside the range as well.
    let del_query = del.join("del_query.fna");
    std::fs::write(
        &del_query,
        b">t ti\x7ftle\nGACATCTAGGTTCAACACAGTGTGTTAAAGGCTGTTTTCGCT\n",
    )
    .unwrap();
    let del_query = del_query.to_str().unwrap().to_string();
    for (query, subject, role, env) in [
        (&utf8_query, &plain_subject, "query", &[][..]),
        (&plain_query, &utf8_subject, "subject", &[][..]),
        (&del_query, &plain_subject, "query (0x7F)", &[][..]),
        // The first batch is completed before the reader error of the second: its queries are
        // written (the archive) by NCBI.
        (
            &input("arch_tbx_bad_u8first.fna"),
            &plain_subject,
            "query of a batch before a reader error",
            &[BATCH][..],
        ),
        // With `-subject` every archive repeats the subjects: the batch before a reader error
        // is written.
        (
            &input("arch_tbx_bad.fna"),
            &utf8_subject,
            "subject with a batch before a reader error",
            &[BATCH][..],
        ),
    ] {
        let args = ["-query", query, "-subject", subject, "-outfmt", "11"];
        let output = run(env, &args, "1");
        assert_ne!(output.status.code(), Some(0), "{role}");
        assert!(output.stdout.is_empty(), "{role}: nothing is written");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("visible ASCII"), "{role}: {stderr}");
        assert!(stderr.contains("-outfmt 11"), "{role}: {stderr}");
        // The same title is fine in outfmt 8 (or reaches the reader error): only the archive
        // refuses it.
        let args = ["-query", query, "-subject", subject, "-outfmt", "8"];
        let output = run(env, &args, "1");
        assert!(
            !String::from_utf8_lossy(&output.stderr).contains("visible ASCII"),
            "{role}: outfmt 8"
        );
    }
    // With `-out` the refusal leaves no file that the run created (as for `tblastn`), and one
    // that was there before stays (emptied, as for every other refusal after the open).
    for (query, subject) in [(&utf8_query, &plain_subject), (&plain_query, &utf8_subject)] {
        let out = del.join("refused.asn");
        let _ = std::fs::remove_file(&out);
        let out_arg = out.to_str().unwrap().to_string();
        let args = [
            "-query", query, "-subject", subject, "-outfmt", "11", "-out", &out_arg,
        ];
        let output = run(&[], &args, "1");
        assert_ne!(output.status.code(), Some(0));
        assert!(!out.exists(), "an empty -out file is left behind");
        std::fs::write(&out, b"old").unwrap();
        let output = run(&[], &args, "1");
        assert_ne!(output.status.code(), Some(0));
        assert!(out.exists(), "a file that was there is kept");
    }
    for database in ["caf\u{e9}", "d\u{7f}b"] {
        let output = run(
            &[],
            &["-query", &plain_query, "-db", database, "-outfmt", "11"],
            "1",
        );
        assert_ne!(output.status.code(), Some(0), "{database:?}");
        assert!(output.stdout.is_empty(), "{database:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("-db name"),
            "{database:?}"
        );
    }
    let _ = std::fs::remove_dir_all(&del);
}

/// A database that does not exist: rc 2 and an error archive whose message is the line on
/// standard error (the message wraps at 78 columns, so the comparison joins the lines).
#[test]
fn a_missing_database_gives_an_error_archive() {
    let query = input("arch_tbx_q.fna");
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

/// A `BATCH_SIZE` NCBI cannot convert stops it with a `CStringException`; LOSAT refuses it
/// before the first byte and says so (outfmt 11 and 6).
#[test]
fn a_batch_size_that_is_not_an_integer_is_refused() {
    let query = input("arch_tbx_q.fna");
    let subject = input("arch_tbx_s.fna");
    for value in ["abc", "", "1.5", "0x10"] {
        for format in ["11", "6"] {
            let output = run(
                &[("BATCH_SIZE", value)],
                &["-query", &query, "-subject", &subject, "-outfmt", format],
                "1",
            );
            assert_ne!(output.status.code(), Some(0), "{value:?} {format}");
            assert!(output.stdout.is_empty(), "{value:?} {format}");
            assert!(
                String::from_utf8_lossy(&output.stderr).contains("BATCH_SIZE"),
                "{value:?} {format}"
            );
        }
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
                    &input("no_such_query.fna"),
                    "-subject",
                    &input("arch_tbx_s.fna"),
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
