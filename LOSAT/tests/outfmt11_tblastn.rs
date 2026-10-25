//! `tblastn -outfmt 11` (the BLAST archive) against the bytes that NCBI BLAST+ 2.17.0 `tblastn
//! -num_threads 1` writes, for the cases that need an environment variable or a failing run (the
//! searches without them are the `sa.arch.tbn.*` rows of `outfmt0_manifest.tsv` and the
//! `db.sa.arch.tbn.*` rows of `db_search_manifest.tsv`): `BATCH_SIZE` (one archive per query
//! batch, `EffectiveSearchSpace` per batch, an unsearched batch repeating the last value, the
//! batch's own SEG and `LCaseMask` masks), a reader error in the second batch (the archive of the
//! first batch, then an error archive, rc 1), `ARCHIVE_FORMAT` (`bin...` in any case: ASN.1
//! binary; any other value that is not empty: nothing is written, rc 0), `-outfmt "11 qseqid"`
//! (the specification is ignored), the errors NCBI turns into an error archive (a `-subject` or
//! `-query` file that does not exist, a `-subject_loc` past the end of a subject, a `-query_loc`
//! past the end of every query, a subject file with no letters), a query file with no letters
//! (rc 3), and an empty query file (`Query is Empty!`, nothing written, rc 0).
//!
//! Frozen bytes: `tests/fixtures/outfmt11/tbn_<case>.out` and `.err`, written by NCBI's `tblastn
//! -num_threads 1 <args> -outfmt 11` run from `LOSAT/` with the case's environment
//! (S/of8/scripts/freeze_outfmt11_tbn_tests.sh of the line OF8 records). Every case runs at 1 and
//! 4 threads, and the multi-batch cases also with each strict-set switch of the search turned off
//! (`LOSAT_X_TBN...=0`). A database that does not exist (rc 2) carries the search path in its
//! message, which depends on where the test runs: only the structure is compared.

use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";
const DB: &str = "tests/db_search/arch_tbn_snucl";

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn frozen(name: &str, ext: &str) -> Vec<u8> {
    std::fs::read(
        engine()
            .join("tests/fixtures/outfmt11")
            .join(format!("tbn_{name}.{ext}")),
    )
    .unwrap_or_default()
}

fn input(name: &str) -> String {
    format!("{INPUTS}/{name}")
}

/// One run of `losat tblastn` from `LOSAT/` in a clean environment plus `env`.
fn run(env: &[(&str, &str)], args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(engine()).arg("tblastn").args(args);
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
    /// `-db tests/db_search/arch_tbn_snucl` in place of `-subject`.
    db: bool,
}

const BATCH: (&str, &str) = ("BATCH_SIZE", "150");
const BIN_FORMAT: (&str, &str) = ("ARCHIVE_FORMAT", "bin");

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
        subject: "arch_tbn_s.fna",
        extra,
        outfmt: "11".to_string(),
        rc,
        db: false,
    }
}

fn cases() -> Vec<Case> {
    vec![
        case("batches_sub", &[BATCH], "arch_tbn_batches.faa", &[], 0),
        case(
            "batches_sub_bin",
            &[BIN_FORMAT, BATCH],
            "arch_tbn_batches.faa",
            &[],
            0,
        ),
        case(
            "batch1_sub",
            &[("BATCH_SIZE", "1")],
            "arch_tbn_batches.faa",
            &[],
            0,
        ),
        case(
            "warn_sub",
            &[("BATCH_SIZE", "60")],
            "arch_tbn_warn.faa",
            &["-max_target_seqs", "2"],
            0,
        ),
        case("bad_sub", &[BATCH], "arch_tbn_bad.faa", &[], 1),
        case(
            "bad_sub_bin",
            &[BIN_FORMAT, BATCH],
            "arch_tbn_bad.faa",
            &[],
            1,
        ),
        case(
            "bad_first_sub",
            &[("BATCH_SIZE", "15000")],
            "arch_tbn_bad.faa",
            &[],
            1,
        ),
        case(
            "invbatch_sub",
            &[("BATCH_SIZE", "40")],
            "arch_tbn_invbatch.faa",
            &[],
            0,
        ),
        case(
            "invfirst_sub",
            &[("BATCH_SIZE", "1")],
            "arch_tbn_invalid.faa",
            &[],
            0,
        ),
        case(
            "invonly_sub",
            &[("BATCH_SIZE", "40")],
            "arch_tbn_invonly.faa",
            &[],
            0,
        ),
        case(
            "lcbatch_sub",
            &[BATCH],
            "arch_tbn_lc.faa",
            &["-lcase_masking"],
            0,
        ),
        case(
            "segbatch_sub",
            &[("BATCH_SIZE", "100")],
            "arch_tbn_lowc.faa",
            &[],
            0,
        ),
        case(
            "qlocbatch_sub",
            &[BATCH],
            "arch_tbn_batches.faa",
            &["-query_loc", "100-120"],
            0,
        ),
        case(
            "zero_sub",
            &[("BATCH_SIZE", "0")],
            "arch_tbn_batches.faa",
            &[],
            3,
        ),
        case(
            "env_bin_upper",
            &[("ARCHIVE_FORMAT", "BIN")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_binary",
            &[("ARCHIVE_FORMAT", "binary")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_empty",
            &[("ARCHIVE_FORMAT", "")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_foo",
            &[("ARCHIVE_FORMAT", "foo")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_spbin",
            &[("ARCHIVE_FORMAT", " bin")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_asn",
            &[("ARCHIVE_FORMAT", "asn")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        case(
            "env_text",
            &[("ARCHIVE_FORMAT", "text")],
            "arch_tbn_q.faa",
            &[],
            0,
        ),
        Case {
            outfmt: "11 qseqid".to_string(),
            ..case("spec", &[], "arch_tbn_q.faa", &[], 0)
        },
        Case {
            subject: "does_not_exist.fa",
            ..case("nosubject", &[], "arch_tbn_q.faa", &[], 1)
        },
        case("noquery_sub", &[], "no_such_query.faa", &[], 1),
        case("empty_query", &[], "arch_tbn_empty.faa", &[], 0),
        case("emptyonly", &[], "arch_tbn_emptyonly.faa", &[], 3),
        Case {
            subject: "arch_tbn_s_emptyonly.fna",
            ..case("emptysubjfile", &[], "arch_tbn_q.faa", &[], 0)
        },
        case(
            "qlocbad",
            &[],
            "arch_tbn_q.faa",
            &["-query_loc", "500-600"],
            3,
        ),
        case(
            "slocbad",
            &[],
            "arch_tbn_q.faa",
            &["-subject_loc", "5000-6000"],
            1,
        ),
        Case {
            db: true,
            ..case("noquery_db", &[], "no_such_query.faa", &[], 1)
        },
    ]
}

/// The runs of the strict-set test: (case, switch, side the default run reaches). `=0` is the
/// reference side of each switch, so a pair is listed only when the default run takes the new
/// side (4 threads, `LOSAT_X_REACH`, OF8b2 review N1); the `=0` run then takes a path other than
/// the default one.
///
/// - `LOSAT_X_TBNSTORE`: new in every `-subject` run (the plans pack per plan); the archive reads
///   the records on the `=0` side whichever the switch says, so one case is enough.
/// - `LOSAT_X_TBNBATCHRULE`: new only where a batch is searched on its own (`batch1_sub`,
///   `warn_sub`).
/// - `LOSAT_X_TBNOWNPOOLS`, `LOSAT_X_TBNFLOWHELP`: new in the multi-batch cases below (ref only
///   in `batch1_sub` and `warn_sub`). `batches_sub_bin` is `batches_sub`'s search and is left out.
/// - `LOSAT_X_TBNKEEPFEW` is not here: with `-subject` and these inputs the default run reaches
///   its ref side only (or not at all), so a `=0` run would repeat the default path, and no
///   fixture reaches the new side (it needs a database of about 30 Mb). The review checked it by
///   hand on a 30.0 Mb database (675 records, `BATCH_SIZE=300`, 8 proteins, 4 threads): the 7
///   archives equal NCBI's with and without `LOSAT_X_TBNKEEPFEW=0`.
const STRICT_SET_RUNS: &[(&str, &str)] = &[
    ("batches_sub", "LOSAT_X_TBNSTORE"),
    ("batch1_sub", "LOSAT_X_TBNBATCHRULE"),
    ("warn_sub", "LOSAT_X_TBNBATCHRULE"),
    ("batches_sub", "LOSAT_X_TBNOWNPOOLS"),
    ("invbatch_sub", "LOSAT_X_TBNOWNPOOLS"),
    ("lcbatch_sub", "LOSAT_X_TBNOWNPOOLS"),
    ("qlocbatch_sub", "LOSAT_X_TBNOWNPOOLS"),
    ("segbatch_sub", "LOSAT_X_TBNOWNPOOLS"),
    ("batches_sub", "LOSAT_X_TBNFLOWHELP"),
    ("invbatch_sub", "LOSAT_X_TBNFLOWHELP"),
    ("lcbatch_sub", "LOSAT_X_TBNFLOWHELP"),
    ("qlocbatch_sub", "LOSAT_X_TBNFLOWHELP"),
    ("segbatch_sub", "LOSAT_X_TBNFLOWHELP"),
];

/// Why `got` (LOSAT's bytes of `what`) is not `want` (the frozen NCBI bytes): the lengths, the first
/// differing byte offset, how many bytes differ in the common length, and about 160 bytes of both
/// sides around the first difference (`<<HERE>>` marks it), so that a run on a machine that
/// cannot be reached shows the difference in its log.
fn byte_difference(what: &str, got: &[u8], want: &[u8]) -> String {
    let at = got
        .iter()
        .zip(want)
        .position(|(a, b)| a != b)
        .unwrap_or(got.len().min(want.len()));
    let differing = got.iter().zip(want).filter(|(a, b)| a != b).count();
    let around = |bytes: &[u8]| {
        let escape = |part: &[u8]| String::from_utf8_lossy(part).escape_debug().to_string();
        let from = at.saturating_sub(64).min(bytes.len());
        let to = (at + 96).min(bytes.len());
        let mid = at.min(bytes.len());
        format!(
            "{}<<HERE>>{}",
            escape(&bytes[from..mid]),
            escape(&bytes[mid..to])
        )
    };
    format!(
        "{what} differs ({} bytes, NCBI wrote {}): first difference at byte {at}, {differing} bytes differ \
         in the common length\n  LOSAT: {}\n  NCBI:  {}",
        got.len(),
        want.len(),
        around(got),
        around(want)
    )
}

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
        "{label}: {}",
        byte_difference("stdout", &output.stdout, &want)
    );
    let want_err = frozen(&case.name, "err");
    assert!(
        String::from_utf8_lossy(&output.stderr) == String::from_utf8_lossy(&want_err),
        "{label}: {}",
        byte_difference("stderr", &output.stderr, &want_err)
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

/// The paths of the strict set (`LOSAT_X_TBNSTORE`, which a run that writes the archive does not
/// take for `-subject`: the archive writes the subjects' letters, `LOSAT_X_TBNOWNPOOLS`,
/// `LOSAT_X_TBNBATCHRULE`, `LOSAT_X_TBNFLOWHELP`) write the archives their `=0` side writes: the
/// runs of `STRICT_SET_RUNS` are compared with the switch turned off.
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
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_tbn_out_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = input("arch_tbn_q.faa");
    let subject = input("arch_tbn_s.fna");
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
    let query = input("arch_tbn_q.faa");
    let subject = input("arch_tbn_s.fna");
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
/// before the first byte, and a bad title in a batch that NCBI completes before a reader error
/// too (the archive of that batch would be written before the error).
#[test]
fn non_ascii_titles_and_database_names_are_refused() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_tbn_title_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let plain_query = engine().join(input("arch_tbn_q.faa"));
    let plain_query = plain_query.to_str().unwrap();
    let plain_subject = engine().join(input("arch_tbn_s.fna"));
    let plain_subject = plain_subject.to_str().unwrap();
    for (title, name) in [(&b"caf\xc3\xa9"[..], "utf8"), (&b"ti\x7ftle"[..], "del")] {
        let mut query = b">t ".to_vec();
        query.extend_from_slice(title);
        query.extend_from_slice(b"\nMWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCW\n");
        let mut subject = b">t ".to_vec();
        subject.extend_from_slice(title);
        subject.extend_from_slice(b"\nATGAAAGTGCTGGCGGCGGGCATTGTGGCGCTGCTGGCGAAAGTGCTGGCG\n");
        let titled_query = dir.join(format!("titled_{name}.faa"));
        let titled_subject = dir.join(format!("titled_{name}.fna"));
        std::fs::write(&titled_query, query).unwrap();
        std::fs::write(&titled_subject, subject).unwrap();
        let titled_query = titled_query.to_str().unwrap();
        let titled_subject = titled_subject.to_str().unwrap();
        for (args, role) in [
            (
                vec!["-query", titled_query, "-subject", plain_subject],
                "query",
            ),
            (
                vec!["-query", plain_query, "-subject", titled_subject],
                "subject",
            ),
        ] {
            let mut full = args.clone();
            full.extend_from_slice(&["-outfmt", "11"]);
            let output = run(&[], &full, "1");
            assert_ne!(output.status.code(), Some(0), "{role}");
            assert!(output.stdout.is_empty(), "{role}: nothing is written");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains("visible ASCII"), "{name} {role}: {stderr}");
            assert!(stderr.contains("-outfmt 11"), "{name} {role}: {stderr}");
            // The same title is fine in outfmt 8: only the archive refuses it.
            let mut plain_format = args.clone();
            plain_format.extend_from_slice(&["-outfmt", "8"]);
            let output = run(&[], &plain_format, "1");
            assert_eq!(output.status.code(), Some(0), "{name} {role}: outfmt 8");
        }
    }
    // A bad title in the first batch and a reader error in the second (`BATCH_SIZE=150`).
    let bad = dir.join("utf_bad.faa");
    std::fs::write(
        &bad,
        b">bq1 caf\xc3\xa9 title\nMWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHP\n\
          CNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG\n\
          >bq2 ok\nSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMA\n\
          >bq3 bad\nMKV==\n",
    )
    .unwrap();
    for format in ["11", "8"] {
        let output = run(
            &[("BATCH_SIZE", "150")],
            &[
                "-query",
                bad.to_str().unwrap(),
                "-subject",
                plain_subject,
                "-outfmt",
                format,
            ],
            "1",
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{format}: {stderr}");
        if format == "11" {
            assert!(output.stdout.is_empty(), "{format}: nothing is written");
            assert!(stderr.contains("visible ASCII"), "{format}: {stderr}");
        } else {
            // Without the archive the run reaches the reader error.
            assert!(stderr.contains("CFastaReader"), "{format}: {stderr}");
        }
    }
    for database in ["caf\u{e9}", "d\u{7f}b"] {
        let output = run(
            &[],
            &["-query", plain_query, "-db", database, "-outfmt", "11"],
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
    let query = input("arch_tbn_q.faa");
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
    let query = input("arch_tbn_q.faa");
    let subject = input("arch_tbn_s.fna");
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
