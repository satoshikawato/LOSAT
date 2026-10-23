//! `blastx -outfmt 11` (the BLAST archive) against the bytes that NCBI BLAST+ 2.17.0 `blastx
//! -num_threads 1` writes, for the cases that need an environment variable or a failing run (the
//! searches without them are the `sa.arch.bx.*` rows of `outfmt0_manifest.tsv` and the
//! `db.sa.arch.bx.*` rows of `db_search_manifest.tsv`): `BATCH_SIZE` (one archive per query
//! batch, `EffectiveSearchSpace` per batch, an unsearched batch repeating the last value, a batch
//! above the default size not split, a batch whose last query crosses the limit), a reader error
//! in the second batch (the archive of the first batch, then an error archive, rc 1),
//! `ARCHIVE_FORMAT` (`bin...` in any case: ASN.1 binary; any other value that is not empty: nothing
//! is written, rc 0), `-outfmt "11 qseqid"` (the specification is ignored), the errors NCBI turns
//! into an error archive (a `-subject` or `-query` file that does not exist, with `-subject` and
//! with `-db`; an option error inside `Run`: `-ungapped` with composition-based statistics; a
//! subject file with no letters), a query file with no letters (rc 3), an empty query file
//! (`Query is Empty!`, nothing written, rc 0), and `BATCH_SIZE` with the other formats (0, 5, 6, 7
//! and 8: the same batches as NCBI's, which `BATCH_SIZE` did not reach before outfmt 11).
//!
//! Frozen bytes: `tests/fixtures/outfmt11/bx_<case>.out` and `.err`, written by NCBI's `blastx
//! -num_threads 1 <args> -outfmt N` run from `LOSAT/` with the case's environment
//! (S/of8/scripts/freeze_outfmt11_bx_tests.sh of the line OF8 records). Every case runs at 1 and 4
//! threads, and the multi-batch cases also with the batches searched in waves instead of the flow
//! (`LOSAT_X_BXFLOW=0`). A database that does not exist (rc 2) carries the search path in its
//! message, which depends on where the test runs: only the structure is compared.

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
            .join(format!("bx_{name}.{ext}")),
    )
    .unwrap_or_default()
}

fn input(name: &str) -> String {
    format!("{INPUTS}/{name}")
}

/// One run of `losat blastx` from `LOSAT/` in a clean environment plus `env`.
fn run(env: &[(&str, &str)], args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(engine()).arg("blastx").args(args);
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
    /// `-db tests/db_search/arch_sprot` in place of `-subject`.
    db: bool,
}

const BATCH: (&str, &str) = ("BATCH_SIZE", "200");
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
        subject: "arch_bx_s.faa",
        extra,
        outfmt: "11".to_string(),
        rc,
        db: false,
    }
}

fn cases() -> Vec<Case> {
    let mut cases = vec![
        case("batches_sub", &[BATCH], "arch_bx_batches.fna", &[], 0),
        case(
            "batches_sub_bin",
            &[BIN_FORMAT, BATCH],
            "arch_bx_batches.fna",
            &[],
            0,
        ),
        case(
            "batch1_sub",
            &[("BATCH_SIZE", "1")],
            "arch_bx_batches.fna",
            &[],
            0,
        ),
        case(
            "batchbig_sub",
            &[("BATCH_SIZE", "30000")],
            "arch_bx_batches.fna",
            &[],
            0,
        ),
        case(
            "batchmix_sub",
            &[("BATCH_SIZE", "400")],
            "arch_bx_mix.fna",
            &[],
            0,
        ),
        case(
            "batchmix_minus_sub",
            &[("BATCH_SIZE", "400")],
            "arch_bx_mix.fna",
            &["-strand", "minus"],
            0,
        ),
        case(
            "warn_sub",
            &[("BATCH_SIZE", "60")],
            "arch_bx_warn.fna",
            &["-max_target_seqs", "2"],
            0,
        ),
        case(
            "bad_sub",
            &[("BATCH_SIZE", "60")],
            "arch_bx_bad.fna",
            &[],
            1,
        ),
        case(
            "bad_sub_bin",
            &[BIN_FORMAT, ("BATCH_SIZE", "60")],
            "arch_bx_bad.fna",
            &[],
            1,
        ),
        case(
            "bad_first_sub",
            &[("BATCH_SIZE", "15000")],
            "arch_bx_bad.fna",
            &[],
            1,
        ),
        case(
            "invbatch_sub",
            &[("BATCH_SIZE", "40")],
            "arch_bx_inv.fna",
            &[],
            0,
        ),
        case(
            "invafter_sub",
            &[("BATCH_SIZE", "40")],
            "arch_bx_invafter.fna",
            &[],
            0,
        ),
        case(
            "invfirst_sub",
            &[("BATCH_SIZE", "1")],
            "arch_bx_invfirst.fna",
            &[],
            0,
        ),
        case(
            "invonly_sub",
            &[("BATCH_SIZE", "1")],
            "arch_bx_invonly.fna",
            &[],
            0,
        ),
        case(
            "lcbatch_sub",
            &[("BATCH_SIZE", "300")],
            "arch_bx_lc.fna",
            &["-lcase_masking"],
            0,
        ),
        case(
            "segbatch_sub",
            &[("BATCH_SIZE", "100")],
            "arch_bx_seg.fna",
            &[],
            0,
        ),
        case(
            "zero_sub",
            &[("BATCH_SIZE", "0")],
            "arch_bx_batches.fna",
            &[],
            3,
        ),
        case(
            "env_bin_upper",
            &[("ARCHIVE_FORMAT", "BIN")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_binary",
            &[("ARCHIVE_FORMAT", "binary")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_empty",
            &[("ARCHIVE_FORMAT", "")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_foo",
            &[("ARCHIVE_FORMAT", "foo")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_spbin",
            &[("ARCHIVE_FORMAT", " bin")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_asn",
            &[("ARCHIVE_FORMAT", "asn")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        case(
            "env_text",
            &[("ARCHIVE_FORMAT", "text")],
            "arch_bx_q.fna",
            &[],
            0,
        ),
        Case {
            outfmt: "11 qseqid".to_string(),
            ..case("spec", &[], "arch_bx_q.fna", &[], 0)
        },
        Case {
            subject: "does_not_exist.fa",
            ..case("nosubject", &[], "arch_bx_q.fna", &[], 1)
        },
        case("noquery_sub", &[], "no_such_query.fna", &[], 1),
        case("empty_query", &[], "arch_bx_empty.fna", &[], 0),
        case("emptyonly", &[], "arch_bx_emptyonly.fna", &[], 3),
        Case {
            subject: "arch_bx_emptyonly.fna",
            ..case("emptysubjfile", &[], "arch_bx_q.fna", &[], 0)
        },
        case(
            "errcbs",
            &[],
            "arch_bx_q.fna",
            &["-ungapped", "-num_descriptions", "5"],
            1,
        ),
        Case {
            db: true,
            ..case("noquery_db", &[], "no_such_query.fna", &[], 1)
        },
        Case {
            db: true,
            ..case("errcbs_db", &[], "arch_bx_q.fna", &["-ungapped"], 1)
        },
    ];
    // `BATCH_SIZE` with the other formats: NCBI's batches in 0, 5, 6, 7 and 8 too.
    for format in ["0", "5", "6", "7", "8"] {
        for (stem, size, query) in [
            ("batches", "200", "arch_bx_batches.fna"),
            ("batchmix", "400", "arch_bx_mix.fna"),
            ("invbatch", "40", "arch_bx_inv.fna"),
            ("invafter", "40", "arch_bx_invafter.fna"),
        ] {
            let env: &'static [(&str, &str)] = Box::leak(Box::new([("BATCH_SIZE", size)]));
            cases.push(Case {
                outfmt: format.to_string(),
                ..case(&format!("{stem}_sub_fmt{format}"), env, query, &[], 0)
            });
        }
    }
    cases.push(Case {
        outfmt: "6".to_string(),
        ..case(
            "bad_sub_fmt6",
            &[("BATCH_SIZE", "60")],
            "arch_bx_bad.fna",
            &[],
            1,
        )
    });
    cases.push(Case {
        outfmt: "6".to_string(),
        ..case(
            "zero_sub_fmt6",
            &[("BATCH_SIZE", "0")],
            "arch_bx_batches.fna",
            &[],
            3,
        )
    });
    cases.push(Case {
        outfmt: "6".to_string(),
        db: true,
        ..case("batches_db_fmt6", &[BATCH], "arch_bx_batches.fna", &[], 0)
    });
    // A reader batch that spans several search units: `arch_bx_units.fna` is 89 queries and
    // 36,608 letters; with `BATCH_SIZE=30000` the first reader batch holds three units (the
    // second is a single query of 10,100 `N`, which has no valid context) and the second batch
    // one. The search keeps units of 10,002 letters, the reports are NCBI's batches'.
    for format in ["11", "6", "5"] {
        let name = if format == "11" {
            "units_sub".to_string()
        } else {
            format!("units_sub_fmt{format}")
        };
        cases.push(Case {
            outfmt: format.to_string(),
            ..case(
                &name,
                &[("BATCH_SIZE", "30000")],
                "arch_bx_units.fna",
                &[],
                0,
            )
        });
    }
    cases
}

/// The cases with more than one batch to a searched end.
const MULTI_BATCH: &[&str] = &[
    "batches_sub",
    "batches_sub_bin",
    "batch1_sub",
    "batchmix_sub",
    "batchmix_minus_sub",
    "warn_sub",
    "invbatch_sub",
    "invafter_sub",
    "invfirst_sub",
    "lcbatch_sub",
    "segbatch_sub",
    "batches_sub_fmt6",
    "batches_sub_fmt5",
    "units_sub",
    "units_sub_fmt6",
];

fn check(case: &Case, extra_env: &[(&str, &str)], threads: &str) {
    let mut env: Vec<(&str, &str)> = case.env.to_vec();
    env.extend_from_slice(extra_env);
    let query = input(case.query);
    let subject = input(case.subject);
    let mut args = vec!["-query", &query];
    if case.db {
        args.extend_from_slice(&["-db", "tests/db_search/arch_sprot"]);
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

/// The flow of batches (`LOSAT_X_BXFLOW`, in the strict set) writes the archives the waves write:
/// the multi-batch cases are compared with `LOSAT_X_BXFLOW=0`.
#[test]
fn the_flow_does_not_change_the_archives() {
    for case in cases()
        .iter()
        .filter(|case| MULTI_BATCH.contains(&case.name.as_str()))
    {
        for threads in ["1", "4"] {
            check(case, &[("LOSAT_X_BXFLOW", "0")], threads);
        }
    }
}

/// `-out FILE` receives the bytes of the standard output; `ARCHIVE_FORMAT` that NCBI does not
/// write still creates the (empty) file.
#[test]
fn the_out_file_receives_the_archive() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_out_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let query = input("arch_bx_q.fna");
    let subject = input("arch_bx_s.faa");
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
    let query = input("arch_bx_q.fna");
    let subject = input("arch_bx_s.faa");
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
    let plain = engine().join(input("arch_bx_s.faa"));
    let plain = plain.to_str().unwrap();
    for (bytes, name) in [
        (
            &b">t caf\xc3\xa9\nATGAAAGTGCTGGCGGCGGGCATTGTGGCGCTGCTG\n"[..],
            "utf8",
        ),
        (
            &b">t ti\x7ftle\nATGAAAGTGCTGGCGGCGGGCATTGTGGCGCTGCTG\n"[..],
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

/// A title the archive cannot write in a batch that NCBI completes before a reader error: the
/// archive of that batch would be written before the error, so LOSAT refuses before the first
/// byte (the prescan checks the titles batch by batch, with the run's batch size). The bad record
/// is in the second batch (`BATCH_SIZE=60`), or after the first 10,002 letters (the default).
#[test]
fn a_non_ascii_title_before_a_reader_error_is_refused_before_the_first_byte() {
    let dir =
        std::env::temp_dir().join(format!("losat_outfmt11_batchtitle_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let subject = engine().join(input("arch_bx_s.faa"));
    let subject = subject.to_str().unwrap();
    let short = dir.join("utf_bad.fna");
    std::fs::write(
        &short,
        b">bq1 caf\xc3\xa9 title\nATGTGGATGGAACACTCGGTGTTGAAAGCAGTATTTGCGGTACGCAGA\n\
          >bq2 ok\nAGTGTATTGAAGGCAGACCAACTCCAGCTAAGAAAAGTCATGGATGAT\n>bq3 bad\nAC1GT%\n",
    )
    .unwrap();
    let long = dir.join("utf_bad_default.fna");
    let mut bytes = b">bq1 caf\xc3\xa9 title\n".to_vec();
    for _ in 0..842 {
        bytes.extend_from_slice(b"ATGAAAGTGCTG");
    }
    bytes.extend_from_slice(
        b"\n>bq2 ok\nAGTGTATTGAAGGCAGACCAACTCCAGCTAAGAAAAGTCATGGATGAT\n>bq3 bad\nAC1GT%\n",
    );
    std::fs::write(&long, bytes).unwrap();
    for (file, batch_size) in [(&short, Some("60")), (&long, None)] {
        let env: Vec<(&str, &str)> = batch_size
            .map(|size| ("BATCH_SIZE", size))
            .into_iter()
            .collect();
        let output = run(
            &env,
            &[
                "-query",
                file.to_str().unwrap(),
                "-subject",
                subject,
                "-outfmt",
                "11",
            ],
            "1",
        );
        let label = format!("{file:?} {batch_size:?}");
        assert_eq!(output.status.code(), Some(1), "{label}");
        assert!(output.stdout.is_empty(), "{label}: nothing is written");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("visible ASCII"), "{label}: {stderr}");
        // The same file without outfmt 11 reaches the reader error (the archive's refusal is
        // not involved).
        let output = run(
            &env,
            &[
                "-query",
                file.to_str().unwrap(),
                "-subject",
                subject,
                "-outfmt",
                "8",
            ],
            "1",
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("CFastaReader"),
            "{label}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// A named FIFO as `-query` with outfmt 11 is refused at once (the prescan would open it a second
/// time and wait for a writer that never comes): the type is decided on the open handle.
#[cfg(unix)]
#[test]
fn a_fifo_query_is_refused_without_waiting() {
    use std::time::{Duration, Instant};
    let dir = std::env::temp_dir().join(format!("losat_outfmt11_fifo_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fifo = dir.join("query_fifo");
    let _ = std::fs::remove_file(&fifo);
    let made = Command::new("mkfifo").arg(&fifo).status().expect("mkfifo");
    assert!(made.success());
    let writer = {
        let fifo = fifo.clone();
        let bytes = std::fs::read(engine().join(input("arch_bx_q.fna"))).unwrap();
        std::thread::spawn(move || {
            use std::io::Write;
            // Opens when the run opens the FIFO, writes, and closes: later opens find no writer.
            let mut file = std::fs::OpenOptions::new().write(true).open(&fifo).unwrap();
            let _ = file.write_all(&bytes);
        })
    };
    let subject = engine().join(input("arch_bx_s.faa"));
    let mut command = Command::new(BIN);
    command
        .current_dir(engine())
        .args(["blastx", "-query", fifo.to_str().unwrap(), "-subject"])
        .arg(&subject)
        .args(["-outfmt", "11", "-num_threads", "1"])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_") || matches!(name.as_str(), "ARCHIVE_FORMAT" | "BATCH_SIZE")
        {
            command.env_remove(name);
        }
    }
    let mut child = command.spawn().expect("spawn losat");
    let started = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().unwrap() {
            break Some(status);
        }
        if started.elapsed() > Duration::from_secs(30) {
            let _ = child.kill();
            let _ = child.wait();
            break None;
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    // The writer finished when the run opened the FIFO; it is not joined when the run never did.
    if writer.is_finished() {
        let _ = writer.join();
    }
    let output = child.wait_with_output().unwrap();
    let status = status.expect("the run must refuse the FIFO, not wait on it");
    assert_eq!(status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("regular file"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// A database that does not exist: rc 2 and an error archive whose message is the line on
/// standard error (the message wraps at 78 columns, so the comparison joins the lines).
#[test]
fn a_missing_database_gives_an_error_archive() {
    let query = input("arch_bx_q.fna");
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

/// A `BATCH_SIZE` NCBI cannot convert stops it with a `CStringException`; LOSAT refuses it
/// before the first byte and says so (outfmt 11 and 6).
#[test]
fn a_batch_size_that_is_not_an_integer_is_refused() {
    let query = input("arch_bx_q.fna");
    let subject = input("arch_bx_s.faa");
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
