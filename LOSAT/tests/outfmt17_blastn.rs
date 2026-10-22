#![allow(warnings, clippy::all)]

//! SAM (outfmt 17) of BLASTN: what the frozen-fixture rows cannot show. NCBI keeps every record until
//! its formatter is destroyed and writes the header last, so a batch that fails to read still has
//! the earlier batches' document written (the formatter is unwound before the error is printed), and
//! a run with several batches writes one document. The bytes of NCBI BLAST+ 2.17.0 are in
//! `tests/fixtures/outfmt17` (`<case>_{sub,db}_{17,17SQSR}.out`, `.err` when stderr was not empty and
//! `.rc`: `blastn -query QUERY (-subject snucl.fna | -db db/snucl) EXTRA -num_threads 1 -outfmt N`, run
//! in that directory with the environment of the case; `S/scripts/probe_of6.sh PART=F`).
//!
//! The one difference from NCBI's bytes is the CL tag of the `@PG` line (LOSAT writes its own command
//! line, DECISIONS 7): both sides are compared without it.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:186-327 (the formatter is a local of
//! the `try` block that holds the batch loop), c++/src/objtools/format/sam_formatter.cpp:313-316
//! (`~CSAM_Formatter` flushes the document).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/outfmt17")
}

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt17_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// LOSAT with no `LOSAT_*`, `NCBI*` or `BLAST*` variable of the caller, an empty home, and `env`.
fn command(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") && text != "LOSAT_X_DEFAULT" || text.starts_with("NCBI") {
            command.env_remove(name);
        }
    }
    for name in ["BATCH_SIZE", "CHUNK_SIZE", "BLASTDB"] {
        command.env_remove(name);
    }
    let home = std::env::temp_dir().join("losat_outfmt17_home");
    fs::create_dir_all(&home).expect("empty home");
    command.env("HOME", &home).env("USERPROFILE", &home);
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt17_worker.log"),
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

/// The SAM document without the CL tag of its @PG line: the bytes from the `\tCL:` that follows
/// `\n@PG\tID:0\tVN:2.17.0+` up to, not including, the first following `\tPN:blastn\n` (the rule
/// of `db_search_fixtures.py strip_sam_cl`).
fn without_cl(document: &[u8]) -> Vec<u8> {
    let find = |needle: &[u8], from: usize| {
        document[from..]
            .windows(needle.len())
            .position(|window| window == needle)
            .map(|at| at + from)
    };
    let head: &[u8] = b"\n@PG\tID:0\tVN:2.17.0+";
    let Some(at) = find(&[head, b"\tCL:"].concat(), 0) else {
        return document.to_vec();
    };
    let start = at + head.len();
    let Some(end) = find(b"\tPN:blastn\n", start + 4) else {
        return document.to_vec();
    };
    [&document[..start], &document[end..]].concat()
}

/// The standard error without the summary lines of the RPTSTREAM shadow.
fn without_shadow(stderr: &[u8]) -> Vec<u8> {
    String::from_utf8_lossy(stderr)
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_SHADOW]"))
        .collect::<String>()
        .into_bytes()
}

const SUBJECT: &str = "snucl.fna";
const DB: &str = "db/snucl";

fn mode_args(mode: &str) -> [&'static str; 2] {
    if mode == "sub" {
        ["-subject", SUBJECT]
    } else {
        ["-db", DB]
    }
}

fn frozen(name: &str, extension: &str) -> Option<Vec<u8>> {
    fs::read(fixtures().join(format!("{name}.{extension}"))).ok()
}

/// LOSAT at `threads` threads writes the stdout (CL tag aside), stderr and exit status that NCBI
/// wrote for `<stem>_<mode>_<tag>`; `query` and `extra` are the case's search, `env` its
/// environment. The environment variable switches of the streaming are tried too.
fn equals_frozen(stem: &str, query: &str, extra: &[&str], env: &[(&str, &str)], rc: i32) {
    let dir = fixtures();
    for (format, tag) in [("17", "17"), ("17 SQ SR", "17SQSR")] {
        for mode in ["sub", "db"] {
            let name = format!("{stem}_{mode}_{tag}");
            let expected =
                frozen(&name, "out").unwrap_or_else(|| panic!("{name}: no frozen bytes"));
            let expected_err = frozen(&name, "err").unwrap_or_default();
            assert_eq!(
                String::from_utf8_lossy(&frozen(&name, "rc").unwrap()).trim(),
                rc.to_string(),
                "{name}: the frozen status"
            );
            let [flag, value] = mode_args(mode);
            for threads in ["1", "4"] {
                for (switch, shadow) in [
                    (None, false),
                    (Some(("LOSAT_X_RPTSTREAM", "0")), false),
                    (Some(("LOSAT_X_RPTSTREAMSHADOW", "1")), true),
                ] {
                    let mut vars = env.to_vec();
                    vars.extend(switch);
                    let mut args = vec!["blastn", "-query", query, flag, value];
                    args.extend_from_slice(extra);
                    args.extend(["-num_threads", threads, "-outfmt", format]);
                    let output = run_in(&dir, &args, &vars);
                    let what = format!("{name} threads {threads} {switch:?}");
                    assert_eq!(output.status.code(), Some(rc), "{what}");
                    assert_eq!(
                        String::from_utf8_lossy(&without_cl(&output.stdout)),
                        String::from_utf8_lossy(&without_cl(&expected)),
                        "{what}: stdout"
                    );
                    let stderr = if shadow {
                        without_shadow(&output.stderr)
                    } else {
                        output.stderr.clone()
                    };
                    assert_eq!(
                        String::from_utf8_lossy(&stderr),
                        String::from_utf8_lossy(&expected_err),
                        "{what}: stderr"
                    );
                }
            }
        }
    }
}

const BATCHES: &[(&str, &str)] = &[("BATCH_SIZE", "100")];

/// NCBI reader error in batch 3 (`BATCH_SIZE=100`: one query per batch; `qnucl_bad.fna`): the
/// document of the two batches before it is written, then the error (exit status 1) --
/// `-subject` and `-db` alike, since SAM holds no query back (unlike outfmt 13-16).
#[test]
fn a_reader_error_in_a_later_batch_writes_the_document_of_the_batches_before() {
    equals_frozen("badbatch", "qnucl_bad.fna", &[], BATCHES, 1);
}

/// The first batch fails to read: no record exists, so NCBI writes no document at all.
#[test]
fn a_reader_error_in_the_first_batch_writes_nothing() {
    equals_frozen("firstbad", "qnucl_bad.fna", &[], &[], 1);
}

/// Three queries in three batches: one document, `@SQ` in the order of the records, `Query_N`
/// continuing across batches.
#[test]
fn several_batches_write_one_document() {
    equals_frozen("batches", "q_multi3.fna", &[], BATCHES, 0);
}

/// A batch of records whose `-query_loc` range starts past their end (here the second query, ten
/// letters long against `-query_loc 20-300`) is a batch without queries: NCBI stops with
/// `Empty CBlastQueryVector` after the document of the first batch.
#[test]
fn a_batch_without_queries_writes_the_document_of_the_batches_before() {
    equals_frozen("skip", "q_skip.fna", &["-query_loc", "20-300"], BATCHES, 3);
}

/// A batch whose only query has no letters (`Sequence contains no data`) stops NCBI after the
/// document of the first batch.
#[test]
fn a_batch_of_queries_without_letters_writes_the_document_of_the_batches_before() {
    equals_frozen("nolet", "q_nolet.fna", &[], BATCHES, 3);
}

/// `-db` whose header file is missing (`db/snucl_nophr`: the data of `snucl` without its `.nhr`): NCBI
/// opens the header file when it reads the final hit lists' definition lines, after the search, and
/// stops there with no document. The path in its message is the database's (`@DIR@` in the frozen
/// bytes is this fixture directory).
#[test]
fn a_header_file_that_cannot_be_read_writes_nothing() {
    let dir = fixtures();
    for (format, tag) in [("17", "17"), ("17 SQ SR", "17SQSR")] {
        let name = format!("nophr_db_{tag}");
        let expected = frozen(&name, "out").unwrap();
        let expected_err = String::from_utf8(frozen(&name, "err").unwrap())
            .unwrap()
            .replace("@DIR@", dir.to_str().unwrap());
        let rc: i32 = String::from_utf8_lossy(&frozen(&name, "rc").unwrap())
            .trim()
            .parse()
            .unwrap();
        for threads in ["1", "4"] {
            let output = run_in(
                &dir,
                &[
                    "blastn",
                    "-query",
                    "qnucl_valid.fna",
                    "-db",
                    "db/snucl_nophr",
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                ],
                &[],
            );
            assert_eq!(output.status.code(), Some(rc), "{name} {threads}");
            assert_eq!(
                without_cl(&output.stdout),
                without_cl(&expected),
                "{name} {threads}"
            );
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                expected_err,
                "{name} {threads}"
            );
        }
    }
}

/// An empty query file is `Query is Empty!` and no bytes (BN:209-212: the formatter does not
/// exist yet).
#[test]
fn an_empty_query_file_writes_nothing() {
    let dir = fixtures();
    let scratch = temp_dir("empty");
    let empty = scratch.join("empty.fna");
    fs::write(&empty, b"").unwrap();
    for format in ["17", "17 SQ SR"] {
        for mode in ["sub", "db"] {
            let [flag, value] = mode_args(mode);
            let output = run_in(
                &dir,
                &[
                    "blastn",
                    "-query",
                    empty.to_str().unwrap(),
                    flag,
                    value,
                    "-num_threads",
                    "1",
                    "-outfmt",
                    format,
                ],
                &[],
            );
            assert!(
                output.stdout.is_empty(),
                "{format} {mode}: {:?}",
                output.stdout
            );
            let name = format!(
                "emptyfile_sub_{}",
                if format == "17" { "17" } else { "17SQSR" }
            );
            let code: i32 = String::from_utf8_lossy(&frozen(&name, "rc").unwrap())
                .trim()
                .parse()
                .unwrap();
            assert_eq!(output.status.code(), Some(code), "{format} {mode}");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&frozen(&name, "err").unwrap_or_default()),
                "{format} {mode}"
            );
        }
    }
}

/// `-out FILE` holds the bytes of standard output, also for the document written by the error
/// path of a failed batch, and the status is the same.
#[test]
fn a_file_sink_holds_the_document_of_standard_output() {
    let dir = fixtures();
    let scratch = temp_dir("out");
    for (query, extra, env, rc) in [
        ("q_multi3.fna", &[][..], BATCHES, 0),
        ("qnucl_bad.fna", &[][..], BATCHES, 1),
    ] {
        for format in ["17", "17 SQ SR"] {
            for mode in ["sub", "db"] {
                let [flag, value] = mode_args(mode);
                let file = scratch.join("report.sam");
                let mut args = vec!["blastn", "-query", query, flag, value];
                args.extend_from_slice(extra);
                args.extend(["-num_threads", "2", "-outfmt", format]);
                let stdout = run_in(&dir, &args, env);
                args.extend(["-out", file.to_str().unwrap()]);
                let to_file = run_in(&dir, &args, env);
                assert_eq!(stdout.status.code(), Some(rc), "{query} {format} {mode}");
                assert_eq!(to_file.status.code(), Some(rc), "{query} {format} {mode}");
                assert!(to_file.stdout.is_empty(), "{query} {format} {mode}");
                assert_eq!(
                    without_cl(&fs::read(&file).unwrap()),
                    without_cl(&stdout.stdout),
                    "{query} {format} {mode}"
                );
                assert_eq!(to_file.stderr, stdout.stderr, "{query} {format} {mode}");
                fs::remove_file(&file).unwrap();
            }
        }
    }
}

/// The CL tag is LOSAT's own command line: every argument followed by one space, the program path
/// first, bytes as the system keeps them; the library has none.
#[test]
fn the_pg_line_shows_the_command_line() {
    let dir = fixtures();
    let output = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "qnucl_valid.fna",
            "-subject",
            SUBJECT,
            "-outfmt",
            "17 SQ",
        ],
        &[],
    );
    assert!(output.status.success());
    let text = String::from_utf8_lossy(&output.stdout);
    let pg = text.lines().find(|line| line.starts_with("@PG")).unwrap();
    let expected = format!(
        "@PG\tID:0\tVN:2.17.0+\tCL:{} blastn -query qnucl_valid.fna -subject {SUBJECT} -outfmt 17 SQ \tPN:blastn",
        env!("CARGO_BIN_EXE_losat")
    );
    assert_eq!(pg, expected);
}
