//! The taxonomy database of a `-subject` search's tabular reports and the `-out` file of the
//! library's `run_to_file`, against NCBI BLAST+ 2.17.0 and the command line (line RG item E1w-2b,
//! `remaining-gaps-20261112/E/SURVEY_E1w2b.md`).
//!
//! The same words go to the command line and, through `examples/library_command_line`
//! (`losat::api::run`, `run_to_file` with `-out`), to the library, in the same directory and
//! environment; a failure of the library ends the example with exit status 1 and the error's
//! message after the call's diagnostics.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104 (cited at
//! `common::database::subject_taxonomy`) and c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
//! (cited at `cli::create_out_file`).
#![cfg(unix)]

use std::ffi::{OsStr, OsString};
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use crate::library_examples::example_binary;

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// Each program's query, subject, and a subject whose reader warns (`w2 hyphens`).
const SEARCHES: [(&str, &str, &str, &str); 6] = [
    (
        "blastn",
        "arch_tbx_q.fna",
        "arch_tbx_s.fna",
        "arch_tbx_warn.fna",
    ),
    (
        "blastp",
        "arch_bp_q.faa",
        "arch_bp_s.faa",
        "arch_bp_warn.faa",
    ),
    (
        "blastx",
        "arch_bx_q.fna",
        "arch_bx_s.faa",
        "arch_bp_warn.faa",
    ),
    (
        "tblastn",
        "arch_tbn_q.faa",
        "arch_tbn_s.fna",
        "arch_tbx_warn.fna",
    ),
    (
        "tblastx",
        "arch_tbx_q.fna",
        "arch_tbx_s.fna",
        "arch_tbx_warn.fna",
    ),
    (
        "psiblast",
        "arch_bp_q.faa",
        "arch_bp_s.faa",
        "arch_bp_warn.faa",
    ),
];

/// A fresh directory with the inputs of `SEARCHES` and the taxonomy database `taxdb/`.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-out-taxdb-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("taxdb")).unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (_, query, subject, warned) in SEARCHES {
        for name in [query, subject, warned] {
            fs::copy(root.join("tests/fasta/outfmt11").join(name), dir.join(name)).unwrap();
        }
    }
    for name in ["taxdb.bti", "taxdb.btd"] {
        fs::copy(
            root.join("tests/seqdb/taxdb/good").join(name),
            dir.join("taxdb").join(name),
        )
        .unwrap();
    }
    dir
}

/// `binary words` in `dir` with an environment of NCBI's probes (`env -i PATH HOME`) and `envs`.
fn run_binary_in(
    binary: &Path,
    dir: &Path,
    words: &[OsString],
    envs: &[(&str, OsString)],
) -> Output {
    let mut command = Command::new(binary);
    command
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .envs(envs.iter().map(|(name, value)| (name, value)))
        .current_dir(dir)
        .args(words)
        .stdin(std::process::Stdio::null());
    // The `--all-features` build's BLASTX diagnostic probe needs its log (as in CI).
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        command.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    command.output().unwrap()
}

fn words(words: &[&str]) -> Vec<OsString> {
    words.iter().map(OsString::from).collect()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The command line's and the library's runs of `words`: the same standard output and error, and
/// exit status 0 where the command line's is 0 (1, the example's, for a failure).
fn both(dir: &Path, words: &[OsString], envs: &[(&str, OsString)]) -> Output {
    let cli = run_binary_in(Path::new(BIN), dir, words, envs);
    let library = run_binary_in(&example_binary("library_command_line"), dir, words, envs);
    let status = cli.status.code().map(|code| i32::from(code != 0));
    assert_eq!(
        (
            library.status.code(),
            text(&library.stdout),
            text(&library.stderr)
        ),
        (status, text(&cli.stdout), text(&cli.stderr)),
        "{words:?} {envs:?}: the library"
    );
    cli
}

/// NCBI BLAST+ 2.17.0 (`tests/fixtures/subject_taxdb`): the plural names warn once per formatter
/// unless `taxdb.bti` is on the search path (outfmt 20 once more, for its line of field names); the
/// singular names look nothing up for subjects without a taxid. A `BLASTDB` setting that is not
/// UTF-8 is a directory of those bytes, where no `taxdb.bti` is: the plural names warn as without
/// a setting (NCBI BLAST+ 2.17.0 with `NCBI_CONFIG__BLAST__BLASTDB` of the bytes `x\xffy`, the
/// six programs: the bytes of `<program>.plur.none.*`; line RG item K10a). The library's searches
/// of records look `taxdb.bti` up as the command line does (they wrote the warning whatever the
/// path held before line RG item E1w-2b).
#[test]
fn subject_taxonomy_is_ncbis_in_the_command_line_and_the_library() {
    let dir = scratch("taxdb");
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/subject_taxdb");
    let taxdb = dir.join("taxdb").into_os_string();
    let unreadable = OsStr::from_bytes(b"x\xffy").to_os_string();
    for (program, query, subject, _) in SEARCHES {
        for (format, outfmt) in [
            (
                "sing",
                "6 qseqid sseqid ssciname scomname sblastname sskingdom",
            ),
            ("plur", "6 qseqid sseqid sscinames"),
            ("csv20", "20 qseqid sscinames"),
        ] {
            for environment in ["none", "taxdb", "badenv"] {
                if (format == "sing") != (environment == "badenv") {
                    continue;
                }
                let envs: Vec<(&str, OsString)> = match environment {
                    "taxdb" => vec![("BLASTDB", taxdb.clone())],
                    "badenv" => vec![("NCBI_CONFIG__BLAST__BLASTDB", unreadable.clone())],
                    _ => Vec::new(),
                };
                let words = words(&[
                    program, "-query", query, "-subject", subject, "-outfmt", outfmt,
                ]);
                let key = format!("{program}.{format}.{environment}");
                let cli = both(&dir, &words, &envs);
                assert_eq!(cli.status.code(), Some(0), "{key}: {cli:?}");
                assert_eq!(
                    (text(&cli.stdout), text(&cli.stderr)),
                    (
                        fs::read_to_string(fixtures.join(format!("{key}.out"))).unwrap(),
                        fs::read_to_string(fixtures.join(format!("{key}.err"))).unwrap()
                    ),
                    "{key}"
                );
            }
        }
        // The plural names with a setting that is not UTF-8: NCBI's output without a setting, in
        // the command line and the library.
        let words = words(&[
            program,
            "-query",
            query,
            "-subject",
            subject,
            "-outfmt",
            "6 qseqid sseqid sscinames",
        ]);
        let envs = [("NCBI_CONFIG__BLAST__BLASTDB", unreadable.clone())];
        let cli = both(&dir, &words, &envs);
        assert_eq!(cli.status.code(), Some(0), "{program}: {cli:?}");
        assert_eq!(
            (text(&cli.stdout), text(&cli.stderr)),
            (
                fs::read_to_string(fixtures.join(format!("{program}.plur.none.out"))).unwrap(),
                fs::read_to_string(fixtures.join(format!("{program}.plur.none.err"))).unwrap()
            ),
            "{program}: plural names, a setting that is not UTF-8"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

/// The library's `run_to_file` makes `-out` where the command line does (after the subjects are
/// read, before the options are processed): the subjects' warnings come before the error of a
/// file that cannot be made (in tblastx's `-logfile` file too), an error of the options after it
/// leaves an empty file, an error of the subjects before it no file; and the parser checks its
/// name as NCBI's `-out` (a name of 256 bytes or more), dry run or not.
#[test]
fn run_to_file_makes_out_where_the_command_line_does() {
    let dir = scratch("out");
    fs::write(dir.join("empty.fa"), "").unwrap();
    let long = "a".repeat(300);
    for (program, query, _, warned) in SEARCHES {
        // The subjects' warnings, then the error of `-out` (a missing directory).
        let mut search = vec![program, "-query", query, "-subject", warned];
        let mut missing = search.clone();
        missing.extend(["-out", "nodir/out.txt"]);
        let cli = both(&dir, &words(&missing), &[]);
        assert_ne!(cli.status.code(), Some(0), "{program}");
        assert!(
            text(&cli.stderr).ends_with(
                "Command line argument error: Argument \"out\". File is not accessible:  `nodir/out.txt'\n"
            ) && text(&cli.stderr).contains("Hyphens are invalid"),
            "{program}: {cli:?}"
        );
        // tblastx's run messages go to `-logfile`: the same file from both.
        if program == "tblastx" {
            let mut logged = missing.clone();
            logged.extend(["-logfile", "log.txt"]);
            let mut logs = Vec::new();
            for binary in [PathBuf::from(BIN), example_binary("library_command_line")] {
                let _ = fs::remove_file(dir.join("log.txt"));
                run_binary_in(&binary, &dir, &words(&logged), &[]);
                logs.push(fs::read_to_string(dir.join("log.txt")).unwrap());
            }
            assert_eq!(logs[1], logs[0], "{program} -logfile");
            assert!(
                logs[0].contains("Hyphens are invalid")
                    && logs[0].contains("File is not accessible")
            );
        }
        // An error of the formatting handler (AIRR, outfmt 19) after `-out` is made: an empty file
        // from both.
        search.extend(["-out", "out.txt", "-outfmt", "19"]);
        for binary in [PathBuf::from(BIN), example_binary("library_command_line")] {
            let _ = fs::remove_file(dir.join("out.txt"));
            let run = run_binary_in(&binary, &dir, &words(&search), &[]);
            assert_ne!(run.status.code(), Some(0), "{program} {binary:?}");
            assert_eq!(
                fs::metadata(dir.join("out.txt"))
                    .map(|meta| meta.len())
                    .ok(),
                Some(0),
                "{program} {binary:?}: outfmt 19"
            );
        }
        // An error of the subjects before it: no file from either.
        let empty = words(&[
            program, "-query", query, "-subject", "empty.fa", "-out", "out.txt",
        ]);
        for binary in [PathBuf::from(BIN), example_binary("library_command_line")] {
            let _ = fs::remove_file(dir.join("out.txt"));
            let run = run_binary_in(&binary, &dir, &empty, &[]);
            assert_ne!(run.status.code(), Some(0), "{program} {binary:?}");
            assert!(
                !dir.join("out.txt").exists(),
                "{program} {binary:?}: empty subjects"
            );
        }
        // A file name of 256 bytes or more: the parser's error, with `-dryrun` too.
        for dry in [false, true] {
            let mut long_out = vec![program, "-query", query, "-subject", warned, "-out", &long];
            if dry {
                long_out.push("-dryrun");
            }
            let library = run_binary_in(
                &example_binary("library_command_line"),
                &dir,
                &words(&long_out),
                &[],
            );
            assert_eq!(library.status.code(), Some(1), "{program} {dry}");
            assert!(
                text(&library.stderr).contains("Illegal value, expected file name length < 256"),
                "{program} {dry}: {library:?}"
            );
        }
        // The BLAST archive (outfmt 11) of an error before `-out` is made (the empty subjects):
        // NCBI's `PrintErrorArchive` makes the file and writes the error archive, also when
        // `ARCHIVE_FORMAT` names no encoding (an empty file then). The same file from both.
        let archive = words(&[
            program, "-query", query, "-subject", "empty.fa", "-out", "out.asn", "-outfmt", "11",
        ]);
        for envs in [vec![], vec![("ARCHIVE_FORMAT", OsString::from("foo"))]] {
            let mut files = Vec::new();
            for binary in [PathBuf::from(BIN), example_binary("library_command_line")] {
                let _ = fs::remove_file(dir.join("out.asn"));
                let run = run_binary_in(&binary, &dir, &archive, &envs);
                assert_ne!(run.status.code(), Some(0), "{program} {binary:?} {envs:?}");
                files.push(fs::read(dir.join("out.asn")).ok());
            }
            // psiblast too (line PS item PS15b-1: its error archive, the subjects' engine error).
            assert!(
                files[0].is_some(),
                "{program} {envs:?}: the command line's file"
            );
            assert_eq!(files[1], files[0], "{program} {envs:?}: outfmt 11");
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

const MULTI_FORMAT_CHILD: &str = "LOSAT_TEST_MULTI_FORMAT_TAXDB_CHILD";

fn run_multi_format_child(case: &str, taxdb_variant: &str, extra_envs: &[(&str, &str)]) {
    let dir = std::env::temp_dir().join(format!("losat-multi-taxdb-{case}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/taxdb")
        .join(taxdb_variant);
    let mut cmd = Command::new(std::env::current_exe().unwrap());
    cmd.env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .env("BLASTDB", &taxdb)
        .env(MULTI_FORMAT_CHILD, case)
        .current_dir(&dir)
        .args([
            "library_out_taxdb::multi_format_taxdb_child",
            "--exact",
            "--nocapture",
            "--test-threads",
            "1",
        ]);
    for &(name, value) in extra_envs {
        cmd.env(name, value);
    }
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        cmd.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    let output = cmd.output().unwrap();
    let _ = fs::remove_dir_all(&dir);
    let stdout = text(&output.stdout);
    assert!(
        output.status.success() && stdout.contains("1 passed; 0 failed"),
        "{case} ({taxdb_variant}):\nstdout:\n{stdout}\nstderr:\n{}",
        text(&output.stderr)
    );
}

/// Multi-format `losat::api::run_collect` with secondary XML2/JSON (`["6", "16"]`, `["5", "16"]`)
/// and primary XML2/JSON (`["16", "6"]`, `["16", "5"]`):
/// - opening-texts taxdb (`bad_magic`): opening texts are written once, not repeated for format 1;
/// - non-collision `offset_last_past_end` (`tm_cache151`): volume offset error is not repeated;
/// - cache-collision `offset_last_past_end` (`tm_cache_gi` / `tm_cache_ngi`): secondary XML2/JSON
///   (`["6", "16"]`, `["5", "16"]`) refuses across `blastx`, `blastn`, and `blastp` (both 1-query
///   and 2-query `BATCH_SIZE=1`, matching `dev`), while primary XML2/JSON (`["16", "6"]`,
///   `["16", "5"]`) succeeds with stdout matching single-format runs and the same 12 diagnostics
///   lines as single `["16"]`.
#[test]
fn multi_format_library_secondary_xml2_shares_open_and_refuses_collision() {
    run_multi_format_child("bad_magic", "bad_magic", &[]);
    run_multi_format_child("tm_cache151", "offset_last_past_end", &[]);
    run_multi_format_child("collision_q1", "offset_last_past_end", &[]);
    run_multi_format_child(
        "collision_q2_b1",
        "offset_last_past_end",
        &[("BATCH_SIZE", "1")],
    );
}

fn first_fasta_record(path: &Path) -> Vec<u8> {
    let raw = fs::read(path).unwrap();
    let second = raw
        .iter()
        .enumerate()
        .skip(1)
        .find(|&(i, &b)| b == b'>' && raw[i - 1] == b'\n')
        .map(|(i, _)| i)
        .unwrap_or(raw.len());
    raw[..second].to_vec()
}

fn run_db_formats(
    program: losat::api::Program,
    db: &Path,
    query: &[u8],
    formats: &[&str],
) -> Result<losat::api::Outcome, losat::api::Error> {
    use losat::api::{read_inputs, run_collect, validate, Inputs};
    let search = validate(program, &["-num_threads", "1"]).unwrap();
    let (queries, _) = read_inputs(&search, query, b"").unwrap();
    run_collect(
        &search,
        Inputs::database(&queries, db.to_str().unwrap()),
        formats,
    )
}

fn assert_pair_reports_match(
    program: losat::api::Program,
    db: &Path,
    query: &[u8],
    first: &str,
    single_16: &losat::api::Outcome,
) -> (losat::api::Outcome, losat::api::Outcome) {
    let single_first = run_db_formats(program, db, query, &[first]).unwrap();
    let both = run_db_formats(program, db, query, &[first, "16"]).unwrap();
    assert_eq!(
        both.reports[0].bytes, single_first.reports[0].bytes,
        "{program:?} {first}"
    );
    assert_eq!(
        both.reports[1].bytes, single_16.reports[0].bytes,
        "{program:?} {first} 16"
    );
    (single_first, both)
}

fn assert_cache_refusal(
    program: losat::api::Program,
    db: &Path,
    query: &[u8],
    formats: &[&str],
    label: &str,
) {
    let err = run_db_formats(program, db, query, formats)
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("256-slot taxonomy cache"),
        "{label} ({program:?} {formats:?}): {err}"
    );
}

#[test]
fn multi_format_taxdb_child() {
    use losat::api::Program;
    let Ok(case) = std::env::var(MULTI_FORMAT_CHILD) else {
        return;
    };
    let db_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let prot_q = fs::read(db_dir.join("tm_cache_q.faa")).unwrap();
    let nucl_q = first_fasta_record(&db_dir.join("tm_cache_ngi.fna"));
    let offset_msg = "Offset error at end of taxdb file.";
    let bad_magic_msg = "Tax database file has wrong magic number.";

    match case.as_str() {
        "bad_magic" => {
            for (program, db_name, query) in [
                (Program::Blastp, "tm_cache_gi", prot_q.as_slice()),
                (Program::Blastn, "tm_cache_ngi", nucl_q.as_slice()),
            ] {
                let db = db_dir.join(db_name);
                let single_16 = run_db_formats(program, &db, query, &["16"]).unwrap();
                for first in ["6", "5"] {
                    let (single_first, both) =
                        assert_pair_reports_match(program, &db, query, first, &single_16);
                    let diag = text(&both.diagnostics);
                    assert_eq!(
                        diag,
                        text(&single_first.diagnostics),
                        "{program:?} [{first}, 16] bad_magic diag vs single:\n{diag}"
                    );
                    assert_eq!(
                        diag.matches(bad_magic_msg).count(),
                        1,
                        "{program:?} [{first}, 16] bad_magic diag:\n{diag}"
                    );
                }
            }
        }
        "tm_cache151" => {
            // `tm_cache151` has 3 sequences: st1 (99990), st2 (151), st3 (99990) — no cache
            // collision. Format 0's prefetch emits 1 volume-load error, format 1 shares primary
            // `volume_lookups` (0 duplicate volume errors) and emits 2 Hit errors (st1, st3) = 3.
            let db = db_dir.join("tm_cache151");
            let single_16 = run_db_formats(Program::Blastp, &db, &prot_q, &["16"]).unwrap();
            for first in ["6", "5"] {
                let (_, both) =
                    assert_pair_reports_match(Program::Blastp, &db, &prot_q, first, &single_16);
                let diag = text(&both.diagnostics);
                assert_eq!(
                    diag,
                    text(&single_16.diagnostics),
                    "blastp [{first}, 16] tm_cache151 diag vs single 16:\n{diag}"
                );
                assert_eq!(
                    diag.matches(offset_msg).count(),
                    3,
                    "blastp [{first}, 16] tm_cache151 diag:\n{diag}"
                );
            }
        }
        "collision_q1" => {
            // Cache-colliding volumes refuse secondary XML2/JSON (`["6","16"]`, `["5","16"]`)
            // across blastx, blastn, and blastp (matching dev), while primary XML2/JSON
            // (`["16","6"]`, `["16","5"]`) succeeds with the same reports and 12 diagnostic lines
            // as single `["16"]`.
            for (program, db_name, query) in [
                (Program::Blastx, "tm_cache_gi", nucl_q.as_slice()),
                (Program::Blastn, "tm_cache_ngi", nucl_q.as_slice()),
                (Program::Blastp, "tm_cache_gi", prot_q.as_slice()),
            ] {
                let db = db_dir.join(db_name);
                let single_16 = run_db_formats(program, &db, query, &["16"]).unwrap();
                let single_16_diag = text(&single_16.diagnostics);
                assert_eq!(
                    single_16_diag.matches(offset_msg).count(),
                    12,
                    "{program:?} [16] {db_name} diag:\n{single_16_diag}"
                );
                for other in ["6", "5"] {
                    assert_cache_refusal(program, &db, query, &[other, "16"], db_name);
                    let single_other = run_db_formats(program, &db, query, &[other]).unwrap();
                    let both = run_db_formats(program, &db, query, &["16", other]).unwrap();
                    assert_eq!(
                        both.reports[0].bytes, single_16.reports[0].bytes,
                        "{program:?} [16, {other}] report 0"
                    );
                    assert_eq!(
                        both.reports[1].bytes, single_other.reports[0].bytes,
                        "{program:?} [16, {other}] report 1"
                    );
                    let diag = text(&both.diagnostics);
                    assert_eq!(
                        diag, single_16_diag,
                        "{program:?} [16, {other}] {db_name} diag vs single 16:\n{diag}"
                    );
                    assert_eq!(
                        diag.matches(offset_msg).count(),
                        12,
                        "{program:?} [16, {other}] {db_name} diag:\n{diag}"
                    );
                }
            }
        }
        "collision_q2_b1" => {
            // 2 identical queries with BATCH_SIZE=1: secondary XML2/JSON on a colliding volume
            // refuses for both `["6","16"]` and `["5","16"]`.
            let mut two_nucl = nucl_q.clone();
            two_nucl.extend_from_slice(b"\n");
            two_nucl.extend_from_slice(&nucl_q);
            for (program, db_name) in [
                (Program::Blastx, "tm_cache_gi"),
                (Program::Blastn, "tm_cache_ngi"),
            ] {
                let db = db_dir.join(db_name);
                for first in ["6", "5"] {
                    assert_cache_refusal(
                        program,
                        &db,
                        &two_nucl,
                        &[first, "16"],
                        "2q BATCH_SIZE=1",
                    );
                }
            }
        }
        other => panic!("unknown case {other}"),
    }
}
