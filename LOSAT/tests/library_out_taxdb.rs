//! The taxonomy database of a `-subject` search's tabular reports and the `-out` file of the
//! library's `run_to_file`, against NCBI BLAST+ 2.17.0 and the command line (line RG item E1w-2b,
//! `remaining-gaps-20261112/E/SURVEY_E1w2b.md`).
//!
//! The same words go to the command line and, through `examples/library_command_line`
//! (`losat::api::run`, `run_to_file` with `-out`), to the library, in the same directory and
//! environment; a failure of the library ends the example with exit status 1 and the error's
//! message after the call's diagnostics.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104 (quoted at
//! `common::database::subject_taxonomy`) and c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
//! (quoted at `cli::create_out_file`).
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
/// singular names look nothing up for subjects without a taxid, so a `BLASTDB` setting that LOSAT
/// cannot read refuses nothing. The library's searches of records look `taxdb.bti` up as the
/// command line does (they wrote the warning whatever the path held before line RG item E1w-2b).
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
        // The plural names with a setting that LOSAT cannot read: the command line's refusal, the
        // library's too.
        let words = words(&[
            program,
            "-query",
            query,
            "-subject",
            subject,
            "-outfmt",
            "6 qseqid sseqid sscinames",
        ]);
        let cli = both(
            &dir,
            &words,
            &[("NCBI_CONFIG__BLAST__BLASTDB", unreadable.clone())],
        );
        assert_ne!(cli.status.code(), Some(0), "{program}: {cli:?}");
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
            ) && text(&cli.stderr).contains("Warning"),
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
            assert!(logs[0].contains("Warning") && logs[0].contains("File is not accessible"));
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
            assert!(!dir.join(&long).exists());
        }
    }
    let _ = fs::remove_dir_all(&dir);
}
