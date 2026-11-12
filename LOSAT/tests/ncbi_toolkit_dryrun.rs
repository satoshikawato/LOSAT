//! The NCBI C++ Toolkit option `-dryrun` of the BLAST programs, makeblastdb and blastdbcmd (line
//! RG item E1b), against NCBI BLAST+ 2.17.0, the official Linux x64 release (probes of 2026-11-13,
//! `remaining-gaps-20261112/E/SURVEY_E1b.md`, run with `env -i` in a fresh directory: the E1
//! probe's `dry*` cases and the sweep of every argument definition of NCBI's `-xmlhelp`).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707, 292-296 (quoted at
//! `cli::main_entry`) and 926-1001 (quoted at `ncbi_environment::preparse`).

#![cfg(unix)]

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-toolkit-dryrun-{name}-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let inputs = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/outfmt11");
    for name in [
        "arch_bp_q.faa",
        "arch_bp_s.faa",
        "arch_bx_q.fna",
        "arch_bx_s.faa",
        "arch_tbn_q.faa",
        "arch_tbn_s.fna",
        "arch_tbx_q.fna",
        "arch_tbx_s.fna",
    ] {
        fs::copy(inputs.join(name), dir.join(name)).unwrap();
    }
    dir
}

/// `losat <words>` in `dir` with an environment of NCBI's probes (`env -i PATH HOME`) and `envs`.
fn run_in(dir: &Path, words: &[&str], envs: &[(&str, &str)]) -> Output {
    let mut command = Command::new(BIN);
    command
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .envs(envs.iter().copied())
        .current_dir(dir)
        .args(words)
        .stdin(std::process::Stdio::null());
    // The `--all-features` build's BLASTX diagnostic probe needs its log (as in CI).
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        command.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    command.output().unwrap()
}

/// The words of each program's run in the probes (rpsblast's database is never opened).
const RUNS: [(&str, &[&str]); 9] = [
    (
        "blastn",
        &["-query", "arch_tbx_q.fna", "-subject", "arch_tbx_s.fna"],
    ),
    (
        "blastp",
        &["-query", "arch_bp_q.faa", "-subject", "arch_bp_s.faa"],
    ),
    (
        "blastx",
        &["-query", "arch_bx_q.fna", "-subject", "arch_bx_s.faa"],
    ),
    (
        "tblastn",
        &["-query", "arch_tbn_q.faa", "-subject", "arch_tbn_s.fna"],
    ),
    (
        "tblastx",
        &["-query", "arch_tbx_q.fna", "-subject", "arch_tbx_s.fna"],
    ),
    (
        "psiblast",
        &["-query", "arch_bp_q.faa", "-subject", "arch_bp_s.faa"],
    ),
    (
        "rpsblast",
        &["-query", "arch_bp_q.faa", "-db", "rpsdb/rps_mini"],
    ),
    // Its `-out` is the database's name: the cases below add one.
    ("makeblastdb", &["-in", "arch_bp_s.faa", "-dbtype", "prot"]),
    ("blastdbcmd", &["-db", "nodb", "-entry", "all"]),
];

const SEARCHES: usize = 7;

fn names(dir: &Path) -> BTreeSet<String> {
    fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// `program words...` and `program` with `words` changed by `edit`.
fn words_of<'a>(
    program: &'a str,
    words: &[&'a str],
    edit: &[&'a str],
    first: bool,
) -> Vec<&'a str> {
    let mut all = vec![program];
    if first {
        all.extend_from_slice(edit);
        all.extend_from_slice(words);
    } else {
        all.extend_from_slice(words);
        all.extend_from_slice(edit);
    }
    all
}

/// NCBI BLAST+ 2.17.0 (probes `dry`, `dry_first`, `dry_twice`, `dry_qerr`, `dry_noquery`,
/// `dry_dbmissing`, `dry_old_out`, `dry_ncbirc`, `dry_env`, `makeblastdb.dry_noin`,
/// `blastdbcmd.dry_nodb`, `blastdbcmd.dry_out`): after the arguments, nothing runs. No input is
/// read (a missing `-query`, `-in` or `-db`, or none), no file is made or truncated (`-out`), a
/// setting that NCBI reads in `Run` (`BL2SEQ_LEGACY`) is not looked at: nothing written, exit
/// status 0. The registry and environment settings that LOSAT refuses in a run are refused with
/// `-dryrun` too (line RG DECISIONS E1b-6: some of them, `DIAG_POST_LEVEL` and the registry's
/// `[DEBUG]`, make NCBI write `DryRun`'s `Info` post).
#[test]
fn a_dry_run_reads_the_arguments_and_runs_nothing() {
    let dir = scratch("nothing");
    let before = names(&dir);
    for (index, (program, words)) in RUNS.iter().enumerate() {
        let mut cases: Vec<(Vec<&str>, Vec<(&str, &str)>)> = vec![
            (
                words_of(program, words, &["-out", "out.txt", "-dryrun"], false),
                vec![],
            ),
            (words_of(program, words, &["-dryrun"], true), vec![]),
            (
                words_of(program, words, &["-dryrun", "-dryrun"], false),
                vec![],
            ),
            (
                words_of(program, words, &["-dryrun"], false),
                vec![("BL2SEQ_LEGACY", "1")],
            ),
        ];
        if index < SEARCHES {
            // A missing query file, and no query (the standard input is not read).
            let mut missing = vec![*program, "-query", "missing.fa"];
            missing.extend_from_slice(&words[2..]);
            missing.push("-dryrun");
            cases.push((missing, vec![]));
            let mut no_query = vec![*program];
            no_query.extend_from_slice(&words[2..]);
            no_query.push("-dryrun");
            cases.push((no_query, vec![]));
        }
        if *program == "makeblastdb" {
            cases.push((
                vec![
                    program,
                    "-in",
                    "missing.fa",
                    "-dbtype",
                    "nucl",
                    "-out",
                    "mdb",
                    "-dryrun",
                ],
                vec![],
            ));
        }
        if *program == "blastn" {
            cases.push((
                vec![
                    program,
                    "-query",
                    "arch_tbx_q.fna",
                    "-db",
                    "nodb",
                    "-outfmt",
                    "6",
                    "-dryrun",
                ],
                vec![],
            ));
        }
        for (words, envs) in &cases {
            let run = run_in(&dir, words, envs);
            assert_eq!(
                (run.status.code(), text(&run.stdout), text(&run.stderr)),
                (Some(0), String::new(), String::new()),
                "{words:?} {envs:?}"
            );
            assert_eq!(names(&dir), before, "{words:?}: a file was made");
        }
        fs::write(dir.join("out.txt"), "OLD OUT\n").unwrap();
        let run = run_in(
            &dir,
            &words_of(program, words, &["-out", "out.txt", "-dryrun"], false),
            &[],
        );
        assert_eq!(run.status.code(), Some(0), "{program}: {run:?}");
        assert_eq!(
            fs::read_to_string(dir.join("out.txt")).unwrap(),
            "OLD OUT\n",
            "{program}"
        );
        fs::remove_file(dir.join("out.txt")).unwrap();
        for (registry, envs) in [
            ("[BLAST]\nLONG_SEQID=1\n", &[][..]),
            ("[DEBUG]\nDIAG_POST_LEVEL=Info\n", &[]),
            ("", &[("DIAG_POST_LEVEL", "Info")]),
        ] {
            if !registry.is_empty() {
                fs::write(dir.join(".ncbirc"), registry).unwrap();
            }
            let plain = run_in(&dir, &words_of(program, words, &["-out", "o"], false), envs);
            let dry = run_in(&dir, &words_of(program, words, &["-dryrun"], false), envs);
            assert_eq!(
                plain.status.code(),
                Some(1),
                "{program} {registry:?} {envs:?}: {plain:?}"
            );
            assert!(
                text(&plain.stderr).contains("not supported by LOSAT"),
                "{plain:?}"
            );
            assert_eq!(
                (dry.status.code(), text(&dry.stderr), dry.stdout.len()),
                (Some(1), text(&plain.stderr), 0),
                "{program} {registry:?} {envs:?}"
            );
            let _ = fs::remove_file(dir.join(".ncbirc"));
            let _ = fs::remove_file(dir.join("o"));
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI BLAST+ 2.17.0 (probe `dry_log`): the `-logfile` file is made as without `-dryrun`, and
/// stays empty.
#[test]
fn a_dry_run_makes_the_logfile() {
    let dir = scratch("log");
    for (program, words) in RUNS {
        let run = run_in(
            &dir,
            &words_of(program, words, &["-dryrun", "-logfile", "log.txt"], false),
            &[],
        );
        assert_eq!(
            (run.status.code(), run.stdout.len(), run.stderr.len()),
            (Some(0), 0, 0),
            "{program}: {run:?}"
        );
        assert_eq!(fs::read(dir.join("log.txt")).unwrap(), b"", "{program}");
        fs::remove_file(dir.join("log.txt")).unwrap();
    }
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI BLAST+ 2.17.0 (probes `dry_argerr`, `dry_unk`, `dry_badtask`, `dry_incompat`,
/// `dry_outval`, `dry_after_dd`, `dry_eq`, `dry_pos`, `makeblastdb.dry_argerr`,
/// `makeblastdb.dry_nodbtype`, `blastdbcmd.dry_argerr`, `blastdbcmd.dry_excl`): an error in the
/// arguments' definitions is NCBI's USAGE and exit status 1 with `-dryrun` as without it. In
/// LOSAT it is the same error as without `-dryrun`: the parser's message and exit status 2
/// (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 1), NCBI's bytes in blastdbcmd (whose reader is
/// NCBI's).
#[test]
fn argument_errors_are_reported_as_without_dryrun() {
    let dir = scratch("argerr");
    for (program, words) in RUNS {
        let mut errors: Vec<&[&str]> =
            vec![&["-bogus"], &["-outfmt", "6", "-outfmt", "6"], &["word"]];
        match program {
            "makeblastdb" => errors.extend([
                &["-blastdb_version", "3"][..],
                &["-title", "a", "-title", "b"],
            ]),
            "blastdbcmd" => errors.extend([&["-line_length", "0"][..], &["-info"]]),
            _ => errors.extend([&["-evalue", "abc"][..], &["-task", "bogus"]]),
        }
        let status = if program == "blastdbcmd" { 1 } else { 2 };
        for error in errors {
            let plain = run_in(&dir, &words_of(program, words, error, false), &[]);
            let mut dry_words = words_of(program, words, error, false);
            dry_words.push("-dryrun");
            let dry = run_in(&dir, &dry_words, &[]);
            assert_eq!(
                plain.status.code(),
                Some(status),
                "{program} {error:?}: {plain:?}"
            );
            assert_eq!(
                (dry.status.code(), &dry.stdout, &dry.stderr),
                (plain.status.code(), &plain.stdout, &plain.stderr),
                "{program} {error:?}"
            );
        }
        // `-out -dryrun`: `-out` has no value; `-- -dryrun` and `-dryrun=1` are not the word.
        for dry_words in [
            words_of(program, words, &["-out", "-dryrun"], false),
            words_of(program, words, &["--", "-dryrun"], false),
            words_of(program, words, &["-dryrun=1"], false),
        ] {
            let dry = run_in(&dir, &dry_words, &[]);
            assert_eq!(dry.status.code(), Some(status), "{dry_words:?}: {dry:?}");
            assert!(
                dry.stdout.is_empty() && !dry.stderr.is_empty(),
                "{dry_words:?}"
            );
        }
    }
    // NCBI's database arguments exclude `-window_masker_taxid` with `-no_taxid_expansion`
    // (blast_args.cpp:2285; probes `wmtaxid_ntx`, `dry_wmtaxid_ntx`, the sweep's
    // `blastn.x.no_taxid_expansion.window_masker_taxid`).
    let words = [
        "blastn",
        "-query",
        "arch_tbx_q.fna",
        "-db",
        "nodb",
        "-window_masker_taxid",
        "1",
        "-no_taxid_expansion",
    ];
    for words in [&words[..], &[&words[..], &["-dryrun"]].concat()] {
        let run = run_in(&dir, words, &[]);
        assert_eq!(run.status.code(), Some(2), "{words:?}: {run:?}");
        assert!(text(&run.stderr).contains("-no_taxid_expansion"), "{run:?}");
    }
    // makeblastdb's `-dbtype` is mandatory.
    let dry = run_in(
        &dir,
        &[
            "makeblastdb",
            "-in",
            "arch_bp_s.faa",
            "-out",
            "mdb",
            "-dryrun",
        ],
        &[],
    );
    assert_eq!(dry.status.code(), Some(2), "{dry:?}");
    let _ = fs::remove_dir_all(&dir);
}

const SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

/// NCBI BLAST+ 2.17.0 (probes `dry_conf_missing`, `dry_conf_bad`, `dry_ini_bad`, `dry_h`,
/// `dry_help`, `dry_ver`): the configuration is loaded before the arguments are read, `-dryrun`
/// or not (a `-conffile` that is not a file aborts with exit status 255, a malformed one or a
/// malformed `<program>.ini` ends with exit status 2), and `-h`, `-help` and `-version` answer
/// as without `-dryrun`.
#[test]
fn the_configuration_help_and_version_come_as_without_dryrun() {
    let dir = scratch("conf");
    for (index, (program, words)) in RUNS.iter().enumerate() {
        let missing = run_in(
            &dir,
            &words_of(
                program,
                words,
                &["-dryrun", "-conffile", "missing.ini"],
                false,
            ),
            &[],
        );
        assert_eq!(
            (missing.status.code(), text(&missing.stderr)),
            (
                Some(255),
                "Critical: (106.23) Registry file \"missing.ini\" not found, aborting.\n"
                    .to_string()
            ),
            "{program}"
        );
        if index < SEARCHES {
            fs::write(dir.join("c.ini"), "[x\n").unwrap();
            let bad = run_in(
                &dir,
                &words_of(program, words, &["-dryrun", "-conffile", "c.ini"], false),
                &[],
            );
            assert_eq!(
                (bad.status.code(), text(&bad.stderr)),
                (Some(2), SECTION_REPORT.repeat(2)),
                "{program} -conffile"
            );
            fs::remove_file(dir.join("c.ini")).unwrap();
            let ini = dir.join(format!("{program}.ini"));
            fs::write(&ini, "[x\n").unwrap();
            let bad = run_in(&dir, &words_of(program, words, &["-dryrun"], false), &[]);
            assert_eq!(
                (bad.status.code(), text(&bad.stderr)),
                (Some(2), SECTION_REPORT.repeat(2)),
                "{program} .ini"
            );
            fs::remove_file(&ini).unwrap();
        }
        for word in ["-h", "-help", "-version"] {
            let plain = run_in(&dir, &[program, word], &[]);
            let dry = run_in(&dir, &[program, "-dryrun", word], &[]);
            assert_eq!(plain.status.code(), Some(0), "{program} {word}");
            assert!(!plain.stdout.is_empty(), "{program} {word}");
            assert_eq!(
                (dry.status.code(), &dry.stdout, &dry.stderr),
                (plain.status.code(), &plain.stdout, &plain.stderr),
                "{program} -dryrun {word}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI BLAST+ 2.17.0 (E1b probes `makeblastdb.taxid_hex` and `makeblastdb.taxid_16`, and the
/// sweep's `makeblastdb.k.taxid.12`): `-taxid` is NCBI's integer, whose hexadecimal form `0x10`
/// is 16 (every file of the database the same bytes as with `-taxid 16`).
#[test]
fn a_hexadecimal_taxid_is_the_same_taxid() {
    let dir = scratch("taxid");
    // One database name in two directories: the LMDB file holds the name.
    for (sub, taxid) in [("hex", "0x10"), ("dec", "16")] {
        fs::create_dir_all(dir.join(sub)).unwrap();
        let run = run_in(
            &dir.join(sub),
            &[
                "makeblastdb",
                "-in",
                "../arch_bp_s.faa",
                "-dbtype",
                "prot",
                "-out",
                "mdb",
                "-taxid",
                taxid,
            ],
            &[],
        );
        assert_eq!(run.status.code(), Some(0), "{taxid}: {run:?}");
    }
    let files = names(&dir.join("hex"));
    assert_eq!(files, names(&dir.join("dec")));
    assert!(files.contains("mdb.ptf"), "{files:?}");
    for file in &files {
        // The index and metadata files hold the creation date.
        if file.ends_with(".pin") || file.ends_with(".pjs") {
            continue;
        }
        assert_eq!(
            fs::read(dir.join("hex").join(file)).unwrap(),
            fs::read(dir.join("dec").join(file)).unwrap(),
            "{file}"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}
