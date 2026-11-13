//! The NCBI C++ Toolkit's file options `-logfile` and `-conffile` of the BLAST programs (line RG
//! item E1), against NCBI BLAST+ 2.17.0, the official Linux x64 release (probes of 2026-11-12,
//! `remaining-gaps-20261112/E/SURVEY_E1.md`, run with `env -i` in a fresh directory).
//!
//! The library API runs the same start in each call (line RG item E1w-2): the same words through
//! `examples/library_command_line` (`LOSAT::api::run`, `run_to_file`) in the same directory and
//! environment give the command line's standard output, standard error and `-logfile` file; a
//! failure's message is the text that the command line writes to its standard error without
//! `-logfile` (the example ends every failure with exit status 1).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-871, 926-1022 (quoted at
//! `ncbi_environment::log_file`, `preparse` and `conf_file_path`), and
//! c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at `ncbi_environment::run_messages_follow_log`).
#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use crate::library_examples::example_binary;

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-toolkit-files-{name}-{}", std::process::id()));
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
    run_binary_in(Path::new(BIN), dir, words, envs)
}

/// `run_in` through the library API: `examples/library_command_line <words>`.
fn library_in(dir: &Path, words: &[&str], envs: &[(&str, &str)]) -> Output {
    run_binary_in(&example_binary("library_command_line"), dir, words, envs)
}

/// `run_in` of `binary` (the command line, or the library's example).
fn run_binary_in(binary: &Path, dir: &Path, words: &[&str], envs: &[(&str, &str)]) -> Output {
    let mut command = Command::new(binary);
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

/// The query and subject of each program's search in the probes.
const SEARCHES: [(&str, [&str; 4]); 6] = [
    (
        "blastp",
        ["-query", "arch_bp_q.faa", "-subject", "arch_bp_s.faa"],
    ),
    (
        "blastx",
        ["-query", "arch_bx_q.fna", "-subject", "arch_bx_s.faa"],
    ),
    (
        "tblastn",
        ["-query", "arch_tbn_q.faa", "-subject", "arch_tbn_s.fna"],
    ),
    (
        "tblastx",
        ["-query", "arch_tbx_q.fna", "-subject", "arch_tbx_s.fna"],
    ),
    (
        "psiblast",
        ["-query", "arch_bp_q.faa", "-subject", "arch_bp_s.faa"],
    ),
    (
        "blastn",
        ["-query", "arch_tbx_q.fna", "-subject", "arch_tbx_s.fna"],
    ),
];

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Runs `words` through the library in `dir` with `envs` (`library_in`) and checks it against
/// `cli`, the command line's run of the same words: the same standard output and standard error,
/// and exit status 0 where the command line's is 0 (1, the example's, for a failure).
fn check_library(dir: &Path, words: &[&str], envs: &[(&str, &str)], cli: &Output) {
    let library = library_in(dir, words, envs);
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
}

/// NCBI BLAST+ 2.17.0 (probes `log_warn`, `log_qerr`, `log_empty`, `log_nodir`, `log_old`): the
/// warning of `-max_target_seqs 1` and the error of a missing `-query` go to the `-logfile` file
/// in tblastx, rpsblast and rpstblastn (its missing `-query` here) and to the standard error in
/// the other programs; the file is made
/// (appended to) in every program, the output is the same as without it.
#[test]
fn run_messages_go_where_ncbi_sends_them() {
    let dir = scratch("route");
    for (program, search) in SEARCHES {
        let to_log = program == "tblastx";
        let warning = format!("Warning: [{program}] Examining 5 or more matches is recommended\n");
        let mut words = vec![program];
        words.extend(search);
        words.extend(["-outfmt", "6", "-max_target_seqs", "1"]);
        let plain = run_in(&dir, &words, &[]);
        assert_eq!(plain.status.code(), Some(0), "{program}: {plain:?}");
        assert_eq!(text(&plain.stderr), warning, "{program}");
        fs::write(dir.join("log.txt"), "OLD LINE\n").unwrap();
        let mut logged_words = words.clone();
        logged_words.extend(["-logfile", "log.txt"]);
        let logged = run_in(&dir, &logged_words, &[]);
        assert_eq!(logged.status.code(), Some(0), "{program}: {logged:?}");
        assert_eq!(logged.stdout, plain.stdout, "{program}");
        let file = fs::read_to_string(dir.join("log.txt")).unwrap();
        if to_log {
            assert_eq!(text(&logged.stderr), "", "{program}");
            assert_eq!(file, format!("OLD LINE\n{warning}"), "{program}");
        } else {
            assert_eq!(text(&logged.stderr), warning, "{program}");
            assert_eq!(file, "OLD LINE\n", "{program}");
        }
        // The library: the same file, the same streams (`run`, and `run_to_file` with `-out`).
        for out in [&[][..], &["-out", "out.txt"]] {
            fs::write(dir.join("log.txt"), "OLD LINE\n").unwrap();
            let mut library_words = logged_words.clone();
            library_words.extend_from_slice(out);
            let library = library_in(&dir, &library_words, &[]);
            assert_eq!(library.status.code(), Some(0), "{program}: {library:?}");
            let stdout = match out {
                [] => library.stdout,
                _ => fs::read(dir.join("out.txt")).unwrap(),
            };
            assert_eq!(
                (
                    stdout,
                    text(&library.stderr),
                    fs::read_to_string(dir.join("log.txt")).unwrap()
                ),
                (logged.stdout.clone(), text(&logged.stderr), file.clone()),
                "{program} library {out:?}"
            );
            let _ = fs::remove_file(dir.join("out.txt"));
        }
        // `-logfile ""` sends the diagnostics nowhere; a file in a missing directory leaves them
        // on the standard error. The library's too, and `-` (the standard error,
        // `ncbi_environment::diag_target`) is the command line's there.
        for (value, quiet) in [
            ("", Some(to_log)),
            ("nodir/log.txt", Some(false)),
            ("-", None),
        ] {
            let mut words = words.clone();
            words.extend(["-logfile", value]);
            let run = run_in(&dir, &words, &[]);
            assert_eq!(run.status.code(), Some(0), "{program} {value:?}");
            assert_eq!(run.stdout, plain.stdout, "{program} {value:?}");
            if let Some(quiet) = quiet {
                let expected = if quiet { "" } else { warning.as_str() };
                assert_eq!(text(&run.stderr), expected, "{program} {value:?}");
            }
            let library = library_in(&dir, &words, &[]);
            assert_eq!(
                (
                    library.status.code(),
                    &library.stdout,
                    text(&library.stderr)
                ),
                (Some(0), &run.stdout, text(&run.stderr)),
                "{program} {value:?} library"
            );
        }
        assert!(!dir.join("nodir").exists());
        assert!(!dir.join("-").exists());
    }
    // The library takes no `-query` file: its caller reads the records (the example does, before
    // the call), so the missing file below is the command line's only.
    // A missing `-query`, rpsblast and rpstblastn too (their error comes before the database is
    // opened; rpstblastn: PS12a's probe `log_qerr` of 2026-11-13, NCBI's official 2.17.0, which
    // writes it into the file as rpsblast does, rpstblastn_app.cpp:111-114).
    let error =
        "Command line argument error: Argument \"query\". File is not accessible:  `missing.fa'\n";
    for program in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpsblast",
        "rpstblastn",
    ] {
        let _ = fs::remove_file(dir.join("log.txt"));
        let other = if matches!(program, "rpsblast" | "rpstblastn") {
            ["-db", "nodb"]
        } else {
            let subject = match program {
                "blastn" | "tblastx" => "arch_tbx_s.fna",
                "tblastn" => "arch_tbn_s.fna",
                _ => "arch_bp_s.faa",
            };
            ["-subject", subject]
        };
        let mut words = vec![program, "-query", "missing.fa"];
        words.extend(other);
        words.extend(["-logfile", "log.txt"]);
        let run = run_in(&dir, &words, &[]);
        assert_eq!(run.status.code(), Some(1), "{program}: {run:?}");
        let file = fs::read_to_string(dir.join("log.txt")).unwrap();
        if matches!(program, "tblastx" | "rpsblast" | "rpstblastn") {
            assert_eq!(
                (text(&run.stderr), file),
                (String::new(), error.to_string()),
                "{program}"
            );
        } else {
            assert_eq!(
                (text(&run.stderr), file),
                (error.to_string(), String::new()),
                "{program}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI's report of a registry file with the line `[x`, which ends the program (exit status 2),
/// and its replay.
const SECTION_REPORT: &str = "Error: (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eSection) Invalid registry section (']' is missing): `[x' (m_Pos = 1)\n";

/// NCBI BLAST+ 2.17.0 (probes `conf_ok`, `conf_eq`, `conf_ini_bad`, `conf_bad`, `conf_dot`,
/// `conf_nolocal_dot`): the file of `-conffile` is read in place of `<program>.ini` (a malformed
/// `<program>.ini` is then not read); a malformed one ends the program as a malformed
/// `<program>.ini` does. A name without a directory is looked for in the registry search path:
/// with `NCBI_DONT_USE_LOCAL_CONFIG` NCBI does not find `c.ini` in the current directory and goes
/// on without it, but reads `./c.ini`.
#[test]
fn the_conffile_is_read_in_place_of_the_programs_ini() {
    let dir = scratch("conf");
    for (program, search) in SEARCHES {
        let mut words = vec![program];
        words.extend(search);
        words.extend(["-outfmt", "6"]);
        let plain = run_in(&dir, &words, &[]);
        assert_eq!(plain.status.code(), Some(0), "{program}: {plain:?}");
        fs::write(dir.join(format!("{program}.ini")), "[x\n").unwrap();
        for conffile in [&["-conffile", "c.ini"][..], &["-conffile=c.ini"]] {
            fs::write(dir.join("c.ini"), "[BLAST]\n").unwrap();
            let mut with = words.clone();
            with.extend_from_slice(conffile);
            let read = run_in(&dir, &with, &[]);
            assert_eq!(
                (read.status.code(), text(&read.stderr), &read.stdout),
                (Some(0), String::new(), &plain.stdout),
                "{program} {conffile:?}"
            );
            check_library(&dir, &with, &[], &read);
            fs::write(dir.join("c.ini"), "[x\n").unwrap();
            let failed = run_in(&dir, &with, &[]);
            assert_eq!(
                (failed.status.code(), text(&failed.stderr)),
                (Some(2), SECTION_REPORT.repeat(2)),
                "{program} {conffile:?}"
            );
            assert!(failed.stdout.is_empty());
            check_library(&dir, &with, &[], &failed);
        }
        // With `-logfile`: the library's file is the command line's, and its error is the text
        // that the command line writes without the file.
        let mut logged = words.clone();
        logged.extend(["-conffile", "c.ini", "-logfile", "log.txt"]);
        let _ = fs::remove_file(dir.join("log.txt"));
        let cli = run_in(&dir, &logged, &[]);
        assert_eq!(cli.status.code(), Some(2), "{program}: {cli:?}");
        let cli_log = fs::read(dir.join("log.txt")).unwrap();
        assert!(!cli_log.is_empty(), "{program}");
        fs::remove_file(dir.join("log.txt")).unwrap();
        let library = library_in(&dir, &logged, &[]);
        assert_eq!(
            (
                library.status.code(),
                text(&library.stderr),
                fs::read(dir.join("log.txt")).unwrap()
            ),
            (Some(1), SECTION_REPORT.repeat(2), cli_log),
            "{program} -logfile: the library"
        );
        fs::remove_file(dir.join("log.txt")).unwrap();
        // A registry that turns NCBI's data loaders off: the library reads its records with them,
        // and they read the same (`api::start::CallStart::settings`).
        fs::write(dir.join("c.ini"), "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        let mut no_loaders = words.clone();
        no_loaders.extend(["-conffile", "c.ini"]);
        let cli = run_in(&dir, &no_loaders, &[]);
        assert_eq!(cli.status.code(), Some(0), "{program}: {cli:?}");
        check_library(&dir, &no_loaders, &[], &cli);
        fs::write(dir.join("c.ini"), "[x\n").unwrap();
        fs::remove_file(dir.join(format!("{program}.ini"))).unwrap();
        let mut bare = words.clone();
        bare.extend(["-conffile", "c.ini"]);
        let not_looked_for = run_in(&dir, &bare, &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")]);
        assert_eq!(
            (not_looked_for.status.code(), &not_looked_for.stdout),
            (Some(0), &plain.stdout),
            "{program}"
        );
        check_library(
            &dir,
            &bare,
            &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")],
            &not_looked_for,
        );
        let mut dotted = words.clone();
        dotted.extend(["-conffile", "./c.ini"]);
        let as_is = run_in(&dir, &dotted, &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")]);
        assert_eq!(as_is.status.code(), Some(2), "{program}: {as_is:?}");
        check_library(
            &dir,
            &dotted,
            &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")],
            &as_is,
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI BLAST+ 2.17.0 (probes `conf_missing`, `conf_dir`, `conf_empty`, `conf_log`, `conf_h`,
/// `conf_ver`, `conf_two_bad`, `makeblastdb.conf_missing`): a `-conffile` that is not a file is
/// NCBI's critical post where the diagnostics go, then its abort (exit status 255), before `-h`
/// and after `-version`.
#[test]
fn a_conffile_that_is_not_a_file_aborts_as_ncbi_aborts() {
    let dir = scratch("abort");
    fs::create_dir_all(dir.join("sub")).unwrap();
    fs::write(dir.join("c.ini"), "[BLAST]\n").unwrap();
    let critical =
        |name: &str| format!("Critical: (106.23) Registry file \"{name}\" not found, aborting.\n");
    for (program, search) in SEARCHES {
        let mut words = vec![program];
        words.extend(search);
        for (conffile, name) in [
            (&["-conffile", "missing.ini"][..], "missing.ini"),
            (&["-conffile", "sub"], "sub"),
            (&["-conffile", ""], ""),
            (&["-conffile="], ""),
            (
                &["-conffile", "c.ini", "-conffile", "missing.ini"],
                "missing.ini",
            ),
        ] {
            let mut with = words.clone();
            with.extend_from_slice(conffile);
            let aborted = run_in(&dir, &with, &[]);
            assert_eq!(
                (
                    aborted.status.code(),
                    text(&aborted.stderr),
                    aborted.stdout.len()
                ),
                (Some(255), critical(name), 0),
                "{program} {conffile:?}"
            );
            // The library's `validate` reads the arguments before a call runs the start: a
            // `-conffile` given twice is the parser's error there (`api::validate`), where the
            // command line's start has aborted before its arguments are read.
            if conffile.len() <= 2 {
                check_library(&dir, &with, &[], &aborted);
            }
        }
        let _ = fs::remove_file(dir.join("log.txt"));
        let mut logged = words.clone();
        logged.extend(["-conffile", "missing.ini", "-logfile", "log.txt"]);
        let aborted = run_in(&dir, &logged, &[]);
        assert_eq!(aborted.status.code(), Some(255), "{program}");
        assert_eq!(
            (
                text(&aborted.stderr),
                fs::read_to_string(dir.join("log.txt")).unwrap()
            ),
            (String::new(), critical("missing.ini")),
            "{program}"
        );
        // The library: the same file; its error is the post.
        fs::remove_file(dir.join("log.txt")).unwrap();
        let library = library_in(&dir, &logged, &[]);
        assert_eq!(
            (
                library.status.code(),
                text(&library.stderr),
                fs::read_to_string(dir.join("log.txt")).unwrap()
            ),
            (Some(1), critical("missing.ini"), critical("missing.ini")),
            "{program}: the library"
        );
        let help = run_in(&dir, &[program, "-h", "-conffile", "missing.ini"], &[]);
        assert_eq!(help.status.code(), Some(255), "{program}");
        let version = run_in(
            &dir,
            &[program, "-version", "-conffile", "missing.ini"],
            &[],
        );
        assert_eq!(
            (version.status.code(), text(&version.stderr)),
            (Some(0), String::new()),
            "{program}"
        );
    }
    let made = run_in(
        &dir,
        &[
            "makeblastdb",
            "-in",
            "arch_bp_s.faa",
            "-dbtype",
            "prot",
            "-out",
            "mdb",
            "-conffile",
            "missing.ini",
        ],
        &[],
    );
    assert_eq!(
        (made.status.code(), text(&made.stderr)),
        (Some(255), critical("missing.ini"))
    );
    assert!(!dir.join("mdb.pin").exists());
    // The library's makeblastdb (`api::run_makeblastdb`) runs the same start.
    let words = [
        "makeblastdb",
        "-in",
        "arch_bp_s.faa",
        "-dbtype",
        "prot",
        "-out",
        "mdb",
        "-conffile",
        "missing.ini",
    ];
    check_library(&dir, &words, &[], &made);
    assert!(!dir.join("mdb.pin").exists());
    let _ = fs::remove_dir_all(&dir);
}

/// NCBI BLAST+ 2.17.0 (probes `conf_unread_dot`, `conf_unread_bare`, `makeblastdb.log_noval`,
/// `makeblastdb.conf_noval`): a `-conffile` with a directory whose file cannot be opened ends the
/// program as a malformed registry file does (exit status 2), one without a directory loads no
/// file; makeblastdb drops a last `-logfile` or `-conffile` without a value.
#[test]
fn an_unreadable_conffile_and_a_last_word_without_value() {
    use std::os::unix::fs::PermissionsExt;
    let dir = scratch("unread");
    fs::write(dir.join("c.ini"), "[BLAST]\n").unwrap();
    fs::set_permissions(dir.join("c.ini"), fs::Permissions::from_mode(0o000)).unwrap();
    // A process that can open the file anyway (root) reads it as NCBI's would.
    let unreadable = fs::File::open(dir.join("c.ini")).is_err();
    let message = "Registry file \"./c.ini\" cannot be opened";
    let report = format!("Error: (CAppException::eNoRegistry) {message}\nError: (106.15) Application's initialization failed (CAppException::eNoRegistry) {message}\n");
    for (program, search) in SEARCHES {
        let mut words = vec![program];
        words.extend(search);
        words.extend(["-outfmt", "6"]);
        let mut dotted = words.clone();
        dotted.extend(["-conffile", "./c.ini"]);
        let run = run_in(&dir, &dotted, &[]);
        if unreadable {
            assert_eq!(
                (run.status.code(), text(&run.stderr)),
                (Some(2), report.repeat(2)),
                "{program}"
            );
        }
        check_library(&dir, &dotted, &[], &run);
        let mut bare = words.clone();
        bare.extend(["-conffile", "c.ini"]);
        let run = run_in(&dir, &bare, &[]);
        assert_eq!(run.status.code(), Some(0), "{program}: {run:?}");
        check_library(&dir, &bare, &[], &run);
    }
    fs::set_permissions(dir.join("c.ini"), fs::Permissions::from_mode(0o644)).unwrap();
    for last in ["-logfile", "-conffile"] {
        let run = run_in(
            &dir,
            &[
                "makeblastdb",
                "-in",
                "arch_bp_s.faa",
                "-dbtype",
                "prot",
                "-out",
                "mdb",
                last,
            ],
            &[],
        );
        assert_eq!(run.status.code(), Some(0), "{last}: {run:?}");
        assert!(dir.join("mdb.pin").exists(), "{last}");
        fs::remove_file(dir.join("mdb.pin")).unwrap();
        let library = library_in(
            &dir,
            &[
                "makeblastdb",
                "-in",
                "arch_bp_s.faa",
                "-dbtype",
                "prot",
                "-out",
                "mdb",
                last,
            ],
            &[],
        );
        // The standard output names the time and the build's duration.
        assert_eq!(
            (library.status.code(), &library.stderr),
            (Some(0), &run.stderr),
            "{last}: the library"
        );
        assert!(dir.join("mdb.pin").exists(), "{last}: the library");
        fs::remove_file(dir.join("mdb.pin")).unwrap();
    }
    let _ = fs::remove_dir_all(&dir);
}

/// No NCBI counterpart: LOSAT's own refusal (here of `BL2SEQ_LEGACY`, which NCBI reads and LOSAT
/// does not) stays on the terminal when tblastx sends the messages of its run to the `-logfile`
/// file (line RG DECISIONS 32).
#[test]
fn losats_refusals_stay_on_the_terminal() {
    let dir = scratch("refusal");
    let (_, search) = SEARCHES[3];
    let mut words = vec!["tblastx"];
    words.extend(search);
    words.extend(["-outfmt", "6", "-logfile", "log.txt"]);
    let refused = run_in(&dir, &words, &[("BL2SEQ_LEGACY", "1")]);
    assert_eq!(refused.status.code(), Some(1), "{refused:?}");
    let stderr = text(&refused.stderr);
    assert!(
        stderr.starts_with("Error: ") && stderr.contains("not supported by LOSAT"),
        "{stderr}"
    );
    assert!(!fs::read_to_string(dir.join("log.txt"))
        .unwrap()
        .contains("not supported by LOSAT"));
    let _ = fs::remove_dir_all(&dir);
}

/// No NCBI counterpart: LOSAT refuses the settings that change NCBI's output (here
/// `BL2SEQ_LEGACY` and `DIAG_POST_LEVEL` in the environment, `[BLAST] LONG_SEQID` in the
/// `-conffile` file) after the start, before the run, dry run or not. A library call with
/// `-logfile`, `-conffile` or `-dryrun` runs the start at its beginning and refuses them alike:
/// its error is the command line's text after `Error: `, and its `-logfile` file is the command
/// line's (line RG item E1w-2).
#[test]
fn the_librarys_start_refuses_what_the_command_line_refuses() {
    let dir = scratch("librefusal");
    fs::write(dir.join("c.ini"), "[BLAST]\nLONG_SEQID=1\n").unwrap();
    let log = dir.join("log.txt");
    for (program, search) in SEARCHES {
        for (extra, envs) in [
            (&["-logfile", "log.txt"][..], &[("BL2SEQ_LEGACY", "1")][..]),
            (&["-logfile", "log.txt"], &[("DIAG_POST_LEVEL", "Info")]),
            (&["-conffile", "c.ini"], &[]),
        ] {
            for dry in [&[][..], &["-dryrun"]] {
                let mut words = vec![program];
                words.extend(search);
                words.extend(["-outfmt", "6"]);
                words.extend_from_slice(extra);
                words.extend_from_slice(dry);
                let _ = fs::remove_file(&log);
                let cli = run_in(&dir, &words, envs);
                let cli_log = fs::read(&log).ok();
                let refused = cli.status.code() == Some(1);
                // Without `-dryrun` each is refused; a dry run refuses the settings that the start
                // reads (`DIAG_POST_LEVEL`, `LONG_SEQID`), not `BL2SEQ_LEGACY`, which the run
                // reads (line RG DECISIONS E1b-6).
                assert!(refused || !dry.is_empty(), "{words:?} {envs:?}: {cli:?}");
                assert!(cli.stdout.is_empty(), "{words:?} {envs:?}: {cli:?}");
                assert!(
                    !refused || text(&cli.stderr).starts_with("Error: "),
                    "{cli:?}"
                );
                let _ = fs::remove_file(&log);
                let library = library_in(&dir, &words, envs);
                let library_stderr = if refused {
                    format!("Error: {}", text(&library.stderr))
                } else {
                    text(&library.stderr)
                };
                assert_eq!(
                    (
                        library.status.code(),
                        library_stderr,
                        text(&library.stdout),
                        fs::read(&log).ok()
                    ),
                    (cli.status.code(), text(&cli.stderr), String::new(), cli_log),
                    "{words:?} {envs:?}: the library"
                );
            }
        }
    }
    let _ = fs::remove_dir_all(&dir);
}
