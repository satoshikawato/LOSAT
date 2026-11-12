//! The NCBI C++ Toolkit's file options `-logfile` and `-conffile` of the BLAST programs (line RG
//! item E1), against NCBI BLAST+ 2.17.0, the official Linux x64 release (probes of 2026-11-12,
//! `remaining-gaps-20261112/E/SURVEY_E1.md`, run with `env -i` in a fresh directory).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-871, 926-1022 (quoted at
//! `ncbi_environment::log_file`, `preparse` and `conf_file_path`), and
//! c++/src/app/blast/tblastx_app.cpp:93-95 (quoted at `ncbi_environment::run_messages_follow_log`).
#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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

/// NCBI BLAST+ 2.17.0 (probes `log_warn`, `log_qerr`, `log_empty`, `log_nodir`, `log_old`): the
/// warning of `-max_target_seqs 1` and the error of a missing `-query` go to the `-logfile` file
/// in tblastx and rpsblast and to the standard error in the other programs; the file is made
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
        // `-logfile ""` sends the diagnostics nowhere; a file in a missing directory leaves them
        // on the standard error.
        for (value, quiet) in [("", to_log), ("nodir/log.txt", false)] {
            let mut words = words.clone();
            words.extend(["-logfile", value]);
            let run = run_in(&dir, &words, &[]);
            assert_eq!(run.status.code(), Some(0), "{program} {value:?}");
            assert_eq!(run.stdout, plain.stdout, "{program} {value:?}");
            let expected = if quiet { "" } else { warning.as_str() };
            assert_eq!(text(&run.stderr), expected, "{program} {value:?}");
        }
        assert!(!dir.join("nodir").exists());
    }
    // A missing `-query`, rpsblast too (its error comes before the database is opened).
    let error =
        "Command line argument error: Argument \"query\". File is not accessible:  `missing.fa'\n";
    for program in [
        "blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast", "rpsblast",
    ] {
        let _ = fs::remove_file(dir.join("log.txt"));
        let other = if program == "rpsblast" {
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
        if matches!(program, "tblastx" | "rpsblast") {
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
            fs::write(dir.join("c.ini"), "[x\n").unwrap();
            let failed = run_in(&dir, &with, &[]);
            assert_eq!(
                (failed.status.code(), text(&failed.stderr)),
                (Some(2), SECTION_REPORT.repeat(2)),
                "{program} {conffile:?}"
            );
            assert!(failed.stdout.is_empty());
        }
        fs::remove_file(dir.join(format!("{program}.ini"))).unwrap();
        let mut bare = words.clone();
        bare.extend(["-conffile", "c.ini"]);
        let not_looked_for = run_in(&dir, &bare, &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")]);
        assert_eq!(
            (not_looked_for.status.code(), &not_looked_for.stdout),
            (Some(0), &plain.stdout),
            "{program}"
        );
        let mut dotted = words.clone();
        dotted.extend(["-conffile", "./c.ini"]);
        let as_is = run_in(&dir, &dotted, &[("NCBI_DONT_USE_LOCAL_CONFIG", "1")]);
        assert_eq!(as_is.status.code(), Some(2), "{program}: {as_is:?}");
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
        let mut bare = words.clone();
        bare.extend(["-conffile", "c.ini"]);
        let run = run_in(&dir, &bare, &[]);
        assert_eq!(run.status.code(), Some(0), "{program}: {run:?}");
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
