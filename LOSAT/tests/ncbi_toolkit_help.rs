//! The NCBI C++ Toolkit option `-help-full` of every program (line RG item E2), against NCBI BLAST+
//! 2.17.0, the official Linux x64 release and bioconda's build for rpstblastn, dustmasker and
//! segmasker, and NCBI's rpsbproc 0.5.1 (probes of 2026-11-13, `remaining-gaps-20261112/E/SURVEY_E2.md`,
//! run with `env -i` in a fresh directory): `-help-full` writes the bytes of `-help` in every place
//! but makeblastdb's, whose hidden `-randomize` LOSAT does not run, and the reader stops at a help
//! flag. LOSAT's help texts are its own (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 1), so the
//! test holds `-help-full` to LOSAT's `-help` in the same place.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3020-3023, 2840-2852 (quoted at
//! `cli::try_parse_from`) and c++/src/corelib/ncbiapp.cpp:685-697
//! ```c++
//!     catch (const CArgHelpException& e) {
//!         x_AddDefaultArgs();
//!         // Print USAGE
//!         if (e.GetErrCode() == CArgHelpException::eHelpXml) {
//!             m_ArgDesc->PrintUsageXml(cout);
//!         } else {
//!             CArgHelpException::TErrCode code =  e.GetErrCode();
//!             string str;
//!             m_ArgDesc->ShowAllArguments(code == CArgHelpException::eHelpShowAll)->PrintUsage
//!                 (str, code == CArgHelpException::eHelpFull || code == CArgHelpException::eHelpShowAll);
//!             cout << str;
//!         }
//!         *exit_code = e.GetErrCode() == CArgHelpException::eHelpErr ? 2 : 0;
//!     }
//! ```

#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn scratch(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-toolkit-help-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// `losat <program> <words>` in `dir` with an environment of NCBI's probes (`env -i PATH HOME`).
fn run_in(dir: &Path, program: &str, words: &[&str]) -> Output {
    let mut command = Command::new(BIN);
    command
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .env("BLAST_USAGE_REPORT", "false")
        .current_dir(dir)
        .arg(program)
        .args(words)
        .stdin(std::process::Stdio::null());
    // The `--all-features` build's BLASTX diagnostic probe needs its log (as in CI).
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        command.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    command.output().unwrap()
}

/// The programs and an option of each that takes a value.
const PROGRAMS: [(&str, &str); 13] = [
    ("blastn", "-out"),
    ("blastp", "-out"),
    ("blastx", "-out"),
    ("tblastn", "-out"),
    ("tblastx", "-out"),
    ("psiblast", "-out"),
    ("rpsblast", "-out"),
    ("rpstblastn", "-out"),
    ("makeblastdb", "-title"),
    ("rpsbproc", "-o"),
    ("blastdbcmd", "-out"),
    ("dustmasker", "-out"),
    ("segmasker", "-out"),
];

#[test]
fn help_full_writes_the_help_in_every_place() {
    let dir = scratch("full");
    for (program, value_option) in PROGRAMS {
        let help = run_in(&dir, program, &["-help"]);
        assert_eq!(help.status.code(), Some(0), "{program}");
        assert!(help.stdout.starts_with(b"USAGE"), "{program}");
        assert!(help.stderr.is_empty(), "{program}");
        // NCBI BLAST+ 2.17.0 probes `hf`, `hf_h`, `hf_unk`, `hf_missing_val2`, `hf_dry`, `hf_xml`,
        // `hf_twice` (`-help-full` and its `-help` twin: the same bytes).
        for words in [
            &["-help-full"][..],
            &["-help-full", "-h"],
            &["-help-full", "-bogus"],
            &["-help", "-bogus"],
            &["-help-full", value_option],
            &["-help", value_option],
            &["-dryrun", "-help-full"],
            &["-help-full", "-xmlhelp"],
            &["-help-full", "-help-full"],
        ] {
            let full = run_in(&dir, program, words);
            assert_eq!(full.status.code(), Some(0), "{program} {words:?}");
            assert_eq!(full.stdout, help.stdout, "{program} {words:?}");
            assert_eq!(full.stderr, help.stderr, "{program} {words:?}");
        }
        // `-h -help-full`, `-h -xmlhelp`: the short help (probes `h_hf`, `h_xml`).
        let short = run_in(&dir, program, &["-h"]);
        for words in [&["-h", "-help-full"][..], &["-h", "-xmlhelp"]] {
            let output = run_in(&dir, program, words);
            assert_eq!(output.status.code(), Some(0), "{program} {words:?}");
            assert_eq!(output.stdout, short.stdout, "{program} {words:?}");
        }
        // The words read before the arguments come first (probes `hf_ver`, `hf_conf_missing`,
        // `hf_log`): the version, the missing registry file (exit status 255), the log file made.
        let version = run_in(&dir, program, &["-help-full", "-version"]);
        assert_eq!(
            version.stdout,
            run_in(&dir, program, &["-version"]).stdout,
            "{program}"
        );
        let conf = run_in(&dir, program, &["-help-full", "-conffile", "missing.ini"]);
        let conf_help = run_in(&dir, program, &["-help", "-conffile", "missing.ini"]);
        assert_eq!(conf.status.code(), conf_help.status.code(), "{program}");
        assert_eq!(conf.stdout, conf_help.stdout, "{program}");
        assert_eq!(conf.stderr, conf_help.stderr, "{program}");
        let log = dir.join("log.txt");
        let _ = fs::remove_file(&log);
        let logged = run_in(&dir, program, &["-help-full", "-logfile", "log.txt"]);
        assert_eq!(logged.stdout, help.stdout, "{program}");
        assert_eq!(fs::read(&log).ok(), Some(Vec::new()), "{program}");
        fs::remove_file(&log).unwrap();
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn xmlhelp_is_refused_by_name_in_every_program() {
    let dir = scratch("xml");
    for (program, _) in PROGRAMS {
        let output = run_in(&dir, program, &["-xmlhelp"]);
        assert_ne!(output.status.code(), Some(0), "{program}");
        assert!(output.stdout.is_empty(), "{program}");
        let text = String::from_utf8_lossy(&output.stderr);
        assert!(
            text.contains("option -xmlhelp is not supported by LOSAT's"),
            "{program}: {text}"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}
