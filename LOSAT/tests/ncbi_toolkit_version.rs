//! The NCBI C++ Toolkit options `-version-full`, `-version-full-xml` and `-version-full-json` of
//! every program (line RG item E3), against NCBI BLAST+ 2.17.0, the official Linux x64 release and
//! bioconda's build for rpstblastn, dustmasker and segmasker, and NCBI's rpsbproc 0.5.1 (probes of
//! 2026-11-13, `remaining-gaps-20261112/E/SURVEY_E3.md`, run with `env -i` in a fresh directory):
//! the application prints the first version word before `--` (in an option's value too, before
//! the configuration is loaded and the arguments are read) and exits 0. The forms are NCBI's, the
//! values LOSAT's (`ncbi_version::print`, RG DECISIONS E3-1): NCBI's first line (`blastn: 2.17.0+`,
//! the same in both builds), LOSAT's package, and none of the lines that NCBI's build writes about
//! itself (they differ between NCBI's two builds).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:980-985 (quoted at `ncbi_version::print`)
//! ```c++
//!             } else if ( NStr::strcmp(argv[i], s_ArgFullVersion) == 0 ) {
//!                 delete[] v;
//!                 // Print full VERSION
//!                 cout << GetFullVersion().Print( appname );
//!                 diag_context.DiscardMessages();
//!                 return 0;
//! ```

#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-toolkit-version-{name}-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// `losat <program> <words>` in `dir` with an environment of NCBI's probes (`env -i PATH HOME`),
/// its standard output on `stdout`.
fn run_to(dir: &Path, program: &str, words: &[&str], stdout: Stdio) -> Output {
    let mut command = Command::new(BIN);
    command
        .env_clear()
        .env("PATH", "/usr/bin:/bin")
        .env("HOME", "/nonexistent")
        .env("BLAST_USAGE_REPORT", "false")
        .current_dir(dir)
        .arg(program)
        .args(words)
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(Stdio::piped());
    // The `--all-features` build's BLASTX diagnostic probe needs its log (as in CI).
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        command.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    command.spawn().unwrap().wait_with_output().unwrap()
}

fn run_in(dir: &Path, program: &str, words: &[&str]) -> Output {
    run_to(dir, program, words, Stdio::piped())
}

/// The programs, NCBI's version of each (the first line of its `-version-full`, probes `vf`, and
/// the numbers of its XML and JSON `version_info`), and an option of each that takes a value.
const PROGRAMS: [(&str, &str, [&str; 3], &str); 13] = [
    ("blastn", "2.17.0+", ["2", "17", "0"], "-out"),
    ("blastp", "2.17.0+", ["2", "17", "0"], "-out"),
    ("blastx", "2.17.0+", ["2", "17", "0"], "-out"),
    ("tblastn", "2.17.0+", ["2", "17", "0"], "-out"),
    ("tblastx", "2.17.0+", ["2", "17", "0"], "-out"),
    ("psiblast", "2.17.0+", ["2", "17", "0"], "-out"),
    ("rpsblast", "2.17.0+", ["2", "17", "0"], "-out"),
    ("rpstblastn", "2.17.0+", ["2", "17", "0"], "-out"),
    ("makeblastdb", "2.17.0+", ["2", "17", "0"], "-title"),
    ("rpsbproc", "0.5.1", ["0", "5", "1"], "-o"),
    ("blastdbcmd", "2.17.0+", ["2", "17", "0"], "-out"),
    ("dustmasker", "1.0.0", ["1", "0", "0"], "-out"),
    ("segmasker", "1.0.0", ["1", "0", "0"], "-out"),
];

/// LOSAT's package: its version, and its numbers in the XML and JSON.
fn losat_numbers() -> (String, String) {
    let (major, minor, patch, pre) = (
        env!("CARGO_PKG_VERSION_MAJOR"),
        env!("CARGO_PKG_VERSION_MINOR"),
        env!("CARGO_PKG_VERSION_PATCH"),
        env!("CARGO_PKG_VERSION_PRE"),
    );
    let (xml_name, json_name) = if pre.is_empty() {
        (String::new(), String::new())
    } else {
        (format!(" name=\"{pre}\""), format!(", \"name\": \"{pre}\""))
    };
    (
        format!("major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"{xml_name}"),
        format!("\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}{json_name}"),
    )
}

/// What LOSAT's `word` prints for `program` (RG DECISIONS E3-1).
fn expected(program: &str, version: &str, numbers: [&str; 3], word: &str) -> String {
    let package = format!(
        " Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    let [major, minor, patch] = numbers;
    let (losat_xml, losat_json) = losat_numbers();
    match word {
        "-version" => format!("{program}: {version}\n{package}"),
        "-version-full" => {
            format!("{program}: {version}\n{package} Package-Config:  \n GI-64bit:  TRUE\n")
        }
        "-version-full-xml" => format!(
            "<?xml version=\"1.0\"?>\n<ncbi_version xmlns=\"ncbi:version\"\n  \
             xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n  \
             xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">\n<appname>{program}</appname>\n\
             <version_info major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"/>\n\
             <package name=\"LOSAT\">\n<version_info {losat_xml}/>\n<build_info>\n</build_info>\n\
             <config></config>\n</package>\n<build_info>\n</build_info>\n</ncbi_version>\n"
        ),
        "-version-full-json" => format!(
            "{{\n  \"ncbi_version\": {{\n    \"appname\": \"{program}\",\n    \
             \"version_info\": {{\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}}},\n    \
             \"component\": [],\n    \"package\": {{\n      \"name\": \"LOSAT\",\n      \
             \"version_info\": {{{losat_json}}},\n      \"build_info\": {{}},\n      \
             \"config\": \"\"}},\n    \"build_info\": {{}}\n  }}\n}}\n"
        ),
        _ => unreachable!("{word}"),
    }
}

#[test]
fn the_full_versions_are_ncbi_s_forms_with_losat_s_values() {
    let dir = scratch("forms");
    for (program, version, numbers, value_option) in PROGRAMS {
        for word in ["-version-full", "-version-full-xml", "-version-full-json"] {
            let text = expected(program, version, numbers, word);
            // NCBI BLAST+ 2.17.0 probes `vf`, `vf_args`, `vf_first`, `vf_val`, `vf_val_args`,
            // `unk_vf`, `vf_unk`, `argerr_vf`, `h_vf`, `vf_h`, `help_vf`, `hf_vf`, `xml_vf`,
            // `vf_dry`, `dry_vf_after`, `vf_twice`, `vf_conf_missing`, `conf_missing_vf`
            // (and their XML and JSON twins): the version word wins wherever it is before `--`.
            for words in [
                &[word][..],
                &[value_option, word],
                &["-bogus", word],
                &[word, "-bogus"],
                &["-h", word],
                &[word, "-h"],
                &["-help", word],
                &["-help-full", word],
                &["-xmlhelp", word],
                &["-dryrun", word],
                &[word, "-dryrun"],
                &[word, word],
                &[word, "-conffile", "missing.ini"],
                &["-conffile", "missing.ini", word],
                &["-conffile=missing.ini", word],
            ] {
                let output = run_in(&dir, program, words);
                assert_eq!(output.status.code(), Some(0), "{program} {words:?}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    text,
                    "{program} {words:?}"
                );
                assert!(output.stderr.is_empty(), "{program} {words:?}");
            }
        }
        // The first version word is printed (probes `ver_vf`, `vf_ver`, `vf_vx`, `vj_vf`).
        for (words, word) in [
            (&["-version", "-version-full"][..], "-version"),
            (&["-version-full", "-version"], "-version-full"),
            (&["-version-full", "-version-full-xml"], "-version-full"),
            (
                &["-version-full-json", "-version-full"],
                "-version-full-json",
            ),
        ] {
            let output = run_in(&dir, program, words);
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                expected(program, version, numbers, word),
                "{program} {words:?}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn the_full_versions_come_before_the_configuration() {
    let dir = scratch("start");
    for (program, version, numbers, _) in PROGRAMS {
        let text = expected(program, version, numbers, "-version-full");
        // A malformed `<program>.ini` and `.ncbirc` are not read (probes `vf_ini_bad`,
        // `vf_ncbirc_bad`); a `-logfile` file is made, empty (`vf_log`, `log_vf`).
        fs::write(dir.join(format!("{program}.ini")), "[x\n").unwrap();
        fs::write(dir.join(".ncbirc"), "[x\n").unwrap();
        let output = run_in(&dir, program, &["-version-full"]);
        assert_eq!(output.status.code(), Some(0), "{program}");
        assert_eq!(String::from_utf8_lossy(&output.stdout), text, "{program}");
        assert!(output.stderr.is_empty(), "{program}");
        fs::remove_file(dir.join(format!("{program}.ini"))).unwrap();
        fs::remove_file(dir.join(".ncbirc")).unwrap();
        for words in [
            &["-version-full", "-logfile", "log.txt"][..],
            &["-logfile", "log.txt", "-version-full"],
        ] {
            let output = run_in(&dir, program, words);
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                text,
                "{program} {words:?}"
            );
            assert_eq!(
                fs::read(dir.join("log.txt")).ok(),
                Some(Vec::new()),
                "{program} {words:?}"
            );
            fs::remove_file(dir.join("log.txt")).unwrap();
        }
        // Not a version word: after `--` (an argument error), a longer word, the value of
        // `-conffile` (NCBI's missing registry file, exit status 255; probes `dd_vf`,
        // `vf_eq`, `vf_suffix`, `vf_confval`).
        for words in [
            &["--", "-version-full"][..],
            &["-version-full=1"],
            &["-version-fullx"],
            &["-conffile", "-version-full"],
        ] {
            let output = run_in(&dir, program, words);
            assert_ne!(output.status.code(), Some(0), "{program} {words:?}");
            assert!(
                !output.stdout.starts_with(format!("{program}: ").as_bytes()),
                "{program} {words:?}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

// NCBI writes the version to `cout` and never checks it: a full disk drops it, exit status 0
// (probes `vf_devfull`, `vx_devfull`, `ver_devfull`). rpsbproc's failed write of its version is
// LOSAT's exception 9 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES): exit status 4 with its message.
// `/dev/full` is Linux's.
#[cfg(target_os = "linux")]
#[test]
fn a_full_disk_drops_the_version_as_ncbi_does() {
    let dir = scratch("full");
    for (program, ..) in PROGRAMS {
        for word in ["-version", "-version-full", "-version-full-xml"] {
            let full = fs::OpenOptions::new()
                .write(true)
                .open("/dev/full")
                .unwrap();
            let output = run_to(&dir, program, &[word], Stdio::from(full));
            if program == "rpsbproc" {
                assert_eq!(output.status.code(), Some(4), "{program} {word}");
                assert!(
                    output
                        .stderr
                        .starts_with(b"rpsbproc failed to write output: "),
                    "{program} {word}"
                );
            } else {
                assert_eq!(output.status.code(), Some(0), "{program} {word}");
                assert!(
                    output.stderr.is_empty(),
                    "{program} {word}: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
            }
        }
    }
    let _ = fs::remove_dir_all(&dir);
}
