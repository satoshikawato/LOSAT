//! The NCBI C++ Toolkit options `-version-full`, `-version-full-xml` and `-version-full-json` of
//! every program (line RG item E3), against NCBI BLAST+ 2.17.0, the official Linux x64 release and
//! bioconda's build for rpstblastn, dustmasker and segmasker, and NCBI's rpsbproc 0.5.1 (probes of
//! 2026-11-13, `remaining-gaps-20261112/E/SURVEY_E3.md`, run with `env -i` in a fresh directory;
//! convert2blastmask, which line RG item K3 added during E3, on the official release's binary):
//! the application prints the first version word before `--` (in an option's value too, before
//! the configuration is loaded and the arguments are read) and exits 0. The forms are NCBI's, the
//! values LOSAT's (`ncbi_version::print`, RG DECISIONS E3-1; the Owner, 2026-11-21): LOSAT's first line
//! (`blastn: LOSAT 0.2.0`; NCBI's is `blastn: 2.17.0+`, the same in both builds), LOSAT's package, and none of the lines that NCBI's build writes about
//! itself (they differ between NCBI's two builds).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:980-985 (cited at `ncbi_version::print`)

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

/// The programs (the first line of the `-version-full` of each is `<program>: LOSAT <version>`:
/// NCBI's is `<program>: 2.17.0+` (rpsbproc
/// `0.5.1`, the maskers `1.0.0`, makembindex `2.17.0`), probes `vf`, and its XML and JSON
/// `version_info` have NCBI's numbers; LOSAT's have LOSAT's), and an option of each that takes a value.
const PROGRAMS: [(&str, &str); 17] = [
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
    ("windowmasker", "-out"),
    ("convert2blastmask", "-out"),
    ("makembindex", "-input"),
    ("blast_formatter", "-out"),
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

/// The components of `program`'s full versions: windowmasker's six (NCBI BLAST+ 2.17.0
/// `windowmasker -version-full*`, probe b1 of line NP M5), each with LOSAT's empty build
/// information; the text, the XML and the JSON list.
fn components(program: &str) -> (String, String, String) {
    if program != "windowmasker" {
        return (String::new(), String::new(), "[]".to_string());
    }
    let list = [
        ("window-masker-algorithm", "", "0"),
        ("windowmasker-statistics-algorithm-version", "", "0"),
        ("windowmasker-statistics-format-version", "ascii ", "1"),
        ("windowmasker-statistics-format-version", "binary ", "1"),
        ("windowmasker-statistics-format-version", "oascii ", "1"),
        ("windowmasker-statistics-format-version", "obinary ", "1"),
    ];
    let text = list
        .iter()
        .map(|(name, prefix, minor)| format!("\n {name}:{prefix}1.{minor}.0\n"))
        .collect();
    let xml = list
        .iter()
        .map(|(name, _, minor)| {
            format!(
                "<component name=\"{name}\">\n\
                 <version_info major=\"1\" minor=\"{minor}\" patch_level=\"0\"/>\n\n\
                 <build_info>\n</build_info>\n</component>\n"
            )
        })
        .collect();
    let json: Vec<String> = list
        .iter()
        .map(|(name, _, minor)| {
            format!(
                "\n      {{ \"name\": \"{name}\", \"version_info\": \
                 {{\"major\": 1, \"minor\": {minor}, \"patch_level\": 0}},\n        \
                 \"build_info\": {{}}}}"
            )
        })
        .collect();
    (text, xml, format!("[{}]", json.join(",")))
}

/// What LOSAT's `word` prints for `program` (RG DECISIONS E3-1).
fn expected(program: &str, word: &str) -> String {
    let version = concat!("LOSAT ", env!("CARGO_PKG_VERSION"));
    let package = format!(
        " Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    let (major, minor, patch) = (
        env!("CARGO_PKG_VERSION_MAJOR"),
        env!("CARGO_PKG_VERSION_MINOR"),
        env!("CARGO_PKG_VERSION_PATCH"),
    );
    let (losat_xml, losat_json) = losat_numbers();
    let (component_text, component_xml, component_json) = components(program);
    match word {
        "-version" => format!("{program}: {version}\n{package}"),
        "-version-full" => {
            format!(
                "{program}: {version}\n{package} Package-Config:  \n GI-64bit:  TRUE\n{component_text}"
            )
        }
        "-version-full-xml" => format!(
            "<?xml version=\"1.0\"?>\n<ncbi_version xmlns=\"ncbi:version\"\n  \
             xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n  \
             xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">\n<appname>{program}</appname>\n\
             <version_info major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"/>\n\
             {component_xml}<package name=\"LOSAT\">\n<version_info {losat_xml}/>\n<build_info>\n</build_info>\n\
             <config></config>\n</package>\n<build_info>\n</build_info>\n</ncbi_version>\n"
        ),
        "-version-full-json" => format!(
            "{{\n  \"ncbi_version\": {{\n    \"appname\": \"{program}\",\n    \
             \"version_info\": {{\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}}},\n    \
             \"component\": {component_json},\n    \"package\": {{\n      \"name\": \"LOSAT\",\n      \
             \"version_info\": {{{losat_json}}},\n      \"build_info\": {{}},\n      \
             \"config\": \"\"}},\n    \"build_info\": {{}}\n  }}\n}}\n"
        ),
        _ => unreachable!("{word}"),
    }
}

#[test]
fn the_full_versions_are_ncbi_s_forms_with_losat_s_values() {
    let dir = scratch("forms");
    for (program, value_option) in PROGRAMS {
        for word in ["-version-full", "-version-full-xml", "-version-full-json"] {
            let text = expected(program, word);
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
                expected(program, word),
                "{program} {words:?}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn the_full_versions_come_before_the_configuration() {
    let dir = scratch("start");
    for (program, _) in PROGRAMS {
        let text = expected(program, "-version-full");
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

/// The library's tools (`api::run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`,
/// `run_segmasker`, `run_windowmasker`, `run_convert2blastmask`, `run_makembindex`, `run_blast_formatter`) run NCBI's
/// application start as the command line does: they print the versions to their output (line RG item E1w-2 for makeblastdb). The search
/// programs of the library return reports, and refuse the words (`cli::version_word_refusal`), which
/// name the Web engine's ABI v2 among the places that print them where it runs the program (line
/// WQA's D14, 2026-11-20).
/// The library runs are `examples/library_command_line`, in the command line's environment.
#[test]
fn the_library_s_tools_print_the_versions_as_the_command_line() {
    let dir = scratch("library");
    let example = crate::library_examples::example_binary("library_command_line");
    let library = |program: &str, words: &[&str]| {
        let mut command = Command::new(&example);
        command
            .env_clear()
            .env("PATH", "/usr/bin:/bin")
            .env("HOME", "/nonexistent")
            .env("BLAST_USAGE_REPORT", "false")
            .current_dir(&dir)
            .arg(program)
            .args(words)
            .stdin(Stdio::null());
        command.output().unwrap()
    };
    for (program, value_option) in PROGRAMS {
        let tool = matches!(
            program,
            "makeblastdb"
                | "blastdbcmd"
                | "dustmasker"
                | "segmasker"
                | "windowmasker"
                | "convert2blastmask"
                | "makembindex"
                | "blast_formatter"
        );
        for word in [
            "-version",
            "-version-full",
            "-version-full-xml",
            "-version-full-json",
        ] {
            let text = expected(program, word);
            for words in [
                &[word][..],
                &[value_option, word],
                &["-dryrun", word],
                &[word, "-conffile", "missing.ini"],
                &["-logfile", "log.txt", word],
            ] {
                // The example takes the search programs' `-out` itself.
                if !tool && words[0] == value_option {
                    continue;
                }
                let output = library(program, words);
                if tool {
                    assert_eq!(
                        (
                            output.status.code(),
                            String::from_utf8_lossy(&output.stdout).into_owned(),
                            output.stderr.len()
                        ),
                        (Some(0), text.clone(), 0),
                        "{program} {words:?}"
                    );
                } else if matches!(program, "rpsblast" | "rpstblastn" | "rpsbproc") {
                    // Not in the library.
                } else {
                    let stderr = String::from_utf8_lossy(&output.stderr);
                    assert_eq!(output.status.code(), Some(1), "{program} {words:?}");
                    let upper = program.to_ascii_uppercase();
                    let refusal = match matches!(program, "blastn" | "blastp" | "blastx" | "tblastn" | "tblastx") {
                        true => format!("the option {word} runs on the command line and in LOSAT Web of LOSAT's {upper} only: the library API does not write NCBI's version text"),
                        false => format!("the option {word} runs on the command line of LOSAT's {upper} only: the library API and the Web engine do not write NCBI's version text"),
                    };
                    assert!(stderr.contains(&refusal), "{program} {words:?}: {stderr}");
                }
            }
        }
        let _ = fs::remove_file(dir.join("log.txt"));
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
