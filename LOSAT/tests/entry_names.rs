//! `losat --install-shims` (`cli::install_shims`) and the executables it links, called by a
//! program's name (`cli::entry_argv`; its rule is checked in `src/utils/entry_names_tests.rs`).
//! No NCBI counterpart: NCBI has one executable per program, and these
//! names stand for them on `PATH`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn losat() -> &'static str {
    env!("CARGO_BIN_EXE_losat")
}

/// A new directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    // In the target directory: the hard link of `the_default_directory_...` needs the binary's
    // file system.
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("losat-entry-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// One run with no `LOSAT_*` variable that chooses a code path, in LOSAT/.
fn run(program: impl AsRef<std::ffi::OsStr>, words: &[&str]) -> Output {
    let mut command = Command::new(program);
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(words)
        .output()
        .unwrap()
}

fn shim(dir: &Path, name: &str) -> PathBuf {
    dir.join(format!("{name}{}", std::env::consts::EXE_SUFFIX))
}

fn install(dir: &Path) -> Output {
    run(losat(), &["--install-shims", dir.to_str().unwrap()])
}

/// The programs that the binary links, in the order it lists them.
const PROGRAMS: [&str; 11] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "makeblastdb",
    "blastdbcmd",
    "dustmasker",
    "rpsbproc",
];

/// The program names that `losat --install-shims` listed on standard output (one path per
/// line): the file names without the executable suffix. The list must hold `PROGRAMS`.
fn installed_names(installed: &Output) -> Vec<String> {
    let names: Vec<String> = String::from_utf8(installed.stdout.clone())
        .unwrap()
        .lines()
        .map(|line| {
            let file = Path::new(line).file_name().unwrap().to_str().unwrap();
            file.strip_suffix(std::env::consts::EXE_SUFFIX)
                .unwrap_or(file)
                .to_string()
        })
        .collect();
    assert!(!names.is_empty(), "{installed:?}");
    for program in PROGRAMS {
        assert!(
            names.iter().any(|name| name == program),
            "{program}: {names:?}"
        );
    }
    names
}

#[test]
fn a_shim_prints_what_the_subcommand_prints() {
    let dir = temp_dir("texts");
    let installed = install(&dir);
    assert!(installed.status.success(), "{installed:?}");
    for name in installed_names(&installed).iter().map(String::as_str) {
        // `blastn -version` is what tools read the version from (`blastn: 2.17.0+`; rpsbproc,
        // which is not part of BLAST+, and dustmasker have their own: `rpsbproc: 0.5.1`,
        // `dustmasker: 1.0.0`).
        let version = run(shim(&dir, name), &["-version"]);
        assert!(version.status.success(), "{name}: {version:?}");
        let first_line = match name {
            "rpsbproc" => "rpsbproc: 0.5.1\n".to_string(),
            "dustmasker" => "dustmasker: 1.0.0\n".to_string(),
            _ => format!("{name}: 2.17.0+\n"),
        };
        assert!(
            version.stdout.starts_with(first_line.as_bytes()),
            "{name}: {version:?}"
        );
        for words in [
            &["-version"][..],
            &["-h"],
            &["-help"],
            &[],
            &["-no_such_option"],
        ] {
            let mut through_losat = vec![name];
            through_losat.extend(words);
            let expected = run(losat(), &through_losat);
            let output = run(shim(&dir, name), words);
            assert_eq!(
                output.status.code(),
                expected.status.code(),
                "{name} {words:?}"
            );
            assert_eq!(output.stdout, expected.stdout, "{name} {words:?}");
            assert_eq!(output.stderr, expected.stderr, "{name} {words:?}");
        }
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_shim_searches_and_builds_as_the_subcommand_does() {
    let dir = temp_dir("runs");
    assert!(install(&dir).status.success());
    let db = dir.join("db");
    let db = db.to_str().unwrap();
    let fasta = "tests/db_search/taxfilter/q.fna";
    let built = run(
        shim(&dir, "makeblastdb"),
        &["-in", fasta, "-dbtype", "nucl", "-out", db],
    );
    assert!(built.status.success(), "{built:?}");
    let searches: [(&str, &[&str]); 3] = [
        ("blastn", &["-query", fasta, "-db", db, "-outfmt", "6"]),
        (
            "tblastx",
            &[
                "-query",
                fasta,
                "-subject",
                fasta,
                "-outfmt",
                "6",
                "-max_target_seqs",
                "1",
            ],
        ),
        ("blastdbcmd", &["-db", db, "-info"]),
    ];
    for (name, words) in searches {
        let mut through_losat = vec![name];
        through_losat.extend(words);
        let expected = run(losat(), &through_losat);
        let output = run(shim(&dir, name), words);
        assert!(expected.status.success(), "{name}: {expected:?}");
        assert!(!expected.stdout.is_empty(), "{name}");
        assert_eq!(output.stdout, expected.stdout, "{name}");
        assert_eq!(output.stderr, expected.stderr, "{name}");
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn installing_twice_is_fine_and_another_file_is_left_alone() {
    let dir = temp_dir("taken");
    let first = install(&dir);
    assert!(first.status.success(), "{first:?}");
    let listed: Vec<PathBuf> = String::from_utf8(first.stdout.clone())
        .unwrap()
        .lines()
        .map(PathBuf::from)
        .collect();
    let dir = fs::canonicalize(&dir).unwrap();
    let names: Vec<PathBuf> = PROGRAMS.iter().map(|name| shim(&dir, name)).collect();
    assert_eq!(listed, names);
    let again = install(&dir);
    assert!(again.status.success(), "{again:?}");
    assert_eq!(again.stdout, first.stdout);

    // NCBI's program in the directory: not replaced, and no link is made beside it.
    let other = temp_dir("taken2");
    fs::write(shim(&other, "blastp"), "NCBI\n").unwrap();
    let refused = install(&other);
    assert_eq!(refused.status.code(), Some(1), "{refused:?}");
    assert!(
        String::from_utf8_lossy(&refused.stderr).contains("blastp"),
        "{refused:?}"
    );
    assert_eq!(fs::read(shim(&other, "blastp")).unwrap(), b"NCBI\n");
    assert!(!shim(&other, "blastn").exists());

    for words in [
        &["--install-shims", "a", "b"][..],
        &["--install-shims", "-h"],
    ] {
        let usage = run(losat(), words);
        assert_eq!(usage.status.code(), Some(1), "{words:?} {usage:?}");
    }
    assert!(!Path::new(env!("CARGO_MANIFEST_DIR")).join("-h").exists());

    // A link left by a `losat` that moved is replaced; the single-dash form of `-help` works.
    let moved = temp_dir("moved");
    #[cfg(unix)]
    std::os::unix::fs::symlink("/no/such/dir/losat", shim(&moved, "blastx")).unwrap();
    #[cfg(not(unix))]
    fs::copy(losat(), shim(&moved, "blastx")).unwrap();
    let replaced = run(losat(), &["-install-shims", moved.to_str().unwrap()]);
    assert!(replaced.status.success(), "{replaced:?}");
    let version = run(shim(&moved, "blastx"), &["-version"]);
    assert!(
        version.stdout.starts_with(b"blastx: 2.17.0+\n"),
        "{version:?}"
    );
    let help = run(losat(), &["-h"]);
    assert!(
        String::from_utf8_lossy(&help.stdout).contains("losat -install-shims [DIR]"),
        "{help:?}"
    );
    fs::remove_dir_all(&moved).unwrap();
    fs::remove_dir_all(&dir).unwrap();
    fs::remove_dir_all(&other).unwrap();
}

#[test]
fn the_default_directory_is_the_executables_and_the_links_are_relative() {
    let dir = temp_dir("default");
    let copy = dir.join(format!("losat{}", std::env::consts::EXE_SUFFIX));
    // A hard link: a copy just written can be busy for the exec of another test's fork.
    fs::hard_link(losat(), &copy).unwrap();
    let installed = run(&copy, &["--install-shims"]);
    assert!(installed.status.success(), "{installed:?}");
    let version = run(shim(&dir, "blastn"), &["-version"]);
    assert!(
        version.stdout.starts_with(b"blastn: 2.17.0+\n"),
        "{version:?}"
    );
    #[cfg(unix)]
    assert_eq!(
        fs::read_link(shim(&dir, "blastn")).unwrap(),
        Path::new("losat")
    );
    fs::remove_dir_all(&dir).unwrap();
}
