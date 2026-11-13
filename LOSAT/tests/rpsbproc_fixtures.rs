//! `losat rpsbproc` against NCBI's rpsbproc 0.5.1 (the binary of the 2025-05-13 package; bioconda's
//! 0.5.1 on the BLAST+ 2.17.0 toolkit gave the same bytes): every case of
//! `tests/fixtures/rpsbproc/CASES.tsv` (README.md there: the archives are NCBI rpsblast 2.17.0
//! `-outfmt 11` of E. coli K-12 MG1655 proteins against CDD v3.22, the data files are CDD v3.22
//! cut to the rows the archives reach). A `compare` case gives NCBI's standard output, standard
//! error and exit status (and the `-o` file `out.txt` and the `-logfile` file `log.txt`); a
//! `refuse` case is refused by LOSAT (exit 1, the
//! last line of standard error `Error: ... is not supported by LOSAT's rpsbproc...`, standard
//! output cut at a line end) with the bytes of `refused/` (LOSAT's own; README.md). NCBI is not
//! run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rpsbproc")
}

/// A directory of its own for one test (the tests of a binary run in parallel) holding a copy
/// of the fixture's `arch/` and `data/` and of the files of `root/` (a copy, not a link:
/// Windows), removed at the end.
struct RunDir(PathBuf);

impl RunDir {
    fn new(stem: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "losat_rpsbproc_{stem}_{}_{nanos}",
            std::process::id()
        ));
        for (from, to) in [
            ("arch", "arch"),
            ("arch/bad", "arch/bad"),
            ("data", "data"),
            ("root", ""),
        ] {
            copy_files(&fixtures().join(from), &dir.join(to));
        }
        RunDir(dir)
    }
}

/// Copies the files of `from` into `to` (made if missing).
fn copy_files(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_file() {
            fs::copy(entry.path(), to.join(entry.file_name())).unwrap();
        }
    }
}

impl Drop for RunDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One run of `losat rpsbproc words...` in `dir` with `stdin` (a file of the fixtures, or
/// nothing): stdout, stderr, exit status. No LOSAT or NCBI variable, no registry file
/// (`NCBI_CONFIG_PATH` names the run directory, which has none).
fn rpsbproc(dir: &Path, words: &[&str], stdin: Option<&Path>) -> (Vec<u8>, Vec<u8>, i32) {
    let stdin = match stdin {
        Some(path) => Stdio::from(fs::File::open(path).unwrap()),
        None => Stdio::null(),
    };
    let output = command(dir).args(words).stdin(stdin).output().unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// `losat rpsbproc` in `dir`, without LOSAT or NCBI variables and registry files.
fn command(dir: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .arg("rpsbproc");
    command
}

/// The rows of CASES.tsv: case, use, stdin, setup, arguments.
fn cases() -> Vec<[String; 5]> {
    let table = fs::read_to_string(fixtures().join("CASES.tsv")).unwrap();
    table
        .lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 5, "{line}");
            [0, 1, 2, 3, 4].map(|i| fields[i].to_string())
        })
        .collect()
}

/// Makes the files of a case's `setup` in `dir`, its parts separated by `,`: `data-with:<a>+<b>`,
/// the directory `<b>/` holding the data files with those of `data_with/<a>/` and then of
/// `data_with/<b>/` in their place; `data-without:<file>`, `data_missing/` without that file;
/// `conf:<file>`, `conf/<file>` as `c.ini` (for `-conffile c.ini`); `ini:<file>`, `conf/<file>`
/// as `rpsbproc.ini`; `ncbirc:<file>`, `conf/<file>` as `.ncbirc`. The files that a run makes
/// and the registry files of the case before are removed first.
fn set_up(dir: &Path, setup: &str) {
    for made in ["out.txt", "log.txt", "c.ini", "rpsbproc.ini", ".ncbirc"] {
        let _ = fs::remove_file(dir.join(made));
    }
    let _ = fs::remove_dir_all(dir.join("data_missing"));
    for part in setup.split(',') {
        set_up_part(dir, part);
    }
}

/// One part of `set_up`.
fn set_up_part(dir: &Path, setup: &str) {
    for (kind, name) in [
        ("conf:", "c.ini"),
        ("ini:", "rpsbproc.ini"),
        ("ncbirc:", ".ncbirc"),
    ] {
        if let Some(file) = setup.strip_prefix(kind) {
            fs::copy(fixtures().join("conf").join(file), dir.join(name)).unwrap();
        }
    }
    if let Some(with) = setup.strip_prefix("data-with:") {
        let layers: Vec<&str> = with.split('+').collect();
        let target = dir.join(layers[layers.len() - 1]);
        copy_files(&dir.join("data"), &target);
        for layer in layers {
            copy_files(&fixtures().join("data_with").join(layer), &target);
        }
    }
    if let Some(without) = setup.strip_prefix("data-without:") {
        fs::create_dir_all(dir.join("data_missing")).unwrap();
        for entry in fs::read_dir(dir.join("data")).unwrap() {
            let entry = entry.unwrap();
            if entry.file_name().to_str() != Some(without) {
                fs::copy(
                    entry.path(),
                    dir.join("data_missing").join(entry.file_name()),
                )
                .unwrap();
            }
        }
    }
}

/// Compares a run's stdout, stderr and exit status with `<dir>/<case>.{out,err,rc}` by bytes
/// (the text is only shown).
fn assert_same(dir: &Path, case: &str, out: &[u8], err: &[u8], status: i32) {
    let expected = |suffix: &str| fs::read(dir.join(format!("{case}.{suffix}"))).unwrap();
    for (name, got, want) in [
        ("stdout", out, expected("out")),
        ("stderr", err, expected("err")),
    ] {
        assert!(
            *got == want,
            "{case}: {name}\nLOSAT:\n{}\nexpected:\n{}",
            String::from_utf8_lossy(got),
            String::from_utf8_lossy(&want)
        );
    }
    let rc: i32 = String::from_utf8(expected("rc"))
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(status, rc, "{case}: exit status");
}

/// Runs the rows of CASES.tsv whose `use` is `usage` and compares them with `oracle/`: the
/// number of rows.
fn compare_rows(stem: &str, usage: &str) -> usize {
    let run_dir = RunDir::new(stem);
    let dir = &run_dir.0;
    let oracle = fixtures().join("oracle");
    let mut compared = 0;
    for [case, row_usage, stdin, setup, arguments] in cases() {
        if row_usage != usage {
            continue;
        }
        set_up(dir, &setup);
        let stdin = (stdin != "-").then(|| fixtures().join(&stdin));
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, stdin.as_deref());
        assert_same(&oracle, &case, &out, &err, status);
        for (suffix, made) in [("outfile", "out.txt"), ("log", "log.txt")] {
            let expected = oracle.join(format!("{case}.{suffix}"));
            let got = fs::read(dir.join(made)).ok();
            assert_eq!(
                got,
                fs::read(&expected).ok(),
                "{case}: the file {made} (expected {})",
                expected.display()
            );
        }
        compared += 1;
    }
    compared
}

#[test]
fn rpsbproc_gives_ncbis_bytes() {
    assert_eq!(compare_rows("compare", "compare"), 183);
}

// The rows of `use` `linux`: `-logfile /dev/stdout` writes the diagnostics into the pipe of the
// standard output, where the report's part of a line that NCBI's `cout` keeps in glibc's buffer
// (of `st_blksize`, 4096 for a pipe on Linux) comes after them (`rpsbproc::NcbiStream::Cout`).
// NCBI's builds are Linux builds; the other platforms' buffers are not NCBI's.
#[cfg(target_os = "linux")]
#[test]
fn rpsbproc_gives_ncbis_bytes_through_a_pipe_on_linux() {
    assert_eq!(compare_rows("linux", "linux"), 1);
}

// The refusals end LOSAT with exit 1 and the refusal as the last line of standard error, and
// leave no part of a line on standard output (the report is cut before the refused site; the
// files of `refused/` hold LOSAT's bytes, NCBI's output up to that line, README.md).
#[test]
fn rpsbproc_refuses_what_it_does_not_give() {
    let run_dir = RunDir::new("refuse");
    let dir = &run_dir.0;
    let expected = fixtures().join("refused");
    let mut refused = 0;
    for [case, usage, stdin, setup, arguments] in cases() {
        if usage != "refuse" {
            continue;
        }
        set_up(dir, &setup);
        let stdin = (stdin != "-").then(|| fixtures().join(&stdin));
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, stdin.as_deref());
        assert_eq!(status, 1, "{case}");
        assert!(out.is_empty() || out.ends_with(b"\n"), "{case}: stdout");
        let last = err
            .strip_suffix(b"\n")
            .and_then(|err| err.split(|&byte| byte == b'\n').next_back())
            .unwrap_or_default();
        assert!(
            last.starts_with(b"Error: ")
                && String::from_utf8_lossy(last).contains(" is not supported by LOSAT's rpsbproc"),
            "{case}: {}",
            String::from_utf8_lossy(&err)
        );
        assert_same(&expected, &case, &out, &err, status);
        refused += 1;
    }
    assert_eq!(refused, 5);
}

// The command line outside NCBI's report (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`): `-version`
// prints rpsbproc's version (NCBI: `rpsbproc: 0.5.1`), so does `-version-full` in NCBI's form
// with LOSAT's values (line RG item E3), `-help` LOSAT's help,
// `-dryrun=1` is an unknown word as in NCBI (`-dryrun` is not an argument), an
// unknown word is named as it was given, and a file name that is not UTF-8 is refused before
// anything is written.
#[test]
fn rpsbproc_command_line() {
    let run_dir = RunDir::new("cli");
    let dir = &run_dir.0;
    let (out, _, status) = rpsbproc(dir, &["-version"], None);
    assert_eq!(status, 0);
    let out = String::from_utf8(out).unwrap();
    assert!(
        out.starts_with("rpsbproc: 0.5.1\n Package: LOSAT "),
        "{out}"
    );
    let (out, _, status) = rpsbproc(dir, &["-help"], None);
    assert_eq!(status, 0);
    assert!(String::from_utf8(out).unwrap().contains("--infile"));
    let (out, err, status) = rpsbproc(
        dir,
        &["-i", "arch/one.asn", "-d", "data", "-version-full"],
        None,
    );
    assert_eq!((status, err), (0, Vec::new()));
    let out = String::from_utf8(out).unwrap();
    assert!(
        out.starts_with("rpsbproc: 0.5.1\n Package: LOSAT ")
            && out.ends_with(" Package-Config:  \n GI-64bit:  TRUE\n"),
        "{out}"
    );
    let (out, err, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-dryrun=1"], None);
    assert!(out.is_empty());
    let err = String::from_utf8(err).unwrap();
    assert_eq!(status, 2, "{err}");
    assert!(
        err.contains("unknown option or argument '-dryrun=1'"),
        "{err}"
    );
    let (_, _, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-m", "xx"], None);
    assert_eq!(status, 2);
    let (_, err, status) = rpsbproc(dir, &["--i", "arch/one.asn"], None);
    let err = String::from_utf8(err).unwrap();
    assert_eq!(status, 2, "{err}");
    assert!(err.contains("unknown option or argument '--i'"), "{err}");
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let name = std::ffi::OsStr::from_bytes(b"x\xff.asn");
        for option in ["-i", "-o", "-d"] {
            let output = command(dir)
                .args(["-t", "doms", option])
                .arg(name)
                .stdin(Stdio::null())
                .output()
                .unwrap();
            let err = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(1), "{option}: {err}");
            assert!(output.stdout.is_empty());
            assert!(
                err.starts_with(&format!("Error: the -{} file name ", &option[1..]))
                    && err.contains("is not UTF-8")
                    && err.ends_with("not supported by LOSAT's rpsbproc\n"),
                "{option}: {err}"
            );
        }
    }
}

// The environment variables that change what NCBI's object stream reads or what its Seq-ids give
// are refused with the values that change it (`rpsbproc::reader_variable`; the gi parameters only
// for an archive with a gi query); the other values run.
#[test]
fn rpsbproc_refuses_reader_variables() {
    let run_dir = RunDir::new("env");
    let dir = &run_dir.0;
    let one = fs::read_to_string(dir.join("arch/one.asn")).unwrap();
    let with_gi = one.replacen(
        "local str \"Query_1\"",
        "local str \"Query_1\",\n gi 16128001",
        1,
    );
    assert_ne!(with_gi, one);
    fs::write(dir.join("gi.asn"), with_gi).unwrap();
    let run = |name: &str, value: &str, input: &str| {
        let output = command(dir)
            .env(name, value)
            .args(["-i", input, "-d", "data", "-t", "doms"])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        (
            output.status.code().unwrap(),
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    for (name, value, input) in [
        ("OBJECTS_SEQ_DESCR_ALLOW_EMPTY", "1", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_READ", "no", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_GET", "DefValue", "arch/one.asn"),
        ("SEQ_ID_AVOID_GI", "1", "gi.asn"),
        ("SEQ_ID_PREFER_ACCESSION_OVER_GI", "true", "gi.asn"),
    ] {
        let (status, err) = run(name, value, input);
        assert_eq!(status, 1, "{name}: {err}");
        assert!(
            err.contains(name) && err.contains("not supported by LOSAT's rpsbproc"),
            "{name}: {err}"
        );
    }
    for (name, value, input) in [
        ("SERIAL_VERIFY_DATA_READ", "yes", "arch/one.asn"),
        ("SERIAL_VERIFY_DATA_GET", "maybe", "arch/one.asn"),
        ("SEQ_ID_AVOID_GI", "1", "arch/one.asn"),
    ] {
        let (status, err) = run(name, value, input);
        assert_eq!(status, 0, "{name}={value}: {err}");
    }
    // A dry run reads no archive (NCBI rpsbproc 0.5.1 passes with these variables: line RG E1c
    // review 1, `r_dry_serial`).
    let output = command(dir)
        .env("SERIAL_VERIFY_DATA_READ", "no")
        .args(["-i", "arch/one.asn", "-d", "data", "-dryrun"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{err}");
    assert!(
        err.ends_with("Dryrun passed, all resources are present and loaded successfully.\n"),
        "{err}"
    );
}

// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 9: a failed write of the report stops LOSAT with
// its message and exit 4 (NCBI's rpsbproc drops the writes and exits 0), with no progress line
// after the failure.
#[cfg(target_os = "linux")]
#[test]
fn rpsbproc_stops_at_a_failed_write() {
    let run_dir = RunDir::new("full");
    let dir = &run_dir.0;
    let (out, err, status) = rpsbproc(
        dir,
        &[
            "-i",
            "arch/three.asn",
            "-d",
            "data",
            "-t",
            "doms",
            "-o",
            "/dev/full",
        ],
        None,
    );
    let err = String::from_utf8(err).unwrap();
    assert!(out.is_empty());
    assert_eq!(status, 4, "{err}");
    assert!(
        err.ends_with("For more information please refer to the README file\nrpsbproc failed to write output: No space left on device (os error 28)\n"),
        "{err}"
    );
    // A standard error whose reader has closed (NCBI: SIGPIPE): exit 4, no panic, also for a
    // refusal and for an error of the command line.
    for words in [
        &["-i", "arch/three.asn", "-d", "data", "-t", "doms"][..],
        &["-i", "arch/three.asn", "-d", "data", "-t", "feats"],
        &["-i", "arch/three.asn", "-x"],
    ] {
        let (reader, writer) = std::io::pipe().unwrap();
        drop(reader);
        let output = command(dir)
            .args(words)
            .stdin(Stdio::null())
            .stderr(Stdio::from(writer))
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(4), "{words:?}");
        assert!(output.stdout.is_empty(), "{words:?}");
    }
}
