//! `losat rpsbproc` against NCBI's rpsbproc 0.5.1 (the binary of the 2025-05-13 package; bioconda's
//! 0.5.1 on the BLAST+ 2.17.0 toolkit gave the same bytes): every case of
//! `tests/fixtures/rpsbproc/CASES.tsv` (README.md there: the archives are NCBI rpsblast 2.17.0
//! `-outfmt 11` of E. coli K-12 MG1655 proteins against CDD v3.22, the data files are CDD v3.22
//! cut to the rows the archives reach). A `compare` case gives NCBI's standard output, standard
//! error and exit status (and the `-o` file); a `refuse` case is refused by LOSAT (the sites and
//! motifs of `-t feats` and `-t both`, item PS16b; the XML input, item PS16d) before any byte of
//! output. NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rpsbproc")
}

/// A directory of its own for one test (the tests of a binary run in parallel) holding a copy
/// of the fixture's `arch/` and `data/` (a copy, not a link: Windows), removed at the end.
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
        for sub in ["arch", "arch/bad", "data"] {
            fs::create_dir_all(dir.join(sub)).unwrap();
            for entry in fs::read_dir(fixtures().join(sub)).unwrap() {
                let entry = entry.unwrap();
                if entry.file_type().unwrap().is_file() {
                    fs::copy(entry.path(), dir.join(sub).join(entry.file_name())).unwrap();
                }
            }
        }
        RunDir(dir)
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") {
            command.env_remove(&name);
        }
    }
    let stdin = match stdin {
        Some(path) => Stdio::from(fs::File::open(path).unwrap()),
        None => Stdio::null(),
    };
    let output = command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .arg("rpsbproc")
        .args(words)
        .stdin(stdin)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
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

#[test]
fn rpsbproc_gives_ncbis_bytes() {
    let run_dir = RunDir::new("compare");
    let dir = &run_dir.0;
    let oracle = fixtures().join("oracle");
    let mut compared = 0;
    for [case, usage, stdin, setup, arguments] in cases() {
        if usage != "compare" {
            continue;
        }
        let _ = fs::remove_file(dir.join("out.txt"));
        let _ = fs::remove_dir_all(dir.join("data_missing"));
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
        let stdin = (stdin != "-").then(|| fixtures().join(&stdin));
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, stdin.as_deref());
        let expected = |suffix: &str| fs::read(oracle.join(format!("{case}.{suffix}"))).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&out),
            String::from_utf8_lossy(&expected("out")),
            "{case}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&err),
            String::from_utf8_lossy(&expected("err")),
            "{case}: stderr"
        );
        let rc: i32 = String::from_utf8(expected("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        assert_eq!(status, rc, "{case}: exit status");
        let outfile = oracle.join(format!("{case}.outfile"));
        if outfile.exists() {
            assert_eq!(
                fs::read(dir.join("out.txt")).unwrap(),
                fs::read(outfile).unwrap(),
                "{case}: -o file"
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 84);
}

#[test]
fn rpsbproc_refuses_what_it_does_not_give() {
    let run_dir = RunDir::new("refuse");
    let dir = &run_dir.0;
    let mut refused = 0;
    for [case, usage, _, _, arguments] in cases() {
        if usage != "refuse" {
            continue;
        }
        let words: Vec<&str> = arguments.split(' ').collect();
        let (out, err, status) = rpsbproc(dir, &words, None);
        let err = String::from_utf8_lossy(&err);
        assert!(
            out.is_empty(),
            "{case}: stdout {}",
            String::from_utf8_lossy(&out)
        );
        assert_ne!(status, 0, "{case}");
        assert!(
            err.contains("not supported by LOSAT's rpsbproc"),
            "{case}: {err}"
        );
        refused += 1;
    }
    assert_eq!(refused, 4);
}

// The command line outside NCBI's report (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`): `-version`
// prints rpsbproc's version (NCBI: `rpsbproc: 0.5.1`), `-help` LOSAT's help, and NCBI's toolkit
// words such as `-dryrun` are refused.
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
    let (out, err, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-d", "data", "-dryrun"], None);
    assert!(out.is_empty());
    assert_eq!(status, 2, "{}", String::from_utf8_lossy(&err));
    let (_, _, status) = rpsbproc(dir, &["-i", "arch/one.asn", "-m", "xx"], None);
    assert_eq!(status, 2);
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
}
