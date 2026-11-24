//! `losat makeblastdb` whose standard output is a pipe that its reader has closed, against
//! makeblastdb 2.17.0 (the official Linux x64 release; probes of 2026-11-24 for line RG item KD,
//! Q-KD-1). NCBI's first line of the log meets the closed pipe and SIGPIPE ends the process (exit
//! status 141, no message) before an existing database of the name is deleted and before any file
//! of the build is made: only the output directory, which the constructor creates first, is
//! there. LOSAT ends at that line with exit status 0 and no message
//! (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10). With `-logfile` the log goes to its file and
//! the build runs as with an open pipe; a `-logfile` that is the closed pipe ends the build as the
//! standard output does.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1165-1177
//! (the constructor: the output directory, the first lines of the log, `DeleteBlastDb`),
//! c++/src/app/blastdb/makeblastdb.cpp:1199-1201 (the log: the file of `-logfile`, or `cout`).

// A closed pipe ends NCBI's process by SIGPIPE outside Windows only; on Windows its write fails
// without a signal and the build goes on (exception 10).
#![cfg(unix)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::SystemTime;

const NUC: &str = ">seq1 first\nACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGC\n>seq2 second\nTTGACCGGTTAACCGGTTAAGGCCTTAAGGCCTTAAACGTACGTACGTAGCTAGCTAGCATCGATCGATCG\n";
// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROT: &str = ">sp|P12345|AAA_HUMAN Alpha protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n";

/// The log of the nucleotide build, its time, output directory and seconds masked (NCBI's
/// `-logfile` log of the probe).
const LOG: &str = "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/d\nNew DB title:  t\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n";

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir = std::fs::canonicalize(std::env::temp_dir())
            .unwrap()
            .join(format!(
                "losat-mkdb-closed-pipe-{name}-{}",
                std::process::id()
            ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A pipe whose reader is closed.
fn closed_pipe() -> Stdio {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// `losat makeblastdb` in `dir` with `words`, its standard output `stdout`.
fn makeblastdb(dir: &Path, words: &[&str], stdout: Stdio) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(dir)
        .env("HOME", dir)
        .arg("makeblastdb")
        .args(words)
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(Stdio::piped());
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy().into_owned();
        if name.starts_with("LOSAT_")
            || name.starts_with("NCBI_CONFIG")
            || name == "BLASTDB"
            || name.starts_with("LMDB_")
        {
            command.env_remove(name);
        }
    }
    command.output().expect("run losat makeblastdb")
}

/// The names of the entries of `dir`, in order.
fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The files of `dir` with their bytes and modification times.
fn snapshot(dir: &Path) -> Vec<(String, Vec<u8>, SystemTime)> {
    names(dir)
        .into_iter()
        .map(|name| {
            let path = dir.join(&name);
            let modified = std::fs::metadata(&path).unwrap().modified().unwrap();
            (name, std::fs::read(&path).unwrap(), modified)
        })
        .collect()
}

/// The time, the output directory and the seconds of a log masked.
fn masked(log: &str, d: &Path) -> String {
    let mut out = log.replace(d.to_str().unwrap(), "<D>");
    for (start, end, fill) in [
        ("current time: ", "\n", "<T>"),
        (" sequences in ", " seconds.", "<S>"),
    ] {
        if let Some(i) = out.find(start) {
            let from = i + start.len();
            if let Some(j) = out[from..].find(end) {
                out.replace_range(from..from + j, fill);
            }
        }
    }
    out
}

#[test]
fn a_closed_pipe_at_the_first_line_makes_no_file() {
    for (name, input, words) in [
        ("nucl", NUC, &["-dbtype", "nucl"][..]),
        ("prot", PROT, &["-dbtype", "prot"][..]),
        (
            "nucl4",
            NUC,
            &["-dbtype", "nucl", "-blastdb_version", "4", "-parse_seqids"][..],
        ),
    ] {
        let tmp = TempDir::new(name);
        std::fs::write(tmp.0.join("in.fa"), input).unwrap();
        let words = [
            &["-in", "in.fa"][..],
            words,
            &["-out", "p/d", "-title", "t"],
        ]
        .concat();
        let output = makeblastdb(&tmp.0, &words, closed_pipe());
        assert_eq!(output.status.code(), Some(0), "{name}");
        assert_eq!(String::from_utf8_lossy(&output.stderr), "", "{name}");
        assert!(names(&tmp.0.join("p")).is_empty(), "{name}");
    }
}

#[test]
fn an_old_database_of_the_name_is_left_as_it_was() {
    let tmp = TempDir::new("old");
    std::fs::write(tmp.0.join("in.fa"), NUC).unwrap();
    let words = ["-in", "in.fa", "-dbtype", "nucl", "-out", "p/d", "-title"];
    let built = makeblastdb(&tmp.0, &[&words[..], &["old"]].concat(), Stdio::null());
    assert_eq!(built.status.code(), Some(0));
    let before = snapshot(&tmp.0.join("p"));
    assert!(before.iter().any(|(name, _, _)| name == "d.nsq"));
    let output = makeblastdb(&tmp.0, &[&words[..], &["t"]].concat(), closed_pipe());
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
    assert!(snapshot(&tmp.0.join("p")) == before);
}

#[test]
fn logfile_builds_with_the_pipe_closed() {
    let tmp = TempDir::new("logfile");
    std::fs::write(tmp.0.join("in.fa"), NUC).unwrap();
    let words = [
        "-in", "in.fa", "-dbtype", "nucl", "-out", "p/d", "-title", "t", "-logfile", "l.log",
    ];
    let output = makeblastdb(&tmp.0, &words, closed_pipe());
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
    assert_eq!(
        names(&tmp.0.join("p")),
        ["d.ndb", "d.nhr", "d.nin", "d.njs", "d.not", "d.nsq", "d.ntf", "d.nto"]
    );
    let log = std::fs::read_to_string(tmp.0.join("l.log")).unwrap();
    assert_eq!(masked(&log, &tmp.0.join("p")), LOG);
}

#[test]
fn a_logfile_that_is_the_closed_pipe_makes_no_file() {
    let tmp = TempDir::new("logfile-stdout");
    std::fs::write(tmp.0.join("in.fa"), NUC).unwrap();
    let words = [
        "-in",
        "in.fa",
        "-dbtype",
        "nucl",
        "-out",
        "p/d",
        "-title",
        "t",
        "-logfile",
        "/dev/stdout",
    ];
    let output = makeblastdb(&tmp.0, &words, closed_pipe());
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
    assert!(names(&tmp.0.join("p")).is_empty());
}
