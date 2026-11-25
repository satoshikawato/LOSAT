//! `losat makembindex` whose standard output (or `-stat` output) is a pipe that its reader has
//! closed, against NCBI BLAST+ 2.17.0 makembindex (probes of 2026-11-25 for line RG item KD).
//!
//! Normal index builds write progress to standard error and write nothing to standard output, so a
//! closed standard-output pipe completes normally in both NCBI and LOSAT. `-show_filters`, `-h`,
//! `-help`, `-help-full`, `-xmlhelp` and `-version*` write to standard output: NCBI 2.17.0 is
//! ended by SIGPIPE (exit status 141, no message) outside Windows, and LOSAT stops with exit
//! status 0 and no message (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 10). `-stat` flushes
//! each non-empty word count line to its file: when that file is a closed pipe (such as `-stat
//! /dev/stdout` with standard output closed, or a pipe whose reader closes mid-run), NCBI 2.17.0
//! is ended by SIGPIPE at that line, leaving the current `.NN.idx` incomplete with no `done`, no
//! `.map`, no `.shd`, and no later volumes; LOSAT stops at the same point with exit status 0 and
//! no message. Non-pipe write failures (`/dev/full`) and volumes with no indexed words finish with
//! exit status 0 in both.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:236-240
//! (`-show_filters` to standard output), c++/src/algo/blast/dbindex/dbindex_factory.cpp:1408-1436
//! (`-stat` flushed per line).

// A closed pipe ends NCBI's process by SIGPIPE outside Windows only; on Windows its write fails
// without a signal and the run goes on (exception 10).
#![cfg(unix)]

use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const SMALL_FA: &str = ">seq1 first\nACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGC\n>seq2 second\nTTGACCGGTTAACCGGTTAAGGCCTTAAGGCCTTAAACGTACGTACGTAGCTAGCTAGCATCGATCGATCG\n";
const TINY_FA: &str = ">t1\nACGTA\n";
const MV_IDMAP_ARGS: [&str; 6] = ["-input", "mvdb", "-nmer", "8", "-idmap", "true"];
const IDX_DONE_ERR: &str = "creating idx.00.idx...done\ncreating idx.01.idx...removed (empty)\n";

fn write_small_fa(dir: &Path) {
    fs::write(dir.join("in.fa"), SMALL_FA).unwrap();
}

fn stderr_text(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir = fs::canonicalize(std::env::temp_dir())
            .unwrap()
            .join(format!(
                "losat-mbi-closed-pipe-{name}-{}",
                std::process::id()
            ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A pipe whose reader is closed.
fn closed_pipe() -> Stdio {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// A clean `losat` command in `dir` with `program` and `words`.
fn losat_cmd(dir: &Path, program: &str, words: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(dir)
        .env("HOME", dir)
        .env("BLAST_USAGE_REPORT", "false")
        .arg(program)
        .args(words)
        .stdin(Stdio::null());
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy().into_owned();
        if (name.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || name.starts_with("NCBI_CONFIG")
            || name == "BLASTDB"
            || name.starts_with("LMDB_")
        {
            command.env_remove(name);
        }
    }
    command
}

fn makembindex(dir: &Path, words: &[&str], stdout: Stdio) -> Output {
    losat_cmd(dir, "makembindex", words)
        .stdout(stdout)
        .stderr(Stdio::piped())
        .output()
        .expect("run losat makembindex")
}

/// Names of index/super-header/idmap files in `dir`, sorted.
fn index_files(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| {
            name.ends_with(".idx") || name.ends_with(".shd") || name.ends_with(".idx.map")
        })
        .collect();
    names.sort();
    names
}

/// Deterministic FASTA with `records` sequences of `len` bases (LCG), enough to populate tens of
/// thousands of distinct 8-mers and exceed the kernel pipe buffer for `-stat`.
fn pseudo_fasta(records: usize, len: usize) -> String {
    let mut state = 1u64;
    let mut out = String::new();
    for r in 0..records {
        out.push_str(&format!(">rec{r}\n"));
        for i in 0..len {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(b"ACGT"[(state >> 62) as usize] as char);
            if (i + 1).is_multiple_of(80) {
                out.push('\n');
            }
        }
        if !len.is_multiple_of(80) {
            out.push('\n');
        }
    }
    out
}

/// Builds a multi-volume nucleotide database `name` in `dir`.
fn build_multivol_db(dir: &Path, name: &str) {
    let mut fasta = String::new();
    let unit = "ACGTACGTAC".repeat(100);
    for i in 0..60 {
        fasta.push_str(&format!(">seq{i} desc {i}\n{unit}\n{unit}\n"));
    }
    let fa = dir.join("mv.fa");
    fs::write(&fa, fasta).unwrap();
    let out = losat_cmd(
        dir,
        "makeblastdb",
        &[
            "-in",
            "mv.fa",
            "-dbtype",
            "nucl",
            "-out",
            name,
            "-max_file_sz",
            "10000B",
        ],
    )
    .stdout(Stdio::null())
    .stderr(Stdio::piped())
    .output()
    .unwrap();
    assert!(out.status.success(), "{:?}", out);
}

#[test]
fn normal_index_builds_complete_with_a_closed_stdout_pipe() {
    let tmp = TempDir::new("normal");
    write_small_fa(&tmp.0);
    let built = losat_cmd(
        &tmp.0,
        "makeblastdb",
        &[
            "-in", "in.fa", "-dbtype", "nucl", "-out", "db", "-title", "t",
        ],
    )
    .stdout(Stdio::null())
    .stderr(Stdio::piped())
    .output()
    .unwrap();
    assert!(built.status.success());

    // New-style single-volume build with -verbosity verbose and -output warning.
    let out = makembindex(
        &tmp.0,
        &[
            "-input",
            "db",
            "-nmer",
            "8",
            "-verbosity",
            "verbose",
            "-output",
            "ignored",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    let err = stderr_text(&out);
    assert!(err.starts_with("Warning: option 'output' is ignored for new style indices\n"));
    assert!(err.contains(".00.idx...done\n"));
    assert!(err.contains(".01.idx...removed (empty)\n"));
    assert_eq!(index_files(&tmp.0), ["db.00.idx", "db.shd"]);

    // Old-style FASTA and BLASTDB builds.
    let old_fa = makembindex(
        &tmp.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "in.fa",
            "-output",
            "old_fa",
            "-nmer",
            "8",
        ],
        closed_pipe(),
    );
    assert_eq!(old_fa.status.code(), Some(0));
    assert!(tmp.0.join("old_fa.00.idx").is_file());

    let old_db = makembindex(
        &tmp.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "blastdb",
            "-input",
            "db",
            "-output",
            "old_db",
            "-nmer",
            "8",
        ],
        closed_pipe(),
    );
    assert_eq!(old_db.status.code(), Some(0));
    assert!(tmp.0.join("old_db.00.idx").is_file());

    // Multi-volume new-style DB with -idmap true.
    let mv = TempDir::new("normal-mv");
    build_multivol_db(&mv.0, "mvdb");
    let out = makembindex(&mv.0, &MV_IDMAP_ARGS, closed_pipe());
    assert_eq!(out.status.code(), Some(0));
    let files = index_files(&mv.0);
    assert!(files.iter().any(|f| f == "mvdb.00.00.idx"), "{files:?}");
    assert!(files.iter().any(|f| f == "mvdb.00.00.idx.map"), "{files:?}");
    assert!(files.iter().any(|f| f == "mvdb.00.shd"), "{files:?}");
    assert!(files.iter().any(|f| f == "mvdb.01.00.idx"), "{files:?}");
    assert!(files.iter().any(|f| f == "mvdb.01.shd"), "{files:?}");
}

#[test]
fn show_filters_help_and_version_on_a_closed_pipe_end_silently() {
    let tmp = TempDir::new("show-help");
    write_small_fa(&tmp.0);
    let mask = losat_cmd(
        &tmp.0,
        "dustmasker",
        &[
            "-in",
            "in.fa",
            "-out",
            "dust.asnb",
            "-outfmt",
            "maskinfo_asn1_bin",
        ],
    )
    .stdout(Stdio::null())
    .stderr(Stdio::piped())
    .output()
    .unwrap();
    assert!(mask.status.success());
    let built = losat_cmd(
        &tmp.0,
        "makeblastdb",
        &[
            "-in",
            "in.fa",
            "-dbtype",
            "nucl",
            "-mask_data",
            "dust.asnb",
            "-out",
            "mdb",
        ],
    )
    .stdout(Stdio::null())
    .stderr(Stdio::piped())
    .output()
    .unwrap();
    assert!(built.status.success());

    for words in [
        &["-input", "mdb", "-show_filters"][..],
        &["-h"],
        &["-help"],
        &["-help-full"],
        &["-xmlhelp"],
        &["-version"],
        &["-version-full"],
    ] {
        let out = makembindex(&tmp.0, words, closed_pipe());
        assert_eq!(out.status.code(), Some(0), "{words:?}");
        assert_eq!(String::from_utf8_lossy(&out.stderr), "", "{words:?}");
    }

    // Earlier warning reaches stderr before -show_filters writes to the closed stdout pipe.
    let warn = makembindex(
        &tmp.0,
        &["-input", "mdb", "-show_filters", "-stride", "5"],
        closed_pipe(),
    );
    assert_eq!(warn.status.code(), Some(0));
    assert_eq!(
        String::from_utf8_lossy(&warn.stderr),
        "Warning: -stride has no effect upon legacy index creation\n"
    );

    // Missing database fails before writing to stdout (exit 3, even on a closed pipe).
    let missing = makembindex(
        &tmp.0,
        &["-input", "nosuch", "-show_filters"],
        closed_pipe(),
    );
    assert_eq!(missing.status.code(), Some(3));
    assert!(
        String::from_utf8_lossy(&missing.stderr).contains("No alias or index file found"),
        "{}",
        String::from_utf8_lossy(&missing.stderr)
    );

    // Non-pipe write failure (/dev/full) finishes silently with exit 0.
    if let Ok(full) = File::options().write(true).open("/dev/full") {
        let out = makembindex(&tmp.0, &["-input", "mdb", "-show_filters"], full.into());
        assert_eq!(out.status.code(), Some(0));
        assert_eq!(String::from_utf8_lossy(&out.stderr), "");
    }
}

#[test]
fn stat_on_a_closed_pipe_stops_at_once() {
    // Multi-volume new-style DB: stops in vol 00's first index volume with no done, no .map,
    // no .shd, and no later DB volumes.
    let mv_ref = TempDir::new("stat-mv-ref");
    build_multivol_db(&mv_ref.0, "mvdb");
    let full_run = makembindex(&mv_ref.0, &MV_IDMAP_ARGS, Stdio::null());
    assert_eq!(full_run.status.code(), Some(0));
    let full_vol0_len = fs::metadata(mv_ref.0.join("mvdb.00.00.idx")).unwrap().len();

    let mv = TempDir::new("stat-mv");
    build_multivol_db(&mv.0, "mvdb");
    let out = makembindex(
        &mv.0,
        &[
            "-input",
            "mvdb",
            "-nmer",
            "8",
            "-idmap",
            "true",
            "-stat",
            "/dev/stdout",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    let err = stderr_text(&out);
    assert!(
        err.starts_with("creating ") && err.ends_with("mvdb.00.00.idx..."),
        "unexpected stderr: {err:?}"
    );
    assert_eq!(index_files(&mv.0), ["mvdb.00.00.idx"]);
    let partial_len = fs::metadata(mv.0.join("mvdb.00.00.idx")).unwrap().len();
    assert!(
        partial_len < full_vol0_len,
        "expected partial volume ({partial_len} < {full_vol0_len})"
    );

    // Old-style multi-volume index (-volsize 1): stops in .00.idx, no .map, no .01.idx.
    let old_mv = TempDir::new("stat-old-mv");
    fs::write(old_mv.0.join("big.fa"), pseudo_fasta(4, 300_000)).unwrap();
    let out = makembindex(
        &old_mv.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "big.fa",
            "-output",
            "idx",
            "-nmer",
            "8",
            "-volsize",
            "1",
            "-idmap",
            "true",
            "-stat",
            "/dev/stdout",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert_eq!(stderr_text(&out), "creating idx.00.idx...");
    assert_eq!(index_files(&old_mv.0), ["idx.00.idx"]);

    // Unwritable -output path still stops at -stat /dev/stdout before done and .01.idx.
    let bad_out = TempDir::new("stat-bad-out");
    write_small_fa(&bad_out.0);
    let out = makembindex(
        &bad_out.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "in.fa",
            "-output",
            "nodir/idx",
            "-nmer",
            "8",
            "-idmap",
            "true",
            "-stat",
            "/dev/stdout",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert_eq!(stderr_text(&out), "creating nodir/idx.00.idx...");
    assert!(index_files(&bad_out.0).is_empty());

    // Volume with zero 8-mers writes nothing to -stat and finishes normally even on a closed pipe.
    let tiny = TempDir::new("stat-tiny");
    fs::write(tiny.0.join("tiny.fa"), TINY_FA).unwrap();
    let out = makembindex(
        &tiny.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "tiny.fa",
            "-output",
            "idx",
            "-nmer",
            "8",
            "-stat",
            "/dev/stdout",
        ],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert_eq!(stderr_text(&out), IDX_DONE_ERR, "tiny 0-nmer");
    assert_eq!(index_files(&tiny.0), ["idx.00.idx"]);

    // Non-pipe failure (-stat /dev/full) is ignored and finishes the full index and .map.
    if File::options().write(true).open("/dev/full").is_ok() {
        let full = TempDir::new("stat-dev-full");
        write_small_fa(&full.0);
        let out = makembindex(
            &full.0,
            &[
                "-old_style_index",
                "true",
                "-iformat",
                "fasta",
                "-input",
                "in.fa",
                "-output",
                "idx",
                "-nmer",
                "8",
                "-idmap",
                "true",
                "-stat",
                "/dev/full",
            ],
            Stdio::null(),
        );
        assert_eq!(out.status.code(), Some(0));
        assert_eq!(stderr_text(&out), IDX_DONE_ERR, "stat /dev/full");
        assert_eq!(index_files(&full.0), ["idx.00.idx", "idx.00.idx.map"]);
    }
}

#[test]
fn stat_pipe_closed_after_first_line_stops_at_once() {
    let tmp = TempDir::new("stat-mid");
    // 250,000 pseudo-random bases cover ~60,000 distinct 8-mers (>800 KiB of -stat output),
    // far larger than the kernel pipe buffer so closing after line 1 reliably hits EPIPE.
    fs::write(tmp.0.join("big.fa"), pseudo_fasta(1, 250_000)).unwrap();

    let full = makembindex(
        &tmp.0,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "big.fa",
            "-output",
            "ref",
            "-nmer",
            "8",
        ],
        Stdio::null(),
    );
    assert_eq!(full.status.code(), Some(0));
    let full_len = fs::metadata(tmp.0.join("ref.00.idx")).unwrap().len();
    fs::remove_file(tmp.0.join("ref.00.idx")).unwrap();

    let mut child = losat_cmd(
        &tmp.0,
        "makembindex",
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "big.fa",
            "-output",
            "idx",
            "-nmer",
            "8",
            "-idmap",
            "true",
            "-stat",
            "/dev/stdout",
        ],
    )
    .stdout(Stdio::piped())
    .stderr(Stdio::piped())
    .spawn()
    .expect("spawn losat makembindex");

    let mut reader = BufReader::new(child.stdout.take().unwrap());
    let mut first = String::new();
    reader.read_line(&mut first).unwrap();
    assert!(!first.is_empty(), "expected at least one -stat line");
    drop(reader);

    let output = child.wait_with_output().unwrap();
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "creating idx.00.idx..."
    );
    assert_eq!(index_files(&tmp.0), ["idx.00.idx"]);
    let partial_len = fs::metadata(tmp.0.join("idx.00.idx")).unwrap().len();
    assert!(
        partial_len < full_len,
        "expected partial volume ({partial_len} < {full_len})"
    );
}
