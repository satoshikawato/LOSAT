//! `losat makembindex` whose standard output, `-stat` output, or index output is a pipe that its
//! reader has closed, against NCBI BLAST+ 2.17.0 makembindex (probes of 2026-11-25 for line RG
//! item KD).
//!
//! Normal index builds write progress to standard error and write nothing to standard output, so a
//! closed standard-output pipe completes normally in both NCBI and LOSAT. `-show_filters`, `-h`,
//! `-help`, `-help-full`, `-xmlhelp` and `-version*` write to standard output: NCBI 2.17.0 is
//! ended by SIGPIPE (exit status 141, no message) outside Windows, and LOSAT stops with exit
//! status 0 and no message (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 10). `-stat` flushes
//! each non-empty word count line to its file: when that file is a closed pipe (such as `-stat
//! /dev/stdout` with standard output closed, or a pipe whose reader closes mid-run), NCBI 2.17.0
//! is ended by SIGPIPE at that line, leaving the current `.NN.idx` incomplete with no `done`, no
//! `.map`, no `.shd`, and no later volumes; LOSAT flushes `-stat` at each 8 KiB `.NN.idx` block
//! boundary, so it stops in the same 8 KiB block with exit status 0 and no message, leaving the
//! same incomplete `.NN.idx` bytes (header plus complete 8 KiB blocks, for both regular files and
//! named pipes) and files as NCBI. Named-pipe index outputs (`.NN.idx`, `.map`, `.shd`) also stop
//! at once with exit status 0 and no message. Non-pipe write failures (`/dev/full`) and volumes
//! with no indexed words finish with exit status 0 in both.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp:236-240
//! (`-show_filters` to standard output), c++/src/algo/blast/dbindex/dbindex_factory.cpp:1408-1436
//! (`-stat` flushed per line), c++/src/algo/blast/dbindex/dbindex_factory.cpp:1604-1634
//! (`SaveHeader` flushed before 8 KiB index stream blocks).

// A closed pipe ends NCBI's process by SIGPIPE outside Windows only; on Windows its write fails
// without a signal and the run goes on (exception 10).
#![cfg(unix)]

use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;

const HEADER_LEN: u64 = 44;
const STREAM_BLOCK: u64 = 8192;
const SMALL_FA: &str = ">seq1 first\nACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGC\n>seq2 second\nTTGACCGGTTAACCGGTTAAGGCCTTAAGGCCTTAAACGTACGTACGTAGCTAGCTAGCATCGATCGATCG\n";
const TINY_FA: &str = ">t1\nACGTA\n";
const MV_IDMAP_ARGS: [&str; 6] = ["-input", "mvdb", "-nmer", "8", "-idmap", "true"];
const N8_NONLEGACY: &[&str] = &["-legacy", "false", "-stride", "1", "-ws_hint", "8"];
const IDX_DONE_ERR: &str = "creating idx.00.idx...done\ncreating idx.01.idx...removed (empty)\n";

fn assert_stopped_vol0(output: &Output, ctx: &str) {
    assert_eq!(output.status.code(), Some(0), "{ctx}");
    assert_eq!(stderr_text(output), "creating idx.00.idx...", "{ctx}");
}

fn make_fifo(path: &Path) {
    assert!(Command::new("mkfifo")
        .arg(path)
        .status()
        .expect("mkfifo")
        .success());
}

fn run_full_old_style_n8(dir: &Path, input: &str) -> Vec<u8> {
    let out = makembindex(
        dir,
        &[
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            input,
            "-output",
            "ref",
            "-nmer",
            "8",
        ],
        Stdio::null(),
    );
    assert_eq!(out.status.code(), Some(0), "full old-style build");
    let bytes = fs::read(dir.join("ref.00.idx")).unwrap();
    fs::remove_file(dir.join("ref.00.idx")).unwrap();
    bytes
}

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
    // no .shd, and no later DB volumes. Its first indexed 8-mer is 0x11b1 (4529), past two 8 KiB
    // blocks of hash-table words, so NCBI leaves 44 + 2 * 8192 = 16,428 bytes.
    let mv_ref = TempDir::new("stat-mv-ref");
    build_multivol_db(&mv_ref.0, "mvdb");
    let full_run = makembindex(&mv_ref.0, &MV_IDMAP_ARGS, Stdio::null());
    assert_eq!(full_run.status.code(), Some(0));
    let full_vol0_bytes = fs::read(mv_ref.0.join("mvdb.00.00.idx")).unwrap();

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
    let partial_bytes = fs::read(mv.0.join("mvdb.00.00.idx")).unwrap();
    assert_eq!(
        partial_bytes.len() as u64,
        HEADER_LEN + 2 * STREAM_BLOCK,
        "mvdb first 8-mer is 0x11b1 (2 complete 8 KiB blocks after header)"
    );
    assert_eq!(partial_bytes, full_vol0_bytes[..partial_bytes.len()]);

    // New-style single-volume DB containing 8-mer 0 (AAAAAAAA): stops right after the 44-byte
    // header, matching NCBI's 44-byte file byte-for-byte.
    let ns = TempDir::new("stat-new-44");
    fs::write(ns.0.join("big.fa"), pseudo_fasta(1, 250_000)).unwrap();
    let built = losat_cmd(
        &ns.0,
        "makeblastdb",
        &["-in", "big.fa", "-dbtype", "nucl", "-out", "db"],
    )
    .stdout(Stdio::null())
    .stderr(Stdio::piped())
    .output()
    .unwrap();
    assert!(built.status.success());
    let full_ns = makembindex(&ns.0, &["-input", "db", "-nmer", "8"], Stdio::null());
    assert_eq!(full_ns.status.code(), Some(0));
    let full_ns_bytes = fs::read(ns.0.join("db.00.idx")).unwrap();
    fs::remove_file(ns.0.join("db.00.idx")).unwrap();
    fs::remove_file(ns.0.join("db.shd")).unwrap();
    let out = makembindex(
        &ns.0,
        &["-input", "db", "-nmer", "8", "-stat", "/dev/stdout"],
        closed_pipe(),
    );
    assert_eq!(out.status.code(), Some(0));
    assert!(stderr_text(&out).ends_with("db.00.idx..."));
    assert_eq!(index_files(&ns.0), ["db.00.idx"]);
    let ns_partial = fs::read(ns.0.join("db.00.idx")).unwrap();
    assert_eq!(ns_partial.len() as u64, HEADER_LEN);
    assert_eq!(ns_partial, full_ns_bytes[..HEADER_LEN as usize]);

    // Old-style multi-volume index (-volsize 1) containing 8-mer 0: stops in .00.idx at 44 bytes,
    // no .map, no .01.idx.
    let old_mv = TempDir::new("stat-old-mv");
    fs::write(old_mv.0.join("big.fa"), pseudo_fasta(4, 300_000)).unwrap();
    let full_old = makembindex(
        &old_mv.0,
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
            "-volsize",
            "1",
        ],
        Stdio::null(),
    );
    assert_eq!(full_old.status.code(), Some(0));
    let full_old_head =
        fs::read(old_mv.0.join("ref.00.idx")).unwrap()[..HEADER_LEN as usize].to_vec();
    for f in index_files(&old_mv.0) {
        fs::remove_file(old_mv.0.join(f)).unwrap();
    }
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
    let old_partial = fs::read(old_mv.0.join("idx.00.idx")).unwrap();
    assert_eq!(old_partial.len() as u64, HEADER_LEN);
    assert_eq!(old_partial, full_old_head);

    // Exact 8 KiB block boundary across hash widths, legacy/non-legacy headers, and the post-loop
    // flush branch: nmer `k` sits in block `(k + 2) / 2048`, matching NCBI 2.17.0 byte-for-byte.
    let boundary_cases: &[(&str, usize, usize, &[&str], &str, u64)] = &[
        (
            "n8-leg0-2045",
            8,
            2045,
            N8_NONLEGACY,
            "TTTTTTTTTT",
            HEADER_LEN,
        ),
        (
            "n8-leg0-2046",
            8,
            2046,
            N8_NONLEGACY,
            "TTTTTTTTTT",
            HEADER_LEN + STREAM_BLOCK,
        ),
        (
            "n8-leg0-2047",
            8,
            2047,
            N8_NONLEGACY,
            "TTTTTTTTTT",
            HEADER_LEN + STREAM_BLOCK,
        ),
        (
            "n8-leg0-6141",
            8,
            6141,
            N8_NONLEGACY,
            "TTTTTTTTTT",
            HEADER_LEN + 2 * STREAM_BLOCK,
        ),
        (
            "n8-leg0-6142",
            8,
            6142,
            N8_NONLEGACY,
            "TTTTTTTTTT",
            HEADER_LEN + 3 * STREAM_BLOCK,
        ),
        (
            "n8-leg1-2046",
            8,
            2046,
            &[],
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            HEADER_LEN + 31 * STREAM_BLOCK,
        ),
        (
            "n9-leg1-100000",
            9,
            100_000,
            &[],
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            HEADER_LEN + 8 * STREAM_BLOCK,
        ),
        (
            "n9-leg0-last",
            9,
            (1usize << 18) - 2,
            &["-legacy", "false", "-stride", "1", "-ws_hint", "9"],
            "",
            HEADER_LEN + 128 * STREAM_BLOCK,
        ),
    ];
    for &(tag, width, first_nmer, extra, tail, expected_len) in boundary_cases {
        let bdir = TempDir::new(&format!("stat-boundary-{tag}"));
        let mut seq = String::with_capacity(width + tail.len());
        for shift in (0..width).rev() {
            seq.push(b"ACGT"[(first_nmer >> (2 * shift)) & 3] as char);
        }
        seq.push_str(tail);
        fs::write(bdir.0.join("in.fa"), format!(">s0\n{seq}\n")).unwrap();
        let width_str = width.to_string();
        let mut common = vec![
            "-old_style_index",
            "true",
            "-iformat",
            "fasta",
            "-input",
            "in.fa",
            "-nmer",
            &width_str,
        ];
        common.extend_from_slice(extra);
        let mut full_args = common.clone();
        full_args.extend(["-output", "ref"]);
        assert_eq!(
            makembindex(&bdir.0, &full_args, Stdio::null())
                .status
                .code(),
            Some(0),
            "{tag}"
        );
        let ref_bytes = fs::read(bdir.0.join("ref.00.idx")).unwrap();
        let mut closed_args = common;
        closed_args.extend(["-output", "idx", "-stat", "/dev/stdout"]);
        let out = makembindex(&bdir.0, &closed_args, closed_pipe());
        assert_stopped_vol0(&out, tag);
        let part_bytes = fs::read(bdir.0.join("idx.00.idx")).unwrap();
        assert_eq!(
            part_bytes.len() as u64,
            expected_len,
            "{tag} expected {expected_len} bytes"
        );
        assert_eq!(
            part_bytes,
            ref_bytes[..expected_len as usize],
            "{tag} bytes mismatch"
        );
    }

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

    let full_bytes = run_full_old_style_n8(&tmp.0, "big.fa");
    let full_len = full_bytes.len() as u64;

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
    let partial_bytes = fs::read(tmp.0.join("idx.00.idx")).unwrap();
    let partial_len = partial_bytes.len() as u64;
    assert!(
        partial_len >= HEADER_LEN && partial_len < full_len,
        "expected partial volume in [{HEADER_LEN}, {full_len}), got {partial_len}"
    );
    assert!(
        (partial_len - HEADER_LEN).is_multiple_of(STREAM_BLOCK),
        "expected header + multiple of 8192 bytes, got {partial_len}"
    );
    assert_eq!(partial_bytes, full_bytes[..partial_bytes.len()]);
}

/// Opens a named FIFO for reading and drops the reader as soon as the writer opens it.
fn spawn_fifo_closer(fifo: &Path) -> (Arc<AtomicBool>, thread::JoinHandle<()>) {
    let opened = Arc::new(AtomicBool::new(false));
    let opened_flag = Arc::clone(&opened);
    let path = fifo.to_path_buf();
    let handle = thread::spawn(move || {
        if let Ok(reader) = File::open(&path) {
            opened_flag.store(true, Ordering::SeqCst);
            drop(reader);
        }
    });
    (opened, handle)
}

fn finish_fifo_closer(fifo: &Path, opened: &AtomicBool, handle: thread::JoinHandle<()>) {
    if !opened.load(Ordering::SeqCst) {
        let _ = OpenOptions::new().write(true).open(fifo);
    }
    handle.join().unwrap();
}

fn run_old_style_fifo_case(tag: &str, fasta: &str, fifo_name: &str, expected_files: &[&str]) {
    let dir = TempDir::new(tag);
    fs::write(dir.0.join("in.fa"), fasta).unwrap();
    let fifo = dir.0.join(fifo_name);
    make_fifo(&fifo);
    let (opened, handle) = spawn_fifo_closer(&fifo);
    let out = makembindex(
        &dir.0,
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
        ],
        Stdio::null(),
    );
    finish_fifo_closer(&fifo, &opened, handle);
    assert!(opened.load(Ordering::SeqCst), "{tag}: FIFO was not opened");
    assert_stopped_vol0(&out, tag);
    assert_eq!(index_files(&dir.0), expected_files, "{tag}");
}

#[test]
fn volume_and_idmap_named_fifos_stop_at_once() {
    // 1. Volume `.00.idx` is a named FIFO (>1 MiB volume exceeds any kernel pipe buffer so
    //    closing the reader on open always hits BrokenPipe, stopping before done, .map, .01.idx).
    run_old_style_fifo_case(
        "fifo-idx",
        &pseudo_fasta(1, 250_000),
        "idx.00.idx",
        &["idx.00.idx"],
    );

    // 2. `.00.idx.map` is a named FIFO (~200 KiB of titles exceeds the kernel pipe buffer so
    //    closing the reader on open always hits BrokenPipe, stopping before done and .01.idx).
    let mut many = String::new();
    let pad = "x".repeat(60);
    for i in 0..2500 {
        many.push_str(&format!(">seq_{i:05}_{pad}\nACGTACGTACGT\n"));
    }
    run_old_style_fifo_case(
        "fifo-map",
        &many,
        "idx.00.idx.map",
        &["idx.00.idx", "idx.00.idx.map"],
    );
}

#[test]
fn volume_fifo_matches_ncbi_partial_bytes_and_stat_progress() {
    let big = pseudo_fasta(1, 250_000);

    // Case A (timing-independent): `.00.idx` is a named FIFO whose reader reads everything to EOF,
    // while `-stat /dev/stdout` is a closed pipe. Because the volume stream writes in 8 KiB blocks
    // and discards the trailing incomplete block on error, the reader receives only the 44-byte
    // header, matching NCBI 2.17.0 byte-for-byte (not 44 + 8188 bytes).
    let dir_a = TempDir::new("fifo-case-a");
    fs::write(dir_a.0.join("big.fa"), &big).unwrap();
    let expected_header = run_full_old_style_n8(&dir_a.0, "big.fa")[..HEADER_LEN as usize].to_vec();

    let fifo_a = dir_a.0.join("idx.00.idx");
    make_fifo(&fifo_a);
    let opened_a = Arc::new(AtomicBool::new(false));
    let opened_a_flag = Arc::clone(&opened_a);
    let fifo_a_path = fifo_a.clone();
    let reader_a = thread::spawn(move || {
        let mut got = Vec::new();
        if let Ok(mut file) = File::open(&fifo_a_path) {
            opened_a_flag.store(true, Ordering::SeqCst);
            file.read_to_end(&mut got).unwrap();
        }
        got
    });
    let out_a = makembindex(
        &dir_a.0,
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
            "-stat",
            "/dev/stdout",
        ],
        closed_pipe(),
    );
    if !opened_a.load(Ordering::SeqCst) {
        let _ = OpenOptions::new().write(true).open(&fifo_a);
    }
    let got_a = reader_a.join().unwrap();
    assert!(opened_a.load(Ordering::SeqCst), "case A FIFO not opened");
    assert_stopped_vol0(&out_a, "case A");
    assert_eq!(index_files(&dir_a.0), ["idx.00.idx"]);
    assert_eq!(got_a.len() as u64, HEADER_LEN);
    assert_eq!(got_a, expected_header);

    // Case B (Linux): `.00.idx` is a named FIFO whose reader reads the 44-byte header, waits
    // briefly so the writer fills the 64 KiB pipe with 8 KiB blocks, then closes, while `-stat`
    // writes to a regular file. Writing `.00.idx` in 8 KiB blocks ensures the failing FIFO write
    // happens before `-stat` runs ahead to the end of the hash table (NCBI stops at nmer 0x47fd).
    #[cfg(target_os = "linux")]
    {
        let dir_b = TempDir::new("fifo-case-b");
        fs::write(dir_b.0.join("big.fa"), &big).unwrap();
        let fifo_b = dir_b.0.join("idx.00.idx");
        make_fifo(&fifo_b);
        let opened_b = Arc::new(AtomicBool::new(false));
        let opened_b_flag = Arc::clone(&opened_b);
        let fifo_b_path = fifo_b.clone();
        let reader_b = thread::spawn(move || {
            if let Ok(mut file) = File::open(&fifo_b_path) {
                opened_b_flag.store(true, Ordering::SeqCst);
                let mut header = [0u8; HEADER_LEN as usize];
                file.read_exact(&mut header).unwrap();
                thread::sleep(std::time::Duration::from_millis(100));
            }
        });
        let out_b = makembindex(
            &dir_b.0,
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
                "-stat",
                "stat.txt",
            ],
            Stdio::null(),
        );
        finish_fifo_closer(&fifo_b, &opened_b, reader_b);
        assert!(opened_b.load(Ordering::SeqCst), "case B FIFO not opened");
        assert_stopped_vol0(&out_b, "case B");
        assert_eq!(index_files(&dir_b.0), ["idx.00.idx"]);
        let stat_text = fs::read_to_string(dir_b.0.join("stat.txt")).unwrap();
        let last_line = stat_text
            .lines()
            .last()
            .expect("expected non-empty stat.txt");
        let last_nmer_hex = last_line
            .split_whitespace()
            .next()
            .expect("expected hex nmer");
        let last_nmer = u32::from_str_radix(last_nmer_hex, 16)
            .unwrap_or_else(|e| panic!("invalid hex nmer {last_nmer_hex:?}: {e}"));
        assert!(
            last_nmer < 0x4800,
            "expected -stat to stop before 0x4800 when .00.idx FIFO closes, got 0x{last_nmer:x} ({} lines)",
            stat_text.lines().count()
        );
    }
}
