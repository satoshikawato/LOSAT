//! How many subjects one entry into the search pool takes in TBLASTN (`LOSAT_X_TBNOIDCHUNK`,
//! `src/algorithm/tblastn/stage_d_pipeline.rs`), on a whole search of many short subjects.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413
//! ```c
//!     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//!
//!     /* iterate over all subject sequences */
//!     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//!            != BLAST_SEQSRC_EOF) {
//! ```
//! Each NCBI search thread takes a hundredth of the subjects at a time. Without the switch one
//! entry into LOSAT's pool takes one subject per thread; with it, one such chunk per thread. The
//! switch is not in the strict set (`src/utils/xswitch.rs`): it is off unless its variable is
//! set. It must change the number of entries and nothing in the report, so the tests compare the
//! report with the one-thread search, which has no pool and searches the subjects in NCBI's
//! loop order, and read the `pool_entries` field of the `[losat-thread-activity]` line
//! (`LOSAT_WASI_THREADS_DEBUG=1`) to check which entries ran.
//!
//! The fixture suites have no search with 200 or more subjects, the least that gives a chunk of
//! two, so the tests cut `tests/fasta/MelaMJNV.fasta` into 240 subjects of 1,000 bases and search
//! them with the first ten proteins of `tests/fasta/PemoMJNVA.faa`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const SUBJECTS: usize = 240;
const SUBJECT_LENGTH: usize = 1_000;
const QUERIES: usize = 10;
const THREADS: usize = 4;

/// The records of a FASTA file as (title line, letters).
fn records(path: &Path) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).expect("fixture");
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        match line.strip_prefix('>') {
            Some(title) => out.push((title.to_string(), String::new())),
            None => out
                .last_mut()
                .expect("title line first")
                .1
                .push_str(line.trim()),
        }
    }
    out
}

/// The query and subject files of the tests, written once: the first `QUERIES` proteins, and
/// the genome cut into `SUBJECTS` pieces of `SUBJECT_LENGTH` bases.
fn inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_oid_chunk");
        std::fs::create_dir_all(&dir).expect("work directory");
        let mut query = String::new();
        for (title, letters) in records(&fasta.join("PemoMJNVA.faa")).iter().take(QUERIES) {
            query.push_str(&format!(">{title}\n{letters}\n"));
        }
        let genome = &records(&fasta.join("MelaMJNV.fasta"))[0].1;
        assert!(
            genome.len() >= SUBJECTS * SUBJECT_LENGTH,
            "genome too short"
        );
        let mut subjects = String::new();
        for index in 0..SUBJECTS {
            let piece = &genome[index * SUBJECT_LENGTH..(index + 1) * SUBJECT_LENGTH];
            subjects.push_str(&format!(">piece{index}\n{piece}\n"));
        }
        let files = (dir.join("query.faa"), dir.join("subjects.fna"));
        std::fs::write(&files.0, query).expect("query file");
        std::fs::write(&files.1, subjects).expect("subject file");
        files
    })
}

/// One TBLASTN run of the inputs (`-outfmt 6` on standard output) with the given variables as
/// its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build with the worker probe
/// (`--all-features`) writes to that file as soon as a pool starts, so the caller's value is
/// kept, or a file of the tests is named. Returns the report and the standard error.
fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let (query, subjects) = inputs();
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// The `pool_entries` field of the one `[losat-thread-activity]` line of the subject stage.
fn pool_entries(stderr: &str) -> usize {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| {
            line.starts_with("[losat-thread-activity] program=tblastn stage=subject_core ")
        })
        .collect();
    assert_eq!(lines.len(), 1, "one query batch, one line:\n{stderr}");
    lines[0]
        .split_whitespace()
        .find_map(|field| field.strip_prefix("pool_entries="))
        .unwrap_or_else(|| panic!("no pool_entries in {}", lines[0]))
        .parse()
        .expect("number of entries")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413
// ```c
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// The one-thread search runs this loop over the subjects in order. The pool's entries of one
// subject per thread and of one chunk per thread must give its report, with the default set and
// with `LOSAT_X_DEFAULT=0`.
#[test]
fn the_report_is_the_one_thread_report_with_and_without_the_switch() {
    let (expected, _) = tblastn(1, &[]);
    let rows = expected
        .split(|&byte| byte == b'\n')
        .filter(|row| !row.is_empty());
    let hit_subjects: std::collections::BTreeSet<_> = rows
        .map(|row| {
            row.split(|&byte| byte == b'\t')
                .nth(1)
                .expect("subject column")
                .to_vec()
        })
        .collect();
    assert!(
        hit_subjects.len() > 1,
        "the search must hit several subjects"
    );
    for variables in [
        &[][..],
        &[("LOSAT_X_TBNOIDCHUNK", "1")],
        &[("LOSAT_X_DEFAULT", "0")],
        &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_TBNOIDCHUNK", "1")],
    ] {
        let (report, _) = tblastn(THREADS, variables);
        assert!(
            report == expected,
            "{variables:?} at {THREADS} threads differs from one thread"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409
// ```c
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
// ```
// 240 subjects give a chunk of two. Without the variable an entry takes `THREADS` subjects (60
// entries); with it, `THREADS` chunks (30 entries).
#[test]
fn the_switch_is_off_unless_set_and_halves_the_entries_here() {
    let debug = ("LOSAT_WASI_THREADS_DEBUG", "1");
    let (plain, plain_err) = tblastn(THREADS, &[debug]);
    let (chunked, chunked_err) = tblastn(THREADS, &[debug, ("LOSAT_X_TBNOIDCHUNK", "1")]);
    assert_eq!(pool_entries(&plain_err), SUBJECTS / THREADS);
    assert_eq!(
        pool_entries(&chunked_err),
        SUBJECTS / (THREADS * (SUBJECTS / 100))
    );
    assert!(plain == chunked, "the reports differ");
}
