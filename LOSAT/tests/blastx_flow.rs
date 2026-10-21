//! BLASTX's query batches searched as a flow of tasks (`LOSAT_X_BXFLOW`,
//! `x_flow` in `src/algorithm/blastx/native.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless `LOSAT_X_BXFLOW=0`); the reference runs here set it
//! to `0`. With two or more threads blastx's path without it queues short query batches (a batch is the queries up to
//! 10,002 nucleotides) and searches them in waves of `threads * 4`; with the switch every batch
//! is a task of its own (spawned once `threads` batches follow it; the last `threads` batches of
//! the input are searched with the pool inside, side by side) and the results, the reader's
//! warnings and the stop at a failing batch are written in the reader's order. So the property
//! under test is that nothing the run writes
//! depends on it: standard output (or the `-out` file), standard error (byte for byte) and the
//! exit status equal those of the same search with the switch off, for two, three and four
//! threads (the bound of `threads * 4` unwritten batches is reached at two and three), outfmt 0,
//! 6 and 7, `-max_target_seqs 1`, a database (`-db`), an input of fewer batches than threads
//! (every batch searched with the pool inside), and an error of the reader in a later batch after
//! the reports of the batches before it. That the bytes equal NCBI's is the business
//! of the BLASTX fixtures.
//!
//! The queries are 55 pieces of AP027131 of 1,500 to 3,900 letters each and, in the middle, one
//! of 30,000 letters (more than two batches, so it takes the "large" path, searched with the
//! pool inside its task): about 13 batches. The subject is the proteins of AP027131, so
//! the pieces have many hits. Later records carry the reader's warnings: hyphens, a title that
//! ends with 20 nucleotide letters, and invalid residues.

#![cfg(feature = "parallel")]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: (&str, &str) = ("LOSAT_X_BXFLOW", "1");
/// The switch off: the waves (the switch is in the strict set, so a run without it flows).
const OFF: (&str, &str) = ("LOSAT_X_BXFLOW", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of a one-record FASTA fixture.
fn genome(name: &str) -> Vec<u8> {
    let text = std::fs::read(fixture(name)).expect("read fixture");
    text.split(|&byte| byte == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A file name in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_flow_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A FASTA file in the temporary directory, removed when dropped.
struct TempFasta(PathBuf);

impl TempFasta {
    fn new(tag: &str, records: &[(String, Vec<u8>)]) -> Self {
        let path = temp_path(&format!("{tag}.fa"));
        let mut text = Vec::new();
        for (title, letters) in records {
            text.extend_from_slice(format!(">{title}\n").as_bytes());
            for line in letters.chunks(70) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        }
        std::fs::write(&path, text).expect("write temporary FASTA");
        TempFasta(path)
    }

    fn path(&self) -> &str {
        self.0.to_str().expect("UTF-8 temporary path")
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// 55 queries cut from AP027131, 1,500 to 3,900 letters each, and the 30,000-letter query
/// `long 30 kb` as record 25.
fn queries() -> Vec<(String, Vec<u8>)> {
    let letters = genome("AP027131.fasta");
    let mut records: Vec<(String, Vec<u8>)> = (0..55)
        .map(|index| {
            let start = index * 10_000;
            let length = 1_500 + (index * 373) % 2_400;
            (
                format!("q{index} piece {start}"),
                letters[start..start + length].to_vec(),
            )
        })
        .collect();
    records.insert(
        25,
        ("long 30 kb".to_string(), letters[300_000..330_000].to_vec()),
    );
    records
}

/// The proteins of AP027131 as the subject: the first `count` records of the fixture.
fn subject(count: usize) -> TempFasta {
    let text = std::fs::read_to_string(fixture("AP027131.faa")).expect("read fixture");
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(title) = line.strip_prefix('>') {
            records.push((title.to_string(), Vec::new()));
        } else if let Some(record) = records.last_mut() {
            record.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    records.truncate(count);
    TempFasta::new("subject", &records)
}

/// The queries with the reader's warnings in later batches: hyphens in a line of record 30, a
/// title ending with 24 nucleotide letters at record 36, and invalid residues (digits and `X`)
/// in lines of record 44.
fn queries_with_warnings() -> Vec<(String, Vec<u8>)> {
    let mut records = queries();
    records[30].1.insert(100, b'-');
    records[30].1.insert(101, b'-');
    records[36].0.push_str(" ACGTACGTTGCATGCAACGTACGT");
    records[44].1.splice(150..150, *b"123");
    records[44].1.splice(300..300, *b"XX");
    records
}

struct Run {
    stdout: Vec<u8>,
    stderr: Vec<u8>,
    code: Option<i32>,
    /// The `-out` file, if the arguments name one.
    file: Option<Vec<u8>>,
}

impl Run {
    fn stderr_text(&self) -> String {
        String::from_utf8_lossy(&self.stderr).into_owned()
    }

    /// The report: standard output, or the `-out` file.
    fn report(&self) -> &[u8] {
        self.file.as_deref().unwrap_or(&self.stdout)
    }
}

/// One `losat blastx` run with `args`, with `env` as its only variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` is a file of its own: a binary built with the `blastx-worker-probe`
/// feature (`cargo test --all-features`, as in CI) needs that file name as soon as a search has
/// more than one thread, and it selects nothing.
fn blastx(args: &[&str], env: &[(&str, &str)]) -> Run {
    let log = temp_path("worker.log");
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastx")
        .args(args)
        .env("LOSAT_BLASTX_WORKER_LOG", &log)
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&log);
    let file = args
        .iter()
        .position(|&arg| arg == "-out")
        .map(|at| std::fs::read(args[at + 1]).unwrap_or_default());
    Run {
        stdout: output.stdout,
        stderr: output.stderr,
        code: output.status.code(),
        file,
    }
}

/// The run with the switch gives what the run without it gives (the report, standard error
/// byte for byte, the exit status); returns the run without it. Unless `fails`, that run must
/// have succeeded and written a report, so that the comparison cannot pass vacuously.
fn same_as_without(what: &str, args: &[&str], fails: bool) -> Run {
    let reference = blastx(args, &[OFF]);
    if fails {
        assert_ne!(
            reference.code,
            Some(0),
            "{what}: the search must fail\n{}",
            reference.stderr_text()
        );
    } else {
        assert_eq!(
            reference.code,
            Some(0),
            "{what}: the reference search failed\n{}",
            reference.stderr_text()
        );
        assert!(
            !reference.report().is_empty(),
            "{what}: the reference search wrote no report"
        );
    }
    let run = blastx(args, &[SWITCH]);
    assert_eq!(
        run.code,
        reference.code,
        "{what}: exit status\n{}",
        run.stderr_text()
    );
    assert!(
        run.stdout == reference.stdout,
        "{what}: standard output differs ({} and {} bytes)",
        run.stdout.len(),
        reference.stdout.len()
    );
    assert!(run.file == reference.file, "{what}: the -out file differs");
    assert_eq!(
        run.stderr_text(),
        reference.stderr_text(),
        "{what}: standard error"
    );
    reference
}

/// The number of distinct query ids (first column) of a tabular report.
fn query_ids(report: &[u8]) -> usize {
    let mut ids: Vec<&[u8]> = report
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty() && !line.starts_with(b"#"))
        .filter_map(|line| line.split(|&byte| byte == b'\t').next())
        .collect();
    ids.sort_unstable();
    ids.dedup();
    ids.len()
}

fn assert_warnings(run: &Run) {
    let stderr = run.stderr_text();
    for expected in [
        "Hyphens are invalid and will be ignored",
        "Title ends with at least 20 valid nucleotide characters",
        "Ignoring invalid residues at position(s)",
    ] {
        assert!(stderr.contains(expected), "no {expected:?}:\n{stderr}");
    }
}

// Two, three and four threads (the bound of `threads * 4` unwritten batches is reached at two and
// three) with each outfmt, on queries with the 30 kb query in the middle and the reader's
// warnings in later batches.
#[test]
fn every_thread_count_and_outfmt_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &queries_with_warnings());
    let subject = subject(400);
    for (threads, outfmt) in [
        ("2", "6"),
        ("3", "6"),
        ("4", "6"),
        ("2", "7"),
        ("4", "7"),
        ("3", "0"),
    ] {
        let args = [
            "-query",
            query.path(),
            "-subject",
            subject.path(),
            "-outfmt",
            outfmt,
            "-num_threads",
            threads,
        ];
        let what = format!("threads {threads} outfmt {outfmt}");
        let run = same_as_without(&what, &args, false);
        assert_warnings(&run);
        if outfmt == "6" {
            assert!(query_ids(&run.stdout) >= 10, "{what}: too few queries hit");
        }
    }
}

// The same queries into an `-out` file (outfmt 0, whose prolog is written with the reports, and
// outfmt 6), at the thread count where the bound is reached.
#[test]
fn an_out_file_gives_the_report_without_the_switch() {
    let query = TempFasta::new("q", &queries_with_warnings());
    let subject = subject(400);
    for outfmt in ["0", "6"] {
        let out = temp_path("out.txt");
        let out_text = out.to_str().expect("UTF-8 path");
        let args = [
            "-query",
            query.path(),
            "-subject",
            subject.path(),
            "-outfmt",
            outfmt,
            "-num_threads",
            "2",
            "-out",
            out_text,
        ];
        let run = same_as_without(&format!("-out outfmt {outfmt}"), &args, false);
        assert!(run.stdout.is_empty());
        assert_warnings(&run);
        let _ = std::fs::remove_file(&out);
    }
}

// `-max_target_seqs 1` (and the default e-value spelt out) against all the proteins.
#[test]
fn one_target_sequence_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &queries_with_warnings());
    let subject = subject(400);
    for threads in ["2", "4"] {
        let args = [
            "-query",
            query.path(),
            "-subject",
            subject.path(),
            "-evalue",
            "1e-5",
            "-max_target_seqs",
            "1",
            "-outfmt",
            "6",
            "-num_threads",
            threads,
        ];
        let what = format!("max_target_seqs 1 threads {threads}");
        let run = same_as_without(&what, &args, false);
        // Queries still hit with one target each.
        assert!(query_ids(&run.stdout) >= 10, "{what}");
    }
}

// An error of the reader in a later batch: record 40 is `@@@@` (a line that is not plausible
// data), so the run ends there after the reports of the batches before it and the warnings
// before the record. Standard output, standard error and the exit status are those without the
// switch (the batches searched past the error are dropped).
#[test]
fn an_error_at_a_later_batch_comes_as_without_the_switch() {
    let subject = subject(400);
    let mut records = queries_with_warnings();
    records.insert(40, ("bad record".to_string(), b"@@@@@@@@".to_vec()));
    let query = TempFasta::new("bad", &records);
    for (threads, outfmt) in [("2", "6"), ("3", "6"), ("4", "0")] {
        let args = [
            "-query",
            query.path(),
            "-subject",
            subject.path(),
            "-outfmt",
            outfmt,
            "-num_threads",
            threads,
        ];
        let what = format!("error threads {threads} outfmt {outfmt}");
        let run = same_as_without(&what, &args, true);
        // The batches before the bad record were reported, with the warning of record 30.
        assert!(!run.stdout.is_empty(), "{what}: nothing reported");
        assert!(
            run.stderr_text().contains("Hyphens are invalid"),
            "{what}:\n{}",
            run.stderr_text()
        );
    }
}

// A protein database made by `losat makedb` as the subject, at three threads.
#[test]
fn a_database_gives_the_reports_without_the_switch() {
    let query = TempFasta::new("q", &queries_with_warnings());
    let dir = temp_path("db");
    std::fs::create_dir_all(&dir).expect("make the database directory");
    let name = dir.join("subject");
    let name = name.to_str().expect("UTF-8 path");
    let proteins = subject(400);
    let made = Command::new(env!("CARGO_BIN_EXE_losat"))
        .args([
            "makedb",
            "-in",
            proteins.path(),
            "-dbtype",
            "prot",
            "-out",
            name,
        ])
        .output()
        .expect("run losat makedb");
    assert!(
        made.status.success(),
        "{}",
        String::from_utf8_lossy(&made.stderr)
    );
    let args = [
        "-query",
        query.path(),
        "-db",
        name,
        "-outfmt",
        "6",
        "-num_threads",
        "3",
    ];
    let run = same_as_without("-db threads 3", &args, false);
    assert_warnings(&run);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn fewer_batches_than_threads_give_the_reports_without_the_switch() {
    // The first three pieces: one batch, and the first twelve: about three batches; all of them
    // are the input's last batches, searched with the pool inside.
    let subject = subject(400);
    for (count, threads) in [(3, "4"), (12, "4"), (12, "2")] {
        let records: Vec<_> = queries().into_iter().take(count).collect();
        let query = TempFasta::new("few", &records);
        let args = [
            "-query",
            query.path(),
            "-subject",
            subject.path(),
            "-outfmt",
            "6",
            "-num_threads",
            threads,
        ];
        same_as_without(&format!("{count} queries, threads {threads}"), &args, false);
    }
}
