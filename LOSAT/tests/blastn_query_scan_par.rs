//! The first reading of LOSAT_X_QSTREAM (each query record's length) read in parts at record starts
//! by the search's threads (`LOSAT_X_QSCANPAR`, `src/algorithm/blastn/blast_engine/x_query_scan_par.rs`)
//! on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless its variable is `0`.
//! It changes who reads the lengths, not the lengths, so the property under test is that nothing a
//! run writes depends on it (exit status, report and standard error), and that the default took its
//! path where its rule says so (`LOSAT_X_REACH`): LOSAT_X_QSTREAM reading the lengths first, two
//! threads or more, a file of `LOSAT_X_QSCANPARMIN` bytes or more (0 here, so that a small file
//! is cut), no `\r` and no reader error. The searches keep LOSAT_X_QSTREAM acting: small batches
//! (`BATCH_SIZE`) and no LOSAT_X_QBATCHPAR slots (whose window would hold the whole file). Cases:
//! 240 pieces of two genomes (lower case, `N` runs, a line with digits that the reader warns
//! about, records without letters and without ids) with dc-megablast and blastn against a subject
//! file (`-subject`) at one, two and four threads and against a database (`-db`) at four; the same
//! file with `\r\n` line ends, with a gap line (`>?`) and with a reader error in its last record,
//! where the default keeps the one-thread reading. That those bytes equal NCBI's is the business of the
//! fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_QSCANPAR";
const SWITCH_OFF: (&str, &str) = ("LOSAT_X_QSCANPAR", "0");
/// Every file is cut, whatever its size (tuning).
const ANY_SIZE: (&str, &str) = ("LOSAT_X_QSCANPARMIN", "0");
/// No slots: LOSAT_X_QSTREAM's rule keeps it off where they would read most of the file ahead.
const NO_SLOTS: (&str, &str) = ("LOSAT_X_QBATCHPAR", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A directory of its own per test and process (the tests of this binary run in parallel).
struct TestDir(PathBuf);

impl TestDir {
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("losat-qscan-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("make the test directory");
        TestDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
}

/// The letters of the first record of a FASTA file of the fixtures.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .take_while(|line| !line.starts_with(b">"))
        .flat_map(|line| line.iter().copied().filter(|b| !b.is_ascii_whitespace()))
        .collect()
}

/// 240 pieces of LC738874 and Sakai of 200 to 2,600 letters in lines of 50 to 80, every fifth with
/// a run of `N` and lower case letters; a line with digits (the reader's warning), a record without
/// letters and one without an id; `tail` appended.
fn query_text(tail: &[u8]) -> Vec<u8> {
    let genomes = [letters("LC738874.fasta"), letters("Sakai.fna")];
    let mut text = Vec::new();
    for i in 0..240usize {
        let genome = &genomes[i % 2];
        let len = 200 + (i * 977) % 2_400;
        let start = (i * 7_919) % (genome.len() - len);
        let mut seq = genome[start..start + len].to_vec();
        if i % 5 == 0 {
            let at = len / 3;
            seq.splice(at..at, std::iter::repeat_n(b'N', 30));
            for byte in &mut seq[at + 40..at + 120] {
                *byte = byte.to_ascii_lowercase();
            }
        }
        match i {
            17 => text.extend_from_slice(b">\n"),
            _ => text.extend_from_slice(format!(">piece{i} words of piece {i}\n").as_bytes()),
        }
        for line in seq.chunks(50 + i % 31) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        match i {
            33 => text.extend_from_slice(b"ACGT12ACGT\n"),
            90 => text.extend_from_slice(b">without letters\n"),
            _ => {}
        }
    }
    text.extend_from_slice(tail);
    text
}

/// One search of `query` with `target` (`-subject FILE` or `-db NAME`) and `task`, `-outfmt 6`,
/// `BATCH_SIZE` 20,000 and the given variables as its only `LOSAT_*` variables.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature (`cargo
/// test --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing (it is set to a file of this test when the environment has none).
/// Returns the exit status, the report and the standard error.
fn search(
    query: &Path,
    target: [&str; 2],
    task: &str,
    threads: u32,
    vars: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, Vec<u8>) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!("losat_qscan_{}_worker.log", std::process::id())),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(query)
        .args(target)
        .args(["-outfmt", "6"])
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "20000")
        .envs([NO_SLOTS, ANY_SIZE])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    (output.status.code(), output.stdout, output.stderr)
}

/// For `query` against `target` with `task`: the run with the switch off at four threads, then by
/// default at `threads`, the same exit status, report and standard error; the default took the
/// switch's path (`own`) at two threads or more, the one-thread reading at one.
fn check(query: &Path, target: [&str; 2], task: &str, threads: &[u32], own: bool) {
    let expected = search(query, target, task, 4, &[SWITCH_OFF]);
    for &n in threads {
        let reach = reach_support::Reach::new("qscan");
        let run = search(query, target, task, n, &[reach.var()]);
        let what = format!("{task} {target:?} threads {n}");
        assert_eq!(run.0, expected.0, "{what}: exit status");
        assert!(run.1 == expected.1, "{what}: the report differs");
        assert_eq!(
            String::from_utf8_lossy(&run.2),
            String::from_utf8_lossy(&expected.2),
            "{what}: standard error"
        );
        // The parts are read with `read_at` of `std::os::unix::fs::FileExt` and need the `parallel`
        // feature: elsewhere (Windows, WebAssembly) the switch never acts and the sides stay
        // `(false, true)` (module doc of `x_query_scan_par.rs`).
        let took = own && n >= 2 && cfg!(all(unix, feature = "parallel"));
        assert_eq!(
            reach.sides(SWITCH),
            (took, !took),
            "{what}: sides of {SWITCH}"
        );
    }
}

#[test]
fn the_parts_give_the_reports_of_the_one_thread_reading() {
    let dir = TestDir::new("subject");
    let query = dir.path("q.fa");
    std::fs::write(&query, query_text(b"")).expect("write the queries");
    let subject = fasta("LC738868.fasta");
    let target = ["-subject", subject.to_str().unwrap()];
    check(&query, target, "dc-megablast", &[1, 2, 4], true);
    check(&query, target, "blastn", &[1, 4], true);
}

#[test]
fn a_database_search_gives_the_same_bytes() {
    let dir = TestDir::new("db");
    let query = dir.path("q.fa");
    std::fs::write(&query, query_text(b"")).expect("write the queries");
    let mut set = Vec::new();
    for name in ["LC738868", "LC738869", "LC738874"] {
        set.extend(std::fs::read(fasta(&format!("{name}.fasta"))).unwrap());
    }
    std::fs::write(dir.path("set.fna"), set).unwrap();
    let db = dir.path("set").to_str().unwrap().to_string();
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(dir.path("set.fna"))
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(output.status.success(), "makedb");
    check(&query, ["-db", &db], "dc-megablast", &[4], true);
}

#[test]
fn a_file_with_cr_a_gap_line_or_a_reader_error_keeps_the_one_thread_reading() {
    let dir = TestDir::new("fallback");
    let subject = fasta("LC738868.fasta");
    let target = ["-subject", subject.to_str().unwrap()];
    let crlf: Vec<u8> = query_text(b"")
        .iter()
        .flat_map(|&byte| {
            if byte == b'\n' {
                vec![b'\r', b'\n']
            } else {
                vec![byte]
            }
        })
        .collect();
    let query = dir.path("crlf.fa");
    std::fs::write(&query, crlf).expect("write the queries");
    check(&query, target, "dc-megablast", &[4], false);
    // A gap line, whose length the reader makes `N`s of.
    let query = dir.path("gap.fa");
    std::fs::write(
        &query,
        query_text(b">gapped\nACGTACGTAC\n>?40\nACGTACGTAC\n"),
    )
    .expect("write the queries");
    check(&query, target, "dc-megablast", &[4], false);
    // A line of mostly invalid letters ends the reading with an error in the last record (an
    // exit that writes no reach file: the bytes alone).
    let query = dir.path("error.fa");
    std::fs::write(&query, query_text(b">bad\n%%%%%%%%%%%%%%%%\n")).expect("write the queries");
    let expected = search(&query, target, "dc-megablast", 4, &[SWITCH_OFF]);
    assert_ne!(
        expected.0,
        Some(0),
        "the error file: {}",
        String::from_utf8_lossy(&expected.2)
    );
    assert_eq!(search(&query, target, "dc-megablast", 4, &[]), expected);
}
