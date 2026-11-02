//! The first reading of LOSAT_X_QSTREAM (each query record's length) read in parts at record starts
//! by the search's threads (`src/algorithm/blastn/blast_engine/x_query_scan_par.rs`) on whole
//! searches.
//!
//! The parts have no switch of their own: wherever LOSAT_X_QSTREAM acts, at two threads or more on
//! unix with the `parallel` feature, the first reading is in parts, and it falls back to one
//! thread on a `\r`, a gap line (`>?`), a reader error and the like. (`LOSAT_X_QSCANPAR` was folded
//! into LOSAT_X_QSTREAM on 2026-11-01, its minimum of 4 MiB removed on 2026-11-02.) It changes who
//! reads the lengths, not the lengths, so the property under test is that nothing a run writes
//! depends on the thread count (exit status, report and standard error), and that the parts were
//! taken where the rule says so: `scan_par=1` in the `[X_STATS] QSTREAM` line (printed with
//! `LOSAT_X_STATS`), `scan_par=0` at one thread and on the fallback files. The searches keep
//! LOSAT_X_QSTREAM acting: small batches (`BATCH_SIZE`; the query files are far above two batches).
//! Cases:
//! 240 pieces of two genomes (lower case, `N` runs, a line with digits that the reader warns
//! about, records without letters and without ids) with dc-megablast and blastn against a subject
//! file (`-subject`) at one, two and four threads and against a database (`-db`) at four; the same
//! file with `\r\n` line ends, with a gap line (`>?`) and with a reader error in its last record,
//! where the one-thread reading is kept. That those bytes equal NCBI's is the business of the
//! fixture suites.

use std::path::{Path, PathBuf};
use std::process::Command;

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
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    (output.status.code(), output.stdout, output.stderr)
}

/// The value of `field` (`scan_par=`) in the one `[X_STATS] QSTREAM` line of `stderr`, and
/// `stderr` without the `[X_STATS]` lines.
fn take_stats(stderr: &[u8], field: &str) -> (String, String) {
    let text = String::from_utf8_lossy(stderr).into_owned();
    let lines: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("[X_STATS] QSTREAM "))
        .collect();
    assert_eq!(lines.len(), 1, "one QSTREAM line:\n{text}");
    assert!(
        lines[0].ends_with(" off=none"),
        "QSTREAM acted: {}",
        lines[0]
    );
    let value = lines[0]
        .split_whitespace()
        .find_map(|word| word.strip_prefix(field))
        .unwrap_or_else(|| panic!("no {field} in: {}", lines[0]))
        .to_string();
    let rest = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_STATS] "))
        .collect();
    (value, rest)
}

/// For `query` against `target` with `task`: the run at one thread, then the runs at `threads`
/// with `LOSAT_X_STATS`, the same exit status, report and standard error (without the stats
/// lines); the first reading was in parts (`scan_par=1`) at two threads or more where `own`, else
/// not.
fn check(query: &Path, target: [&str; 2], task: &str, threads: &[u32], own: bool) {
    let expected = search(query, target, task, 1, &[]);
    for &n in threads {
        let run = search(query, target, task, n, &[("LOSAT_X_STATS", "1")]);
        let what = format!("{task} {target:?} threads {n}");
        let (scan_par, stderr) = take_stats(&run.2, "scan_par=");
        assert_eq!(run.0, expected.0, "{what}: exit status");
        assert!(run.1 == expected.1, "{what}: the report differs");
        assert_eq!(
            stderr,
            String::from_utf8_lossy(&expected.2),
            "{what}: standard error"
        );
        // The parts are read with `read_at` of `std::os::unix::fs::FileExt` and need the `parallel`
        // feature: elsewhere (Windows, WebAssembly) the first reading stays on one thread
        // (module doc of `x_query_scan_par.rs`).
        let took = own && n >= 2 && cfg!(all(unix, feature = "parallel"));
        assert_eq!(scan_par, if took { "1" } else { "0" }, "{what}: scan_par");
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
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.path("set.fna"))
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(output.status.success(), "makeblastdb");
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
    let expected = search(&query, target, "dc-megablast", 1, &[]);
    assert_ne!(
        expected.0,
        Some(0),
        "the error file: {}",
        String::from_utf8_lossy(&expected.2)
    );
    assert_eq!(search(&query, target, "dc-megablast", 4, &[]), expected);
}
