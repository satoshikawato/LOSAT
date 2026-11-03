//! LOSAT_X_BNQPARTS on whole searches: the store of LOSAT_X_BNQPACK held in parts
//! (`src/algorithm/blastn/blast_engine/x_query_drop.rs`, `Part`), read on the search's threads where
//! the query file can be cut at record starts (`x_query_scan_par.rs`, `store_parts`), else by the
//! search's reader in parts of about 1 MiB (`read_parts`), and the batch's records made again from
//! their parts when the loop takes them.
//!
//! The switch changes where the records wait, not the records, so the property under test is that
//! nothing a run writes depends on it or on the thread count (exit status, report and standard
//! error, against `LOSAT_X_BNQPARTS=0` at one thread: the store as before), and that the parts were
//! taken where they should be: the `[X_STATS] BNQPARTS records= parts= in_parts=` line (printed with
//! `LOSAT_X_STATS`) says `in_parts=1` at two threads or more (unix, `parallel` feature) and
//! `in_parts=0` at one thread and where the reader reads the file (a reader message after the first
//! part, `-parse_deflines`). `LOSAT_X_QSCANPARMIN=0` (tuning) lets the parts cut a small file;
//! `BATCH_SIZE` 20,000 makes batches that cross the parts. Cases: megablast of 900 pieces of two
//! genomes (lower case, `N` runs, a record without letters, one without an id; about 1.3 MB) against
//! a subject file and a database, outfmt 6
//! and 7, with `-lcase_masking`, with `-parse_deflines` (the ids the reader parsed are kept beside
//! the parts), and the same file with a hyphen in its last record (a reader message after the first
//! part). `LOSAT_X_QSTREAM=0` keeps the records read before the search (the file is longer than
//! twice its batch). That those bytes equal NCBI's is the business of the fixture suites.

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
        let dir = std::env::temp_dir().join(format!("losat-qparts-{test}-{}", std::process::id()));
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

/// 900 pieces of LC738874 and Sakai of 200 to 2,600 letters in lines of 50 to 80, every fifth with
/// a run of `N` and lower case letters; a record without an id and one without letters; `tail`
/// appended.
fn query_text(tail: &[u8]) -> Vec<u8> {
    let genomes = [letters("LC738874.fasta"), letters("Sakai.fna")];
    let mut text = Vec::new();
    for i in 0..900usize {
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
        if i == 90 {
            text.extend_from_slice(b">without letters\n");
        }
    }
    text.extend_from_slice(tail);
    text
}

/// One megablast search of `query` with `target` (`-subject FILE` or `-db NAME`) and `args`,
/// `BATCH_SIZE` 20,000, `LOSAT_X_QSCANPARMIN=0`, `LOSAT_X_QSTREAM=0` and the given variables as its
/// only `LOSAT_*` variables. `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the
/// `blastx-worker-probe` feature (`cargo test --all-features`, as in CI) needs that file name as
/// soon as a search has more than one thread, and it selects nothing. Returns the exit status, the
/// report and the standard error.
fn search(
    query: &Path,
    target: [&str; 2],
    args: &[&str],
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
            std::env::temp_dir().join(format!("losat_qparts_{}_worker.log", std::process::id())),
        );
    }
    let output = command
        .arg("blastn")
        .args(["-task", "megablast"])
        .arg("-query")
        .arg(query)
        .args(target)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", "20000")
        .env("LOSAT_X_QSCANPARMIN", "0")
        .env("LOSAT_X_QSTREAM", "0")
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    (output.status.code(), output.stdout, output.stderr)
}

/// The `[X_STATS] BNQPARTS` line of `stderr` (its `records=`, `parts=` and `in_parts=`), and
/// `stderr` without the `[X_STATS]` lines.
fn take_stats(stderr: &[u8]) -> ([usize; 3], String) {
    let text = String::from_utf8_lossy(stderr).into_owned();
    let lines: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with("[X_STATS] BNQPARTS "))
        .collect();
    assert_eq!(lines.len(), 1, "one BNQPARTS line:\n{text}");
    let value = |field: &str| -> usize {
        lines[0]
            .split_whitespace()
            .find_map(|word| word.strip_prefix(field))
            .unwrap_or_else(|| panic!("no {field} in: {}", lines[0]))
            .parse()
            .expect("a number")
    };
    let rest = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_STATS] "))
        .collect();
    (
        [value("records="), value("parts="), value("in_parts=")],
        rest,
    )
}

/// For `query` against `target` with `args`: the run at one thread with `LOSAT_X_BNQPARTS=0`, then
/// the runs with the switch at `threads` with `LOSAT_X_STATS`, the same exit status, report and
/// standard error (without the stats lines); the store held `records` records, in parts read on
/// the search's threads at two threads or more where `cut`, else read by the search's reader.
fn check(query: &Path, target: [&str; 2], args: &[&str], records: usize, cut: bool) {
    let expected = search(query, target, args, 1, &[("LOSAT_X_BNQPARTS", "0")]);
    assert_eq!(
        expected.0,
        Some(0),
        "{}",
        String::from_utf8_lossy(&expected.2)
    );
    assert!(!expected.1.is_empty(), "the reference report has hits");
    for n in [1, 2, 4] {
        let run = search(query, target, args, n, &[("LOSAT_X_STATS", "1")]);
        let what = format!("{args:?} {target:?} threads {n}");
        let ([held, parts, in_parts], stderr) = take_stats(&run.2);
        assert_eq!(run.0, expected.0, "{what}: exit status");
        assert!(run.1 == expected.1, "{what}: the report differs");
        assert_eq!(
            stderr,
            String::from_utf8_lossy(&expected.2),
            "{what}: standard error"
        );
        assert_eq!(held, records, "{what}: records");
        // The parts are read with `read_at` of `std::os::unix::fs::FileExt` and need the `parallel`
        // feature: elsewhere (Windows, WebAssembly) the search's reader reads them.
        let took = cut && n >= 2 && cfg!(all(unix, feature = "parallel"));
        assert_eq!(in_parts, usize::from(took), "{what}: in_parts");
        assert!(parts >= if took { 2 } else { 1 }, "{what}: {parts} parts");
    }
}

#[test]
fn the_parts_give_the_reports_of_the_store() {
    let dir = TestDir::new("subject");
    let query = dir.path("q.fa");
    std::fs::write(&query, query_text(b"")).expect("write the queries");
    let subject = fasta("LC738868.fasta");
    let target = ["-subject", subject.to_str().unwrap()];
    check(&query, target, &["-outfmt", "6"], 901, true);
    check(
        &query,
        target,
        &["-outfmt", "7", "-lcase_masking"],
        901,
        true,
    );
    // The ids the reader parsed wait beside the parts; the search's reader reads the file.
    check(
        &query,
        target,
        &["-outfmt", "6", "-parse_deflines"],
        901,
        false,
    );
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
    check(&query, ["-db", &db], &["-outfmt", "6"], 901, true);
}

#[test]
fn a_reader_message_after_the_first_part_leaves_the_file_to_the_reader() {
    let dir = TestDir::new("message");
    let query = dir.path("hyphen.fa");
    std::fs::write(&query, query_text(b">hyphen\nACGTACGT-ACGTACGTAC\n")).expect("write");
    let subject = fasta("LC738868.fasta");
    let target = ["-subject", subject.to_str().unwrap()];
    check(&query, target, &["-outfmt", "6"], 902, false);
    let run = search(&query, target, &["-outfmt", "6"], 4, &[]);
    assert!(
        String::from_utf8_lossy(&run.2).contains("Hyphens are invalid"),
        "the reader's message: {}",
        String::from_utf8_lossy(&run.2)
    );
}
