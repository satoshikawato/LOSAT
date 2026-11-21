//! BLASTN's search of a database (`blastn -db NAME`, stage DB-4 of the database search
//! campaign) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the report of a database search is the same
//! at every thread count (the database is read in blocks searched in parallel and merged in
//! the order of the ordinal ids, `src/algorithm/blastn/blast_engine/database.rs`), and for
//! a set of sequences without ambiguity codes the tabular rows of `-db` are those of
//! `-subject` with the same records (NCBI scans a subject set as a database; the two differ
//! only where makeblastdb has replaced an ambiguity code by a random base). That the bytes
//! equal NCBI BLAST+ `blastn -db` is the business of the campaign's comparison runs.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// A temporary directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastn-db-search-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The letters of the first record of a FASTA file.
fn letters(path: &Path) -> String {
    let text = std::fs::read_to_string(path).expect("read the FASTA file");
    text.lines()
        .skip_while(|line| !line.starts_with('>'))
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect()
}

/// The subject set (seven genomes without ambiguity codes, 2.1 Mb) as one FASTA file, and
/// its database made by `losat makeblastdb`.
fn subject_set(dir: &Path) -> (PathBuf, String) {
    let set = dir.join("set.fna");
    let mut bytes = Vec::new();
    for name in [
        "LC738868", "LC738869", "LC738870", "LC738871", "LC738873", "LC738874", "LC738875",
    ] {
        bytes.extend(std::fs::read(fasta(&format!("{name}.fasta"))).unwrap());
    }
    std::fs::write(&set, bytes).unwrap();
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (set, db)
}

/// One search on standard output; returns the report.
fn blastn(subjects: &[&str], task: &str, outfmt: &str, threads: u32) -> Vec<u8> {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .args(subjects)
        .args(["-task", task, "-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{subjects:?} {task} {outfmt} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

#[test]
fn a_database_search_gives_the_same_report_at_every_thread_count() {
    let dir = temp_dir("threads");
    let (set, db) = subject_set(&dir);
    for (task, outfmt) in [
        ("megablast", "6"),
        ("megablast", "0"),
        ("megablast", "7"),
        ("dc-megablast", "6"),
    ] {
        let one = blastn(&["-db", &db], task, outfmt, 1);
        assert!(!one.is_empty(), "{task} {outfmt}");
        for threads in [2, 4] {
            assert!(
                blastn(&["-db", &db], task, outfmt, threads) == one,
                "{task} outfmt {outfmt}: {threads} threads differ from 1"
            );
        }
        if outfmt == "6" {
            // No ambiguity code: the rows of the subject set.
            let subject = blastn(&["-subject", set.to_str().unwrap()], task, outfmt, 1);
            assert!(subject == one, "{task}: -db and -subject rows differ");
        }
    }
    // The reports name the database: outfmt 7 by the `-db` value, outfmt 0 by its title
    // (makeblastdb's default title: the input file's name) and its numbers.
    let comments = String::from_utf8(blastn(&["-db", &db], "megablast", "7", 1)).unwrap();
    assert!(
        comments.contains(&format!("\n# Database: {db}\n")),
        "{comments}"
    );
    let pairwise = String::from_utf8(blastn(&["-db", &db], "megablast", "0", 1)).unwrap();
    assert!(
        pairwise.contains("Database: ") && pairwise.contains("7 sequences;"),
        "{pairwise}"
    );
    assert!(!pairwise.contains("Posted date:  Unknown"), "{pairwise}");
    std::fs::remove_dir_all(&dir).unwrap();
}

// LOSAT_X_BNDECODE (`src/seqdb/volume.rs` `x_blastna_by_byte`): a database sequence's BLASTNA
// letters made four bases per stored byte give the reports of the decode one base at a time
// (`LOSAT_X_BNDECODE=0`), for the whole decode of the preliminary search and, with
// `LOSAT_X_BNTBRANGE=0`, of the traceback, and for the traceback's decode of ranges (the default
// for a sequence of 1 MiB or more). The database is EDL933 (5.5 Mb, 2,466 ambiguity runs of 11
// codes, as makeblastdb stores them); the queries are its letters around 24 of those runs, the
// runs written as `A`, so the alignments cross them and `sseq` prints their codes.
#[test]
fn the_by_byte_decode_gives_the_reports_of_the_decode_by_base() {
    let dir = temp_dir("decode");
    let genome = fasta("EDL933.fna");
    let db = dir.join("edl933").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(&genome)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let letters = letters(&genome);
    let runs: Vec<usize> = letters
        .bytes()
        .enumerate()
        .filter(|&(_, letter)| !b"ACGT".contains(&letter.to_ascii_uppercase()))
        .map(|(at, _)| at)
        .collect();
    let mut queries = String::new();
    let mut last = 0usize;
    for &at in &runs {
        if queries.matches('>').count() == 24 {
            break;
        }
        if at < last + 20_000 || at < 400 || at + 400 > letters.len() {
            continue;
        }
        last = at;
        let window: String = letters[at - 397..at + 403]
            .chars()
            .map(|c| {
                if "ACGT".contains(c.to_ascii_uppercase()) {
                    c
                } else {
                    'A'
                }
            })
            .collect();
        queries.push_str(&format!(">q{at}\n{window}\n"));
    }
    let query = dir.join("queries.fna");
    std::fs::write(&query, queries).unwrap();
    let search = |vars: &[(&str, &str)]| -> Vec<u8> {
        let output = losat()
            .arg("blastn")
            .arg("-query")
            .arg(&query)
            .args(["-db", &db, "-outfmt", "6 std sseq", "-num_threads", "1"])
            .envs(vars.iter().copied())
            .output()
            .expect("run LOSAT");
        assert!(
            output.status.success() && output.stderr.is_empty(),
            "{vars:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    };
    for whole in [&[][..], &[("LOSAT_X_BNTBRANGE", "0")][..]] {
        let by_byte = search(whole);
        let rows = String::from_utf8(by_byte.clone()).unwrap();
        assert!(rows.lines().count() >= 24, "{whole:?}: {rows}");
        let sseq: String = rows
            .lines()
            .map(|row| row.rsplit('\t').next().unwrap())
            .collect();
        assert!(
            sseq.bytes().any(|letter| !b"ACGT-".contains(&letter)),
            "{whole:?}: no ambiguity code in sseq"
        );
        let mut off = whole.to_vec();
        off.push(("LOSAT_X_BNDECODE", "0"));
        assert!(
            search(&off) == by_byte,
            "{whole:?}: LOSAT_X_BNDECODE=0 differs"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_args.cpp:2364-2381: `-db` excludes `-subject` and `-subject_loc` (LOSAT's parser
// reports it); blast_app_util.hpp:185-203: a database that is not found is NCBI's
// `LOSAT database error` (NCBI's `BLAST Database error`), exit status 2.
#[test]
fn database_arguments_are_checked_as_ncbi_checks_them() {
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .args(["-db", "x", "-subject"])
        .arg(fasta("LC738874.fasta"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let dir = temp_dir("missing");
    let missing = dir.join("no_such_db");
    let output = losat()
        .arg("blastn")
        .arg("-query")
        .arg(fasta("LC738874.fasta"))
        .arg("-db")
        .arg(&missing)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with(
            "LOSAT database error: No alias or index file found for nucleotide database ["
        ),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
