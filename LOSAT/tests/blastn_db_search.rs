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

/// The subject set (seven genomes without ambiguity codes, 2.1 Mb) as one FASTA file, and
/// its database made by `losat makedb`.
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
        .arg("makedb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
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
    // (makedb's default title: the input file's name) and its numbers.
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

// blast_args.cpp:2364-2381: `-db` excludes `-subject` and `-subject_loc` (LOSAT's parser
// reports it); blast_app_util.hpp:185-203: a database that is not found is NCBI's
// `BLAST Database error`, exit status 2.
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
            "BLAST Database error: No alias or index file found for nucleotide database ["
        ),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
