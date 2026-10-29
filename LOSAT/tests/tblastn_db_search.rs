//! TBLASTN's search of a database (`tblastn -db NAME`, stage DB-5 of the database search
//! campaign) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the report of a database search is the same
//! at every thread count; for a subject set without ambiguity codes whose shortest record is
//! 10 letters or fewer, the tabular rows of `-db` are those of `-subject` with the same
//! records (NCBI scans a subject set as a database; the statistics differ only where a set's
//! shortest record is longer than the 10 letters a database answers, and the bases only
//! where makeblastdb has replaced an ambiguity code by a random base); and a database keeps
//! the weak hits that the shortest-subject cutoff of a FASTA set of long records drops. That
//! the bytes equal NCBI BLAST+ `tblastn -db` is the business of the campaign's comparison runs.

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
        "losat-tblastn-db-search-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The database of `fasta_file` made by `losat makeblastdb` in `dir` under `name`.
fn makedb(dir: &Path, fasta_file: &Path, name: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

/// One search of `query` on standard output; returns the report.
fn tblastn(query: &Path, subjects: &[&str], outfmt: &str, threads: u32) -> Vec<u8> {
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .args(subjects)
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{subjects:?} {outfmt} threads {threads}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn rows(report: &[u8]) -> usize {
    report.iter().filter(|&&byte| byte == b'\n').count()
}

#[test]
fn a_database_search_gives_the_same_report_at_every_thread_count() {
    let dir = temp_dir("threads");
    // The first 30 proteins of PemoMJNVA against the MelaMJNV genome (no ambiguity code) and a
    // record of 10 letters, so that the set's shortest record is a database's.
    let proteins = std::fs::read_to_string(fasta("PemoMJNVA.faa")).unwrap();
    let records: Vec<&str> = proteins.split_inclusive('>').collect();
    let query = dir.join("q30.faa");
    std::fs::write(
        &query,
        format!(">{}", records[1..=30].concat().trim_end_matches('>')),
    )
    .unwrap();
    let set = dir.join("set.fna");
    let mut bytes = std::fs::read(fasta("MelaMJNV.fasta")).unwrap();
    bytes.extend_from_slice(b">ten ten bases\nACGTTGCAAC\n");
    std::fs::write(&set, bytes).unwrap();
    let db = makedb(&dir, &set, "set");
    for outfmt in ["6", "0", "7"] {
        let one = tblastn(&query, &["-db", &db], outfmt, 1);
        assert!(!one.is_empty(), "outfmt {outfmt}");
        for threads in [2, 4] {
            assert!(
                tblastn(&query, &["-db", &db], outfmt, threads) == one,
                "outfmt {outfmt}: {threads} threads differ from 1"
            );
        }
        if outfmt == "6" {
            assert!(rows(&one) > 30, "{}", rows(&one));
            let subject = tblastn(&query, &["-subject", set.to_str().unwrap()], outfmt, 1);
            assert!(subject == one, "-db and -subject rows differ");
        }
    }
    // The reports name the database: outfmt 7 by the `-db` value, outfmt 0 by its title
    // (makedb's default title: the input file's name) and its numbers.
    let comments = String::from_utf8(tblastn(&query, &["-db", &db], "7", 1)).unwrap();
    assert!(
        comments.contains(&format!("\n# Database: {db}\n")),
        "{comments}"
    );
    let pairwise = String::from_utf8(tblastn(&query, &["-db", &db], "0", 1)).unwrap();
    assert!(
        pairwise.contains("Database: ") && pairwise.contains("2 sequences;"),
        "{pairwise}"
    );
    assert!(
        !pairwise.contains("User specified sequence set"),
        "{pairwise}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_setup.c:968-975 with seqdbimpl.cpp:131-136: a database answers 10 for its shortest
// sequence, so the hit-saving cutoff of `-db` is that of a 3-residue subject; a FASTA set
// answers its shortest record (the MelaMJNV genome). bench5's tblastn pair: NCBI BLAST+ 2.17.0
// keeps 548 rows with `-db` and 505 with `-subject` (dbsearch-20261012 DESIGN.md section 1.2).
#[test]
fn a_database_keeps_the_weak_hits_of_its_shortest_length() {
    let dir = temp_dir("d1");
    let db = makedb(&dir, &fasta("MelaMJNV.fasta"), "mela");
    let query = fasta("PemoMJNVA.faa");
    let database = tblastn(&query, &["-db", &db], "6", 1);
    let subject = tblastn(
        &query,
        &["-subject", fasta("MelaMJNV.fasta").to_str().unwrap()],
        "6",
        1,
    );
    assert_eq!((rows(&database), rows(&subject)), (548, 505));
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one. A database without
// its `.nsq` answers a query of X letters (no valid context) as NCBI BLAST+ 2.17.0 does: a
// report without hits and the query's warning, exit status 0; with a second query in a batch
// of its own (`BATCH_SIZE=1`), that batch fails after the report of the first one.
#[test]
fn the_sequence_file_is_read_by_the_first_batch_with_a_valid_context() {
    let dir = temp_dir("nsq");
    let db = makedb(&dir, &fasta("MelaMJNV.fasta"), "mela");
    std::fs::remove_file(format!("{db}.nsq")).unwrap();
    let proteins = std::fs::read_to_string(fasta("PemoMJNVA.faa")).unwrap();
    let first = proteins.split_inclusive('>').nth(1).unwrap();
    let all_x = format!(">xx all X\n{}\n", "X".repeat(32));
    let query = dir.join("xx.faa");
    std::fs::write(&query, &all_x).unwrap();
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(&query)
        .args(["-db", &db])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("Query= xx all X") && report.contains("***** No hits found *****"),
        "{report}"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Could not calculate ungapped Karlin-Altschul parameters"),
        "{stderr}"
    );
    std::fs::write(&query, format!("{all_x}>{}", first.trim_end_matches('>'))).unwrap();
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(&query)
        .args(["-db", &db])
        .env("BATCH_SIZE", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(3));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("Query= xx all X") && report.matches("Query= ").count() == 1,
        "{report}"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.ends_with("BLAST engine error: Database memory map file error\n"),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_args.cpp:2364-2381: `-db` excludes `-subject` and `-subject_loc` (LOSAT's parser
// reports it); blast_app_util.hpp:185-203: a database that is not found is NCBI's
// `BLAST Database error`, exit status 2.
#[test]
fn database_arguments_are_checked_as_ncbi_checks_them() {
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(fasta("PemoMJNVA.faa"))
        .args(["-db", "x", "-subject"])
        .arg(fasta("MelaMJNV.fasta"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let dir = temp_dir("missing");
    let missing = dir.join("no_such_db");
    let output = losat()
        .arg("tblastn")
        .arg("-query")
        .arg(fasta("PemoMJNVA.faa"))
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
