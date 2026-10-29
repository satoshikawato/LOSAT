//! BLASTX's search of a database (`blastx -db NAME`, stage DB-5 of the database search
//! campaign) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the report of a database search is the same
//! at every thread count; for a subject set whose shortest record is 10 letters or fewer,
//! the tabular rows of `-db` are those of `-subject` with the same records (NCBI scans a
//! subject set as a database; the statistics differ only where a set's shortest record is
//! longer than the 10 letters a database answers); and a database keeps the weak hits that
//! the shortest-subject cutoff of a FASTA set of long records drops. That the bytes equal
//! NCBI BLAST+ `blastx -db` is the business of the campaign's comparison runs.

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
        "losat-blastx-db-search-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The protein database of `fasta_file` made by `losat makeblastdb` in `dir` under `name`.
fn makedb(dir: &Path, fasta_file: &Path, name: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "prot", "-out", &db])
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
fn blastx(query: &Path, subjects: &[&str], outfmt: &str, threads: u32) -> Vec<u8> {
    let output = losat()
        .arg("blastx")
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

/// The first `bases` bases of the MelaMJNV genome as one query.
fn genome_slice(dir: &Path, bases: usize) -> PathBuf {
    let genome = std::fs::read_to_string(fasta("MelaMJNV.fasta")).unwrap();
    let letters: String = genome
        .lines()
        .skip(1)
        .flat_map(str::chars)
        .take(bases)
        .collect();
    let query = dir.join("slice.fna");
    std::fs::write(&query, format!(">slice MelaMJNV slice\n{letters}\n")).unwrap();
    query
}

#[test]
fn a_database_search_gives_the_same_report_at_every_thread_count() {
    let dir = temp_dir("threads");
    // A 30,000-base slice of the MelaMJNV genome (three query batches' worth of chunks)
    // against the PemoMJNVA proteins and a record of 10 letters, so that the set's shortest
    // record is a database's.
    let query = genome_slice(&dir, 30_000);
    let set = dir.join("set.faa");
    let mut bytes = std::fs::read(fasta("PemoMJNVA.faa")).unwrap();
    bytes.extend_from_slice(b">ten ten residues\nMKVLAAGIVW\n");
    std::fs::write(&set, bytes).unwrap();
    let db = makedb(&dir, &set, "set");
    for outfmt in ["6", "0", "7"] {
        let one = blastx(&query, &["-db", &db], outfmt, 1);
        assert!(!one.is_empty(), "outfmt {outfmt}");
        for threads in [2, 4] {
            assert!(
                blastx(&query, &["-db", &db], outfmt, threads) == one,
                "outfmt {outfmt}: {threads} threads differ from 1"
            );
        }
        if outfmt == "6" {
            assert!(rows(&one) > 10, "{}", rows(&one));
            let subject = blastx(&query, &["-subject", set.to_str().unwrap()], outfmt, 1);
            assert!(subject == one, "-db and -subject rows differ");
        }
    }
    // The reports name the database: outfmt 7 by the `-db` value, outfmt 0 by its title
    // (makedb's default title: the input file's name) and its numbers.
    let comments = String::from_utf8(blastx(&query, &["-db", &db], "7", 1)).unwrap();
    assert!(
        comments.contains(&format!("\n# Database: {db}\n")),
        "{comments}"
    );
    let pairwise = String::from_utf8(blastx(&query, &["-db", &db], "0", 1)).unwrap();
    assert!(pairwise.contains("Database: "), "{pairwise}");
    assert!(
        !pairwise.contains("User specified sequence set"),
        "{pairwise}"
    );
    assert!(pairwise.contains("\n>"), "{pairwise}");
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_setup.c:968-975 with seqdbimpl.cpp:131-136: a database answers 10 for its shortest
// sequence (BLASTX's subject is not translated, so the value is not divided), a FASTA set
// its shortest record. bench5's blastx pair: NCBI BLAST+ 2.17.0 keeps 229 rows with `-db`
// and 226 with `-subject` (dbsearch-20261012 DESIGN.md section 1.2).
#[test]
fn a_database_keeps_the_weak_hits_of_its_shortest_length() {
    let dir = temp_dir("d1");
    let db = makedb(&dir, &fasta("PemoMJNVA.faa"), "pemo");
    let query = fasta("MelaMJNV.fasta");
    let database = blastx(&query, &["-db", &db], "6", 1);
    let subject = blastx(
        &query,
        &["-subject", fasta("PemoMJNVA.faa").to_str().unwrap()],
        "6",
        1,
    );
    assert_eq!((rows(&database), rows(&subject)), (229, 226));
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one. On a database
// without its `.psq`, NCBI BLAST+ 2.17.0 reports a first batch of one query of N bases (no
// valid context: a batch of its own, longer than the 10,002 bases of a BLASTX batch) and
// fails the second batch with `BLAST engine error: Database memory map file error`.
#[test]
fn the_sequence_file_is_read_by_the_first_batch_with_a_valid_context() {
    let dir = temp_dir("psq");
    let db = makedb(&dir, &fasta("PemoMJNVA.faa"), "pemo");
    std::fs::remove_file(format!("{db}.psq")).unwrap();
    let all_n = format!(">nn all N\n{}\n", "N".repeat(10_005));
    let query = dir.join("nn.fna");
    std::fs::write(&query, &all_n).unwrap();
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(&query)
        .args(["-db", &db])
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(report.contains("Query= nn all N"), "{report}");
    let slice = genome_slice(&dir, 3_000);
    let mut both = all_n.into_bytes();
    both.extend_from_slice(&std::fs::read(&slice).unwrap());
    std::fs::write(&query, both).unwrap();
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(&query)
        .args(["-db", &db])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(3));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("Query= nn all N") && report.matches("Query= ").count() == 1,
        "{report}"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.ends_with("BLAST engine error: Database memory map file error\n"),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_args.cpp:2364-2381: `-db` excludes `-subject` (LOSAT's parser reports it);
// blast_args.cpp:2558-2562: neither is NCBI's options error; blast_app_util.hpp:185-203: a
// database that is not found is NCBI's `BLAST Database error`, exit status 2.
#[test]
fn database_arguments_are_checked_as_ncbi_checks_them() {
    let query = fasta("MelaMJNV.fasta");
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(&query)
        .args(["-db", "x", "-subject"])
        .arg(fasta("PemoMJNVA.faa"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(&query)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Either a BLAST database or subject sequence(s) must be specified"),
        "{stderr}"
    );
    let dir = temp_dir("missing");
    let missing = dir.join("no_such_db");
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(&query)
        .arg("-db")
        .arg(&missing)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with(
            "BLAST Database error: No alias or index file found for protein database ["
        ),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
