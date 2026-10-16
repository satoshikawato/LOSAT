//! TBLASTX's search of a database (`tblastx -db NAME`, stage DB-5 of the database search
//! campaign) on a whole run of the binary.
//!
//! The properties under test are LOSAT's own: the report of a database search is the same
//! at every thread count; for a subject set of A, C, G and T only, the tabular rows of `-db`
//! are those of `-subject` with the same records (NCBI scans a subject set as a database, and
//! TBLASTX has no Gumbel block, so the statistics are the same); `-dbsize` gives the same
//! statistics to both; and the database's errors come where NCBI raises them. That the bytes
//! equal NCBI BLAST+ `tblastx -db` is the business of the campaign's comparison runs.

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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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
        "losat-tblastx-db-search-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The letters of a one-record FASTA file.
fn letters(name: &str) -> String {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

/// A FASTA record of `letters` in lines of 70.
fn record(defline: &str, letters: &str) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.as_bytes().chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The first 6 kb of MjeNMV, and the MelaMJNV and PemoMJNVA genomes as a set of two records
/// (bench5's tblastx inputs, without ambiguity codes).
fn inputs(dir: &Path) -> (PathBuf, PathBuf) {
    let query = dir.join("q.fna");
    std::fs::write(
        &query,
        record(
            "mje6k first 6 kb of MjeNMV",
            &letters("MjeNMV.fasta")[..6000],
        ),
    )
    .unwrap();
    let set = dir.join("set.fna");
    std::fs::write(
        &set,
        record("mela MelaMJNV", &letters("MelaMJNV.fasta"))
            + &record("pemo PemoMJNVA", &letters("PemoMJNVA.fasta")),
    )
    .unwrap();
    (query, set)
}

/// The database of `fasta_file` made by `LOSAT makedb` in `dir` under `name`.
fn makedb(dir: &Path, fasta_file: &Path, name: &str) -> String {
    let db = dir.join(name).to_str().unwrap().to_string();
    let output = losat()
        .arg("makedb")
        .arg("-in")
        .arg(fasta_file)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    db
}

/// One search of `query` on standard output; returns the report.
fn tblastx(query: &Path, subjects: &[&str], outfmt: &str, threads: u32) -> Vec<u8> {
    let output = losat()
        .arg("tblastx")
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
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set, "set");
    for outfmt in ["6", "0", "7"] {
        let one = tblastx(&query, &["-db", &db], outfmt, 1);
        assert!(!one.is_empty(), "outfmt {outfmt}");
        for threads in [2, 4] {
            assert!(
                tblastx(&query, &["-db", &db], outfmt, threads) == one,
                "outfmt {outfmt}: {threads} threads differ from 1"
            );
        }
        if outfmt == "6" {
            assert!(rows(&one) > 100, "{}", rows(&one));
            let subject = tblastx(&query, &["-subject", set.to_str().unwrap()], outfmt, 1);
            assert!(subject == one, "-db and -subject rows differ");
        }
    }
    // The reports name the database: outfmt 7 by the `-db` value, outfmt 0 by its title
    // (makedb's default title: the input file's name) and its numbers.
    let comments = String::from_utf8(tblastx(&query, &["-db", &db], "7", 1)).unwrap();
    assert!(
        comments.contains(&format!("\n# Database: {db}\n")),
        "{comments}"
    );
    let pairwise = String::from_utf8(tblastx(&query, &["-db", &db], "0", 1)).unwrap();
    assert!(
        pairwise.contains("Database: ") && pairwise.contains("2 sequences;"),
        "{pairwise}"
    );
    assert!(
        !pairwise.contains("User specified sequence set"),
        "{pairwise}"
    );
    // A database sequence's heading is `>` and its stored title.
    assert!(pairwise.contains("\n>mela MelaMJNV\n"), "{pairwise}");
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_setup.c:716-735: `-dbsize` replaces the database length of the effective lengths, with
// `-db` and with `-subject` alike (the hit-saving average and the link cutoffs keep the real
// total: blast_setup.c:976-978, blast_engine.c:1407). A database length whose search space
// overflows NCBI's Int8 is LOSAT's refusal.
#[test]
fn the_database_size_gives_both_paths_the_same_statistics() {
    let dir = temp_dir("dbsize");
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set, "set");
    let subject = set.to_str().unwrap();
    let given_size = tblastx(&query, &["-db", &db, "-dbsize", "100000000"], "6", 1);
    assert!(
        tblastx(
            &query,
            &["-subject", subject, "-dbsize", "100000000"],
            "6",
            1
        ) == given_size,
        "-db and -subject differ with -dbsize"
    );
    let default_size = tblastx(&query, &["-db", &db], "6", 1);
    assert!(given_size != default_size && rows(&given_size) <= rows(&default_size));
    for value in ["0", "-5"] {
        assert!(
            tblastx(&query, &["-db", &db, "-dbsize", value], "6", 1) == default_size,
            "-dbsize {value}"
        );
    }
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(&query)
        .args(["-db", &db, "-dbsize", "9000000000000000000", "-outfmt", "6"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's TBLASTX"
        ),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// local_blast.cpp:177-224 (`CheckInternalData`): NCBI does not search a query batch without a
// valid context, so the sequence file is read by the first batch with one. A database without
// its `.nsq` answers a query of N letters (no valid context) as NCBI BLAST+ 2.17.0 does: a
// report without hits, exit status 0; with a second query in a batch of its own
// (`BATCH_SIZE=1`), that batch fails after the report of the first one.
#[test]
fn the_sequence_file_is_read_by_the_first_batch_with_a_valid_context() {
    let dir = temp_dir("nsq");
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set, "set");
    std::fs::remove_file(format!("{db}.nsq")).unwrap();
    let all_n = record("nn all N", &"N".repeat(60));
    let unsearched = dir.join("nn.fna");
    std::fs::write(&unsearched, &all_n).unwrap();
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(&unsearched)
        .args(["-db", &db])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(0));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("Query= nn all N") && report.contains("***** No hits found *****"),
        "{report}"
    );
    let two = dir.join("nn_q.fna");
    std::fs::write(&two, all_n + &std::fs::read_to_string(&query).unwrap()).unwrap();
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(&two)
        .args(["-db", &db])
        .env("BATCH_SIZE", "1")
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

// tblastx_app.cpp:203-209: NCBI writes each query batch's report after its search, and reads
// the titles of the batch's reported sequences there. In batches of one query (`BATCH_SIZE=1`)
// two queries that report the same sequences give the rows of `-subject` (the reported
// sequences are named by rank over all batches); without the `.nhr`, the second batch, the first
// with hits, fails with NCBI's database error after the first batch's report.
#[test]
fn query_batches_read_the_titles_of_their_reported_sequences() {
    let dir = temp_dir("batches");
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set, "set");
    let genome = letters("MjeNMV.fasta");
    let two = dir.join("two.fna");
    std::fs::write(
        &two,
        record("q1 MjeNMV 1-3000", &genome[..3000])
            + &record("q2 MjeNMV 3001-6000", &genome[3000..6000]),
    )
    .unwrap();
    let search = |subjects: &[&str], query: &Path| {
        let output = losat()
            .arg("tblastx")
            .arg("-query")
            .arg(query)
            .args(subjects)
            .args(["-outfmt", "6"])
            .env("BATCH_SIZE", "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    };
    let rows = search(&["-db", &db], &two);
    let text = String::from_utf8_lossy(&rows);
    for (query_id, subject_id) in [
        ("q1", "mela"),
        ("q1", "pemo"),
        ("q2", "mela"),
        ("q2", "pemo"),
    ] {
        assert!(
            text.lines()
                .any(|row| row.starts_with(&format!("{query_id}\t{subject_id}\t"))),
            "{query_id} x {subject_id}: {text}"
        );
    }
    assert!(rows == search(&["-subject", set.to_str().unwrap()], &two));
    std::fs::remove_file(format!("{db}.nhr")).unwrap();
    let first_unsearched = dir.join("nn_q.fna");
    std::fs::write(
        &first_unsearched,
        record("nn all N", &"N".repeat(60)) + &std::fs::read_to_string(&query).unwrap(),
    )
    .unwrap();
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(&first_unsearched)
        .args(["-db", &db])
        .env("BATCH_SIZE", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(
        report.contains("Query= nn all N") && report.matches("Query= ").count() == 1,
        "{report}"
    );
    // The first batch's warning (no valid context), then the error.
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.starts_with("Warning: [tblastx] Query_1 nn all N: ")
            && stderr
                .lines()
                .last()
                .is_some_and(|line| line.starts_with("BLAST Database error: ")),
        "{stderr}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// A sequence file shorter than its index says is LOSAT's refusal where NCBI reads past its
// mapping: it is raised by the first batch with a valid context, as NCBI's read errors are, and
// stops the run there without the reports of the batches before (as TBLASTN's does).
#[test]
fn a_truncated_sequence_file_stops_the_run() {
    let dir = temp_dir("cut");
    let (query, set) = inputs(&dir);
    let db = makedb(&dir, &set, "set");
    let nsq = format!("{db}.nsq");
    let bytes = std::fs::read(&nsq).unwrap();
    std::fs::write(&nsq, &bytes[..bytes.len() / 2]).unwrap();
    let two = dir.join("nn_q.fna");
    std::fs::write(
        &two,
        record("nn all N", &"N".repeat(60)) + &std::fs::read_to_string(&query).unwrap(),
    )
    .unwrap();
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(&two)
        .args(["-db", &db])
        .env("BATCH_SIZE", "1")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    let report = String::from_utf8_lossy(&output.stdout);
    assert!(!report.contains("Query= "), "{report}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("passes the end of the file"), "{stderr}");
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_args.cpp:2364-2381: `-db` excludes `-subject` and `-subject_loc` (LOSAT's parser
// reports it); blast_app_util.hpp:185-203: a database that is not found is NCBI's
// `BLAST Database error`, exit status 2.
#[test]
fn database_arguments_are_checked_as_ncbi_checks_them() {
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(fasta("MjeNMV.fasta"))
        .args(["-db", "x", "-subject"])
        .arg(fasta("MelaMJNV.fasta"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let dir = temp_dir("missing");
    let missing = dir.join("no_such_db");
    let output = losat()
        .arg("tblastx")
        .arg("-query")
        .arg(fasta("MjeNMV.fasta"))
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
