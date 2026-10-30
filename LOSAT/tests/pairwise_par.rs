//! The pairwise reports (outfmt 0) formatted on the search's pool (`LOSAT_X_PAIRWISEPAR`,
//! `src/report/x_pairwise_par.rs`) on whole searches: BLASTN (the blastn and megablast tasks,
//! `-subject` and `-db`), TBLASTX (`-subject` and `-db`) and BLASTP (`-subject` and `-db`, in one
//! search and in LOSAT_X_BPQBLOCK's query blocks).
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): it is on unless `LOSAT_X_PAIRWISEPAR`
//! is `0`. With it on and two or more threads, each
//! query's head, runs of its subjects' alignments and its footer are formatted into buffers by the
//! pool's threads and written in order.
//!
//! The property under test is that the reports, the messages and the exit status do not depend on
//! the switch or on the thread count: each search is compared with the same LOSAT search with the
//! switch off at one thread. `LOSAT_X_REACH` shows that the pool's writer ran at two threads or
//! more (new side) and not at one thread. The inputs have a record without letters between queries
//! (a warning before its report), queries with hits and without, several subjects, a query whose
//! alignments are cut into several pieces, and (BLASTN at two threads) more than one window of
//! queries.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_PAIRWISEPAR", "0");
/// The thread counts of `check`.
const THREADS: [u32; 4] = [1, 2, 4, 16];

fn manifest(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
}

fn letters(name: &str) -> Vec<u8> {
    std::fs::read(manifest("tests/fasta").join(name))
        .expect("test FASTA")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// A FASTA file of `records` (title, letters) in this test's temporary directory.
fn write_fasta(file: &str, records: &[(String, &[u8])]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-pairwisepar-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let path = dir.join(file);
    let mut text = Vec::new();
    for (title, seq) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in seq.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).expect("write FASTA");
    path
}

/// BLASTN queries: 30 parts of 4 kb of MelaMJNV (titles of several words), a record without
/// letters after the fifth, and a part of AP027131 that has no hits in the subjects; written
/// once (the tests run in parallel).
fn nucleotide_queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let mela = letters("MelaMJNV.fasta");
        let other = letters("AP027131.fasta");
        let mut records: Vec<(String, &[u8])> = Vec::new();
        for part in 0..30 {
            records.push((
                format!("mela part {part}"),
                &mela[part * 4_000..part * 4_000 + 4_000],
            ));
            if part == 4 {
                records.push(("empty record".to_string(), b""));
            }
        }
        records.push(("no hits".to_string(), &other[..3_000]));
        write_fasta("nucleotide_queries.fa", &records)
    })
}

/// Two BLASTN subjects (MejoMJNV and MeenMJNV), written once.
fn nucleotide_subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let mejo = letters("MejoMJNV.fasta");
        let meen = letters("MeenMJNV.fasta");
        write_fasta(
            "nucleotide_subjects.fa",
            &[
                ("mejo genome".to_string(), &mejo[..]),
                ("meen genome".to_string(), &meen[..]),
            ],
        )
    })
}

/// TBLASTX queries: a part of AP027131 with more than 2,048 HSPs (several pieces), a record
/// without letters and a second part; written once.
fn tblastx_queries() -> &'static Path {
    static QUERIES: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    QUERIES.get_or_init(|| {
        let a = letters("AP027131.fasta");
        write_fasta(
            "tblastx_queries.fa",
            &[
                ("a".to_string(), &a[..120_000]),
                ("empty record".to_string(), b""),
                ("a second part".to_string(), &a[400_000..430_000]),
            ],
        )
    })
}

/// A TBLASTX subject from AP027133, written once.
fn tblastx_subjects() -> &'static Path {
    static SUBJECTS: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    SUBJECTS.get_or_init(|| {
        let b = letters("AP027133.fasta");
        write_fasta(
            "tblastx_subjects.fa",
            &[("b first part".to_string(), &b[..120_000])],
        )
    })
}

/// A search with the given variables as its only `LOSAT_*` variables that choose a code path
/// (`BATCH_SIZE` of the test's environment is dropped too; `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature needs it
/// as soon as a search has more than one thread, and it selects nothing).
fn losat(args: &[&str], threads: u32, switches: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG") || text == "BATCH_SIZE"
        {
            command.env_remove(name);
        }
    }
    command
        .envs(switches.iter().copied())
        .args(args)
        .args(["-outfmt", "0", "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// Checks the property for one search at the thread counts `THREADS`.
fn check(name: &str, args: &[&str]) {
    check_with(name, args, &[], &THREADS)
}

/// `check` with more variables in every run, at `thread_counts`.
fn check_with(name: &str, args: &[&str], extra: &[(&str, &str)], thread_counts: &[u32]) {
    assert!(
        thread_counts.iter().any(|&threads| threads > 1),
        "{name}: a thread count above one"
    );
    let mut off_vars = extra.to_vec();
    off_vars.push(OFF);
    let off = losat(args, 1, &off_vars);
    assert!(off.status.success(), "{name}: {}", stderr(&off));
    assert!(!off.stdout.is_empty(), "{name}");
    for &threads in thread_counts {
        let reach = reach_support::Reach::new("pairwisepar");
        let mut on_vars = extra.to_vec();
        on_vars.push(reach.var());
        let on = losat(args, threads, &on_vars);
        assert!(
            on.status.success(),
            "{name}, threads {threads}: {}",
            stderr(&on)
        );
        assert_eq!(on.stdout, off.stdout, "{name}, threads {threads}");
        assert_eq!(stderr(&on), stderr(&off), "{name}, threads {threads}");
        // The pool's writer needs two threads: at one thread the serial writer runs, which marks
        // the reference side whenever the switch is on.
        let sides = reach.sides("LOSAT_X_PAIRWISEPAR");
        if threads == 1 {
            assert_eq!(
                sides,
                (false, true),
                "{name}, threads {threads}: the serial writer"
            );
        } else {
            assert!(
                sides.0,
                "{name}, threads {threads}: the pool's writer: {sides:?}"
            );
        }
    }
    let mut off_at_four_vars = extra.to_vec();
    off_at_four_vars.push(OFF);
    let off_at_four = losat(args, 4, &off_at_four_vars);
    assert!(
        off_at_four.status.success(),
        "{name}: {}",
        stderr(&off_at_four)
    );
    assert_eq!(
        off_at_four.stdout, off.stdout,
        "{name}, LOSAT_X_PAIRWISEPAR=0 at 4 threads"
    );
}

fn nucleotide(task: &str) -> Vec<String> {
    vec![
        "blastn".to_string(),
        "-task".to_string(),
        task.to_string(),
        "-query".to_string(),
        nucleotide_queries().display().to_string(),
        "-subject".to_string(),
        nucleotide_subjects().display().to_string(),
    ]
}

fn strings(args: &[String]) -> Vec<&str> {
    args.iter().map(String::as_str).collect()
}

#[test]
fn the_blastn_pairwise_report_does_not_depend_on_the_switch() {
    let args = nucleotide("blastn");
    check("blastn", &strings(&args));
}

#[test]
fn the_megablast_pairwise_report_does_not_depend_on_the_switch() {
    let args = nucleotide("megablast");
    check("megablast", &strings(&args));
}

#[test]
fn fewer_alignments_than_subjects_do_not_depend_on_the_switch() {
    let mut args = nucleotide("blastn");
    args.extend(["-max_target_seqs".to_string(), "1".to_string()]);
    check("blastn -max_target_seqs 1", &strings(&args));
}

#[test]
fn the_blastn_database_report_does_not_depend_on_the_switch() {
    let db = manifest("tests/db_search/nslice");
    let query = manifest("tests/db_search/nslice_q.fna");
    for task in ["blastn", "megablast"] {
        let args = [
            "blastn",
            "-task",
            task,
            "-query",
            query.to_str().unwrap(),
            "-db",
            db.to_str().unwrap(),
        ];
        check(&format!("blastn -db {task}"), &args);
    }
}

#[test]
fn the_tblastx_pairwise_report_does_not_depend_on_the_switch() {
    let queries = tblastx_queries().display().to_string();
    let subjects = tblastx_subjects().display().to_string();
    let args = [
        "tblastx",
        "-query",
        &queries,
        "-subject",
        &subjects,
        "-query_gencode",
        "4",
        "-db_gencode",
        "4",
    ];
    check("tblastx", &args);
}

#[test]
fn the_tblastx_database_report_does_not_depend_on_the_switch() {
    let db = manifest("tests/db_search/nslice");
    let query = manifest("tests/db_search/nslice_q.fna");
    let args = [
        "tblastx",
        "-query",
        query.to_str().unwrap(),
        "-db",
        db.to_str().unwrap(),
    ];
    check("tblastx -db", &args);
}

#[test]
fn a_failed_write_ends_the_search_as_without_the_switch() {
    if !Path::new("/dev/full").exists() {
        return;
    }
    let mut args = nucleotide("megablast");
    args.extend(["-out".to_string(), "/dev/full".to_string()]);
    let args = strings(&args);
    let off = losat(&args, 4, &[OFF]);
    assert!(!off.status.success(), "{}", stderr(&off));
    let on = losat(&args, 4, &[]);
    assert_eq!(on.status.code(), off.status.code(), "{}", stderr(&on));
    assert_eq!(stderr(&on), stderr(&off));
    assert_eq!(on.stdout, off.stdout);
}

/// A protein database of MejoMJNV (losat makeblastdb), made once.
fn protein_database() -> &'static Path {
    static DATABASE: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
    DATABASE.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("losat-pairwisepar-db-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temporary directory");
        let name = dir.join("MejoMJNV");
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .arg("makeblastdb")
            .arg("-in")
            .arg(manifest("tests/fasta/MejoMJNV.faa"))
            .args(["-dbtype", "prot"])
            .arg("-out")
            .arg(&name)
            .output()
            .expect("run losat makeblastdb");
        assert!(output.status.success(), "makeblastdb: {}", stderr(&output));
        name
    })
}

#[test]
fn the_blastp_pairwise_report_does_not_depend_on_the_switch() {
    let query = manifest("tests/fasta/MeenMJNV.faa");
    let fasta = manifest("tests/fasta/MejoMJNV.faa");
    let db = protein_database();
    // (subject option, its value, variables, thread counts). Each search of the fixture
    // pair costs CI seconds, and in the writer the thread count above one sets only how the pool
    // cuts the report (the unit tests of `src/report/x_pairwise_par.rs` write the serial report in
    // any chunks), so each case has two thread counts:
    // - one search of all queries with `-subject`, at one thread (the serial writer: the pool
    //   needs two threads) and four;
    // - `-db` in LOSAT_X_BPQBLOCK's blocks of 50,000 letters (two reports) and in one search of
    //   all queries, at two and four threads: outfmt 0 with LOSAT_X_BPDBCUT's cut, which acts at
    //   two threads or more (the search's own cuts: `blastp_db_cut.rs`).
    // `-subject` in blocks: `the_blastp_block_reports_do_not_depend_on_the_switch` below.
    let one_search: &[(&str, &str)] = &[("LOSAT_X_BPQBLOCK", "0")];
    let cases = [
        ("-subject", fasta.as_path(), one_search, [1, 4]),
        ("-db", db, &[][..], [2, 4]),
        ("-db", db, one_search, [2, 4]),
    ];
    for (subject, value, extra, thread_counts) in cases {
        let args = [
            "blastp",
            "-query",
            query.to_str().unwrap(),
            subject,
            value.to_str().unwrap(),
            "-evalue",
            "0.001",
        ];
        let name = format!("blastp {subject} {extra:?}");
        check_with(&name, &args, extra, &thread_counts);
    }
}

#[test]
fn the_blastp_block_reports_do_not_depend_on_the_switch() {
    // Small NCBI query batches and a block per batch (LOSAT_X_BPQBLOCK=1): one report per block.
    let query = manifest("tests/fasta/MeenMJNV.faa");
    let fasta = manifest("tests/fasta/MejoMJNV.faa");
    let args = [
        "blastp",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        fasta.to_str().unwrap(),
        "-evalue",
        "0.001",
    ];
    let blocks = [("BATCH_SIZE", "3000"), ("LOSAT_X_BPQBLOCK", "1")];
    let mut off_vars = blocks.to_vec();
    off_vars.push(OFF);
    let off = losat(&args, 1, &off_vars);
    assert!(off.status.success(), "{}", stderr(&off));
    for threads in [2, 4] {
        let reach = reach_support::Reach::new("pairwisepar-blocks");
        let mut on_vars: Vec<(&str, &str)> = blocks.to_vec();
        on_vars.push(reach.var());
        let on = losat(&args, threads, &on_vars);
        assert!(on.status.success(), "threads {threads}: {}", stderr(&on));
        assert_eq!(on.stdout, off.stdout, "threads {threads}");
        assert_eq!(stderr(&on), stderr(&off), "threads {threads}");
        let (new, _) = reach.sides("LOSAT_X_PAIRWISEPAR");
        assert!(new, "threads {threads}: the pool's writer");
    }
}
