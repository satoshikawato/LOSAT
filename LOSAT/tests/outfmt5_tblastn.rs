//! BLAST XML (outfmt 5) of TBLASTN: the one title that NCBI BLAST+ 2.17.0 cannot make (punctuation
//! alone, such as `, ,`, on which `x_CleanAndCompress` reads past the end of the string) is refused
//! before the first byte of stdout or `-out` at 1, 4, 8 and 16 threads, also in a later batch, and the
//! document across several query batches (one header, global Iteration numbers, one trailer: TBLASTN
//! renders once per run). The titles NCBI does write (a byte >= 0x80, a control byte, a title whose
//! encoding it cannot guess, the `-db` name) are the rows `db.xb.tblastn.*` of
//! `tests/db_search_manifest/tblastn.tsv` and `tests/xml_dbname.rs`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-733

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use losat::api::{read_fasta, run_batches, validate, Inputs, Program, Role};

mod temp_support;
use temp_support::{Finished, TempFile};

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt5_tblastn_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn run_str(args: &[&str]) -> Finished<Output, TempFile> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    // An --all-features binary at 2+ threads keeps this variable (see the other tests).
    let log = TempFile::keep_on_failure("losat_outfmt5_worker", "worker.log");
    command.env("LOSAT_BLASTX_WORKER_LOG", &log);
    let output = command.args(args).output().expect("run LOSAT");
    Finished::new(output, log)
}

const PROTEIN: &str = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";

/// One codon per amino acid.
fn codon(letter: char) -> &'static str {
    match letter {
        'A' => "GCT",
        'C' => "TGC",
        'D' => "GAT",
        'E' => "GAA",
        'F' => "TTC",
        'G' => "GGT",
        'H' => "CAC",
        'I' => "ATC",
        'K' => "AAA",
        'L' => "CTG",
        'M' => "ATG",
        'N' => "AAC",
        'P' => "CCT",
        'Q' => "CAG",
        'R' => "CGT",
        'S' => "TCT",
        'T' => "ACC",
        'V' => "GTG",
        'W' => "TGG",
        'Y' => "TAC",
        other => panic!("no codon for {other}"),
    }
}

fn back_translate(protein: &str) -> String {
    protein.chars().map(codon).collect()
}

fn fasta(path: &Path, records: &[(&[u8], &[u8])]) {
    let mut text = Vec::new();
    for (title, sequence) in records {
        text.extend_from_slice(b">");
        text.extend_from_slice(title);
        text.push(b'\n');
        text.extend_from_slice(sequence);
        text.push(b'\n');
    }
    fs::write(path, text).expect("write FASTA");
}

fn refused(output: &Output, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line is not supported by LOSAT's TBLASTN"
        )),
        "{stderr}"
    );
}

/// The title of punctuation that NCBI's `x_CleanAndCompress` reads past the end of.
const PUNCTUATION: &[u8] = b", ,";

/// Queries in the file of `many_queries`: 150 of `PROTEIN.len()` (140) letters are 21,000
/// letters, so the first batch (TBLASTN's query batch takes queries until 20,000 letters are
/// read) holds 143 and the second the other 7.
const QUERIES: usize = 150;

/// The query batches of `many_queries` as the search makes them: the ranges of query numbers that
/// `losat::api::run_batches` gives for the same queries against one subject.
fn query_batches() -> Vec<std::ops::Range<usize>> {
    let mut text = Vec::new();
    for n in 0..QUERIES {
        text.extend(format!(">q{n} plain\n{PROTEIN}\n").bytes());
    }
    let queries = read_fasta(Program::Tblastn, Role::Query, &text).unwrap();
    let subject = format!(">s1 subject\n{}\n", back_translate(PROTEIN));
    let subjects = read_fasta(Program::Tblastn, Role::Subject, subject.as_bytes()).unwrap();
    let search = validate(Program::Tblastn, &[]).unwrap();
    run_batches(&search, Inputs::new(&queries, &subjects))
        .unwrap()
        .map(|batch| batch.unwrap().query_indices)
        .collect()
}

/// 150 queries of `PROTEIN.len()` letters: more than a query batch holds, so two batches
/// (`query_batches`, asserted here).
fn many_queries(path: &Path, bad_title: Option<(usize, &[u8])>) {
    let batches = query_batches();
    assert!(batches.len() >= 2, "{batches:?}");
    if let Some((index, _)) = bad_title {
        assert!(index >= batches[0].end, "the title is in a later batch");
    }
    let titles: Vec<String> = (0..QUERIES).map(|n| format!("q{n} plain")).collect();
    let mut records: Vec<(&[u8], &[u8])> = Vec::new();
    for n in 0..QUERIES {
        records.push((titles[n].as_bytes(), PROTEIN.as_bytes()));
    }
    if let Some((index, title)) = bad_title {
        records[index].0 = title;
    }
    fasta(path, &records);
}

#[test]
fn several_batches_make_one_document_with_global_iteration_numbers() {
    let dir = temp_dir("batches");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    many_queries(&query, None);
    fasta(
        &subject,
        &[(b"s1 subject", back_translate(PROTEIN).as_bytes())],
    );
    for threads in ["1", "4", "8", "16"] {
        let output = run_str(&[
            "tblastn",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "5",
            "-num_threads",
            threads,
        ]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{stderr}");
        let text = String::from_utf8(output.stdout.clone()).expect("ASCII");
        assert_eq!(text.matches("<?xml").count(), 1);
        assert_eq!(text.matches("<BlastOutput_program>").count(), 1);
        assert_eq!(text.matches("<Iteration>").count(), QUERIES);
        let numbers: Vec<u32> = text
            .split("<Iteration_iter-num>")
            .skip(1)
            .map(|rest| rest.split('<').next().unwrap().parse().unwrap())
            .collect();
        assert_eq!(numbers, (1..=QUERIES as u32).collect::<Vec<u32>>());
        assert!(text.ends_with("</Iteration>\n</BlastOutput_iterations>\n</BlastOutput>\n\n"));
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_refusal_in_a_later_batch_writes_nothing() {
    // TBLASTN formats once, after every batch is searched: the document is checked as a whole
    // before its first byte.
    let dir = temp_dir("later");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    many_queries(&query, Some((148, PUNCTUATION)));
    fasta(
        &subject,
        &[(b"s1 subject", back_translate(PROTEIN).as_bytes())],
    );
    for threads in ["1", "4", "8", "16"] {
        let output = run_str(&[
            "tblastn",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "5",
            "-num_threads",
            threads,
        ]);
        refused(&output, "a query title");
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn an_empty_seekable_query_input_writes_nothing() {
    // `-query -` with standard input on /dev/null (seekable): NCBI's tblastn stops at "Query is
    // Empty!" before the formatter writes anything (no header, no epilog): 0 bytes and exit
    // status 0. An empty pipe is another path (one newline): `outfmt5_no_query.rs`.
    let dir = temp_dir("empty");
    let subject = dir.join("s.fa");
    fasta(
        &subject,
        &[(b"s1 subject", back_translate(PROTEIN).as_bytes())],
    );
    for threads in ["1", "4", "8", "16"] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            if name.to_string_lossy().starts_with("LOSAT_") {
                command.env_remove(name);
            }
        }
        let log = TempFile::keep_on_failure("losat_outfmt5_worker", "worker.log");
        command.env("LOSAT_BLASTX_WORKER_LOG", &log);
        let output = command
            .args([
                "tblastn",
                "-query",
                "-",
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                "5",
                "-num_threads",
                threads,
            ])
            .stdin(std::process::Stdio::null())
            .output()
            .expect("run LOSAT");
        assert_eq!(output.status.code(), Some(0));
        assert!(output.stdout.is_empty());
        assert_eq!(output.stderr, b"Warning: [tblastn] Query is Empty!\n");
    }
    fs::remove_dir_all(dir).ok();
}
