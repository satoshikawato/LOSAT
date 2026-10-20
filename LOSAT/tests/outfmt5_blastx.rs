#![allow(warnings, clippy::all)]

//! BLAST XML (outfmt 5) of BLASTX: the refusals of what NCBI BLAST+ 2.17.0 cannot write (a
//! byte >= 0x80 in a title or in the `-db` name, DECISIONS 19, Q2) and the document across
//! several query batches (one header, global Iteration numbers, one trailer).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-733
//! ```c++
//!     if( incremental_output && incremental_struct->m_IterationNum == 0) {
//!         ...
//!         *out_stream << serial_xml_start << "\n<BlastOutput_iterations>" ;
//!         incremental_struct->m_SerialXmlEnd = "\n</BlastOutput_iterations>" + serial_xml_end;
//!     }
//! ```

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt5_blastx_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn run_str(args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_X_DEFAULT" {
            command.env_remove(name);
        }
    }
    // An --all-features binary at 2+ threads keeps this variable (see the other tests).
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt5_worker.log"),
    );
    command.args(args).output().expect("run LOSAT")
}

const PROTEIN: &str = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
const OTHER_PROTEIN: &str = "WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM";

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
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTX"
        )),
        "{stderr}"
    );
}

#[test]
fn a_query_title_outside_us_ascii_is_refused() {
    let dir = temp_dir("query");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let nucleotide = back_translate(PROTEIN);
    for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
        for position in 0..2 {
            let mut records: Vec<(&[u8], &[u8])> = vec![
                (b"q1 ok", nucleotide.as_bytes()),
                (b"q2 ok", nucleotide.as_bytes()),
            ];
            records[position] = (title, nucleotide.as_bytes());
            fasta(&query, &records);
            fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
            let output = run_str(&[
                "blastx",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                "5",
            ]);
            if title == b"caf\xe9" {
                // BLASTX's reader refuses a title that is not UTF-8 for every format.
                assert_eq!(output.status.code(), Some(1));
                assert!(output.stdout.is_empty());
                assert!(
                    String::from_utf8_lossy(&output.stderr).contains("invalid UTF-8 FASTA title")
                );
            } else {
                refused(&output, "a query title");
            }
        }
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_subject_title_outside_us_ascii_is_refused_where_it_is_written() {
    let dir = temp_dir("subject");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let nucleotide = back_translate(PROTEIN);
    fasta(&query, &[(b"q1 ok", nucleotide.as_bytes())]);
    // The subject with the bad title has a hit: NCBI cuts the document after the Iteration's
    // hits start; LOSAT refuses.
    fasta(
        &subject,
        &[
            (b"s1 ok", OTHER_PROTEIN.as_bytes()),
            ("s2 caf\u{e9}".as_bytes(), PROTEIN.as_bytes()),
        ],
    );
    let args = [
        "blastx",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "5",
    ];
    refused(&run_str(&args), "a subject title");
    // Without a hit the subject is never written, and the search succeeds (an all-N query
    // has none).
    fasta(&query, &[(b"q1 ok", &[b'N'; 60])]);
    fasta(&subject, &[("s1 caf\u{e9}".as_bytes(), PROTEIN.as_bytes())]);
    let output = run_str(&args);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output
        .stdout
        .starts_with(b"<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n"));
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_database_name_outside_us_ascii_is_refused() {
    let dir = temp_dir("db");
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("bslice_caf\u{e9}");
    for extension in ["phr", "pin", "psq"] {
        fs::copy(
            tests.join(format!("bslice.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    let output = run_str(&[
        "blastx",
        "-query",
        tests.join("bslice_q.fna").to_str().unwrap(),
        "-db",
        copy.to_str().unwrap(),
        "-outfmt",
        "5",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(output.stdout.is_empty(), "nothing is written");
    assert!(
        stderr.contains("a database name with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTX"),
        "{stderr}"
    );
    fs::remove_dir_all(dir).ok();
}

/// 30 queries of about 1,500 letters (the batch holds 10,002 letters, so several batches),
/// each the back-translation of the protein, on alternate strands.
fn many_queries(path: &Path, bad_title: Option<(usize, &[u8])>) {
    let mut forward = back_translate(PROTEIN);
    forward.push_str(&"ACGTTGCAAGCTTGCA".repeat(60));
    let reverse: String = forward
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let titles: Vec<String> = (0..30).map(|n| format!("q{n} plain")).collect();
    let mut records: Vec<(&[u8], &[u8])> = Vec::new();
    for n in 0..30 {
        let sequence = if n % 2 == 0 { &forward } else { &reverse };
        records.push((titles[n].as_bytes(), sequence.as_bytes()));
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
    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
    for threads in ["1", "4"] {
        let output = run_str(&[
            "blastx",
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
        let text = String::from_utf8(output.stdout).expect("ASCII");
        assert_eq!(text.matches("<?xml").count(), 1);
        assert_eq!(text.matches("<BlastOutput_program>").count(), 1);
        assert_eq!(text.matches("<Iteration>").count(), 30);
        for n in 1..=30 {
            assert!(
                text.contains(&format!("<Iteration_iter-num>{n}</Iteration_iter-num>")),
                "{n}"
            );
        }
        // the numbers rise in the order of the document
        let numbers: Vec<u32> = text
            .split("<Iteration_iter-num>")
            .skip(1)
            .map(|rest| rest.split('<').next().unwrap().parse().unwrap())
            .collect();
        assert_eq!(numbers, (1..=30).collect::<Vec<u32>>());
        assert!(text.ends_with("</Iteration>\n</BlastOutput_iterations>\n</BlastOutput>\n\n"));
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_refusal_in_a_later_batch_ends_the_run_without_the_trailer() {
    // A query is checked when its batch is reached (the reader delivers the queries batch by
    // batch), so the batches before it are already written; the run stops with the refusal
    // and no trailer.
    let dir = temp_dir("later");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    many_queries(&query, Some((25, "q25 caf\u{e9}".as_bytes())));
    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
    let output = run_str(&[
        "blastx",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "5",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("a query title with a character outside the US-ASCII range"));
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(!text.contains("</BlastOutput>"));
    fs::remove_dir_all(dir).ok();
}
