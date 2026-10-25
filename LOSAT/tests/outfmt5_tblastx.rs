#![allow(warnings, clippy::all)]

//! BLAST XML (outfmt 5) of TBLASTX: the refusals of what NCBI BLAST+ 2.17.0 cannot write (a
//! byte >= 0x80 in a title or in the `-db` name, DECISIONS 19, Q2), made before the first byte
//! of stdout or `-out` at 1 and 4 threads, and the document across several query batches (one
//! header, global Iteration numbers, one trailer: TBLASTX renders once per run, after every
//! batch is searched), and the run that stops at a failing batch (the earlier batches' Iterations,
//! no trailer).
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
        "losat_outfmt5_tblastx_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn run_str(args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
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
const OTHER_NUCLEOTIDE: &[u8] = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";

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
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's TBLASTX"
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
            fasta(&subject, &[(b"s1 subject", nucleotide.as_bytes())]);
            for threads in ["1", "4"] {
                let out = dir.join("out.xml");
                let output = run_str(&[
                    "tblastx",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    "5",
                    "-num_threads",
                    threads,
                    "-out",
                    out.to_str().unwrap(),
                ]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{stderr}");
                assert!(stderr.contains("a query title with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's TBLASTX"));
                // NCBI creates the -out file before it checks anything (blast_args.cpp:3456-3481),
                // so the file exists, empty.
                assert_eq!(
                    fs::metadata(&out).map(|m| m.len()).unwrap_or(0),
                    0,
                    "-out stays empty"
                );
                // to stdout
                let output = run_str(&[
                    "tblastx",
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
    // hits start; LOSAT refuses before the first byte.
    fasta(
        &subject,
        &[
            (b"s1 ok", OTHER_NUCLEOTIDE),
            ("s2 caf\u{e9}".as_bytes(), nucleotide.as_bytes()),
        ],
    );
    for threads in ["1", "4"] {
        let args = [
            "tblastx",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "5",
            "-num_threads",
            threads,
        ];
        refused(&run_str(&args), "a subject title");
        // Without a hit the subject is never written, and the search succeeds (an all-N query
        // has none).
        fasta(&query, &[(b"q1 ok", &[b'N'; 90])]);
        fasta(
            &subject,
            &[("s1 caf\u{e9}".as_bytes(), nucleotide.as_bytes())],
        );
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
        fasta(&query, &[(b"q1 ok", nucleotide.as_bytes())]);
        fasta(
            &subject,
            &[
                (b"s1 ok", OTHER_NUCLEOTIDE),
                ("s2 caf\u{e9}".as_bytes(), nucleotide.as_bytes()),
            ],
        );
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_database_name_outside_us_ascii_is_refused() {
    let dir = temp_dir("db");
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("xslice_caf\u{e9}");
    for extension in ["nhr", "nin", "nsq"] {
        fs::copy(
            tests.join(format!("xslice.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for threads in ["1", "4"] {
        let output = run_str(&[
            "tblastx",
            "-query",
            tests.join("xslice_q.fna").to_str().unwrap(),
            "-db",
            copy.to_str().unwrap(),
            "-outfmt",
            "5",
            "-num_threads",
            threads,
        ]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{stderr}");
        assert!(output.stdout.is_empty(), "nothing is written");
        assert!(
            stderr.contains("a database name with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's TBLASTX"),
            "{stderr}"
        );
    }
    fs::remove_dir_all(dir).ok();
}

/// 120 queries of 390 letters (a query batch holds 10,002 letters, so several batches).
fn many_queries(path: &Path, bad_title: Option<(usize, &[u8])>) {
    let titles: Vec<String> = (0..120).map(|n| format!("q{n} plain")).collect();
    let nucleotide = back_translate(PROTEIN);
    let mut records: Vec<(&[u8], &[u8])> = Vec::new();
    for n in 0..120 {
        records.push((titles[n].as_bytes(), nucleotide.as_bytes()));
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
    for threads in ["1", "4"] {
        let output = run_str(&[
            "tblastx",
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
        assert_eq!(text.matches("<Iteration>").count(), 120);
        let numbers: Vec<u32> = text
            .split("<Iteration_iter-num>")
            .skip(1)
            .map(|rest| rest.split('<').next().unwrap().parse().unwrap())
            .collect();
        assert_eq!(numbers, (1..=120).collect::<Vec<u32>>());
        assert!(text.ends_with("</Iteration>\n</BlastOutput_iterations>\n</BlastOutput>\n\n"));
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_refusal_in_a_later_batch_writes_nothing() {
    // TBLASTX formats once, after every batch is searched: the document is checked as a whole
    // before its first byte.
    let dir = temp_dir("later");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    many_queries(&query, Some((110, "q110 caf\u{e9}".as_bytes())));
    fasta(
        &subject,
        &[(b"s1 subject", back_translate(PROTEIN).as_bytes())],
    );
    for threads in ["1", "4"] {
        let output = run_str(&[
            "tblastx",
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
    // `-query -` with standard input on /dev/null (seekable): NCBI's tblastx stops at "Query is
    // Empty!" before the formatter writes anything (no header, no epilog): 0 bytes and exit
    // status 0. An empty pipe is another path (one newline): `outfmt5_no_query.rs`.
    let dir = temp_dir("empty");
    let subject = dir.join("s.fa");
    fasta(
        &subject,
        &[(b"s1 subject", back_translate(PROTEIN).as_bytes())],
    );
    for threads in ["1", "4"] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            if name.to_string_lossy().starts_with("LOSAT_") {
                command.env_remove(name);
            }
        }
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat_outfmt5_worker.log"),
        );
        let output = command
            .args([
                "tblastx",
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
        assert_eq!(output.stderr, b"Warning: [tblastx] Query is Empty!\n");
    }
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_failing_batch_ends_the_run_after_the_earlier_iterations_without_the_trailer() {
    // A reader error in the second query batch (BATCH_SIZE 500 letters: q1 and q2 fill the
    // first batch) comes after the first batch's report, and PrintEpilog is never reached: the
    // document ends after the second Iteration, as NCBI's does (checked against NCBI BLAST+
    // 2.17.0 by hand, RESULTS4.md). A reader error in the first batch writes nothing.
    let dir = temp_dir("failing");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let nucleotide = back_translate(PROTEIN);
    fasta(&subject, &[(b"s1 subject", nucleotide.as_bytes())]);
    let mut text = Vec::new();
    for title in ["q1 plain", "q2 plain"] {
        text.extend_from_slice(format!(">{title}\n{nucleotide}\n").as_bytes());
    }
    let good = text.clone();
    text.extend_from_slice(b">q3\nACGT%%%%\n");
    for threads in ["1", "4"] {
        let run = |text: &[u8]| {
            fs::write(&query, text).expect("write query");
            let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
            for (name, _) in std::env::vars_os() {
                if name.to_string_lossy().starts_with("LOSAT_") {
                    command.env_remove(name);
                }
            }
            command.env(
                "LOSAT_BLASTX_WORKER_LOG",
                std::env::temp_dir().join("losat_outfmt5_worker.log"),
            );
            command
                .env("BATCH_SIZE", "500")
                .args([
                    "tblastx",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    "5",
                    "-num_threads",
                    threads,
                ])
                .output()
                .expect("run LOSAT")
        };
        let complete = run(&good);
        assert_eq!(complete.status.code(), Some(0));
        let complete = String::from_utf8(complete.stdout).expect("ASCII");
        assert_eq!(complete.matches("<Iteration>").count(), 2);
        assert!(complete.ends_with("</Iteration>\n</BlastOutput_iterations>\n</BlastOutput>\n\n"));
        let failed = run(&text);
        let stderr = String::from_utf8_lossy(&failed.stderr);
        assert_eq!(failed.status.code(), Some(1), "{stderr}");
        assert!(
            stderr.starts_with("BLAST query error: CFastaReader: "),
            "{stderr}"
        );
        let partial = String::from_utf8(failed.stdout).expect("ASCII");
        // The same document without its trailer.
        let trailer = "\n</BlastOutput_iterations>\n</BlastOutput>\n\n";
        assert_eq!(partial, complete.strip_suffix(trailer).expect("trailer"));
        // A reader error in the first batch (BATCH_SIZE 5000: q1, q2 and the bad record share
        // the batch) writes nothing.
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        for (name, _) in std::env::vars_os() {
            if name.to_string_lossy().starts_with("LOSAT_") {
                command.env_remove(name);
            }
        }
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat_outfmt5_worker.log"),
        );
        let first = command
            .env("BATCH_SIZE", "5000")
            .args([
                "tblastx",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                "5",
                "-num_threads",
                threads,
            ])
            .output()
            .expect("run LOSAT");
        assert_eq!(first.status.code(), Some(1));
        assert!(first.stdout.is_empty());
    }
    fs::remove_dir_all(dir).ok();
}
