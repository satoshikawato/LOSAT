#![allow(warnings, clippy::all)]

//! BLAST XML (outfmt 5) of BLASTX: the refusals of what NCBI BLAST+ 2.17.0 cannot write (a
//! byte >= 0x80 in a title or in the `-db` name, DECISIONS 19, Q2) and the document across
//! several query batches (one header, global Iteration numbers, one trailer), and the refusals
//! of a later batch (OF4b part 2): every refusal comes before the first byte of the output.
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
    run_env(args, &[])
}

fn run_env(args: &[&str], env: &[(&str, &str)]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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
    command.envs(env.iter().copied());
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
    refused_as(output, what, what)
}

fn refused_as(output: &Output, what: &str, label: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{label}: {stderr}");
    assert!(output.stdout.is_empty(), "{label}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTX"
        )),
        "{label}: {stderr}"
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

/// The run of `blastx -outfmt 5` on `query` and `subject` at 1 thread, at 4 threads, at 4
/// threads with the waves of queued batches instead of the flow, and with `-out`: each has
/// the same result. Returns the first.
fn runs(
    dir: &Path,
    query: &Path,
    subject: &Path,
    extra: &[&str],
) -> Vec<(String, Output, Vec<u8>)> {
    let mut results = Vec::new();
    let variants: [(&str, &str, &[(&str, &str)]); 3] = [
        ("1 thread", "1", &[]),
        ("4 threads", "4", &[]),
        ("4 threads, waves", "4", &[("LOSAT_X_BXFLOW", "0")]),
    ];
    for (what, threads, env) in variants {
        let out = dir.join(format!("out_{threads}_{}.xml", env.len()));
        let mut args = vec![
            "blastx",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "5",
            "-num_threads",
            threads,
        ];
        args.extend_from_slice(extra);
        let mut with_out = args.clone();
        with_out.extend_from_slice(&["-out", out.to_str().unwrap()]);
        let output = run_env(&with_out, env);
        let file = fs::read(&out).expect("the -out file is made before the search");
        results.push((what.to_string(), output, file));
        let output = run_env(&args, env);
        results.push((format!("{what}, stdout"), output, Vec::new()));
    }
    results
}

/// `refused` for every run of `runs`: rc 1, the message, nothing on standard output, and an
/// empty `-out` file.
fn all_refused(results: &[(String, Output, Vec<u8>)], what: &str) {
    for (variant, output, file) in results {
        refused_as(output, what, &format!("{what} ({variant})"));
        assert!(
            file.is_empty(),
            "{what} ({variant}): the -out file is empty"
        );
    }
}

#[test]
fn a_bad_query_title_in_a_later_batch_writes_nothing() {
    let dir = temp_dir("later_query");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    many_queries(&query, Some((25, "q25 caf\u{e9}".as_bytes())));
    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
    all_refused(&runs(&dir, &query, &subject, &[]), "a query title");
    // The last query of the file, and the first of a later batch.
    for index in [29, 12] {
        many_queries(&query, Some((index, "caf\u{e9}".as_bytes())));
        all_refused(&runs(&dir, &query, &subject, &[]), "a query title");
    }
    fs::remove_dir_all(dir).ok();
}

/// 30 queries: the first 25 are the back-translation of `PROTEIN` (hits for the subject `s1`),
/// the others of `OTHER_PROTEIN` (a hit for the subject `s2`, which comes only in the batches
/// after the first). The titles of the subjects are those given.
fn late_hit_files(dir: &Path, s2_title: &[u8]) -> (PathBuf, PathBuf) {
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let padding = "ACGTTGCAAGCTTGCA".repeat(60);
    let first = format!("{}{padding}", back_translate(PROTEIN));
    let late = format!("{}{padding}", back_translate(OTHER_PROTEIN_LONG));
    let titles: Vec<String> = (0..30).map(|n| format!("q{n} plain")).collect();
    let mut records: Vec<(&[u8], &[u8])> = Vec::new();
    for n in 0..30 {
        let sequence = if n < 25 { &first } else { &late };
        records.push((titles[n].as_bytes(), sequence.as_bytes()));
    }
    fasta(&query, &records);
    fasta(
        &subject,
        &[
            (b"s1 subject", PROTEIN.as_bytes()),
            (s2_title, OTHER_PROTEIN_LONG.as_bytes()),
        ],
    );
    (query, subject)
}

/// Only the exact matches are hits, so a query hits the subject it was made from alone.
const STRICT: &[&str] = &["-evalue", "1e-10"];

const OTHER_PROTEIN_LONG: &str =
    "MKVLAAGIVGLLSQPRTEWDNHYCFPSMKQDLERAIHTGWNVYFKCMLPQSTDERAGHIKNVWYFLMPQRSTE";

#[test]
fn a_bad_subject_title_whose_only_hit_is_in_a_later_batch_writes_nothing() {
    let dir = temp_dir("later_subject");
    // Control: the same searches with a plain title; s2 is hit by the late queries alone.
    let (query, subject) = late_hit_files(&dir, b"s2 plain");
    let control = run_str(&[
        "blastx",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "5",
        "-evalue",
        "1e-10",
    ]);
    assert_eq!(control.status.code(), Some(0));
    let text = String::from_utf8(control.stdout).expect("ASCII");
    assert_eq!(
        text.matches("<Hit_def>s2 plain</Hit_def>").count(),
        5,
        "{text}"
    );
    assert!(text.matches("<Hit_def>s1 subject</Hit_def>").count() >= 25);
    let (query, subject) = late_hit_files(&dir, "s2 caf\u{e9}".as_bytes());
    all_refused(&runs(&dir, &query, &subject, STRICT), "a subject title");
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_bad_subject_title_without_a_hit_gives_the_whole_document() {
    let dir = temp_dir("subject_no_hit");
    // No query hits s2: its title is never written, so the document is NCBI's (it equals the
    // document of the same searches with a plain title, which has been compared with
    // `blastx -outfmt 5 -num_threads 1` of NCBI BLAST+ 2.17.0).
    let write = |s2_title: &[u8]| {
        let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
        many_queries(&query, None);
        fasta(
            &subject,
            &[
                (b"s1 subject", PROTEIN.as_bytes()),
                (s2_title, OTHER_PROTEIN_LONG.as_bytes()),
            ],
        );
        (query, subject)
    };
    let (query, subject) = write(b"s2 plain");
    let expected = runs(&dir, &query, &subject, STRICT);
    let (query, subject) = write("s2 caf\u{e9}".as_bytes());
    let results = runs(&dir, &query, &subject, STRICT);
    let reference = &expected[1].1.stdout;
    assert_eq!(expected[1].1.status.code(), Some(0));
    let text = String::from_utf8_lossy(reference);
    assert_eq!(text.matches("<Iteration>").count(), 30);
    assert!(!text.contains("s2 plain"), "no query hits s2");
    assert!(text.ends_with("</Iteration>\n</BlastOutput_iterations>\n</BlastOutput>\n\n"));
    for (variant, output, file) in &results {
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{variant}: {stderr}");
        if variant.ends_with(", stdout") {
            assert_eq!(&output.stdout, reference, "{variant}");
        } else {
            assert!(output.stdout.is_empty(), "{variant}");
            assert_eq!(file, reference, "{variant}: the -out file");
        }
    }
    fs::remove_dir_all(dir).ok();
}

// A query that is not a regular file (OF4b review finding 1, DECISIONS 22): the pre-scan of the
// query titles (DECISIONS 21) reads the file a second time, which a pipe cannot give, so
// outfmt 5 refuses such a query input before the search, with nothing written. The other
// formats read the input once and keep working.

/// What the refusal says (`native.rs`).
#[cfg(unix)]
const NOT_REGULAR: &str = "unsupported BLASTX query input for -outfmt 5: use a regular file";

#[cfg(unix)]
fn one_query_files(dir: &Path) -> (Vec<u8>, PathBuf) {
    let subject = dir.join("s.fa");
    let query = format!(
        "{}{}",
        back_translate(PROTEIN),
        "ACGTTGCAAGCTTGCA".repeat(20)
    );
    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
    (format!(">q1 plain\n{query}\n").into_bytes(), subject)
}

/// `blastx` reading the query from a named pipe that a writer thread fills.
#[cfg(unix)]
fn run_fifo(dir: &Path, text: &[u8], subject: &Path, format: &str, out: Option<&Path>) -> Output {
    use std::io::Write;
    let fifo = dir.join(format!("query_{format}.fifo"));
    let made = Command::new("mkfifo").arg(&fifo).status().expect("mkfifo");
    assert!(made.success());
    let writer = {
        let (fifo, text) = (fifo.clone(), text.to_vec());
        // The open waits for the reader; the write may fail once a refusing run has gone.
        std::thread::spawn(move || {
            if let Ok(mut file) = fs::OpenOptions::new().write(true).open(fifo) {
                let _ = file.write_all(&text);
            }
        })
    };
    let mut args = vec![
        "blastx",
        "-query",
        fifo.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        format,
        "-evalue",
        "1e-10",
    ];
    if let Some(out) = out {
        args.extend_from_slice(&["-out", out.to_str().unwrap()]);
    }
    let output = run_str(&args);
    writer.join().expect("writer thread");
    fs::remove_file(&fifo).ok();
    output
}

/// `blastx -query /dev/stdin` with a pipe on standard input.
#[cfg(unix)]
fn run_dev_stdin(text: &[u8], subject: &Path, format: &str) -> Output {
    use std::io::Write;
    use std::process::Stdio;
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_X_DEFAULT" {
            command.env_remove(name);
        }
    }
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt5_worker.log"),
    );
    let mut child = command
        .args([
            "blastx",
            "-query",
            "/dev/stdin",
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            format,
            "-evalue",
            "1e-10",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    let mut stdin = child.stdin.take().expect("stdin");
    let text = text.to_vec();
    let writer = std::thread::spawn(move || {
        let _ = stdin.write_all(&text);
    });
    let output = child.wait_with_output().expect("wait");
    writer.join().expect("writer thread");
    output
}

#[cfg(unix)]
fn refused_not_regular(output: &Output, label: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{label}: {stderr}");
    assert!(output.stdout.is_empty(), "{label}: nothing is written");
    assert!(stderr.contains(NOT_REGULAR), "{label}: {stderr}");
}

#[cfg(unix)]
#[test]
fn a_named_pipe_query_is_refused_for_outfmt_5_and_works_for_outfmt_6() {
    let dir = temp_dir("fifo_query");
    let (text, subject) = one_query_files(&dir);
    refused_not_regular(&run_fifo(&dir, &text, &subject, "5", None), "stdout");
    let out = dir.join("out.xml");
    let output = run_fifo(&dir, &text, &subject, "5", Some(&out));
    refused_not_regular(&output, "-out");
    assert_eq!(fs::read(&out).expect("the -out file is made"), b"", "-out");
    // The other formats read the pipe once, as NCBI does.
    let tabular = run_fifo(&dir, &text, &subject, "6", None);
    assert_eq!(tabular.status.code(), Some(0));
    let rows = String::from_utf8(tabular.stdout).unwrap();
    assert!(rows.starts_with("q1\ts1\t"), "{rows}");
    // The same query as a regular file gives the document.
    let query = dir.join("q.fa");
    fs::write(&query, &text).unwrap();
    let document = run_str(&[
        "blastx",
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "5",
        "-evalue",
        "1e-10",
    ]);
    assert_eq!(document.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&document.stdout).contains("<Hit_def>s1 subject</Hit_def>"));
    fs::remove_dir_all(dir).ok();
}

#[cfg(unix)]
#[test]
fn dev_stdin_from_a_pipe_is_refused_for_outfmt_5_and_works_for_outfmt_6() {
    let dir = temp_dir("dev_stdin_query");
    let (text, subject) = one_query_files(&dir);
    refused_not_regular(&run_dev_stdin(&text, &subject, "5"), "/dev/stdin");
    let tabular = run_dev_stdin(&text, &subject, "6");
    assert_eq!(tabular.status.code(), Some(0));
    let rows = String::from_utf8(tabular.stdout).unwrap();
    assert!(rows.starts_with("q1\ts1\t"), "{rows}");
    fs::remove_dir_all(dir).ok();
}

#[test]
fn a_held_document_is_written_when_the_run_ends_with_a_reader_error() {
    // OF4b review finding 11: the hold (DECISIONS 21) ends a run that fails for any reason but
    // a refusal as the run without a hold does: the Iterations of the batches before the failed
    // one are written, the document has no trailer, rc 1 and the same message. The hold is on
    // when some `-subject` title could still be refused (s2, which no query hits).
    let dir = temp_dir("held_reader_error");
    let query = dir.join("q.fa");
    many_queries(&query, None);
    let mut text = fs::read(&query).unwrap();
    // 30 queries of about 1,400 letters: eight to a batch. The bad record is the 31st.
    text.extend_from_slice(b">q31 bad\nACGT%%%%\n");
    fs::write(&query, &text).unwrap();
    let subject = dir.join("s.fa");
    let mut results: Vec<(String, Vec<u8>, Vec<u8>)> = Vec::new();
    for (label, s2_title) in [
        ("held", "s2 caf\u{e9}".as_bytes()),
        ("not held", b"s2 plain".as_slice()),
    ] {
        fasta(
            &subject,
            &[
                (b"s1 subject", PROTEIN.as_bytes()),
                (s2_title, OTHER_PROTEIN_LONG.as_bytes()),
            ],
        );
        for threads in ["1", "4"] {
            for to_file in [false, true] {
                let out = dir.join(format!("out_{label}_{threads}.xml"));
                let mut args = vec![
                    "blastx",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    "5",
                    "-num_threads",
                    threads,
                    "-evalue",
                    "1e-10",
                ];
                if to_file {
                    args.extend_from_slice(&["-out", out.to_str().unwrap()]);
                }
                let output = run_str(&args);
                let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
                assert_eq!(output.status.code(), Some(1), "{label} {threads}: {stderr}");
                assert!(
                    stderr.contains("CFastaReader"),
                    "{label} {threads}: {stderr}"
                );
                let written = if to_file {
                    assert!(output.stdout.is_empty(), "{label} {threads}");
                    fs::read(&out).unwrap()
                } else {
                    output.stdout
                };
                results.push((
                    format!(
                        "{label}, {threads} threads, {}",
                        if to_file { "-out" } else { "stdout" }
                    ),
                    written,
                    output.stderr,
                ));
            }
        }
    }
    let (held, plain) = results.split_at(results.len() / 2);
    for ((label, written, stderr), (_, expected, expected_stderr)) in held.iter().zip(plain) {
        let text = String::from_utf8_lossy(written);
        assert!(text.matches("<Iteration>").count() >= 1, "{label}: {text}");
        assert!(!text.contains("</BlastOutput>"), "{label}: no trailer");
        assert_eq!(
            written, expected,
            "{label}: what the run without a hold writes"
        );
        assert_eq!(stderr, expected_stderr, "{label}");
    }
    fs::remove_dir_all(dir).ok();
}
