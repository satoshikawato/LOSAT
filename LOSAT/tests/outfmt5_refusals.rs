#![allow(warnings, clippy::all)]

//! BLAST XML (outfmt 5) of BLASTP and BLASTN: the searches whose XML NCBI BLAST+ 2.17.0 cannot
//! write (a byte >= 0x80 in a title or in the `-db` name: an exception, rc 255 with no output
//! for the first query, a document cut after the failed element for a later one) are refused
//! (DECISIONS 19, Q2) instead of reproducing the failure.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6900-6903
//! ```c++
//!     ... "Failed to convert symbol to requested encoding" ...
//! ```
//! (`CUtf8::SymbolToChar`, reached from `CObjectOStreamXml::WriteEncodedChar` of the US-ASCII
//! stream of blastxml_format.cpp:786.)

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
        "losat_outfmt5_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn run(args: &[&std::ffi::OsStr]) -> Output {
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

fn run_str(args: &[&str]) -> Output {
    let args: Vec<&std::ffi::OsStr> = args.iter().map(std::ffi::OsStr::new).collect();
    run(&args)
}

const PROTEIN: &[u8] = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
const NUCLEOTIDE: &[u8] = b"ATGGCTAGCTAGGATCCGATCGATTTAGGCTAGCTAGGCTAATCGGATCGATCGGGATCGATTAGCCTAGGATCGATCGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATGGCTAGCTAGGATCCGATCGATTTAGGCTAGCTAGGCTAATCGGATCGATCGGGATCGATTAGCCTAGGATCGATCGAAATTTCCCGGGAAATTTCCCGGG";
const OTHER_PROTEIN: &[u8] = b"WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM";
const OTHER_NUCLEOTIDE: &[u8] = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";

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

fn refused(output: &Output, program: &str, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's {program}"
        )),
        "{stderr}"
    );
}

fn programs() -> [(&'static str, &'static str, &'static [u8], &'static [u8]); 2] {
    [
        ("blastp", "BLASTP", PROTEIN, OTHER_PROTEIN),
        ("blastn", "BLASTN", NUCLEOTIDE, OTHER_NUCLEOTIDE),
    ]
}

#[test]
fn a_query_title_outside_us_ascii_is_refused() {
    for (program, name, hit, _) in programs() {
        let dir = temp_dir("query");
        let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
        for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
            // the first query, and a later one (NCBI catches the later one and goes on)
            for position in 0..2 {
                let mut records: Vec<(&[u8], &[u8])> = vec![(b"q1 ok", hit), (b"q2 ok", hit)];
                records[position] = (title, hit);
                fasta(&query, &records);
                fasta(&subject, &[(b"s1 subject", hit)]);
                let output = run_str(&[
                    program,
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    "5",
                ]);
                refused(&output, name, "a query title");
            }
        }
        fs::remove_dir_all(dir).ok();
    }
}

#[test]
fn a_subject_title_outside_us_ascii_is_refused_where_it_is_written() {
    for (program, name, hit, miss) in programs() {
        let dir = temp_dir("subject");
        let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
        fasta(&query, &[(b"q1 ok", hit)]);
        // the subject with the bad title has a hit: NCBI cuts the document after the Iteration's
        // hits start; LOSAT refuses
        fasta(
            &subject,
            &[(b"s1 ok", miss), ("s2 caf\u{e9}".as_bytes(), hit)],
        );
        let args = [
            program,
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "5",
        ];
        refused(&run_str(&args), name, "a subject title");
        // without a hit the subject is never written, and the search succeeds (an invalid query
        // has none)
        let invalid: Vec<u8> = vec![if program == "blastp" { b'X' } else { b'N' }; 40];
        fasta(&query, &[(b"q1 ok", &invalid)]);
        fasta(&subject, &[("s1 caf\u{e9}".as_bytes(), hit)]);
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
}

#[test]
fn a_database_name_outside_us_ascii_is_refused() {
    let dir = temp_dir("db");
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    for (program, name, db, query, extensions) in [
        (
            "blastn",
            "BLASTN",
            "nslice",
            "nslice_q.fna",
            ["nhr", "nin", "nsq"],
        ),
        (
            "blastp",
            "BLASTP",
            "residues",
            "residues_q.faa",
            ["phr", "pin", "psq"],
        ),
    ] {
        let copy = dir.join(format!("{db}_caf\u{e9}"));
        for extension in extensions {
            fs::copy(
                tests.join(format!("{db}.{extension}")),
                format!("{}.{extension}", copy.display()),
            )
            .expect("copy database");
        }
        let output = run(&[
            std::ffi::OsStr::new(program),
            std::ffi::OsStr::new("-query"),
            tests.join(query).as_os_str(),
            std::ffi::OsStr::new("-db"),
            copy.as_os_str(),
            std::ffi::OsStr::new("-outfmt"),
            std::ffi::OsStr::new("5"),
        ]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        // the database may be unsupported for another reason (its title); the name is the point
        if output.status.code() == Some(1) && stderr.contains("a database name") {
            assert!(output.stdout.is_empty());
            assert!(
                stderr.contains(&format!("a database name with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's {name}")),
                "{stderr}"
            );
        } else {
            panic!(
                "{program}: expected the database name to be refused: rc {:?} {stderr}",
                output.status.code()
            );
        }
    }
    fs::remove_dir_all(dir).ok();
}
