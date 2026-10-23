#![allow(warnings, clippy::all)]

//! BLAST XML (outfmt 5) of BLASTP and BLASTN when no query is formatted (an empty standard
//! input as the query): NCBI writes the epilog's end of line alone, one newline.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2242-2243
//! ```c++
//!     if (m_FormatType == CFormattingArgs::eXml) {
//!         m_Outfile << m_BlastXMLIncremental->m_SerialXmlEnd << endl;
//! ```
//! (`m_SerialXmlEnd` is empty until a query is formatted; frozen NCBI 2.17.0 output: `\n`.)

use std::path::Path;
use std::process::{Command, Stdio};

const PROTEIN_SUBJECT: &[u8] =
    b">s1 subject\nMWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPK\n";
const NUCLEOTIDE_SUBJECT: &[u8] =
    b">s1 subject\nATGGCTAGCTAGGATCCGATCGATTTAGGCTAGCTAGGCTAATCGGATCGATCGGGATCGATTAGCCTAGG\n";

fn run_with_empty_stdin(args: &[&str]) -> std::process::Output {
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
    let mut child = command
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    drop(child.stdin.take());
    child.wait_with_output().expect("wait for LOSAT")
}

#[test]
fn an_empty_query_input_writes_one_newline() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt5_noquery_{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    for (program, subject, db) in [
        ("blastp", PROTEIN_SUBJECT, "residues"),
        ("blastn", NUCLEOTIDE_SUBJECT, "nslice"),
    ] {
        let subject_path = dir.join(format!("{program}.fa"));
        std::fs::write(&subject_path, subject).expect("write subject");
        let db_path = tests.join(db);
        for subject_args in [
            ["-subject", subject_path.to_str().unwrap()],
            ["-db", db_path.to_str().unwrap()],
        ] {
            for threads in ["1", "4"] {
                let mut args = vec![
                    program,
                    "-query",
                    "-",
                    "-outfmt",
                    "5",
                    "-num_threads",
                    threads,
                ];
                args.extend(subject_args);
                let output = run_with_empty_stdin(&args);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{args:?}: {stderr}");
                assert_eq!(output.stdout, b"\n", "{args:?}");
            }
        }
    }
    std::fs::remove_dir_all(dir).ok();
}

/// TBLASTN and TBLASTX (OF4b): the same newline, also with a matrix name TBLASTN runs and the
/// XML formatting table does not know (`-matrix blosum62`, NCBI's lookup of the run ignores
/// case, the XML writer's `strcmp` does not; with no formatted query NCBI never reads the
/// table).
#[test]
fn an_empty_query_pipe_writes_one_newline_for_tblastn_and_tblastx() {
    let dir = std::env::temp_dir().join(format!("losat_outfmt5_noquery_tx_{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let subject_path = dir.join("subject.fa");
    std::fs::write(&subject_path, NUCLEOTIDE_SUBJECT).expect("write subject");
    let cases: [(&str, &[&str]); 4] = [
        ("tblastn", &[]),
        ("tblastn", &["-matrix", "blosum62"]),
        (
            "tblastn",
            &[
                "-matrix",
                "BLOSUM45",
                "-word_size",
                "2",
                "-comp_based_stats",
                "0",
            ],
        ),
        ("tblastx", &[]),
    ];
    for (program, extra) in cases {
        for threads in ["1", "4"] {
            let mut args = vec![
                program,
                "-query",
                "-",
                "-subject",
                subject_path.to_str().unwrap(),
                "-outfmt",
                "5",
                "-num_threads",
                threads,
            ];
            args.extend_from_slice(extra);
            let output = run_with_empty_stdin(&args);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(0), "{args:?}: {stderr}");
            assert_eq!(output.stdout, b"\n", "{args:?}");
        }
    }
    std::fs::remove_dir_all(dir).ok();
}
