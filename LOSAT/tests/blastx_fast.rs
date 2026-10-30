//! BLASTX `-task blastx-fast` on the call ABRicate makes (bin/abricate:129-141): the query comes on
//! stdin and there is no `-query`. The expected bytes are the frozen NCBI BLAST+ 2.17.0 row
//! `db.bxfast.abr.6` of `db_search_manifest.tsv`, the same search with `-query`.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Output, Stdio};

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn run_with_stdin(args: &[&str], stdin: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root())
        .args(args)
        .env_remove("BLASTDB")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    // A run that refuses its options exits without reading: the write may find the pipe closed.
    let _ = child.stdin.take().unwrap().write_all(stdin);
    child.wait_with_output().unwrap()
}

const ABRICATE: &[&str] = &[
    "blastx",
    "-task",
    "blastx-fast",
    "-seg",
    "no",
    "-db",
    "tests/db_search/bslice",
    "-outfmt",
    "6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen \
     stitle",
    "-num_threads",
    "1",
    "-evalue",
    "1E-20",
    "-culling_limit",
    "1",
    "-max_target_seqs",
    "10000",
];

#[test]
fn abricate_protein_database_call_with_the_query_on_stdin() {
    let query = std::fs::read(root().join("tests/db_search/bslice_q.fna")).unwrap();
    let expected =
        std::fs::read(root().join("tests/fixtures/db_search/db.bxfast.abr.6.out")).unwrap();
    assert!(!expected.is_empty(), "the search has hits");
    let omitted = run_with_stdin(ABRICATE, &query);
    let mut args = ABRICATE.to_vec();
    args.extend(["-query", "-"]);
    let named = run_with_stdin(&args, &query);
    for output in [omitted, named] {
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, expected);
        assert!(output.stderr.is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-859: a pipe has no position,
// so an empty one is not "Query is Empty!": NCBI prints nothing and exits 0.
#[test]
fn an_empty_pipe_is_no_query_and_no_warning() {
    let output = run_with_stdin(ABRICATE, b"");
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
}

// The formats that read the query file once before the search (outfmt 5, 11, 13 to 16) refuse
// standard input as they refuse any other query that cannot be read twice.
#[test]
fn formats_that_read_the_query_twice_refuse_standard_input() {
    for format in ["5", "11"] {
        let output = run_with_stdin(
            &[
                "blastx",
                "-subject",
                "tests/db_search/bslice.faa",
                "-outfmt",
                format,
            ],
            b">q\nACGT\n",
        );
        assert_eq!(output.status.code(), Some(1));
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr)
            .contains("unsupported BLASTX query input for -outfmt"));
    }
}

#[test]
fn a_subject_on_standard_input_stays_refused() {
    let output = run_with_stdin(&["blastx", "-subject", "-", "-outfmt", "6"], b">q\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("unsupported BLASTX stdin subject"));
}
