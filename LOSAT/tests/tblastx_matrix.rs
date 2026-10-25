#![allow(warnings, clippy::all)]

//! TBLASTX `-matrix` (series SCT, SC7): the runs that stop, which the frozen-fixture rows
//! (`score.tblastx.*` of `outfmt0_manifest.tsv`, exit status 0) cannot hold. The bytes of NCBI
//! BLAST+ 2.17.0 are in `tests/fixtures/tblastx_matrix` (`<case>.out`, `<case>.err`: `tblastx
//! ARGS -num_threads 1`, run from LOSAT/ with an empty pipe as standard input). LOSAT prints them
//! at 1 and 4 threads.
//!
//! NCBI reference: c++/src/algo/blast/core/blast_message.c:242-247 (a name without a table:
//! `Blast_ScoreBlkMatrixFill` returns -1, one message per query of the batch)
//! ```c
//!     default:
//!         {
//!             char buf[512];
//!             snprintf(buf, sizeof(buf) - 1, "Unknown error code %d", error_code);
//! ```
//! and c++/src/algo/blast/format/data4xmlformat.cpp:211-218 (the XML reports, when the first
//! report is made: a name that is not spelled as the table's)
//! ```c++
//!     else {
//!         string prog_name = Blast_ProgramNameFromType(
//!                                            m_Options.GetProgramType());
//!         if (prog_name != "blastn" && prog_name != "megablast") {
//!             NCBI_THROW(blast::CBlastException, eInvalidArgument,
//!                         "unsupported score matrix");
//! ```

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const QUERY: &str = "tests/fasta/outfmt0/tblastx_multi_query.fasta";
const SUBJECT: &str = "tests/fasta/outfmt0/tblastx_multi_subject.fasta";
const DB_QUERY: &str = "tests/db_search/xslice_q.fna";
const DB: &str = "tests/db_search/xslice";

fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_matrix")
        .join(name);
    fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn run(args: &[&str], threads: &str) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") || name_text.starts_with("NCBI") {
            command.env_remove(&name);
        }
    }
    for name in ["BATCH_SIZE", "CHUNK_SIZE", "BLASTDB", "BLASTMAT"] {
        command.env_remove(name);
    }
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_tblastx_matrix_worker.log"),
    );
    let mut child = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("tblastx")
        .args(args)
        .args(["-num_threads", threads])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    drop(child.stdin.take());
    child.wait_with_output().expect("LOSAT output")
}

/// LOSAT prints NCBI's bytes of `case` and exits with `status`, at 1 and 4 threads.
fn equals_ncbi(case: &str, args: &[&str], status: i32) {
    let (out, err) = (
        fixture(&format!("{case}.out")),
        fixture(&format!("{case}.err")),
    );
    for threads in ["1", "4"] {
        let output = run(args, threads);
        assert_eq!(output.status.code(), Some(status), "{case} n{threads}");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&out),
            "{case} n{threads} stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            String::from_utf8_lossy(&err),
            "{case} n{threads} stderr"
        );
    }
}

/// A name without a table (`FOO`) stops the first batch with exit status 3 and one message per
/// query of the batch (two of the four queries): after the prolog (outfmt 0: the header;
/// 16: the XML2 head; 15: `{\n"BlastOutput2": [\n`), with nothing on the other formats.
#[test]
fn an_unknown_matrix_stops_the_first_batch() {
    for format in ["0", "6", "7", "5", "16", "15"] {
        equals_ncbi(
            &format!("unknown_{format}"),
            &[
                "-query", QUERY, "-subject", SUBJECT, "-matrix", "FOO", "-outfmt", format,
            ],
            3,
        );
    }
    equals_ncbi(
        "unknown_db_0",
        &["-query", DB_QUERY, "-db", DB, "-matrix", "FOO"],
        3,
    );
}

/// The XML reports read the name as typed: `pam30` searches (outfmt 0 and the tabular formats
/// are `score.tblastx.lc_pam30.0`), and the first report stops the run with exit status 3.
#[test]
fn the_xml_reports_refuse_a_name_in_lower_case() {
    for format in ["5", "16", "15"] {
        equals_ncbi(
            &format!("lower_{format}"),
            &[
                "-query", QUERY, "-subject", SUBJECT, "-matrix", "pam30", "-outfmt", format,
            ],
            3,
        );
    }
    equals_ncbi(
        "lower_db_16",
        &[
            "-query", DB_QUERY, "-db", DB, "-matrix", "pam30", "-outfmt", "16",
        ],
        3,
    );
}

/// The options are checked before the matrix (exit status 1), and a run without a query reads no
/// matrix: `-query -` from an empty pipe prints the prolog and the epilog (`Matrix: FOO`) with
/// exit status 0. (From `/dev/null` NCBI warns `Query is Empty!` and prints nothing; LOSAT does
/// the same, probe `s_FOO_o0` of series SCT.)
#[test]
fn an_unknown_matrix_comes_after_the_options_and_the_queries() {
    equals_ncbi(
        "unknown_threshold0",
        &[
            "-query",
            QUERY,
            "-subject",
            SUBJECT,
            "-matrix",
            "FOO",
            "-threshold",
            "0",
        ],
        1,
    );
    equals_ncbi(
        "empty_unknown",
        &["-query", "-", "-subject", SUBJECT, "-matrix", "FOO"],
        0,
    );
}

// No NCBI counterpart: NCBI 2.17.0 reads a name that starts with `identity` (any case) as
// IDENTITY's table without IDENTITY's statistics and aborts (`free(): invalid size`, exit status
// 134); LOSAT refuses it before the search, as its other refusals of NCBI defects.
#[test]
fn a_name_read_as_identity_without_its_statistics_is_refused() {
    for name in ["IDENTITYX", "identityx", "Identity1"] {
        let output = run(
            &["-query", QUERY, "-subject", SUBJECT, "-matrix", name],
            "1",
        );
        assert_eq!(output.status.code(), Some(1), "{name}");
        assert!(output.stdout.is_empty(), "{name}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(&format!("-matrix {name} is read as IDENTITY"))
                && stderr.contains("not supported by LOSAT"),
            "{name}: {stderr}"
        );
    }
}
