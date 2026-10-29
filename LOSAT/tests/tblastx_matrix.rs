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
use std::path::Path;
use std::process::{Command, Output, Stdio};

const QUERY: &str = "tests/fasta/outfmt0/tblastx_multi_query.fasta";
const SUBJECT: &str = "tests/fasta/outfmt0/tblastx_multi_subject.fasta";
const DB_QUERY: &str = "tests/db_search/xslice_q.fna";
const DB: &str = "tests/db_search/xslice";
/// 180 bases of the first query of `QUERY` (the queries of `-query -` with a record without
/// letters).
const SEQUENCE: &str = "AGACTAAAGATTTAAGATTAGAAAAGAAAAAATGTCTAAATAATCAATTACCTTCATTCTATACATGTATGACGAATAATCCATCTTCTTGTATGAATACGACACTATTTGACGGTATATATAACAATGATTACATTTATAGTGATGAAAATGATGTTGACATGAGTGATGCTGACTCTG";

fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_matrix")
        .join(name);
    fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn run(args: &[&str], threads: &str) -> Output {
    run_with_stdin(args, threads, b"")
}

/// LOSAT with `stdin` on a pipe (written, then closed).
fn run_with_stdin(args: &[&str], threads: &str, stdin: &[u8]) -> Output {
    run_with_env(args, threads, stdin, &[])
}

/// `run_with_stdin` with the variables `envs` set for the run.
fn run_with_env(args: &[&str], threads: &str, stdin: &[u8], envs: &[(&str, &Path)]) -> Output {
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
    for (name, value) in envs {
        command.env(name, value);
    }
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
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(stdin)
        .expect("write stdin");
    child.wait_with_output().expect("LOSAT output")
}

/// LOSAT prints NCBI's bytes of `case` and exits with `status`, at 1 and 4 threads.
fn equals_ncbi(case: &str, args: &[&str], status: i32) {
    equals_ncbi_with_stdin(case, args, b"", status)
}

/// `equals_ncbi` with `stdin` on the pipe (the standard input of NCBI's run).
fn equals_ncbi_with_stdin(case: &str, args: &[&str], stdin: &[u8], status: i32) {
    let (out, err) = (
        fixture(&format!("{case}.out")),
        fixture(&format!("{case}.err")),
    );
    for threads in ["1", "4"] {
        let output = run_with_stdin(args, threads, stdin);
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

/// A query without letters keeps its own message, after the error of its list (the error has
/// the id -1, the message the query's index; `RemoveDuplicates` sorts each list).
#[test]
fn an_unknown_matrix_lists_the_messages_of_queries_without_letters() {
    let first = format!(">empty\n>q1\n{SEQUENCE}\n");
    equals_ncbi_with_stdin(
        "norecord_unknown",
        &["-query", "-", "-subject", SUBJECT, "-matrix", "FOO"],
        first.as_bytes(),
        3,
    );
    let last = format!(">q1\n{SEQUENCE}\n>empty\n");
    equals_ncbi_with_stdin(
        "norecord_unknown_last",
        &[
            "-query", "-", "-subject", SUBJECT, "-matrix", "FOO", "-outfmt", "6",
        ],
        last.as_bytes(),
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
    // Without queries NCBI makes no report data: the document's header and end, exit status 0.
    for format in ["5", "16"] {
        equals_ncbi(
            &format!("empty_lower_{format}"),
            &[
                "-query", "-", "-subject", SUBJECT, "-matrix", "pam30", "-outfmt", format,
            ],
            0,
        );
    }
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
    equals_ncbi(
        "empty_unknown_15",
        &[
            "-query", "-", "-subject", SUBJECT, "-matrix", "FOO", "-outfmt", "15",
        ],
        0,
    );
}

// No NCBI counterpart: NCBI 2.17.0 reads a name that starts with `identity` (any case) as
// IDENTITY's table without IDENTITY's statistics and aborts (`free(): invalid size`, exit status
// 134) when it sets up the first batch with letters; LOSAT refuses that batch before anything is
// written, as its other refusals of NCBI defects. Without queries NCBI reads no matrix
// (`empty_identityx`: the prolog and the epilog with `Matrix: IDENTITYX`, exit status 0).
#[test]
fn a_name_read_as_identity_without_its_statistics_is_refused() {
    for (name, subject) in [
        ("IDENTITYX", ["-subject", SUBJECT]),
        ("identityx", ["-subject", SUBJECT]),
        ("Identity1", ["-db", DB]),
    ] {
        let query = if subject[0] == "-db" { DB_QUERY } else { QUERY };
        let output = run(
            &["-query", query, subject[0], subject[1], "-matrix", name],
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
    equals_ncbi(
        "empty_identityx",
        &["-query", "-", "-subject", SUBJECT, "-matrix", "IDENTITYX"],
        0,
    );
}

/// The runs of LOSAT_X_PAIRPAR's own side that `LOSAT_X_REACH` recorded in `dir` (the second
/// column of its rows).
fn units_reached(dir: &Path) -> u64 {
    let mut runs = 0;
    for entry in fs::read_dir(dir).expect("reach directory") {
        let text = fs::read_to_string(entry.expect("reach file").path()).expect("reach file");
        for row in text.lines() {
            let fields: Vec<&str> = row.split('\t').collect();
            if fields.first() == Some(&"LOSAT_X_PAIRPAR") {
                runs += fields[1].parse::<u64>().expect("reach count");
            }
        }
    }
    runs
}

/// Series SCT, SC7b: another matrix's search takes the unit path (PAIRLUT's group tables with its
/// neighbour words, PAIRPAR's units, the table walks of TBXUNGAP with its byte table) and prints
/// what the reference frame loop prints. A query of 3 kb (6,000 letters, above the one-thread gate
/// of 5,000) against 20 kb of the related genome: at one thread with the units (`LOSAT_X_REACH`
/// records LOSAT_X_PAIRPAR's own side) and without them (`LOSAT_X_PAIRPARN1=0`, the reference
/// frame loop), and at four threads (the units of the one subject).
#[test]
fn another_matrix_takes_the_unit_path_with_the_reference_output() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir =
        std::env::temp_dir().join(format!("losat_tblastx_matrix_units_{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temporary directory");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let head =
        |path: &str, bytes: usize| fs::read(root.join(path)).expect("FASTA")[..bytes].to_vec();
    fs::write(
        &query,
        head("tests/fasta/truncated/AP027131_50k.fasta", 3100),
    )
    .expect("query");
    fs::write(
        &subject,
        head("tests/fasta/truncated/AP027133_50k.fasta", 20400),
    )
    .expect("subject");
    let (query, subject) = (query.to_str().unwrap(), subject.to_str().unwrap());
    let off = Path::new("0");
    for matrix in [
        "BLOSUM45", "BLOSUM50", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250", "IDENTITY",
    ] {
        let args = [
            "-query", query, "-subject", subject, "-matrix", matrix, "-outfmt", "6",
        ];
        let reach = dir.join(format!("reach_{matrix}"));
        fs::create_dir_all(&reach).expect("reach directory");
        let units = run_with_env(&args, "1", b"", &[("LOSAT_X_REACH", &reach)]);
        assert!(units_reached(&reach) > 0, "{matrix}: the units did not run");
        let reach_off = dir.join(format!("reach_{matrix}_off"));
        fs::create_dir_all(&reach_off).expect("reach directory");
        let reference = run_with_env(
            &args,
            "1",
            b"",
            &[("LOSAT_X_PAIRPARN1", off), ("LOSAT_X_REACH", &reach_off)],
        );
        assert_eq!(units_reached(&reach_off), 0, "{matrix}: the units ran");
        let four = run(&args, "4");
        for output in [&units, &reference, &four] {
            assert_eq!(output.status.code(), Some(0), "{matrix}");
            assert!(output.stderr.is_empty(), "{matrix}");
        }
        assert!(!reference.stdout.is_empty(), "{matrix}: no hits");
        assert_eq!(
            String::from_utf8_lossy(&units.stdout),
            String::from_utf8_lossy(&reference.stdout),
            "{matrix} n1"
        );
        assert_eq!(
            String::from_utf8_lossy(&four.stdout),
            String::from_utf8_lossy(&reference.stdout),
            "{matrix} n4"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}
