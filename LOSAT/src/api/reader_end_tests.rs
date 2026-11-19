//! A query input whose reader fails on a later record, or that holds blank and comment lines only
//! (line WQA's D10, 2026-11-19): the library's search of the records read before the error
//! ([`read_input_records`]) writes the reports of the query batches before the batch that reads
//! the record, then fails with the reader's error, as the command line run of the same files
//! does; [`read_inputs`] fails before any report. The command line half runs in this process,
//! through each program's own `run`, and writes its report to `-out`.

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;
use crate::api::{run, run_batches, validate, Program};
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::{try_parse_from, Cli, Commands};

/// The command line's message of the bad line of `with_bad_record` (its line number follows).
const BAD_LINE: &str = "BLAST query error: CFastaReader: Near line ";

/// A path of its own under the test directory (the counter tells apart two of one test).
fn temp_path(name: &str) -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    crate::utils::test_dir::temp_dir().join(format!(
        "losat_reader_end_{}_{}_{name}",
        crate::utils::test_dir::process_id(),
        CALLS.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Letters of `alphabet`, the same for the same `seed`.
fn letters(alphabet: &[u8], length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            alphabet[((state >> 33) % alphabet.len() as u64) as usize]
        })
        .collect()
}

/// A nucleotide sequence whose first frame reads `protein` (one codon per residue).
fn coding(protein: &[u8]) -> Vec<u8> {
    protein
        .iter()
        .flat_map(|residue| match residue {
            b'A' => *b"GCT",
            b'C' => *b"TGT",
            b'D' => *b"GAT",
            b'E' => *b"GAA",
            b'F' => *b"TTT",
            b'G' => *b"GGT",
            b'H' => *b"CAT",
            b'I' => *b"ATT",
            b'K' => *b"AAA",
            b'L' => *b"CTG",
            b'M' => *b"ATG",
            b'N' => *b"AAT",
            b'P' => *b"CCG",
            b'Q' => *b"CAG",
            b'R' => *b"CGT",
            b'S' => *b"AGC",
            b'T' => *b"ACC",
            b'V' => *b"GTG",
            b'W' => *b"TGG",
            _ => *b"TAT",
        })
        .collect()
}

/// FASTA of `count` records `q1`, `q2`, ... of `sequence`, 60 letters a line.
fn fasta(count: usize, sequence: &[u8]) -> Vec<u8> {
    let mut text = Vec::new();
    for index in 1..=count {
        text.extend_from_slice(format!(">q{index}\n").as_bytes());
        for line in sequence.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    text
}

/// `text` and then a record whose first line of letters NCBI's reader refuses (`CheckDataLine`).
fn with_bad_record(mut text: Vec<u8>) -> Vec<u8> {
    text.extend_from_slice(b">bad\nACGTAC%%%%\n");
    text
}

/// One program's inputs: the query and subject sequences, and the number of query records that
/// fill more than NCBI's first query batch of that program (blastp 10000 letters, tblastn 20000,
/// blastx and tblastx 10002; blastn's first batch, `CBatchSizeMixer`'s, is smaller than 40 records
/// of 5000 letters).
struct Case {
    program: Program,
    query: Vec<u8>,
    subject: Vec<u8>,
    records: usize,
}

fn cases() -> Vec<Case> {
    let protein = |length| letters(b"ACDEFGHIKLMNPQRSTVWY", length, 7);
    vec![
        Case {
            program: Program::Blastn,
            query: letters(b"ACGT", 5000, 3),
            subject: letters(b"ACGT", 5000, 3),
            records: 40,
        },
        Case {
            program: Program::Blastp,
            query: protein(1100),
            subject: protein(1100),
            records: 12,
        },
        Case {
            program: Program::Tblastn,
            query: protein(2100),
            subject: coding(&protein(2100)),
            records: 12,
        },
        Case {
            program: Program::Blastx,
            query: coding(&protein(1100)),
            subject: protein(1100),
            records: 5,
        },
        Case {
            program: Program::Tblastx,
            query: coding(&protein(1100)),
            subject: coding(&protein(1100)),
            records: 5,
        },
    ]
}

/// The command line's run of `query` against `subject` in outfmt 6, in this process: its report
/// (`-out`) and its error.
fn command_line(program: Program, query: &[u8], subject: &[u8]) -> (Vec<u8>, Option<String>) {
    let (query_path, subject_path, out) = (temp_path("q.fa"), temp_path("s.fa"), temp_path("out"));
    std::fs::write(&query_path, query).expect("write the query");
    std::fs::write(&subject_path, subject).expect("write the subject");
    let argv = [
        "losat",
        program.name(),
        "-query",
        query_path.to_str().unwrap(),
        "-subject",
        subject_path.to_str().unwrap(),
        "-outfmt",
        "6",
        "-out",
        out.to_str().unwrap(),
    ];
    let cli: Cli = try_parse_from(argv).expect("the argv parses");
    let settings = ApplicationSettings::default();
    let result = match cli.command {
        Commands::Blastn(args) => crate::algorithm::blastn::run(args, settings),
        Commands::Blastp(args) => crate::algorithm::blastp::run(args, settings),
        Commands::Tblastn(args) => crate::algorithm::tblastn::TblastnArgs::run(args, settings),
        Commands::Blastx(args) => crate::algorithm::blastx::BlastxArgs::run(args, settings),
        Commands::Tblastx(args) => crate::algorithm::tblastx::run(args, settings),
        _ => unreachable!("a search program"),
    };
    let report = std::fs::read(&out).unwrap_or_default();
    for path in [query_path, subject_path, out] {
        let _ = std::fs::remove_file(path);
    }
    let error = result
        .err()
        .map(|error| format!("{error:#}").trim_end().to_string());
    (report, error)
}

/// The library's search of the records of `query` and `subject` that `read_input_records` reads,
/// in outfmt 6: its report and its error.
fn library(program: Program, query: &[u8], subject: &[u8]) -> (Vec<u8>, Option<String>) {
    let search = validate(program, &[]).expect("no options");
    let records = read_input_records(&search, query, subject).expect("the subject reads");
    let mut report = Vec::new();
    let mut diagnostics = Vec::new();
    let result = run(
        &search,
        records.inputs(),
        &mut [("6", &mut report)],
        &mut diagnostics,
    );
    let error = result
        .err()
        .map(|error| error.message().trim_end().to_string());
    (report, error)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:216-222
// ```c++
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
//         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
// ...
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ```
#[test]
fn a_reader_error_ends_the_search_in_its_batch_as_on_the_command_line() {
    for case in cases() {
        let name = case.program.name();
        let records = fasta(case.records, &case.query);
        let query = with_bad_record(records.clone());
        let (cli_report, cli_error) = command_line(case.program, &query, &case.subject);
        let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: the command line fails"));
        assert!(cli_error.starts_with(BAD_LINE), "{name}: {cli_error}");
        assert!(!cli_report.is_empty(), "{name}: a batch before the error");
        let (report, error) = library(case.program, &query, &case.subject);
        assert_eq!(error.as_deref(), Some(cli_error.as_str()), "{name}");
        assert!(
            report == cli_report,
            "{name}: the reports of the batches before"
        );
        // The batch of the record is not searched: the whole input has more rows.
        let (whole, whole_error) = library(case.program, &records, &case.subject);
        assert_eq!(whole_error, None, "{name}");
        assert!(
            whole.len() > report.len(),
            "{name}: the last batch was searched"
        );
        assert!(whole.starts_with(&report), "{name}");
        // `read_inputs` fails before any report.
        let search = validate(case.program, &[]).unwrap();
        let early = read_inputs(&search, &query, &case.subject).unwrap_err();
        assert_eq!(early, Error::Input(cli_error), "{name}");
    }
}

// The batches of `run_batches` end alike: the batches before, then the error in the batch that
// reads the record; a reader error in the first batch is the first batch's error.
#[test]
fn run_batches_ends_with_the_reader_error_in_its_batch() {
    let case = cases().remove(1);
    let records = fasta(case.records, &case.query);
    let query = with_bad_record(records.clone());
    let search = validate(case.program, &[]).unwrap();
    let read = read_input_records(&search, &query, &case.subject).unwrap();
    assert_eq!(read.query.len(), case.records);
    assert!(matches!(read.query_end, QueryEnd::Error(_)));
    let batches: Vec<_> = run_batches(&search, read.inputs()).unwrap().collect();
    let [first, second] = &batches[..] else {
        panic!("two batches, the second the error: {}", batches.len())
    };
    let first = first.as_ref().expect("the first batch");
    assert_eq!(first.query_indices, 0..10);
    assert!(!first.hsps.is_empty());
    let error = second.as_ref().expect_err("the batch of the record");
    assert!(error.message().starts_with(BAD_LINE), "{error}");
    let (_, cli_error) = command_line(case.program, &query, &case.subject);
    assert_eq!(error.message().trim_end(), cli_error.unwrap());

    // An error in the first batch: no batch is searched.
    let query = with_bad_record(fasta(2, &case.query));
    let read = read_input_records(&search, &query, &case.subject).unwrap();
    let batches: Vec<_> = run_batches(&search, read.inputs()).unwrap().collect();
    assert_eq!(batches.len(), 1);
    assert!(batches[0].is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
// ```c++
// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
//     : m_QueryVector(& queries)
// {
//     if (queries.Empty()) {
//         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
//     }
// ```
// A query of blank and comment lines only is not NCBI's `Query is Empty!` (white space only): its
// first batch has no query, which every program refuses.
#[test]
fn a_query_of_comment_lines_only_is_a_batch_without_a_query() {
    for case in cases() {
        let name = case.program.name();
        let query = b"\n; a comment\n\n";
        let (cli_report, cli_error) = command_line(case.program, query, &case.subject);
        let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: the command line fails"));
        assert!(
            cli_error.contains("Empty CBlastQueryVector"),
            "{name}: {cli_error}"
        );
        let search = validate(case.program, &[]).unwrap();
        let read = read_input_records(&search, query, &case.subject).unwrap();
        assert!(read.query.is_empty());
        assert_eq!(read.query_end, QueryEnd::BlankLines, "{name}");
        let (report, error) = library(case.program, query, &case.subject);
        assert_eq!(error.as_deref(), Some(cli_error.as_str()), "{name}");
        assert_eq!(report, cli_report, "{name}");
    }
}
