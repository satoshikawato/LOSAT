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

/// The command line's message of the bad line of `with_bad_record` (its line number follows), in
/// LOSAT's words (`crate::wording::QUERY_ERROR`, line NW2).
fn bad_line() -> String {
    format!("{} CFastaReader: Near line ", crate::wording::QUERY_ERROR)
}

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
    command_line_with(program, query, subject, &[], "6")
}

/// `command_line` with the options `options` and the format `outfmt`.
fn command_line_with(
    program: Program,
    query: &[u8],
    subject: &[u8],
    options: &[&str],
    outfmt: &str,
) -> (Vec<u8>, Option<String>) {
    let (query_path, subject_path, out) = (temp_path("q.fa"), temp_path("s.fa"), temp_path("out"));
    std::fs::write(&query_path, query).expect("write the query");
    std::fs::write(&subject_path, subject).expect("write the subject");
    let mut argv = vec![
        "losat",
        program.name(),
        "-query",
        query_path.to_str().unwrap(),
        "-subject",
        subject_path.to_str().unwrap(),
        "-outfmt",
        outfmt,
        "-out",
        out.to_str().unwrap(),
    ];
    argv.extend_from_slice(options);
    let cli: Cli = try_parse_from(argv).expect("the argv parses");
    let settings = ApplicationSettings::default();
    let result = match cli.command {
        Commands::Blastn(args) => crate::algorithm::blastn::run(args, settings),
        Commands::Blastp(args) => crate::algorithm::blastp::run(args, settings),
        Commands::Tblastn(args) => crate::algorithm::tblastn::TblastnArgs::run(args, settings),
        Commands::Blastx(args) => crate::algorithm::blastx::BlastxArgs::run(args, settings),
        Commands::Tblastx(args) => crate::algorithm::tblastx::run(args, settings),
        Commands::Psiblast(args) => crate::algorithm::psiblast::run(args, settings),
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
    library_with(program, query, subject, &[], "6")
}

/// `library` with the options `options` and the format `outfmt`.
fn library_with(
    program: Program,
    query: &[u8],
    subject: &[u8],
    options: &[&str],
    outfmt: &str,
) -> (Vec<u8>, Option<String>) {
    let search = validate(program, options).expect("the options");
    let records = read_input_records(&search, query, subject).expect("the subject reads");
    let mut report = Vec::new();
    let mut diagnostics = Vec::new();
    let result = run(
        &search,
        records.inputs(),
        &mut [(outfmt, &mut report)],
        &mut diagnostics,
    );
    let error = result.err().map(|error| {
        String::from_utf8_lossy(error.message())
            .trim_end()
            .to_string()
    });
    (report, error)
}

/// The rows of `run_batches`'s search of `query` and `subject` with `options` (outfmt 6's, one an
/// HSP) before its error, and the error: of the call, or of the batch that ends the batches.
fn batches_with(
    program: Program,
    query: &[u8],
    subject: &[u8],
    options: &[&str],
) -> (usize, Option<String>) {
    let search = validate(program, options).expect("the options");
    let records = read_input_records(&search, query, subject).expect("the subject reads");
    let batches = match run_batches(&search, records.inputs()) {
        Ok(batches) => batches,
        Err(error) => {
            return (
                0,
                Some(
                    String::from_utf8_lossy(error.message())
                        .trim_end()
                        .to_string(),
                ),
            )
        }
    };
    let mut rows = 0;
    for batch in batches {
        match batch {
            Ok(batch) => rows += batch.hsps.len(),
            Err(error) => {
                return (
                    rows,
                    Some(
                        String::from_utf8_lossy(error.message())
                            .trim_end()
                            .to_string(),
                    ),
                )
            }
        }
    }
    (rows, None)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:216-222
// With `-parse_deflines` (the batches' Seq-id scopes, `query_batch::QueryScopes`) and with
// `-query_loc` (a cut of each record, the batch counting the whole record) the run ends alike
// (review L3, 2026-11-20).
#[test]
fn a_reader_error_ends_the_search_in_its_batch_as_on_the_command_line() {
    for options in [&[][..], &["-parse_deflines"], &["-query_loc", "10-200"]] {
        for case in cases() {
            let name = format!("{} {options:?}", case.program.name());
            // A subject with a defline, which `-parse_deflines` reads.
            let subject = fasta_named(&["s1".to_string()], &case.subject);
            let records = fasta(case.records, &case.query);
            let query = with_bad_record(records.clone());
            let (cli_report, cli_error) =
                command_line_with(case.program, &query, &subject, options, "6");
            let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: the command line fails"));
            assert!(cli_error.starts_with(&bad_line()), "{name}: {cli_error}");
            assert!(!cli_report.is_empty(), "{name}: a batch before the error");
            let (report, error) = library_with(case.program, &query, &subject, options, "6");
            assert_eq!(error.as_deref(), Some(cli_error.as_str()), "{name}");
            assert!(
                report == cli_report,
                "{name}: the reports of the batches before"
            );
            // The batch of the record is not searched: the whole input has more rows.
            let (whole, whole_error) = library_with(case.program, &records, &subject, options, "6");
            assert_eq!(whole_error, None, "{name}");
            assert!(
                whole.len() > report.len(),
                "{name}: the last batch was searched"
            );
            assert!(whole.starts_with(&report), "{name}");
            // `run_batches`: the rows of the batches before, then the error.
            let (rows, batches_error) = batches_with(case.program, &query, &subject, options);
            assert_eq!(
                rows,
                cli_report.split(|&byte| byte == b'\n').count() - 1,
                "{name}"
            );
            assert_eq!(batches_error.as_deref(), Some(cli_error.as_str()), "{name}");
            // `read_inputs` fails before any report.
            let search = validate(case.program, options).unwrap();
            let early = read_inputs(&search, &query, &subject).unwrap_err();
            assert_eq!(early, Error::Input(cli_error.into_bytes()), "{name}");
        }
    }
}

/// FASTA of records named `names`, each of `sequence`, 60 letters a line.
fn fasta_named(names: &[String], sequence: &[u8]) -> Vec<u8> {
    let mut text = Vec::new();
    for name in names {
        text.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in sequence.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    text
}

// Under `-parse_deflines` with `-subject` LOSAT refuses two queries of the run with one id
// (`query_batch::QueryScopes`), and a scope holds the batches that the reader delivers: a batch
// that a reader error ends adds none, so its error is the run's (review M3, 2026-11-20). The
// library's run and `run_batches` refuse where the command line does, in outfmt 6 and 5: after the
// reports of the batches before the one with the second id (line WQA D10 round 4, 2026-11-20).
#[test]
fn a_scope_holds_the_batches_that_the_reader_delivers() {
    let options = ["-parse_deflines"];
    for case in cases() {
        let subject = fasta_named(&["s1".to_string()], &case.subject);
        // One record more than `cases`, so that the last batch holds two.
        let count = case.records + 1;
        let mut names: Vec<String> = (1..=count).map(|index| format!("q{index}")).collect();
        for (variant, other) in [("in its batch", count - 2), ("in the first batch", 0)] {
            names[count - 1] = names[other].clone();
            let duplicate = fasta_named(&names, &case.query);
            let name = format!("{} {variant}", case.program.name());
            for outfmt in ["6", "5"] {
                // The second id alone: LOSAT's refusal.
                let (cli_report, cli_error) =
                    command_line_with(case.program, &duplicate, &subject, &options, outfmt);
                let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: a refusal"));
                assert!(!cli_error.starts_with(&bad_line()), "{name}: {cli_error}");
                let (report, error) =
                    library_with(case.program, &duplicate, &subject, &options, outfmt);
                assert_eq!(
                    error.as_deref(),
                    Some(cli_error.as_str()),
                    "{name} {outfmt}"
                );
                assert!(report == cli_report, "{name} {outfmt}: the reports");
                // Then a record that the reader fails on, in the batch of the second id.
                let failing = with_bad_record(duplicate.clone());
                let (cli_report, cli_error) =
                    command_line_with(case.program, &failing, &subject, &options, outfmt);
                let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: the reader's error"));
                assert!(
                    cli_error.starts_with(&bad_line()),
                    "{name} {outfmt}: {cli_error}"
                );
                assert!(!cli_report.is_empty(), "{name} {outfmt}: the first batch");
                let (report, error) =
                    library_with(case.program, &failing, &subject, &options, outfmt);
                assert_eq!(
                    error.as_deref(),
                    Some(cli_error.as_str()),
                    "{name} {outfmt}"
                );
                assert!(report == cli_report, "{name} {outfmt}: the reports");
            }
            for query in [duplicate.clone(), with_bad_record(duplicate)] {
                let (cli_report, cli_error) =
                    command_line_with(case.program, &query, &subject, &options, "6");
                let (rows, error) = batches_with(case.program, &query, &subject, &options);
                assert_eq!(error, cli_error, "{name}: run_batches");
                assert_eq!(
                    rows,
                    cli_report.split(|&byte| byte == b'\n').count() - 1,
                    "{name}: run_batches"
                );
            }
        }
    }
}

/// The command line's run of `query` against `subject` with `options` in outfmt 6 rows HSPs and
/// then fails with an error that holds `error`; the library's run and `run_batches` give the same
/// rows and error (`name` tells the case).
fn ends_as_the_command_line(
    name: &str,
    program: Program,
    query: &[u8],
    subject: &[u8],
    options: &[&str],
    error: &str,
) {
    let (cli_report, cli_error) = command_line_with(program, query, subject, options, "6");
    let cli_error = cli_error.unwrap_or_else(|| panic!("{name}: the command line fails"));
    assert!(cli_error.contains(error), "{name}: {cli_error}");
    let cli_rows = cli_report.split(|&byte| byte == b'\n').count() - 1;
    assert!(cli_rows > 0, "{name}: the batches before");
    let (report, library_error) = library_with(program, query, subject, options, "6");
    assert_eq!(library_error.as_deref(), Some(cli_error.as_str()), "{name}");
    assert!(report == cli_report, "{name}: the reports");
    let (rows, batches_error) = batches_with(program, query, subject, options);
    assert_eq!(
        (rows, batches_error.as_deref()),
        (cli_rows, Some(cli_error.as_str())),
        "{name}: run_batches"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652 (a batch whose
// queries all have no data stops the run, cited at `query_batch::no_data_batch_error`)
// A batch that NCBI fails is never searched and adds nothing to the scope of the run's queries,
// so a query with the id of an earlier one in that batch is not LOSAT's refusal: the batches
// before are reported and the run ends with NCBI's error, in the library's run and in
// `run_batches` as on the command line (review M-b of line WQA's D10, 2026-11-20).
#[test]
fn a_batch_that_ncbi_fails_adds_nothing_to_the_scope() {
    let protein = |length| letters(b"ACDEFGHIKLMNPQRSTVWY", length, 7);
    let options = ["-parse_deflines"];
    let named = |name: &str, sequence: &[u8]| fasta_named(&[name.to_string()], sequence);
    // A query `q1` longer than the first batch, then a second `q1` without residues, alone in the
    // next batch: `Sequence contains no data` after the first batch's report.
    for (program, query, subject) in [
        (Program::Blastp, protein(10500), protein(2000)),
        (Program::Tblastn, protein(20500), coding(&protein(2000))),
        (
            Program::Tblastx,
            coding(&protein(3500)),
            coding(&protein(700)),
        ),
    ] {
        let mut text = named("q1", &query);
        text.extend_from_slice(b">q1\n");
        ends_as_the_command_line(
            program.name(),
            program,
            &text,
            &named("s1", &subject),
            &options,
            "Sequence contains no data",
        );
    }
    // TBLASTX stops a batch whose queries all have fewer than three letters when one has none,
    // which BLASTP's rule (every query without letters) does not: the last batch, with the second
    // `q1`, fails on its search, after the first batch's report.
    let mut text = named("q1", &coding(&protein(3500)));
    text.extend_from_slice(b">e\n>t\nAC\n>q1\nA\n");
    ends_as_the_command_line(
        "tblastx: a batch of short queries",
        Program::Tblastx,
        &text,
        &named("s1", &coding(&protein(700))),
        &options,
        "Sequence contains no data",
    );
    // psiblast's batches hold one query with letters and the queries without letters before it:
    // `q1`, then `q2` (no residues) with `q3`, which `CPsiBlast` refuses, before the second `q1`.
    let mut text = named("q1", &protein(300));
    text.extend_from_slice(b">q2\n");
    text.extend_from_slice(&named("q3", &protein(300)));
    text.extend_from_slice(&named("q1", &protein(300)));
    ends_as_the_command_line(
        "psiblast",
        Program::Psiblast,
        &text,
        &named("s1", &protein(300)),
        &options,
        "PSI-BLAST only accepts one protein sequence as query",
    );
}

// TBLASTX refuses two queries with one id at the batch of the second, after the reports of the
// batches before (before line WQA D10 round 4, 2026-11-20: without them); NCBI's error of that
// batch's set-up for its subjects comes first (review M-b, 2026-11-20; since round 5 from the
// set-up of the batch's records with new ids, before its search, `query_batch::refused_batch_error`),
// in `run_batches` as on the command line.
#[test]
fn tblastx_refuses_two_queries_of_one_id_after_their_batch() {
    let protein = |length| letters(b"ACDEFGHIKLMNPQRSTVWY", length, 7);
    let options = ["-parse_deflines"];
    let query = coding(&protein(1100));
    let names: Vec<String> = ["q1", "q2", "q3", "q4", "q5", "q1"]
        .iter()
        .map(|name| name.to_string())
        .collect();
    let text = fasta_named(&names, &query);
    // Subjects with letters: the reports of the first batch (`q1` to `q4`), then the refusal.
    let subject = fasta_named(&["s1".to_string()], &coding(&protein(1100)));
    ends_as_the_command_line(
        "tblastx: the refusal",
        Program::Tblastx,
        &text,
        &subject,
        &options,
        "two query records have the Seq-id lcl|q1",
    );
    // A subject without letters: NCBI's error of the first batch's search comes first.
    let empty = b">s1\n".to_vec();
    let (cli_report, cli_error) = command_line_with(Program::Tblastx, &text, &empty, &options, "6");
    let cli_error = cli_error.expect("the search's error");
    assert!(cli_error.contains("average subject length"), "{cli_error}");
    assert!(cli_report.is_empty());
    let (rows, error) = batches_with(Program::Tblastx, &text, &empty, &options);
    assert_eq!((rows, error.as_deref()), (0, Some(cli_error.as_str())));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2783-2799 (a query whose letters
// give no Karlin-Altschul parameters, cited at `report::query_warnings::INVALID_QUERY_MESSAGE`)
// An error of NCBI's in the search of a batch comes before LOSAT's refusal of two queries with one
// id in a later batch, on the command line, in the library's run and in `run_batches` (line WQA D10
// round 4, 2026-11-20; NCBI BLAST+ 2.17.0, probes `fix-reader-error/r3-psiblast-probe` and
// `r4-probe`): psiblast's first query of only `X` stops the run there with no report (exit 3), as
// without the second `q1`, and so do subjects too short for `-searchsp` in BLASTP and in BLASTX's
// BLAST XML; before, LOSAT refused the second `q1`. Without such an error the run writes the
// reports of the batches before the refusal.
#[test]
fn an_error_of_ncbi_before_the_refusal_comes_first() {
    let options = ["-parse_deflines"];
    let protein = b"WNQQGYFFYHSSPALDGYTNTQRWQFAGQPMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
    let named = |name: &str, sequence: &[u8]| fasta_named(&[name.to_string()], sequence);
    let subject = named("s1", protein);
    let mut text = named("q1", &[b'X'; 30]);
    text.extend_from_slice(&named("q2", protein));
    text.extend_from_slice(&named("q1", protein));
    // NCBI's stdout (empty) and stderr of `psiblast -query dup.fa -subject s.fa -parse_deflines
    // -outfmt 6` (exit 3).
    // NCBI's text in LOSAT's words (line NW2).
    let ncbi = crate::wording::losat_text("BLAST engine error: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options").into_owned();
    let ncbi = ncbi.as_str();
    let (cli_report, cli_error) =
        command_line_with(Program::Psiblast, &text, &subject, &options, "6");
    assert_eq!(cli_error.as_deref(), Some(ncbi));
    assert!(cli_report.is_empty(), "NCBI writes no report");
    let (report, error) = library_with(Program::Psiblast, &text, &subject, &options, "6");
    assert_eq!((report.is_empty(), error.as_deref()), (true, Some(ncbi)));
    let (rows, error) = batches_with(Program::Psiblast, &text, &subject, &options);
    assert_eq!((rows, error.as_deref()), (0, Some(ncbi)));
    // `q1` with letters: the report of `q1` and `q2`, then the refusal of the second `q1`.
    let mut text = named("q1", protein);
    text.extend_from_slice(&named("q2", protein));
    text.extend_from_slice(&named("q1", protein));
    ends_as_the_command_line(
        "psiblast: the refusal after the batches before",
        Program::Psiblast,
        &text,
        &subject,
        &options,
        "two query records have the Seq-id lcl|q1",
    );
    // NCBI's error of the first batch's search (subjects of 2, 0 and 0 letters with `-searchsp
    // 1000`: `The average subject length is too short`, exit 3, no report; NCBI BLAST+ 2.17.0, probe
    // `fix-reader-error/r4-probe` `bp_avg`, `nuc_avg`): BLASTP before the second `q1` of the next
    // batch, and BLASTX in BLAST XML, whose read before the search met the refusal first.
    // NCBI's text in LOSAT's words (line NW2).
    let ncbi = crate::wording::losat_text("BLAST engine error: The average subject length is too short").into_owned();
    let ncbi = ncbi.as_str();
    let options = ["-parse_deflines", "-searchsp", "1000"];
    let protein = |length| letters(b"ACDEFGHIKLMNPQRSTVWY", length, 7);
    let subjects = b">s1\nMK\n>s2\n>s3\n";
    for (program, first, other, outfmt) in [
        (Program::Blastp, protein(10500), protein(300), "6"),
        (
            Program::Blastx,
            coding(&protein(3500)),
            coding(&protein(300)),
            "5",
        ),
    ] {
        let name = program.name();
        let mut text = named("q1", &first);
        text.extend_from_slice(&named("q2", &other));
        text.extend_from_slice(&named("q1", &other));
        let (cli_report, cli_error) = command_line_with(program, &text, subjects, &options, outfmt);
        assert_eq!(cli_error.as_deref(), Some(ncbi), "{name}");
        assert!(cli_report.is_empty(), "{name}: NCBI writes no report");
        let (report, error) = library_with(program, &text, subjects, &options, outfmt);
        assert_eq!(
            (report.is_empty(), error.as_deref()),
            (true, Some(ncbi)),
            "{name}"
        );
        let (rows, error) = batches_with(program, &text, subjects, &options);
        assert_eq!((rows, error.as_deref()), (0, Some(ncbi)), "{name}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 (`BLAST_GapAlignSetUp`,
// cited at `query_batch::refused_batch_error`)
// The batch of a repeated id stops at one point in every program, on the command line, in the
// library's run and in `run_batches` (line WQA D10 round 5, 2026-11-20; review 3 of D10, M1): where
// the batch's records with new ids have a valid context, NCBI's error of the batch's set-up for
// its subjects comes in place of LOSAT's refusal, as NCBI BLAST+ 2.17.0 gives it on that very input
// (`q1`, `q2`, `q1` in one batch; subjects of 2, 0 and 0 letters with `-searchsp 1000`: exit 3,
// rows `db.pd.q1c.*.avgdup.subj.6`); before, BLASTP and TBLASTN refused. Where they have none
// (a first `q1` of only `N`), NCBI sets up no valid context and exits 0 with the repeat read as the
// first's letters, so LOSAT refuses; before, BLASTX gave the subjects' error.
#[test]
fn the_batch_of_a_repeated_id_stops_at_one_point_in_every_program() {
    // NCBI's text in LOSAT's words (line NW2).
    let ncbi = crate::wording::losat_text("BLAST engine error: The average subject length is too short").into_owned();
    let ncbi = ncbi.as_str();
    let refusal = "two query records have the Seq-id lcl|q1";
    let options = ["-parse_deflines", "-searchsp", "1000"];
    let protein = |length, seed| letters(b"ACDEFGHIKLMNPQRSTVWY", length, seed);
    let nucleotide = |length, seed| letters(b"ACGT", length, seed);
    let named = |name: &str, sequence: &[u8]| fasta_named(&[name.to_string()], sequence);
    let (protein_subjects, nucleotide_subjects) = (b">s1\nMK\n>s2\n>s3\n", b">s1\nAC\n>s2\n>s3\n");
    for program in [
        Program::Blastp,
        Program::Tblastn,
        Program::Blastn,
        Program::Blastx,
        Program::Tblastx,
    ] {
        let name = program.name();
        let protein_query = matches!(program, Program::Blastp | Program::Tblastn);
        let record = |seed| {
            if protein_query {
                protein(300, seed)
            } else {
                nucleotide(900, seed)
            }
        };
        let subjects: &[u8] = match program {
            Program::Blastp | Program::Blastx => protein_subjects,
            _ => nucleotide_subjects,
        };
        let mut text = named("q1", &record(11));
        text.extend_from_slice(&named("q2", &record(12)));
        text.extend_from_slice(&named("q1", &record(13)));
        let (cli_report, cli_error) = command_line_with(program, &text, subjects, &options, "6");
        assert_eq!(cli_error.as_deref(), Some(ncbi), "{name}");
        assert!(cli_report.is_empty(), "{name}: NCBI writes no report");
        let (report, error) = library_with(program, &text, subjects, &options, "6");
        assert_eq!(
            (report.is_empty(), error.as_deref()),
            (true, Some(ncbi)),
            "{name}"
        );
        let (rows, error) = batches_with(program, &text, subjects, &options);
        assert_eq!((rows, error.as_deref()), (0, Some(ncbi)), "{name}");
    }
    // BLASTX: a first `q1` of only `N` and a `q2` of only `N` (no valid context; NCBI exits 0).
    let mut text = named("q1", &[b'N'; 90]);
    text.extend_from_slice(&named("q2", &[b'N'; 120]));
    text.extend_from_slice(&named("q1", &nucleotide(60, 14)));
    let (_, cli_error) = command_line_with(Program::Blastx, &text, protein_subjects, &options, "6");
    let cli_error = cli_error.expect("the refusal");
    assert!(cli_error.contains(refusal), "{cli_error}");
    let (_, error) = library_with(Program::Blastx, &text, protein_subjects, &options, "6");
    assert_eq!(error.as_deref(), Some(cli_error.as_str()));
    let (rows, error) = batches_with(Program::Blastx, &text, protein_subjects, &options);
    assert_eq!((rows, error.as_deref()), (0, Some(cli_error.as_str())));
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
    assert!(error.message().starts_with(bad_line().as_bytes()), "{error}");
    let (_, cli_error) = command_line(case.program, &query, &case.subject);
    assert_eq!(
        String::from_utf8_lossy(error.message()).trim_end(),
        cli_error.unwrap()
    );

    // An error in the first batch: no batch is searched.
    let query = with_bad_record(fasta(2, &case.query));
    let read = read_input_records(&search, &query, &case.subject).unwrap();
    let batches: Vec<_> = run_batches(&search, read.inputs()).unwrap().collect();
    assert_eq!(batches.len(), 1);
    assert!(batches[0].is_err());
}

// `run_batches` cuts its batches as NCBI's reader counts the records (`batches::batch_count`,
// `Batches::next_range`; review 5 of line WQA D10, L1): with `-parse_deflines` and `-subject`, a
// repeated id counts its first record's letters. BLASTP's batches hold 10000 letters: a `q1` of
// 9950 and a second `q1` of 10 (g2a's form) fill the first batch, which LOSAT refuses; a `q1` of
// 100 and a second `q1` of 9950 (g2b's) count 200, so the batch reads the line that the reader
// refuses, the reader's error. Counted by their own letters, the two would swap.
#[test]
fn run_batches_counts_a_repeated_id_by_its_first_records_letters() {
    let protein = |length, seed| letters(b"ACDEFGHIKLMNPQRSTVWY", length, seed);
    let named = |name: &str, sequence: &[u8]| fasta_named(&[name.to_string()], sequence);
    let subject = named("s1", &protein(200, 21));
    let options = ["-parse_deflines"];
    for (first, second, refused) in [(9950, 10, true), (100, 9950, false)] {
        let mut text = named("q1", &protein(first, 22));
        text.extend_from_slice(&named("q1", &protein(second, 23)));
        let text = with_bad_record(text);
        let name = format!("q1 of {first}, q1 of {second}");
        let (_, cli_error) = command_line_with(Program::Blastp, &text, &subject, &options, "6");
        let cli_error = cli_error.expect("an error");
        if refused {
            assert!(
                cli_error.contains("two query records have the Seq-id lcl|q1"),
                "{name}: {cli_error}"
            );
        } else {
            assert!(cli_error.starts_with(&bad_line()), "{name}: {cli_error}");
        }
        let (rows, error) = batches_with(Program::Blastp, &text, &subject, &options);
        assert_eq!(
            (rows, error.as_deref()),
            (0, Some(cli_error.as_str())),
            "{name}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
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
