//! The library API (`losat::api`, docs/losat_library_surface_plan.md §6): a search of records in
//! memory writes the bytes of the command line, one report per `-outfmt`, with the command line's
//! warnings; the typed values are the columns of `-outfmt 6` and the `Statistics` of outfmt 5; two
//! searches at once give the reports of each alone; the example writes nothing to standard error;
//! the options the API owns are refused.
//!
//! The command line is the oracle here: its reports are the frozen NCBI bytes of the regression
//! fixtures (`fast output regressions`), so a report equal to the command line's is NCBI's.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::library_examples::example_binary;

mod temp_support;
use temp_support::TempPath;

use losat::api::{
    describe, read_fasta, read_input_file, read_inputs, run, run_batches, run_blastdbcmd,
    run_collect, run_convert2blastmask, run_dustmasker, run_makeblastdb, run_makembindex,
    run_segmasker, run_to_file, run_windowmasker, validate, Error, FastaRecord, HspRecord, Inputs,
    Outcome, Program, Role, Strand, ValueKind, API_OWNED,
};

/// One search: a program, two fixtures of `tests/fasta`, options and the formats to compare.
struct Case {
    program: Program,
    query: &'static str,
    subject: &'static str,
    options: &'static [&'static str],
    formats: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        program: Program::Blastp,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_subject.faa",
        options: &["-evalue", "1e-5"],
        formats: &["6", "0", "7"],
    },
    // Two hit list sizes (`-num_alignments` without `-max_target_seqs`): two searches, each
    // format as the command line writes it alone.
    Case {
        program: Program::Blastp,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_subject.faa",
        options: &["-num_alignments", "1"],
        formats: &["0", "6"],
    },
    Case {
        program: Program::Blastn,
        query: "small_test.fasta",
        subject: "small_test.fasta",
        options: &["-task", "megablast"],
        formats: &["6", "0", "7"],
    },
    Case {
        program: Program::Tblastx,
        query: "outfmt0/tblastx_many_query.fasta",
        subject: "outfmt0/tblastx_many_subject.fasta",
        options: &["-query_gencode", "11", "-db_gencode", "11"],
        formats: &["6", "0", "7"],
    },
    Case {
        program: Program::Tblastn,
        query: "outfmt0/e2d_t_pq.faa",
        subject: "outfmt0/e2e_tblastn_subject.fna",
        options: &[],
        formats: &["6", "0", "7"],
    },
    Case {
        program: Program::Blastx,
        query: "outfmt0/sml_bx_split_q.fna",
        subject: "outfmt0/sml_bp_split_s.faa",
        options: &["-num_threads", "2"],
        formats: &["6", "0", "7"],
    },
    // A protein file read as a nucleotide query: the reader's warning is the command line's
    // standard error, and the API's diagnostics.
    Case {
        program: Program::Blastx,
        query: "outfmt0/sml_bx_tie_s.faa",
        subject: "outfmt0/sml_bx_tie_s.faa",
        options: &[],
        formats: &["6", "0"],
    },
    // psiblast's first round: blastp's search, a batch of one query, psiblast's report.
    Case {
        program: Program::Psiblast,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_subject.faa",
        options: &["-evalue", "1e-5"],
        formats: &["6", "0", "7"],
    },
    // psiblast's rounds (line RG item Q5): NCBI 2.17.0 finds 300 subjects in round 1 and 299 in
    // round 2, then writes `Search has CONVERGED!`.
    Case {
        program: Program::Psiblast,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_subject.faa",
        options: &["-evalue", "1e-5", "-num_iterations", "3"],
        formats: &["6", "0", "7"],
    },
];

/// The psiblast rounds case above.
const ROUNDS_CASE: usize = 8;

/// The warnings case above.
const WARNINGS_CASE: usize = 6;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

fn temp_path(name: &str) -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_library_api_{}_{}_{name}",
        std::process::id(),
        CALLS.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A command of this crate's binaries without LOSAT's own variables, which change LOSAT's behavior
/// (`LOSAT_BLASTX_WORKER_LOG`, which the `blastx-worker-probe` feature requires, stays). NCBI's
/// variables (`NCBI_*`, `BLASTDB`, ...) stay: the library's calls in this process read them at
/// their start, as the command line does, so that both see the same environment.
fn command(binary: &str) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(&name);
        }
    }
    command
}

/// The command line's report for one `-outfmt`, and its standard error.
fn cli(case: &Case, query: &Path, subject: &Path, outfmt: &str) -> (Vec<u8>, Vec<u8>) {
    let out = temp_path("cli.out");
    let output = command(env!("CARGO_BIN_EXE_losat"))
        .arg(case.program.name())
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(case.options)
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-out")
        .arg(&out)
        .output()
        .expect("run the command line");
    assert!(
        output.status.success(),
        "{} -outfmt {outfmt} failed: {}",
        case.program,
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = std::fs::read(&out).expect("the command line's report");
    let _ = std::fs::remove_file(&out);
    (bytes, output.stderr)
}

/// The library's error for the case with `outfmt`.
fn api_error(case: &Case, query: &Path, subject: &Path, outfmt: &str) -> Error {
    let queries = read_fasta(case.program, Role::Query, &std::fs::read(query).unwrap()).unwrap();
    let subjects = read_fasta(
        case.program,
        Role::Subject,
        &std::fs::read(subject).unwrap(),
    )
    .unwrap();
    let search = validate(case.program, case.options).unwrap();
    run_collect(&search, Inputs::new(&queries, &subjects), &[outfmt]).unwrap_err()
}

/// The library's outcome for the case, with the fixtures' paths as the inputs' names (what the
/// command line shows for its `-query` and `-subject` arguments).
fn api(case: &Case, query: &Path, subject: &Path, extra: &[&str]) -> Outcome {
    let query_bytes = std::fs::read(query).expect("query fixture");
    let subject_bytes = std::fs::read(subject).expect("subject fixture");
    let query_records = read_fasta(case.program, Role::Query, &query_bytes).expect("query records");
    let subject_records =
        read_fasta(case.program, Role::Subject, &subject_bytes).expect("subject records");
    let options: Vec<&str> = case
        .options
        .iter()
        .copied()
        .chain(extra.iter().copied())
        .collect();
    let search = validate(case.program, &options).expect("the options of the case");
    let inputs = Inputs::new(&query_records, &subject_records)
        .named(query.to_str().unwrap(), subject.to_str().unwrap());
    run_collect(&search, inputs, case.formats)
        .unwrap_or_else(|error| panic!("{} {:?}: {error}", case.program, case.options))
}

/// The line ranges of a tabular report: one per row (not the comments of outfmt 7, nor the empty
/// line and `Search has CONVERGED!` that psiblast writes after a query's rounds converged).
fn rows(report: &[u8]) -> Vec<std::ops::Range<usize>> {
    let mut rows = Vec::new();
    let mut start = 0;
    for (at, &byte) in report.iter().enumerate() {
        let line = &report[start..at];
        if byte == b'\n'
            && !line.starts_with(b"#")
            && !line.is_empty()
            && line != b"Search has CONVERGED!"
        {
            rows.push(start..at + 1);
        }
        if byte == b'\n' {
            start = at + 1;
        }
    }
    rows
}

#[test]
fn the_reports_and_the_diagnostics_are_the_command_lines() {
    let (mut shown_links, mut shown_methods) = (0, 0);
    for case in CASES {
        let (query, subject) = (fixture(case.query), fixture(case.subject));
        let outcome = api(case, &query, &subject, &[]);
        assert_eq!(outcome.reports.len(), case.formats.len());
        let context = format!("{} {:?}", case.program, case.options);
        // The diagnostics are the first search's: with `-num_alignments` or `-num_descriptions`
        // the formats split by hit list size and the tabular formats are searched first
        // (`search_groups`), so they are the command line's standard error for the first tabular
        // format; otherwise for the first format.
        let split = case
            .options
            .iter()
            .any(|option| matches!(*option, "-num_alignments" | "-num_descriptions"));
        let diagnostics_format = case
            .formats
            .iter()
            .position(|outfmt| !split || outfmt.parse::<u32>().is_ok_and(|number| number > 4))
            .unwrap();
        let mut first_stderr = None;
        for (position, outfmt) in case.formats.iter().enumerate() {
            let (bytes, stderr) = cli(case, &query, &subject, outfmt);
            assert_eq!(outcome.reports[position].outfmt, *outfmt);
            assert!(
                outcome.reports[position].bytes == bytes,
                "{context} -outfmt {outfmt}: the report differs from the command line's\n--- api\n{}\n--- cli\n{}",
                String::from_utf8_lossy(&outcome.reports[position].bytes),
                String::from_utf8_lossy(&bytes)
            );
            if position == diagnostics_format {
                first_stderr = Some(stderr);
            }
        }
        assert_eq!(
            outcome.diagnostics,
            first_stderr.unwrap_or_default(),
            "{context}: diagnostics"
        );
        // Every HSP has its row in the tabular report of the records' search, and every row is
        // one HSP's: the rows of outfmt 6, in order, are the records' spans.
        let tabular = case
            .formats
            .iter()
            .position(|outfmt| *outfmt == "6")
            .unwrap();
        let report = &outcome.reports[tabular].bytes;
        let spans: Vec<_> = outcome.reports[tabular]
            .spans
            .iter()
            .map(|span| span.clone().expect("a row in outfmt 6"))
            .collect();
        assert_eq!(spans, rows(report), "{context}: the rows of outfmt 6");
        for report in &outcome.reports {
            assert_eq!(report.spans.len(), outcome.hsps.len(), "{context}");
            assert_eq!(report.subject_headings.len(), outcome.hsps.len());
        }
        for (index, hsp) in outcome.hsps.iter().enumerate() {
            assert_eq!(hsp.index, index);
            let row = &report[spans[index].clone()];
            let fields: Vec<&[u8]> = row[..row.len() - 1].split(|&b| b == b'\t').collect();
            assert_eq!(fields.len(), 12, "{context}: row {index}");
            // The printed coordinates are the record's.
            let number = |field: &[u8]| {
                std::str::from_utf8(field)
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            };
            assert_eq!(
                (
                    number(fields[6]),
                    number(fields[7]),
                    number(fields[8]),
                    number(fields[9])
                ),
                (hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end),
                "{context}: coordinates of row {index}"
            );
            // The pairwise sections and headings, where format 0 is in the records' search: an
            // HSP has a section when the pairwise report shows its subject (the first
            // `-num_alignments` subjects, 250 by default), and then its subject's heading.
            if let Some(pairwise) = case.formats.iter().position(|outfmt| *outfmt == "0") {
                let report = &outcome.reports[pairwise];
                if case.options.contains(&"-num_alignments") {
                    assert!(
                        report.spans[index].is_none(),
                        "{context}: another search's report"
                    );
                } else {
                    let bytes = &report.bytes;
                    match (&report.spans[index], &report.subject_headings[index]) {
                        (Some(section), Some(heading)) => {
                            assert!(
                                bytes[section.clone()].starts_with(b" Score ="),
                                "{context}: section {index}"
                            );
                            assert!(
                                bytes[heading.clone()].starts_with(b">"),
                                "{context}: heading {index}"
                            );
                            assert!(heading.end <= section.start, "{context}: heading {index}");
                            let text = String::from_utf8_lossy(&bytes[section.clone()]);
                            // The score line shows the linked HSPs and the composition method.
                            match hsp.sum_n.filter(|linked| *linked > 1) {
                                Some(linked) => {
                                    assert!(
                                        text.contains(&format!("Expect({linked})")),
                                        "{context}: {text}"
                                    );
                                    shown_links += 1;
                                }
                                None => assert!(!text.contains("Expect("), "{context}: {text}"),
                            }
                            let method = match hsp.comp_adjust_method {
                                Some(1) => ", Method: Composition-based stats.",
                                Some(2) => ", Method: Compositional matrix adjust.",
                                _ => "",
                            };
                            if method.is_empty() {
                                assert!(!text.contains("Method:"), "{context}: {text}");
                            } else {
                                assert!(text.contains(method), "{context}: {text}");
                                shown_methods += 1;
                            }
                        }
                        (None, None) => {}
                        other => panic!("{context}: HSP {index} has a section without a heading or the reverse: {other:?}"),
                    }
                }
            }
        }
        // Every section of the pairwise report belongs to one HSP of the records' search.
        if let Some(pairwise) = case.formats.iter().position(|outfmt| *outfmt == "0") {
            if !case.options.contains(&"-num_alignments") {
                let bytes = &outcome.reports[pairwise].bytes;
                let shown = bytes
                    .split(|&b| b == b'\n')
                    .filter(|line| line.starts_with(b" Score ="))
                    .count();
                let sections = outcome.reports[pairwise].spans.iter().flatten().count();
                assert_eq!(sections, shown, "{context}: sections of outfmt 0");
                assert_eq!(
                    sections > 0,
                    !outcome.hsps.is_empty(),
                    "{context}: a case with hits has sections"
                );
            }
        }
        if !outcome.hsps.is_empty() {
            assert!(outcome.hsps.iter().any(|hsp| hsp.rank == 0));
        }
    }
    assert!(shown_links > 0, "a case with linked HSPs");
    assert!(
        shown_methods > 0,
        "a case with composition-based statistics"
    );
    // The warnings case has warnings, so the comparison above checked them.
    let case = &CASES[WARNINGS_CASE];
    let outcome = api(case, &fixture(case.query), &fixture(case.subject), &[]);
    assert!(
        String::from_utf8_lossy(&outcome.diagnostics).contains("Ignoring invalid residues"),
        "{}",
        String::from_utf8_lossy(&outcome.diagnostics)
    );
}

/// Searches with the coordinates and strands that the plain cases do not reach: intervals of the
/// records (`-query_loc`, `-subject_loc`) and the minus strand of a nucleotide query.
const RANGE_CASES: &[Case] = &[
    Case {
        program: Program::Blastn,
        query: "small_test.fasta",
        subject: "small_test.fasta",
        options: &["-task", "blastn", "-strand", "minus", "-evalue", "1"],
        formats: &[],
    },
    Case {
        program: Program::Blastn,
        query: "small_test.fasta",
        subject: "small_test.fasta",
        options: &[
            "-task",
            "megablast",
            "-query_loc",
            "101-3000",
            "-subject_loc",
            "201-5000",
        ],
        formats: &[],
    },
    Case {
        program: Program::Blastp,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_subject.faa",
        options: &["-evalue", "1e-5", "-query_loc", "11-150"],
        formats: &[],
    },
    Case {
        program: Program::Tblastx,
        query: "outfmt0/tblastx_many_query.fasta",
        subject: "outfmt0/tblastx_many_subject.fasta",
        options: &[
            "-query_gencode",
            "11",
            "-db_gencode",
            "11",
            "-query_loc",
            "11-170",
            "-subject_loc",
            "5-170",
        ],
        formats: &[],
    },
    Case {
        program: Program::Tblastn,
        query: "outfmt0/e2d_t_pq.faa",
        subject: "outfmt0/e2e_tblastn_subject.fna",
        options: &["-query_loc", "11-200", "-subject_loc", "101-6000"],
        formats: &[],
    },
    Case {
        program: Program::Blastx,
        query: "outfmt0/sml_bx_split_q.fna",
        subject: "outfmt0/sml_bp_split_s.faa",
        options: &["-query_loc", "12001-15000"],
        formats: &[],
    },
];

/// Every typed column of the HSP records, in the order of `check_columns`.
const COLUMNS: &str = "6 length nident positive mismatch gapopen gaps qlen slen qstart qend \
     sstart send qframe sframe sstrand score bitscore evalue qaccver saccver pident";

/// One row of `COLUMNS` against the record of its HSP.
fn check_columns(context: &str, hsp: &HspRecord, row: &[u8]) {
    let row = std::str::from_utf8(row).unwrap().trim_end();
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields.len(), 21, "{context}: {row}");
    assert_eq!(
        (hsp.query_id.as_bytes(), hsp.subject_id.as_bytes()),
        (fields[18].as_bytes(), fields[19].as_bytes()),
        "{context}: the ids of {row}"
    );
    assert_eq!(
        format!("{:.3}", hsp.percent_identity()),
        fields[20],
        "{context}: pident of {row}"
    );
    let number = |at: usize| fields[at].parse::<i64>().unwrap();
    assert_eq!(
        (
            hsp.length as i64,
            hsp.num_ident as i64,
            hsp.num_positives as i64,
            hsp.mismatch,
            hsp.gap_open as i64,
            hsp.gaps as i64,
        ),
        (
            number(0),
            number(1),
            number(2),
            number(3),
            number(4),
            number(5)
        ),
        "{context}: the counts of {row}"
    );
    assert_eq!(
        (hsp.query_length as i64, hsp.subject_length as i64),
        (number(6), number(7)),
        "{context}: the lengths of {row}"
    );
    assert_eq!(
        (
            hsp.q_start as i64,
            hsp.q_end as i64,
            hsp.s_start as i64,
            hsp.s_end as i64
        ),
        (number(8), number(9), number(10), number(11)),
        "{context}: the coordinates of {row}"
    );
    // A frame exists where the program translates; the columns print a number for the others.
    if let Some(frame) = hsp.query_frame {
        assert_eq!(i64::from(frame), number(12), "{context}: qframe of {row}");
    }
    if let Some(frame) = hsp.subject_frame {
        assert_eq!(i64::from(frame), number(13), "{context}: sframe of {row}");
    }
    let strand = match hsp.subject_strand {
        Some(Strand::Plus) => "plus",
        Some(Strand::Minus) => "minus",
        None => "N/A",
    };
    assert_eq!(strand, fields[14], "{context}: sstrand of {row}");
    assert_eq!(
        i64::from(hsp.raw_score),
        number(15),
        "{context}: score of {row}"
    );
    // The reports cut these to a few digits (`report::tabular`); a bit score in (99.9, 100) is
    // NCBI's ` 99` (`report::outfmt6::format_bitscore_ncbi`).
    let bits: f64 = fields[16].trim_start().parse().unwrap();
    assert!(
        (hsp.bit_score - bits).abs() <= 1.0 + hsp.bit_score / 1000.0,
        "{context}: bitscore of {row}: {}",
        hsp.bit_score
    );
    let evalue: f64 = fields[17].parse().unwrap();
    if evalue == 0.0 {
        assert!(hsp.e_value < 1e-179, "{context}: evalue of {row}");
    } else {
        let ratio = hsp.e_value / evalue;
        assert!(
            (0.5..=1.5).contains(&ratio),
            "{context}: evalue of {row}: {}",
            hsp.e_value
        );
    }
}

/// The letters of a row of an alignment, without its gaps, in upper case.
fn letters(row: &str) -> Vec<u8> {
    row.bytes()
        .filter(|&byte| byte != b'-')
        .map(|byte| byte.to_ascii_uppercase())
        .collect()
}

/// The reverse complement of nucleotide letters of the test fixtures, in upper case.
fn reverse_complement(letters: &[u8]) -> Vec<u8> {
    letters
        .iter()
        .rev()
        .map(|&letter| match letter.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// The intervals of a record against the letters of its rows: `query_range` and `subject_range`
/// hold the letters of the rows on the record's plus strand, and three letters of a translated
/// sequence hold one of its row.
fn check_ranges(
    context: &str,
    program: Program,
    hsp: &HspRecord,
    query: &FastaRecord,
    subject: &FastaRecord,
) {
    let (query_row, subject_row) = (
        letters(hsp.query_aligned.as_deref().expect("a query row")),
        letters(hsp.subject_aligned.as_deref().expect("a subject row")),
    );
    let upper = |record: &FastaRecord, range: std::ops::Range<usize>| -> Vec<u8> {
        record.seq()[range].to_ascii_uppercase()
    };
    let (query_range, subject_range) = (hsp.query_range(), hsp.subject_range());
    assert_eq!(query_range.end, hsp.q_start.max(hsp.q_end), "{context}");
    assert_eq!(subject_range.start + 1, hsp.s_start.min(hsp.s_end));
    match hsp.query_frame {
        Some(_) => assert_eq!(query_range.len(), 3 * query_row.len(), "{context}"),
        None => assert_eq!(upper(query, query_range), query_row, "{context}"),
    }
    match (hsp.subject_frame, hsp.subject_strand) {
        (Some(_), _) => assert_eq!(subject_range.len(), 3 * subject_row.len(), "{context}"),
        (None, Some(Strand::Minus)) => assert_eq!(
            reverse_complement(&upper(subject, subject_range)),
            subject_row,
            "{context}"
        ),
        (None, _) => assert_eq!(upper(subject, subject_range), subject_row, "{context}"),
    }
    // Which rows have frames and strands.
    let translated_query = matches!(program, Program::Blastx | Program::Tblastx);
    let translated_subject = matches!(program, Program::Tblastn | Program::Tblastx);
    assert_eq!(hsp.query_frame.is_some(), translated_query, "{context}");
    assert_eq!(hsp.subject_frame.is_some(), translated_subject, "{context}");
    let protein_query = matches!(
        program,
        Program::Blastp | Program::Psiblast | Program::Tblastn
    );
    let protein_subject = matches!(
        program,
        Program::Blastp | Program::Psiblast | Program::Blastx
    );
    assert_eq!(hsp.query_strand.is_none(), protein_query, "{context}");
    assert_eq!(hsp.subject_strand.is_none(), protein_subject, "{context}");
    if let Some(strand) = hsp.query_strand {
        assert_eq!(
            strand == Strand::Minus,
            hsp.q_start > hsp.q_end,
            "{context}"
        );
    }
}

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
/// A protein, and a nucleotide sequence that translates to it.
const PROTEIN: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
const NUCLEOTIDE: &str = "ATGAAAACCGCGTATATTGCGAAACAGCGCCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGCCAGCTG\
                          GAAGAACGCCTGGGCCTGATTGAAGTGCAGGCGCCGATTCTGAGCCGCGTGGGCGATGGCACCCAGGAT\
                          AACCTGAGCGGCGCGGAAAAAGCGGTGCAGGTGAAAGTGAAAGCGCTGCCGGATGCGCAG";

/// The letters that `program` reads for the input of `role`.
fn letters_of(program: Program, role: Role) -> &'static str {
    match (program, role) {
        (Program::Blastp | Program::Psiblast, _)
        | (Program::Blastx, Role::Subject)
        | (Program::Tblastn, Role::Query) => PROTEIN,
        _ => NUCLEOTIDE,
    }
}

/// The search of `case` with its options as written, `formats` asked for, from `query_bytes`.
fn search_case(case: &Case, query_bytes: &[u8], subject_bytes: &[u8], formats: &[&str]) -> Outcome {
    let queries = read_fasta(case.program, Role::Query, query_bytes).expect("query records");
    let subjects = read_fasta(case.program, Role::Subject, subject_bytes).expect("subject records");
    let search = validate(case.program, case.options).expect("the options of the case");
    run_collect(&search, Inputs::new(&queries, &subjects), formats)
        .unwrap_or_else(|error| panic!("{} {:?}: {error}", case.program, case.options))
}

#[test]
fn the_typed_values_are_the_columns_of_the_tabular_report() {
    let mut programs = std::collections::HashSet::new();
    for case in CASES.iter().chain(RANGE_CASES) {
        let context = format!("{} {:?}", case.program, case.options);
        let query_bytes = std::fs::read(fixture(case.query)).expect("query fixture");
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        let outcome = search_case(case, &query_bytes, &subject_bytes, &[COLUMNS]);
        // Every case but the warnings case, whose query is no nucleotide sequence.
        assert_eq!(
            outcome.hsps.is_empty(),
            std::ptr::eq(case, &CASES[WARNINGS_CASE]),
            "{context}: the case's hits"
        );
        let report = &outcome.reports[0];
        assert_eq!(
            rows(&report.bytes).len(),
            outcome.hsps.len(),
            "{context}: rows"
        );
        // The records with their records' letters: the whole records, as the reports count them.
        let queries = read_fasta(case.program, Role::Query, &query_bytes).unwrap();
        let subjects = read_fasta(case.program, Role::Subject, &subject_bytes).unwrap();
        for hsp in &outcome.hsps {
            let row = report.spans[hsp.index]
                .clone()
                .unwrap_or_else(|| panic!("{context}: HSP {} has no row", hsp.index));
            check_columns(&context, hsp, &report.bytes[row]);
            check_ranges(
                &context,
                case.program,
                hsp,
                &queries[hsp.q_idx],
                &subjects[hsp.s_idx],
            );
        }
        programs.insert(case.program);
    }
    assert_eq!(
        programs.len(),
        Program::ALL.len(),
        "every program is a case"
    );
}

/// `<Statistics_NAME>value` of an `Iteration` of outfmt 5.
fn statistic(iteration: &str, name: &str) -> f64 {
    let open = format!("<Statistics_{name}>");
    let start = iteration.find(&open).expect("the statistic") + open.len();
    let end = start + iteration[start..].find('<').unwrap();
    iteration[start..end].parse().unwrap()
}

/// Each query's statistics against the `Statistics` of its `Iteration`, the first `iterations`
/// of the report being those of `searched`, in order; a query without a valid context prints -1.
fn check_statistics(context: &str, outcome: &Outcome, xml: &str, searched: &[usize]) {
    let iterations: Vec<&str> = xml.split("<Iteration>").skip(1).collect();
    assert_eq!(iterations.len(), searched.len(), "{context}: iterations");
    for (iteration, &query) in iterations.iter().zip(searched) {
        let close = |value: f64, name: &str| {
            let printed = statistic(iteration, name);
            assert!(
                (value - printed).abs() <= value.abs() * 1e-4,
                "{context}: query {query}: {name}: {value} and {printed}"
            );
        };
        match outcome.query_stats[query] {
            Some(stats) => {
                assert_eq!(
                    f64::from(stats.hsp_len),
                    statistic(iteration, "hsp-len"),
                    "{context}: query {query}"
                );
                close(stats.eff_space as f64, "eff-space");
                close(stats.lambda, "lambda");
                close(stats.kappa, "kappa");
                close(stats.entropy, "entropy");
            }
            None => assert_eq!(
                statistic(iteration, "kappa"),
                -1.0,
                "{context}: query {query}"
            ),
        }
    }
}

/// The command line's outfmt 5 report of `case`.
fn xml_report(case: &Case, query: &Path, subject: &Path) -> String {
    String::from_utf8(cli(case, query, subject, "5").0).unwrap()
}

#[test]
fn the_statistics_are_those_of_the_xml_report() {
    // LOSAT's psiblast writes no outfmt 5.
    let cases = CASES.iter().chain(RANGE_CASES);
    for case in cases.filter(|case| case.program != Program::Psiblast) {
        let context = format!("{} {:?}", case.program, case.options);
        let query_bytes = std::fs::read(fixture(case.query)).expect("query fixture");
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        let queries = read_fasta(case.program, Role::Query, &query_bytes).unwrap();
        let typed = search_case(case, &query_bytes, &subject_bytes, &["6"]);
        assert_eq!(typed.query_stats.len(), queries.len(), "{context}");
        let xml = xml_report(case, &fixture(case.query), &fixture(case.subject));
        let searched: Vec<usize> = (0..queries.len()).collect();
        check_statistics(&context, &typed, &xml, &searched);
        assert!(
            typed.query_stats.iter().any(Option::is_some),
            "{context}: a query with statistics"
        );
    }
}

#[test]
fn a_run_of_batches_numbers_the_statistics_by_query() {
    // Copies of a query file past NCBI's batch size: each program searches the batches apart, and
    // every copy of a query has the numbers of the first.
    let blastn = Case {
        program: Program::Blastn,
        query: "small_test.fasta",
        subject: "small_test.fasta",
        options: &["-task", "blastn"],
        formats: &[],
    };
    let blastp = Case {
        program: Program::Blastp,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_query.faa",
        options: &[],
        formats: &[],
    };
    for (case, copies) in [
        (&blastn, 20),
        // Past two of BLASTP's query blocks (`LOSAT_X_BPQBLOCK`, 50,000 letters).
        (&blastp, 600),
        (&CASES[3], 200),
        (&CASES[4], 100),
        (&CASES[5], 2),
    ] {
        let context = format!("{} {:?}", case.program, case.options);
        let once = std::fs::read(fixture(case.query)).expect("query fixture");
        let copied = once.repeat(copies);
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        let queries = read_fasta(case.program, Role::Query, &copied).unwrap();
        let typed = search_case(case, &copied, &subject_bytes, &["6"]);
        let path = temp_path("copies.fa");
        std::fs::write(&path, &copied).unwrap();
        let xml = xml_report(case, &path, &fixture(case.subject));
        let _ = std::fs::remove_file(&path);
        let searched: Vec<usize> = (0..queries.len()).collect();
        check_statistics(&context, &typed, &xml, &searched);
        let per_copy = queries.len() / copies;
        for (query, stats) in typed.query_stats.iter().enumerate() {
            assert_eq!(*stats, typed.query_stats[query % per_copy], "{context}");
        }
        assert!(typed.query_stats[0].is_some(), "{context}");
    }
}

#[test]
fn a_query_without_a_valid_context_or_a_search_has_no_statistics() {
    let good = "ACGTTGCAAGGCTTAACCGGATATCGCGTACGATCGATTAGCTAGGATCTTAGGCTAACGTGCCATATTAGC\
                GGCATTCGATCGGATTACAGGCATCGATCGGCTAAGCTTAGGCATTCAGGCATTAGCGCATATCGATTAGGC";
    let fasta = format!(">tiny\nACG\n>good\n{good}\n>all_n\n{}\n", "N".repeat(150));
    let queries = read_fasta(Program::Blastn, Role::Query, fasta.as_bytes()).unwrap();
    let subject = format!(">s\n{good}\n");
    let subjects = read_fasta(Program::Blastn, Role::Subject, subject.as_bytes()).unwrap();
    let search = validate(Program::Blastn, &["-task", "blastn"]).unwrap();
    let outcome = run_collect(&search, Inputs::new(&queries, &subjects), &["6"]).unwrap();
    assert!(outcome.query_stats[1].is_some());
    assert_eq!(outcome.query_stats[2], None);
    // The query before the `-query_loc` interval is not searched and keeps its place.
    let search = validate(Program::Blastn, &["-task", "blastn", "-query_loc", "5-120"]).unwrap();
    let outcome = run_collect(&search, Inputs::new(&queries, &subjects), &["6"]).unwrap();
    assert_eq!(outcome.query_stats.len(), 3);
    assert_eq!(outcome.query_stats[0], None);
    assert!(outcome.query_stats[1].is_some());
    assert_eq!(outcome.query_stats[2], None);
    assert!(!outcome.hsps.is_empty());
    assert!(outcome.hsps.iter().all(|hsp| hsp.q_idx == 1));
}

#[test]
fn the_records_and_the_statistics_change_nothing_at_the_edges_of_the_inputs() {
    // `run_collect` asks each program for its hit list and the query statistics, `run` for
    // neither: without records, without residues, or with one of them missing, both give the
    // same reports and diagnostics, or the same error.
    for program in Program::ALL {
        let inputs = |role: Role| -> Vec<Vec<FastaRecord>> {
            let letters = letters_of(program, role);
            [
                String::new(),
                ">empty\n".to_string(),
                format!(">good\n{letters}\n"),
                format!(">empty\n>good\n{letters}\n>empty2\n"),
            ]
            .iter()
            .map(|fasta| read_fasta(program, role, fasta.as_bytes()).unwrap())
            .collect()
        };
        let search = validate(program, &[]).unwrap();
        let formats = ["6", "0"];
        for query in &inputs(Role::Query) {
            for subject in &inputs(Role::Subject) {
                let inputs = Inputs::new(query, subject);
                let context = format!("{program}: {} {}", query.len(), subject.len());
                let collected = run_collect(&search, inputs, &formats);
                let (mut tabular, mut pairwise, mut diagnostics) =
                    (Vec::new(), Vec::new(), Vec::new());
                let written = run(
                    &search,
                    inputs,
                    &mut [("6", &mut tabular), ("0", &mut pairwise)],
                    &mut diagnostics,
                );
                match (collected, written) {
                    (Ok(outcome), Ok(())) => {
                        assert!(outcome.reports[0].bytes == tabular, "{context}");
                        assert!(outcome.reports[1].bytes == pairwise, "{context}");
                        assert_eq!(outcome.diagnostics, diagnostics, "{context}");
                        assert_eq!(outcome.query_stats.len(), query.len(), "{context}");
                    }
                    (Err(collected), Err(written)) => assert_eq!(collected, written, "{context}"),
                    other => panic!("{context}: {other:?}"),
                }
            }
        }
    }
}

/// The batches of `case` one after the other against the whole search, and how many there are.
fn check_batches(case: &Case, query_bytes: &[u8], subject_bytes: &[u8]) -> usize {
    let context = format!("{} {:?}", case.program, case.options);
    let whole = search_case(case, query_bytes, subject_bytes, &["6"]);
    let queries = read_fasta(case.program, Role::Query, query_bytes).unwrap();
    let subjects = read_fasta(case.program, Role::Subject, subject_bytes).unwrap();
    let search = validate(case.program, case.options).unwrap();
    let batches = run_batches(&search, Inputs::new(&queries, &subjects)).unwrap();
    assert_eq!(batches.query_count(), queries.len(), "{context}");
    let (mut hsps, mut stats, mut rounds, mut diagnostics) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut count = 0;
    for batch in batches {
        let batch = batch.unwrap_or_else(|error| panic!("{context}: {error}"));
        assert_eq!(batch.query_indices.start, stats.len(), "{context}");
        assert_eq!(
            batch.query_stats.len(),
            batch.query_indices.len(),
            "{context}"
        );
        assert_eq!(batch.rounds.len(), batch.query_indices.len(), "{context}");
        // The views of a batch number the queries as the inputs do.
        for (query, view) in batch.query_indices.clone().zip(batch.queries()) {
            assert_eq!(view.q_idx, query, "{context}");
            assert!(view.hsps.iter().all(|hsp| hsp.q_idx == query), "{context}");
        }
        hsps.extend(batch.hsps);
        stats.extend(batch.query_stats);
        rounds.extend(batch.rounds);
        diagnostics.extend(batch.diagnostics);
        count += 1;
    }
    assert_eq!(stats.len(), queries.len(), "{context}: every query");
    assert!(hsps == whole.hsps, "{context}: the records");
    assert_eq!(stats, whole.query_stats, "{context}: the statistics");
    assert_eq!(rounds, whole.rounds, "{context}: the rounds");
    assert_eq!(
        String::from_utf8_lossy(&diagnostics),
        String::from_utf8_lossy(&whole.diagnostics),
        "{context}: the diagnostics"
    );
    count
}

#[test]
fn the_batches_one_after_the_other_are_the_whole_search() {
    for case in CASES.iter().chain(RANGE_CASES) {
        let query_bytes = std::fs::read(fixture(case.query)).expect("query fixture");
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        check_batches(case, &query_bytes, &subject_bytes);
    }
    // Copies of a query file past NCBI's batch size.
    let blastn = Case {
        program: Program::Blastn,
        query: "small_test.fasta",
        subject: "small_test.fasta",
        options: &["-task", "blastn"],
        formats: &[],
    };
    let blastp = Case {
        program: Program::Blastp,
        query: "outfmt0/e2e_many_query.faa",
        subject: "outfmt0/e2e_many_query.faa",
        options: &[],
        formats: &[],
    };
    for (case, copies) in [
        (&blastn, 20),
        (&blastp, 60),
        (&CASES[3], 200),
        (&CASES[4], 100),
        (&CASES[5], 2),
        // Two queries of two rounds each: a batch per query, its rounds in it.
        (&CASES[ROUNDS_CASE], 2),
    ] {
        let once = std::fs::read(fixture(case.query)).expect("query fixture");
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        let count = check_batches(case, &once.repeat(copies), &subject_bytes);
        assert!(count > 1, "{} batches of {}", count, case.program);
    }
}

#[test]
fn the_messages_of_the_run_come_once_with_the_first_batch() {
    // `-max_target_seqs 1` has a warning of the run, a subject without residues another.
    for program in Program::ALL {
        let options: &'static [&'static str] = &["-max_target_seqs", "1"];
        let case = Case {
            program,
            query: "",
            subject: "",
            options,
            formats: &[],
        };
        // A query longer than each program's batch size, twice, and one without residues
        // (with psiblast, that one and the next are a batch that `CPsiBlast` does not take).
        let letters = letters_of(program, Role::Query);
        let long = letters.repeat(25_000 / letters.len());
        let empty = if program == Program::Psiblast {
            ""
        } else {
            ">q3\n"
        };
        let query = format!(">q1\n{long}\n>q2\n{long}\n{empty}>q4\n{letters}\n");
        let letters = letters_of(program, Role::Subject);
        let subject = format!(">s1\n{letters}\n>empty\n");
        let batches = check_batches(&case, query.as_bytes(), subject.as_bytes());
        assert!(batches > 1, "{program}: {batches} batches");
        // The check above compared them with the search of all queries; they are there.
        let whole = search_case(&case, query.as_bytes(), subject.as_bytes(), &["6"]);
        let text = String::from_utf8_lossy(&whole.diagnostics).into_owned();
        assert_eq!(
            text.matches("Subject sequence contains no data").count(),
            1,
            "{program}: {text}"
        );
    }
}

#[test]
fn a_caller_stops_between_batches() {
    let case = &CASES[3];
    let once = std::fs::read(fixture(case.query)).expect("query fixture");
    let copied = once.repeat(200);
    let queries = read_fasta(case.program, Role::Query, &copied).unwrap();
    let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
    let subjects = read_fasta(case.program, Role::Subject, &subject_bytes).unwrap();
    let search = validate(case.program, case.options).unwrap();
    let mut batches = run_batches(&search, Inputs::new(&queries, &subjects)).unwrap();
    let first = batches.next().unwrap().unwrap();
    assert!(first.query_indices.end < batches.query_count());
    drop(batches);
    let whole = search_case(case, &copied, &subject_bytes, &["6"]);
    assert_eq!(first.hsps[..], whole.hsps[..first.hsps.len()]);
}

#[test]
fn the_views_hold_each_query_and_each_pair_together() {
    let mut pairs_seen = 0;
    for case in CASES.iter().chain(RANGE_CASES) {
        let context = format!("{} {:?}", case.program, case.options);
        let query_bytes = std::fs::read(fixture(case.query)).expect("query fixture");
        let subject_bytes = std::fs::read(fixture(case.subject)).expect("subject fixture");
        let outcome = search_case(case, &query_bytes, &subject_bytes, &["6"]);
        assert_eq!(outcome.queries().count(), outcome.query_stats.len());
        let mut seen = 0;
        let mut pairs = std::collections::HashSet::new();
        let mut first_rounds = std::collections::HashSet::new();
        for (q_idx, query) in outcome.queries().enumerate() {
            assert_eq!(query.q_idx, q_idx, "{context}");
            assert_eq!(
                query.stats,
                outcome.query_stats[q_idx].as_ref(),
                "{context}"
            );
            assert_eq!(query.round, None, "{context}");
            assert_eq!(outcome.query(q_idx), Some(query), "{context}");
            for subject in query.subjects() {
                // A pair is one run of the hit list in each round: it comes once a round.
                assert!(
                    pairs.insert((q_idx, subject.s_idx, subject.round)),
                    "{context}"
                );
                if first_rounds.insert((q_idx, subject.s_idx)) {
                    assert_eq!(query.subject(subject.s_idx), Some(subject), "{context}");
                }
                for hsp in subject.hsps {
                    assert_eq!(
                        (hsp.q_idx, hsp.s_idx, hsp.round),
                        (q_idx, subject.s_idx, subject.round),
                        "{context}"
                    );
                    assert_eq!(hsp.index, seen, "{context}: in the order of the hit list");
                    seen += 1;
                }
            }
            // The rounds of the query: its HSPs round by round, each with its statistics.
            let rounds = &outcome.rounds[q_idx];
            assert_eq!(rounds.stats[0], outcome.query_stats[q_idx], "{context}");
            let views: Vec<_> = query.rounds().collect();
            assert_eq!(views.len(), rounds.stats.len(), "{context}");
            let mut held = Vec::new();
            for (offset, view) in views.iter().enumerate() {
                let round = offset as u32 + 1;
                assert_eq!(view.round, Some(round), "{context}");
                assert_eq!(view.stats, rounds.stats[offset].as_ref(), "{context}");
                assert!(view.hsps.iter().all(|hsp| hsp.round == round), "{context}");
                assert_eq!(view.rounds().collect::<Vec<_>>(), [*view], "{context}");
                held.extend_from_slice(view.hsps);
            }
            assert!(held == query.hsps, "{context}: the rounds hold every HSP");
        }
        assert_eq!(seen, outcome.hsps.len(), "{context}: every HSP");
        assert_eq!(outcome.query(outcome.query_stats.len()), None);
        pairs_seen += pairs.len();
    }
    assert!(pairs_seen > 0);
}

#[test]
fn parsed_ids_show_as_on_the_command_line() {
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
        Program::Psiblast,
    ] {
        let query = format!(
            ">gb|AB000001.1| a query\n{0}\n>lcl|mine a second query\n{0}\n",
            letters_of(program, Role::Query)
        );
        let subject = format!(
            ">sp|P69905.2|HBA_HUMAN a subject\n{0}\n>gnl|db|42 a second subject\n{0}\n",
            letters_of(program, Role::Subject)
        );
        let (query_path, subject_path) = (temp_path("ids_q.fa"), temp_path("ids_s.fa"));
        std::fs::write(&query_path, &query).unwrap();
        std::fs::write(&subject_path, &subject).unwrap();
        let options: &'static [&'static str] = &["-parse_deflines"];
        let case = Case {
            program,
            query: "",
            subject: "",
            options,
            formats: &[],
        };
        let search = validate(program, options).unwrap();
        let (queries, subjects) =
            read_inputs(&search, query.as_bytes(), subject.as_bytes()).unwrap();
        let inputs = Inputs::new(&queries, &subjects)
            .named(query_path.to_str().unwrap(), subject_path.to_str().unwrap());
        let formats = ["6", "0", "7"];
        let outcome = run_collect(&search, inputs, &formats).unwrap();
        // The ids of the records are those of the default report's columns.
        let (ids, _) = cli(&case, &query_path, &subject_path, "6 qaccver saccver");
        let ids: Vec<String> = String::from_utf8(ids)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect();
        let shown: Vec<String> = outcome
            .hsps
            .iter()
            .map(|hsp| {
                let (query, subject) = (hsp.query_id.to_str(), hsp.subject_id.to_str());
                format!("{}\t{}", query.unwrap(), subject.unwrap())
            })
            .collect();
        assert_eq!(shown, ids, "{program:?}");
        for (position, outfmt) in formats.iter().enumerate() {
            let (bytes, _) = cli(&case, &query_path, &subject_path, outfmt);
            assert!(
                outcome.reports[position].bytes == bytes,
                "{program} -outfmt {outfmt}\n--- api\n{}\n--- cli\n{}",
                String::from_utf8_lossy(&outcome.reports[position].bytes),
                String::from_utf8_lossy(&bytes)
            );
        }
        assert!(!outcome.hsps.is_empty(), "{program}");
        // The records of `read_fasta` have no ids to show.
        let (plain_queries, plain_subjects) = (
            read_fasta(program, Role::Query, query.as_bytes()).unwrap(),
            read_fasta(program, Role::Subject, subject.as_bytes()).unwrap(),
        );
        let plain = Inputs::new(&plain_queries, &plain_subjects);
        let error = run_collect(&search, plain, &["6"]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{error}");
        assert!(run_batches(&search, plain).is_err());
        let _ = std::fs::remove_file(&query_path);
        let _ = std::fs::remove_file(&subject_path);
    }
}

/// The search of `program` with `options` on the FASTA bytes `query` and `subject`: its
/// `-outfmt 6` report is the command line's, and each HSP's ids are the bytes of its row's first
/// two columns (`Report::spans`) and those of `query_ids[q_idx]` and `subject_ids[s_idx]`. Returns
/// the query and subject numbers of the HSPs.
fn ids_of_the_rows(
    program: Program,
    options: &'static [&'static str],
    (query, subject): (&[u8], &[u8]),
    (query_ids, subject_ids): (&[&[u8]], &[&[u8]]),
) -> std::collections::BTreeSet<(usize, usize)> {
    let context = format!("{program} {options:?}");
    let (query_path, subject_path) = (temp_path("bytes_q.fa"), temp_path("bytes_s.fa"));
    std::fs::write(&query_path, query).unwrap();
    std::fs::write(&subject_path, subject).unwrap();
    let case = Case {
        program,
        query: "",
        subject: "",
        options,
        formats: &[],
    };
    let search = validate(program, options).unwrap();
    let (queries, subjects) = read_inputs(&search, query, subject).unwrap();
    let outcome = run_collect(&search, Inputs::new(&queries, &subjects), &["6"]).unwrap();
    let report = &outcome.reports[0];
    let (cli, _) = cli(&case, &query_path, &subject_path, "6");
    assert!(
        report.bytes == cli,
        "{context}: the report is the command line's"
    );
    let mut pairs = std::collections::BTreeSet::new();
    for hsp in &outcome.hsps {
        let row = &report.bytes[report.spans[hsp.index].clone().unwrap()];
        let columns: Vec<&[u8]> = row.split(|&byte| byte == b'\t').collect();
        assert_eq!(hsp.query_id, columns[0], "{context}");
        assert_eq!(hsp.subject_id, columns[1], "{context}");
        assert_eq!(hsp.query_id, query_ids[hsp.q_idx], "{context}");
        assert_eq!(hsp.subject_id, subject_ids[hsp.s_idx], "{context}");
        pairs.insert((hsp.q_idx, hsp.s_idx));
    }
    let _ = std::fs::remove_file(&query_path);
    let _ = std::fs::remove_file(&subject_path);
    pairs
}

/// FASTA records of `ids`, each with the title `title` and the letters that `program` reads for
/// `role`.
fn records_of(program: Program, role: Role, ids: &[&[u8]], title: &[u8]) -> Vec<u8> {
    let mut fasta = Vec::new();
    for id in ids {
        fasta.extend_from_slice(b">");
        fasta.extend_from_slice(id);
        fasta.extend_from_slice(title);
        fasta.extend_from_slice(b"\n");
        fasta.extend_from_slice(letters_of(program, role).as_bytes());
        fasta.extend_from_slice(b"\n");
    }
    fasta
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:214-217
// The ids are bytes (`std::string`) that the report writes as they are.
/// Titles whose first word holds bytes that are not UTF-8 (Latin-1's `é`, 0xE9; a continuation
/// byte alone, 0x80; 0xFF) or UTF-8's own letters of two and three bytes (`é`, `€`): in every
/// program, each record's id is the bytes of its column in the `-outfmt 6` row of the same run
/// (`Report::spans`), which is the command line's row.
#[test]
fn the_ids_are_the_bytes_of_the_tabular_report() {
    let query_ids: [&[u8]; 2] = [b"q\xE9", "q\u{e9}".as_bytes()];
    let subject_ids: [&[u8]; 2] = [b"s\x80\xFF", "s\u{e9}\u{20ac}".as_bytes()];
    for program in Program::ALL {
        let query = records_of(program, Role::Query, &query_ids, b" a title");
        let subject = records_of(program, Role::Subject, &subject_ids, b" a title");
        let pairs = ids_of_the_rows(program, &[], (&query, &subject), (&query_ids, &subject_ids));
        assert_eq!(pairs.len(), 4, "{program}: every query hits every subject");
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:77-86
/// With `-parse_deflines`, an id with a byte that is not ASCII (`lcl|q\xE9`, `gb|AB\xE9001.1|`)
/// stops before it, and the rest of the word goes to the title: NCBI BLAST+ 2.17.0 writes `q` and
/// `AB` in `qaccver` for these deflines (blastp, probed on 2026-11-20), and so does LOSAT in every
/// program, in the records too.
#[test]
fn a_parsed_id_stops_before_a_byte_that_is_not_ascii() {
    let deflines: [&[u8]; 2] = [b"lcl|q\xE9", b"gb|AB\xE9001.1|"];
    let (query_ids, subject_ids): ([&[u8]; 2], [&[u8]; 1]) = ([b"q", b"AB"], [b"s1"]);
    for program in Program::ALL {
        let query = records_of(program, Role::Query, &deflines, b" t");
        let subject = records_of(program, Role::Subject, &[b"s1"], b" t");
        let pairs = ids_of_the_rows(
            program,
            &["-parse_deflines"],
            (&query, &subject),
            (&query_ids, &subject_ids),
        );
        assert_eq!(pairs.len(), 2, "{program}: both queries hit the subject");
    }
}

#[test]
fn the_value_kinds_read_the_defaults_of_the_options() {
    let kind_of = |program: Program, flag: &str| {
        describe(program)
            .parameters
            .iter()
            .find(|parameter| parameter.flag == flag)
            .unwrap_or_else(|| panic!("{program} {flag}"))
            .value_kind
    };
    assert_eq!(kind_of(Program::Blastn, "-evalue"), ValueKind::Real);
    assert_eq!(kind_of(Program::Blastp, "-num_threads"), ValueKind::Integer);
    assert_eq!(kind_of(Program::Blastn, "-ungapped"), ValueKind::Flag);
    assert_eq!(kind_of(Program::Blastp, "-task"), ValueKind::Text);
    assert_eq!(
        kind_of(Program::Tblastn, "-soft_masking"),
        ValueKind::Boolean
    );
    for program in Program::ALL {
        for parameter in describe(program).parameters {
            let flag = format!("{program} {}", parameter.flag);
            assert_eq!(
                parameter.takes_value,
                parameter.value_kind != ValueKind::Flag,
                "{flag}"
            );
            let Some(default) = parameter.default.as_deref() else {
                continue;
            };
            match parameter.value_kind {
                ValueKind::Integer => assert!(default.parse::<i64>().is_ok(), "{flag}: {default}"),
                ValueKind::Real => assert!(default.parse::<f64>().is_ok(), "{flag}: {default}"),
                ValueKind::Boolean => assert!(
                    matches!(default, "true" | "false" | "t" | "f" | "T" | "F"),
                    "{flag}: {default}"
                ),
                _ => {}
            }
        }
    }
}

#[test]
fn every_structured_format_is_the_command_lines() {
    // The formats that hold the whole run in one document, from the search of each program.
    for case in [&CASES[0], &CASES[2], &CASES[3], &CASES[4], &CASES[5]] {
        const FORMATS: &[&str] = &["5", "8", "9", "10", "11", "12", "15", "16", "20"];
        let (query, subject) = (fixture(case.query), fixture(case.subject));
        let outcome = api(
            &Case {
                program: case.program,
                query: case.query,
                subject: case.subject,
                options: case.options,
                formats: FORMATS,
            },
            &query,
            &subject,
            &[],
        );
        for (report, outfmt) in outcome.reports.iter().zip(FORMATS) {
            let (bytes, _) = cli(case, &query, &subject, outfmt);
            assert!(
                report.bytes == bytes,
                "{} -outfmt {outfmt}: the report differs from the command line's",
                case.program
            );
        }
    }
}

/// The file of `outfmt` and the side files next to it, by name, removed.
fn take_files(directory: &Path) -> Vec<(String, Vec<u8>)> {
    let mut files: Vec<(String, Vec<u8>)> = std::fs::read_dir(directory)
        .unwrap()
        .map(|entry| {
            let path = entry.unwrap().path();
            let bytes = std::fs::read(&path).unwrap();
            std::fs::remove_file(&path).unwrap();
            (
                path.file_name().unwrap().to_string_lossy().into_owned(),
                bytes,
            )
        })
        .collect();
    files.sort();
    files
}

#[test]
fn a_report_in_a_file_has_its_side_files_as_on_the_command_line() {
    let directory = temp_path("files");
    std::fs::create_dir_all(&directory).unwrap();
    // Both write to one name: the master file names its side files by it.
    let out = directory.join("report.out");
    for case in [&CASES[0], &CASES[2], &CASES[3], &CASES[4], &CASES[5]] {
        for outfmt in ["13", "14", "6"] {
            let context = format!("{} -outfmt {outfmt}", case.program);
            let (query, subject) = (fixture(case.query), fixture(case.subject));
            let queries =
                read_fasta(case.program, Role::Query, &std::fs::read(&query).unwrap()).unwrap();
            let subjects = read_fasta(
                case.program,
                Role::Subject,
                &std::fs::read(&subject).unwrap(),
            )
            .unwrap();
            let search = validate(case.program, case.options).unwrap();
            let inputs = Inputs::new(&queries, &subjects)
                .named(query.to_str().unwrap(), subject.to_str().unwrap());
            let mut diagnostics = Vec::new();
            run_to_file(&search, inputs, outfmt, &out, &mut diagnostics).unwrap();
            let written = take_files(&directory);
            let status = command(env!("CARGO_BIN_EXE_losat"))
                .arg(case.program.name())
                .arg("-query")
                .arg(&query)
                .arg("-subject")
                .arg(&subject)
                .args(case.options)
                .args(["-outfmt", outfmt, "-out"])
                .arg(&out)
                .output()
                .unwrap();
            assert!(status.status.success(), "{context}");
            let frozen = take_files(&directory);
            assert_eq!(
                written.iter().map(|file| &file.0).collect::<Vec<_>>(),
                frozen.iter().map(|file| &file.0).collect::<Vec<_>>(),
                "{context}"
            );
            assert!(written == frozen, "{context}: the bytes");
            assert_eq!(diagnostics, status.stderr, "{context}: the warnings");
            assert_eq!(
                written.len() > 1,
                outfmt != "6",
                "{context}: the side files"
            );
        }
    }
    // Without a file there is no place for the side files: the command line's message.
    let (query, subject) = (fixture(CASES[0].query), fixture(CASES[0].subject));
    let error = api_error(&CASES[0], &query, &subject, "14");
    assert!(error.to_string().contains("provide a file name"), "{error}");
    // `-` is the command line's standard output, which the library does not write.
    let search = validate(Program::Blastp, &[]).unwrap();
    let error = run_to_file(
        &search,
        Inputs::new(&[], &[]),
        "6",
        Path::new("-"),
        &mut std::io::sink(),
    )
    .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    std::fs::remove_dir_all(&directory).unwrap();
}

// Series Q item Q6a: psiblast's PSSM files are written where the command line writes them, with
// its names: `-save_each_pssm` numbers them over the run, also over `run_batches`' batches of one
// query each (NCBI's `CAutoOutputFileReset`, made once per run, blast_args.cpp:1304-1318). Two
// queries (each makes a PSSM), with rounds and with the PSSM after the last round.
#[test]
fn the_pssm_files_are_the_command_lines() {
    let directory = temp_path("pssm");
    std::fs::create_dir_all(&directory).unwrap();
    let (query, subject) = (fixture("psiblast/q4c_tk_pol.faa"), fixture("PajaWSV.faa"));
    let queries = read_fasta(
        Program::Psiblast,
        Role::Query,
        &std::fs::read(&query).unwrap(),
    )
    .unwrap();
    let subjects = read_fasta(
        Program::Psiblast,
        Role::Subject,
        &std::fs::read(&subject).unwrap(),
    )
    .unwrap();
    assert_eq!(queries.len(), 2);
    let inputs =
        Inputs::new(&queries, &subjects).named(query.to_str().unwrap(), subject.to_str().unwrap());
    let (asn, ascii) = (directory.join("p.asn"), directory.join("p.txt"));
    let files = [
        "-out_pssm",
        asn.to_str().unwrap(),
        "-out_ascii_pssm",
        ascii.to_str().unwrap(),
    ];
    for options in [
        &["-num_iterations", "2", "-save_each_pssm"][..],
        &["-num_iterations", "2"],
        &["-save_pssm_after_last_round", "-save_each_pssm"],
        &[
            "-num_iterations",
            "3",
            "-save_pssm_after_last_round",
            "-save_each_pssm",
        ],
    ] {
        let options: Vec<&str> = options.iter().chain(&files).copied().collect();
        let context = options[..options.len() - 4].join(" ");
        let out = directory.join("report.out");
        let status = command(env!("CARGO_BIN_EXE_losat"))
            .arg("psiblast")
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .args(&options)
            .args(["-outfmt", "6", "-out"])
            .arg(&out)
            .output()
            .unwrap();
        assert!(status.status.success(), "{context}");
        let report = std::fs::read(&out).unwrap();
        std::fs::remove_file(&out).unwrap();
        let frozen = take_files(&directory);
        let names: Vec<&str> = frozen.iter().map(|file| file.0.as_str()).collect();
        if options.contains(&"-save_each_pssm") {
            assert!(
                names.contains(&"p.asn.2") && names.contains(&"p.txt.2"),
                "{context}: {names:?}"
            );
        } else {
            assert_eq!(names, ["p.asn", "p.txt"], "{context}");
        }
        let search = validate(Program::Psiblast, &options).unwrap();
        let collected = run_collect(&search, inputs, &["6"]).unwrap();
        assert_eq!(collected.reports[0].bytes, report, "{context}: run_collect");
        assert!(
            take_files(&directory) == frozen,
            "{context}: run_collect's files"
        );
        let mut written = Vec::new();
        run(
            &search,
            inputs,
            &mut [("6", &mut written)],
            &mut std::io::sink(),
        )
        .unwrap();
        assert_eq!(written, report, "{context}: run");
        assert!(take_files(&directory) == frozen, "{context}: run's files");
        run_to_file(&search, inputs, "6", &out, &mut std::io::sink()).unwrap();
        assert_eq!(
            std::fs::read(&out).unwrap(),
            report,
            "{context}: run_to_file"
        );
        std::fs::remove_file(&out).unwrap();
        assert!(
            take_files(&directory) == frozen,
            "{context}: run_to_file's files"
        );
        let batches: Vec<_> = run_batches(&search, inputs)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(batches.len(), 2, "{context}: one query per batch");
        assert_eq!(
            batches.iter().map(|batch| batch.hsps.len()).sum::<usize>(),
            collected.hsps.len(),
            "{context}"
        );
        assert!(
            take_files(&directory) == frozen,
            "{context}: the batches' files"
        );
    }
    // A call whose formats take two searches (`-num_descriptions` with a report and a table, or
    // psiblast's outfmt 5 or 11 with another format: their hit lists are cut apart) would write the
    // files of two searches: refused before any search.
    let search = validate(
        Program::Psiblast,
        &[
            &["-num_iterations", "2", "-num_descriptions", "5"][..],
            &files,
        ]
        .concat(),
    )
    .unwrap();
    let refused = |error: &Error| {
        matches!(error, Error::Argv(_))
            && error.to_string().starts_with(
                "-out_pssm is not accepted by the library API with formats that take two searches",
            )
    };
    let error = run_collect(&search, inputs, &["0", "6"]).unwrap_err();
    assert!(refused(&error), "{error}");
    assert!(take_files(&directory).is_empty());
    let other_cuts = validate(
        Program::Psiblast,
        &[&["-num_iterations", "2"][..], &files].concat(),
    )
    .unwrap();
    for formats in [["11", "6"], ["5", "6"]] {
        let error = run_collect(&other_cuts, inputs, &formats).unwrap_err();
        assert!(refused(&error), "{formats:?}: {error}");
        assert!(take_files(&directory).is_empty(), "{formats:?}");
    }
    // One format, or formats of one hit list size, search once.
    run_collect(&search, inputs, &["6", "7"]).unwrap();
    assert!(!take_files(&directory).is_empty());
    // A query title that the ASN.1 PSSM cannot write as it is (here the third query's `é`, two
    // bytes outside the visible ASCII range): NCBI writes each such byte as `#` and posts error
    // 803.7 for it as it writes the file, once per file written, and ends with status 0
    // (objostrasn.cpp:413-444; psiblast_app.cpp:385-407 `SavePssmToFile`; NCBI BLAST+ 2.17.0's
    // standard error for this command, line NP M4a-1b2). Nothing is refused: the library's calls,
    // and its batches of one query each, write the command line's files and give its posts where
    // the library gives its diagnostics.
    let mut records = std::fs::read(&query).unwrap();
    records.extend_from_slice(b">third query \xc3\xa9\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n");
    let titled = temp_path("titled.faa");
    std::fs::write(&titled, &records).unwrap();
    let queries = read_fasta(Program::Psiblast, Role::Query, &records).unwrap();
    assert_eq!(queries.len(), 3);
    let inputs =
        Inputs::new(&queries, &subjects).named(titled.to_str().unwrap(), subject.to_str().unwrap());
    let options = [&["-save_pssm_after_last_round"][..], &files].concat();
    let written = command(env!("CARGO_BIN_EXE_losat"))
        .arg("psiblast")
        .arg("-query")
        .arg(&titled)
        .arg("-subject")
        .arg(&subject)
        .args(&options)
        .args(["-outfmt", "6"])
        .output()
        .unwrap();
    std::fs::remove_file(&titled).unwrap();
    assert!(written.status.success(), "{written:?}");
    let posts: Vec<u8> = [0xc3u8, 0xa9]
        .iter()
        .flat_map(|byte| {
            let mut post = format!(
                "Error: (803.7) [psiblast] PssmWithParameters.pssm.query.seq.descr.E.title\nBad char [0x{byte:02X}] in string at line 12\n"
            )
            .into_bytes();
            post.extend_from_slice(b"third query \xc3\xa9\n");
            post
        })
        .collect();
    assert_eq!(
        String::from_utf8_lossy(&written.stderr),
        String::from_utf8_lossy(&posts),
        "the command line's posts"
    );
    // Without `-save_each_pssm` each query's PSSM replaces the file of the query before: the
    // files left are the third query's.
    let frozen = take_files(&directory);
    let names: Vec<&str> = frozen.iter().map(|file| file.0.as_str()).collect();
    assert_eq!(names, ["p.asn", "p.txt"]);
    let asn = String::from_utf8(frozen[0].1.clone()).unwrap();
    assert!(asn.contains("        title \"third query ##\""), "{asn}");
    let search = validate(Program::Psiblast, &options).unwrap();
    let collected = run_collect(&search, inputs, &["6"]).unwrap();
    assert_eq!(collected.reports[0].bytes, written.stdout, "run_collect");
    assert_eq!(collected.diagnostics, posts, "run_collect's posts");
    assert!(take_files(&directory) == frozen, "run_collect's files");
    let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
    run(&search, inputs, &mut [("6", &mut report)], &mut diagnostics).unwrap();
    assert_eq!(report, written.stdout, "run");
    assert_eq!(diagnostics, posts, "run's posts");
    assert!(take_files(&directory) == frozen, "run's files");
    let out = temp_path("titled.out");
    let mut diagnostics = Vec::new();
    run_to_file(&search, inputs, "6", &out, &mut diagnostics).unwrap();
    assert_eq!(std::fs::read(&out).unwrap(), written.stdout, "run_to_file");
    std::fs::remove_file(&out).unwrap();
    assert_eq!(diagnostics, posts, "run_to_file's posts");
    assert!(take_files(&directory) == frozen, "run_to_file's files");
    // The posts come with the batch of the query whose file NCBI writes then.
    let batches: Vec<_> = run_batches(&search, inputs)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(batches.len(), 3, "one query per batch");
    let batch_posts: Vec<&[u8]> = batches
        .iter()
        .map(|batch| batch.diagnostics.as_slice())
        .collect();
    let none: &[u8] = &[];
    assert_eq!(
        batch_posts,
        [none, none, posts.as_slice()],
        "the batches' posts"
    );
    assert_eq!(
        batches.iter().map(|batch| batch.hsps.len()).sum::<usize>(),
        collected.hsps.len()
    );
    assert!(take_files(&directory) == frozen, "the batches' files");
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn the_writers_of_a_run_receive_the_reports_of_the_collected_run() {
    for case in [&CASES[0], &CASES[1], &CASES[WARNINGS_CASE]] {
        let context = format!("{} {:?}", case.program, case.options);
        let (query, subject) = (fixture(case.query), fixture(case.subject));
        let collected = api(case, &query, &subject, &[]);
        let query_bytes = std::fs::read(&query).unwrap();
        let subject_bytes = std::fs::read(&subject).unwrap();
        let queries = read_fasta(case.program, Role::Query, &query_bytes).unwrap();
        let subjects = read_fasta(case.program, Role::Subject, &subject_bytes).unwrap();
        let search = validate(case.program, case.options).unwrap();
        let inputs = Inputs::new(&queries, &subjects)
            .named(query.to_str().unwrap(), subject.to_str().unwrap());
        let mut writers: Vec<Vec<u8>> = vec![Vec::new(); case.formats.len()];
        let mut diagnostics = Vec::new();
        {
            let mut reports: Vec<(&str, &mut (dyn std::io::Write + Send))> = case
                .formats
                .iter()
                .zip(writers.iter_mut())
                .map(|(outfmt, writer)| (*outfmt, writer as &mut (dyn std::io::Write + Send)))
                .collect();
            run(&search, inputs, &mut reports, &mut diagnostics)
                .unwrap_or_else(|error| panic!("{context}: {error}"));
        }
        for (written, report) in writers.iter().zip(&collected.reports) {
            assert!(written == &report.bytes, "{context}: {}", report.outfmt);
        }
        assert_eq!(diagnostics, collected.diagnostics, "{context}: diagnostics");
    }
}

#[test]
fn two_searches_at_once_give_the_reports_of_each_alone() {
    let cases = [&CASES[0], &CASES[3]];
    let alone: Vec<Outcome> = cases
        .iter()
        .map(|case| {
            api(
                case,
                &fixture(case.query),
                &fixture(case.subject),
                &["-num_threads", "2"],
            )
        })
        .collect();
    let together: Vec<Outcome> = std::thread::scope(|scope| {
        let handles: Vec<_> = cases
            .iter()
            .map(|case| {
                scope.spawn(move || {
                    api(
                        case,
                        &fixture(case.query),
                        &fixture(case.subject),
                        &["-num_threads", "2"],
                    )
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("a search thread"))
            .collect()
    });
    for ((case, alone), together) in cases.iter().zip(&alone).zip(&together) {
        assert_eq!(alone, together, "{}: at once versus alone", case.program);
        assert!(
            !alone.hsps.is_empty(),
            "{}: the case has hits",
            case.program
        );
    }
}

/// One run of `run_batches`, split in two: the start (`run_batches`), then its batches.
struct BatchRun {
    search: losat::api::Search,
    queries: Vec<FastaRecord>,
    subjects: Vec<FastaRecord>,
    pssm: Option<(PathBuf, String)>,
}

impl BatchRun {
    /// psiblast `-in_pssm` of the PSSM whose query's annotations make NCBI's object manager fail
    /// its main thread's diagnostic stream (`Failed to parse location of`, 1204.6), after which
    /// its later posts print their headers alone (the fixture `db.ps.inpssm.d3.om.k01_failed_then_warn`).
    fn muting(directory: &Path, threads: &str) -> BatchRun {
        let case = PssmCase {
            program: Program::Psiblast,
            pssm: "fixtures/pssm/d3om_k01_failed_then_alwarn.asn.gz",
            subject: "db_search/paja",
            database: true,
            options: &[],
        };
        BatchRun {
            search: validate(
                Program::Psiblast,
                &["-comp_based_stats", "0", "-num_threads", threads],
            )
            .unwrap(),
            queries: Vec::new(),
            subjects: Vec::new(),
            pssm: Some((case.file(directory), engine_tests(case.subject).to_string())),
        }
    }

    /// blastp of a query that NCBI warns about (`Warning: [blastp] Query_1 p1: Could not
    /// calculate ungapped Karlin-Altschul parameters ...`), a post of the main thread.
    fn warning(threads: &str) -> BatchRun {
        let subjects = read_fasta(
            Program::Blastp,
            Role::Subject,
            &std::fs::read(fixture("outfmt0/e2e_many_subject.faa")).unwrap(),
        )
        .unwrap();
        BatchRun {
            search: validate(Program::Blastp, &["-num_threads", threads]).unwrap(),
            queries: read_fasta(
                Program::Blastp,
                Role::Query,
                b">p1\nXXXXXXXXXXXXXXXXXXXXXXXXX\n",
            )
            .unwrap(),
            subjects,
            pssm: None,
        }
    }

    fn inputs(&self) -> Inputs<'_> {
        match &self.pssm {
            Some((pssm, database)) => Inputs::pssm_database(pssm, database),
            None => Inputs::new(&self.queries, &self.subjects),
        }
    }

    /// The diagnostics of every batch.
    fn diagnostics(batches: losat::api::Batches<'_>) -> Vec<u8> {
        batches
            .flat_map(|batch| batch.unwrap().diagnostics)
            .collect()
    }

    fn alone(&self) -> Vec<u8> {
        BatchRun::diagnostics(run_batches(&self.search, self.inputs()).unwrap())
    }
}

// Line RG D3-1 review 4 B1: NCBI's object manager fails the main thread's diagnostic stream of
// one run (one application); in the library each call is a run, and a call's failed stream mutes
// no other call's posts, nor does another call's start unmute them: two calls whose batches one
// thread takes in turn, and two calls at once on two threads, write the diagnostics of each alone.
#[test]
fn a_failed_diagnostic_stream_is_its_own_runs() {
    let directory = temp_path("muted_streams");
    std::fs::create_dir_all(&directory).unwrap();
    let muting = BatchRun::muting(&directory, "1");
    let warning = BatchRun::warning("1");
    let (muted_alone, warning_alone) = (muting.alone(), warning.alone());
    assert!(
        muted_alone.starts_with(b"Error: (1204.6) [psiblast] Failed to parse location of \n"),
        "{}",
        String::from_utf8_lossy(&muted_alone)
    );
    assert!(
        warning_alone.starts_with(b"Warning: [blastp] Query_1 p1: Could not calculate"),
        "{}",
        String::from_utf8_lossy(&warning_alone)
    );
    // One thread: both calls start, then the muting call's batch, then the other's; and the other
    // way round.
    for muting_first in [true, false] {
        let mut muted = run_batches(&muting.search, muting.inputs()).unwrap();
        let mut warned = run_batches(&warning.search, warning.inputs()).unwrap();
        let (mut muted_text, mut warned_text) = (Vec::new(), Vec::new());
        if muting_first {
            muted_text.extend(muted.next().unwrap().unwrap().diagnostics);
            warned_text.extend(BatchRun::diagnostics(warned));
            muted_text.extend(BatchRun::diagnostics(muted));
        } else {
            warned_text.extend(warned.next().unwrap().unwrap().diagnostics);
            muted_text.extend(BatchRun::diagnostics(muted));
            warned_text.extend(BatchRun::diagnostics(warned));
        }
        assert_eq!(muted_text, muted_alone, "muting first {muting_first}");
        assert_eq!(warned_text, warning_alone, "muting first {muting_first}");
    }
    // Two threads: both calls start, then each searches on its thread (the muting call first),
    // with pools of two threads.
    let (muting, warning) = (BatchRun::muting(&directory, "2"), BatchRun::warning("2"));
    let started = std::sync::Barrier::new(2);
    let muted_done = std::sync::Barrier::new(2);
    let (muted_text, warned_text) = std::thread::scope(|scope| {
        let muted = scope.spawn(|| {
            let batches = run_batches(&muting.search, muting.inputs()).unwrap();
            started.wait();
            let text = BatchRun::diagnostics(batches);
            muted_done.wait();
            text
        });
        let warned = scope.spawn(|| {
            let batches = run_batches(&warning.search, warning.inputs()).unwrap();
            started.wait();
            muted_done.wait();
            BatchRun::diagnostics(batches)
        });
        (muted.join().unwrap(), warned.join().unwrap())
    });
    assert_eq!(muted_text, muting.alone(), "two threads");
    assert_eq!(warned_text, warning.alone(), "two threads");
    assert_eq!(muted_text, muted_alone, "two threads: as with one");
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn the_example_prints_the_command_lines_report_and_writes_only_its_warnings_to_stderr() {
    let example = example_binary("library_search");
    for (case, warnings) in [(&CASES[0], false), (&CASES[WARNINGS_CASE], true)] {
        let (query, subject) = (fixture(case.query), fixture(case.subject));
        let (bytes, stderr) = cli(case, &query, &subject, "6");
        let output = command(example.to_str().unwrap())
            .arg(case.program.name())
            .arg(&query)
            .arg(&subject)
            .args(case.options)
            .output()
            .expect("run the example");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stdout == bytes,
            "{}: the example's report differs",
            case.program
        );
        assert_eq!(output.stderr, stderr, "{}: standard error", case.program);
        assert_eq!(
            warnings,
            !stderr.is_empty(),
            "{}: the case's warnings",
            case.program
        );
    }
}

/// One database search: a program, a query and a database under `LOSAT/tests`, and options. The
/// whole set of database searches through the API is `tests/db_search_fixtures.py check --library`.
struct DatabaseCase {
    program: Program,
    query: &'static str,
    db: &'static str,
    options: &'static [&'static str],
}

const DATABASE_CASES: &[DatabaseCase] = &[
    DatabaseCase {
        program: Program::Blastn,
        query: "seqdb/nucl.fna",
        db: "seqdb/db/nucl_v5",
        options: &["-task", "blastn"],
    },
    DatabaseCase {
        program: Program::Blastp,
        query: "seqdb/prot.faa",
        db: "seqdb/db/prot_v5",
        options: &[],
    },
    // A sequence with parsed Seq-ids shows the ids of its definition lines.
    DatabaseCase {
        program: Program::Blastp,
        query: "db_search/ids/ids_q.faa",
        db: "db_search/ids/prot_v5",
        options: &[],
    },
    DatabaseCase {
        program: Program::Blastx,
        query: "db_search/bslice_q.fna",
        db: "db_search/bslice",
        options: &[],
    },
    DatabaseCase {
        program: Program::Blastx,
        query: "db_search/ids/ids_q.fna",
        db: "db_search/ids/prot_v5",
        options: &[],
    },
    DatabaseCase {
        program: Program::Tblastn,
        query: "db_search/tslice_q.faa",
        db: "db_search/tslice",
        options: &[],
    },
    DatabaseCase {
        program: Program::Tblastx,
        query: "db_search/xslice_q.fna",
        db: "db_search/xslice",
        options: &[],
    },
    DatabaseCase {
        program: Program::Psiblast,
        query: "fasta/psiblast/tk.faa",
        db: "db_search/paja",
        options: &[],
    },
    // psiblast's rounds: 9 subjects in round 1, 5 in round 2, converged (the fixture
    // `db.ps.it.it3.tk.7`); each round's hit list ranks its own subjects.
    DatabaseCase {
        program: Program::Psiblast,
        query: "fasta/psiblast/tk.faa",
        db: "db_search/paja",
        options: &["-num_iterations", "3"],
    },
    // psiblast's rounds of the second query record: `-query_loc` starts past the end of the first
    // (10 letters), which NCBI's reader skips (the fixture `db.ps.it.qloc.skip.6`).
    DatabaseCase {
        program: Program::Psiblast,
        query: "fasta/psiblast/q4a_skip.faa",
        db: "db_search/paja",
        options: &["-num_iterations", "2", "-query_loc", "20-300"],
    },
];

fn engine_tests(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join(name)
        .to_str()
        .unwrap()
        .to_string()
}

impl DatabaseCase {
    /// The query file's records and the search, as `-db` runs them.
    fn read(&self) -> (Vec<u8>, Vec<FastaRecord>, losat::api::Search) {
        let bytes = std::fs::read(engine_tests(self.query)).unwrap();
        let search = validate(self.program, self.options).unwrap();
        let (records, none) = read_inputs(&search, &bytes, b"").unwrap();
        assert!(none.is_empty());
        (bytes, records, search)
    }

    /// The command line's report for one `-outfmt`.
    fn cli(&self, outfmt: &str) -> Vec<u8> {
        let output = command(env!("CARGO_BIN_EXE_losat"))
            .arg(self.program.name())
            .arg("-query")
            .arg(engine_tests(self.query))
            .arg("-db")
            .arg(engine_tests(self.db))
            .args(self.options)
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .expect("run the command line");
        assert!(
            output.status.success(),
            "{} {} -outfmt {outfmt}: {}",
            self.program,
            self.db,
            String::from_utf8_lossy(&output.stderr)
        );
        output.stdout
    }
}

#[test]
fn a_database_search_is_the_command_lines() {
    for case in DATABASE_CASES {
        let context = format!("{} {}", case.program, case.db);
        let (_, records, search) = case.read();
        let db = engine_tests(case.db);
        let query_name = engine_tests(case.query);
        let inputs = Inputs::database(&records, &db).named(&query_name, "ignored");
        // psiblast writes the first three and the columns (`PsiblastArgs::check_outfmt`).
        let formats: &[&str] = match case.program {
            Program::Psiblast => &["6", "0", "7", COLUMNS],
            _ => &["6", "0", "7", COLUMNS, "5", "11"],
        };
        let outcome = run_collect(&search, inputs, formats).unwrap();
        for (report, outfmt) in outcome.reports.iter().zip(formats) {
            assert!(
                report.bytes == case.cli(outfmt),
                "{context} -outfmt {outfmt}"
            );
        }
        assert!(!outcome.hsps.is_empty(), "{context}");
        // The typed values are the columns, each row's span in every report that has rows.
        let table = &outcome.reports[3].bytes;
        let rows = rows(table);
        assert_eq!(rows.len(), outcome.hsps.len(), "{context}");
        for (hsp, row) in outcome.hsps.iter().zip(rows) {
            check_columns(&context, hsp, &table[row]);
            assert!(outcome.reports[0].spans[hsp.index].is_some(), "{context}");
            assert!(outcome.reports[1].spans[hsp.index].is_some(), "{context}");
            assert!(
                outcome.reports[1].subject_headings[hsp.index].is_some(),
                "{context}"
            );
        }
        // One sequence has one ordinal id, whichever query hit it.
        let mut ids = std::collections::HashMap::new();
        for hsp in &outcome.hsps {
            let id = ids
                .entry(hsp.s_idx)
                .or_insert_with(|| hsp.subject_id.clone());
            assert_eq!(*id, hsp.subject_id, "{context}");
        }
        // The statistics are those of the XML report.
        if let Some(xml) = outcome.reports.get(4) {
            let xml = String::from_utf8_lossy(&xml.bytes);
            let searched: Vec<usize> = (0..records.len()).collect();
            check_statistics(&context, &outcome, &xml, &searched);
        }
    }
}

/// A PSSM query (series Q item Q6b): `Inputs::pssm` and `Inputs::pssm_database`, the command
/// line's `-in_pssm` of psiblast and tblastn (PSI-TBLASTN), against records and a database. The
/// whole set of PSSM queries through the API is `tests/db_search_fixtures.py check --library`.
struct PssmCase {
    program: Program,
    /// Under `tests/`; `.gz` is inflated into a file of its own first (the harness does the same).
    pssm: &'static str,
    /// Under `tests/`: a FASTA subject set (`-subject`), or a database (`-db`).
    subject: &'static str,
    database: bool,
    options: &'static [&'static str],
}

const PSSM_CASES: &[PssmCase] = &[
    PssmCase {
        program: Program::Psiblast,
        pssm: "fixtures/pssm/ps4_tk_paja.asn.gz",
        subject: "fasta/PajaWSV.faa",
        database: false,
        options: &[],
    },
    PssmCase {
        program: Program::Psiblast,
        pssm: "fixtures/pssm/ps4_tk_paja.asn.gz",
        subject: "db_search/paja",
        database: true,
        options: &[],
    },
    // The rounds of a PSSM query (`-in_pssm` with `-num_iterations`, item PS5b of line PS).
    PssmCase {
        program: Program::Psiblast,
        pssm: "fixtures/pssm/ps4_tk_paja.asn.gz",
        subject: "db_search/paja",
        database: true,
        options: &["-num_iterations", "3"],
    },
    PssmCase {
        program: Program::Tblastn,
        pssm: "db_search/of1/tnq1.pssm",
        subject: "db_search/of1/tns.fa",
        database: false,
        options: &[],
    },
    PssmCase {
        program: Program::Tblastn,
        pssm: "fixtures/pssm/ps7_ts1_tslice.asn.gz",
        subject: "db_search/tslice",
        database: true,
        options: &[],
    },
];

impl PssmCase {
    /// The PSSM file as the command line reads it: the fixture, or its inflated copy in `directory`.
    fn file(&self, directory: &Path) -> PathBuf {
        let path = PathBuf::from(engine_tests(self.pssm));
        if !self.pssm.ends_with(".gz") {
            return path;
        }
        let mut bytes = Vec::new();
        std::io::Read::read_to_end(
            &mut flate2::read::GzDecoder::new(std::fs::read(&path).unwrap().as_slice()),
            &mut bytes,
        )
        .unwrap();
        let inflated = directory.join(path.file_stem().unwrap());
        std::fs::write(&inflated, bytes).unwrap();
        inflated
    }

    /// The command line's report for one `-outfmt`, and its standard error.
    fn cli(&self, pssm: &Path, extra: &[&str], outfmt: &str) -> (Vec<u8>, Vec<u8>) {
        let output = command(env!("CARGO_BIN_EXE_losat"))
            .arg(self.program.name())
            .arg("-in_pssm")
            .arg(pssm)
            .arg(if self.database { "-db" } else { "-subject" })
            .arg(engine_tests(self.subject))
            .args(self.options)
            .args(extra)
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .expect("run the command line");
        assert!(
            output.status.success(),
            "{} {} -outfmt {outfmt}: {}",
            self.program,
            self.pssm,
            String::from_utf8_lossy(&output.stderr)
        );
        (output.stdout, output.stderr)
    }
}

#[test]
fn a_pssm_query_is_the_command_lines() {
    let directory = temp_path("pssm_query");
    std::fs::create_dir_all(&directory).unwrap();
    for case in PSSM_CASES {
        let context = format!(
            "{} {} {} {:?}",
            case.program, case.pssm, case.subject, case.options
        );
        let pssm = case.file(&directory);
        let search = validate(case.program, case.options).unwrap();
        let subject = engine_tests(case.subject);
        let subjects = match case.database {
            true => Vec::new(),
            false => read_fasta(
                case.program,
                Role::Subject,
                &std::fs::read(&subject).unwrap(),
            )
            .unwrap(),
        };
        let inputs = match case.database {
            true => Inputs::pssm_database(&pssm, &subject),
            false => Inputs::pssm(&pssm, &subjects).named("ignored", &subject),
        };
        let formats = ["6", "0", "7", COLUMNS];
        let outcome = run_collect(&search, inputs, &formats)
            .unwrap_or_else(|error| panic!("{context}: {error}"));
        for (report, outfmt) in outcome.reports.iter().zip(formats) {
            let (stdout, stderr) = case.cli(&pssm, &[], outfmt);
            assert!(report.bytes == stdout, "{context} -outfmt {outfmt}");
            if outfmt == "6" {
                assert_eq!(
                    outcome.diagnostics, stderr,
                    "{context}: the diagnostics are the command line's"
                );
            }
        }
        // One query, the PSSM's: query 0, with its statistics and rounds.
        assert!(!outcome.hsps.is_empty(), "{context}");
        assert!(outcome.hsps.iter().all(|hsp| hsp.q_idx == 0), "{context}");
        assert_eq!(outcome.query_stats.len(), 1, "{context}");
        assert!(outcome.query_stats[0].is_some(), "{context}");
        assert_eq!(outcome.rounds.len(), 1, "{context}");
        let rounds = outcome.hsps.iter().map(|hsp| hsp.round).max().unwrap();
        assert_eq!(outcome.rounds[0].stats.len() as u32, rounds, "{context}");
        if case.options.contains(&"-num_iterations") {
            assert!(rounds > 1, "{context}");
        }
        let table = &outcome.reports[3].bytes;
        let rows = rows(table);
        assert_eq!(rows.len(), outcome.hsps.len(), "{context}");
        for (hsp, row) in outcome.hsps.iter().zip(rows) {
            check_columns(&context, hsp, &table[row]);
            assert!(outcome.reports[0].spans[hsp.index].is_some(), "{context}");
        }
        // One batch, the same records: NCBI reads no query batches for a PSSM query.
        let batches = run_batches(&search, inputs).unwrap();
        assert_eq!(batches.query_count(), 1, "{context}");
        let batches: Vec<_> = batches.map(Result::unwrap).collect();
        assert_eq!(batches.len(), 1, "{context}");
        assert_eq!(batches[0].query_indices, 0..1, "{context}");
        assert_eq!(batches[0].hsps, outcome.hsps, "{context}");
        assert_eq!(batches[0].query_stats, outcome.query_stats, "{context}");
        assert_eq!(batches[0].rounds, outcome.rounds, "{context}");
        // The writer of `run` and the file of `run_to_file` hold the same bytes.
        let mut written = Vec::new();
        let mut diagnostics = Vec::new();
        run(
            &search,
            inputs,
            &mut [("6", &mut written)],
            &mut diagnostics,
        )
        .unwrap();
        assert!(written == outcome.reports[0].bytes, "{context}");
        let out = directory.join("out.6");
        run_to_file(&search, inputs, "6", &out, &mut Vec::new()).unwrap();
        assert!(
            std::fs::read(&out).unwrap() == outcome.reports[0].bytes,
            "{context}"
        );
    }
    // The rounds of a PSSM query write the command line's PSSM files.
    let case = &PSSM_CASES[2];
    let pssm = case.file(&directory);
    let (api_file, cli_file) = (directory.join("api.txt"), directory.join("cli.txt"));
    let options = [
        case.options,
        &["-out_ascii_pssm", api_file.to_str().unwrap()],
    ]
    .concat();
    let search = validate(case.program, &options).unwrap();
    let database = engine_tests(case.subject);
    let outcome = run_collect(&search, Inputs::pssm_database(&pssm, &database), &["6"]).unwrap();
    let (stdout, _) = case.cli(&pssm, &["-out_ascii_pssm", cli_file.to_str().unwrap()], "6");
    assert!(outcome.reports[0].bytes == stdout);
    assert!(std::fs::read(&api_file).unwrap() == std::fs::read(&cli_file).unwrap());
    // Each format is checked as the command line checks it alone: psiblast's `-in_pssm` writes
    // outfmt 0, 5 to 11 and 20 (`in_pssm::check_in_pssm_options`; the others are item PS15b-2's),
    // whatever the order of the formats.
    let search = validate(Program::Psiblast, &[]).unwrap();
    for formats in [&["6", "1"][..], &["1", "6"], &["0", "15"]] {
        let error =
            run_collect(&search, Inputs::pssm_database(&pssm, &database), formats).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("is not supported by LOSAT's PSIBLAST yet"),
            "{formats:?}: {error}"
        );
    }
    let outcome = run_collect(
        &search,
        Inputs::pssm_database(&pssm, &database),
        &["6", "10", "20"],
    )
    .unwrap();
    assert_eq!(outcome.reports.len(), 3);
    // A PSSM file the command line cannot read is its error, and so is `-query_loc` with it.
    let missing = directory.join("none.asn");
    let search = validate(Program::Psiblast, &[]).unwrap();
    let subjects = read_fasta(
        Program::Psiblast,
        Role::Subject,
        &std::fs::read(fixture("PajaWSV.faa")).unwrap(),
    )
    .unwrap();
    let error = run_collect(&search, Inputs::pssm(&missing, &subjects), &["6"]).unwrap_err();
    let output = command(env!("CARGO_BIN_EXE_losat"))
        .args(["psiblast", "-in_pssm"])
        .arg(&missing)
        .arg("-subject")
        .arg(fixture("PajaWSV.faa"))
        .args(["-outfmt", "6"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(
        error.message().trim_ascii_end(),
        output.stderr.trim_ascii_end(),
        "{error}"
    );
    let ranged = validate(Program::Tblastn, &["-query_loc", "1-10"]).unwrap();
    assert!(matches!(
        run_collect(&ranged, Inputs::pssm(&missing, &[]), &["6"]),
        Err(Error::Argv(_))
    ));
    // Under `-parse_deflines` a subject record with the PSSM query's id is refused, as on the
    // command line (the fixtures `db.ps.q1b.rej.subj.pd`, `db.tbn.inpssm.q1b.rej.subj.pd`): the
    // records of `read_inputs` were read in a scope of their own.
    for (program, subject) in [
        (Program::Psiblast, "fasta/psiblast/q1_subj_ids.faa"),
        (Program::Tblastn, "fasta/psiblast/q1b_gi_nt.fna"),
    ] {
        let search = validate(program, &["-parse_deflines"]).unwrap();
        let (_, subjects) =
            read_inputs(&search, b"", &std::fs::read(engine_tests(subject)).unwrap()).unwrap();
        let pssm = PathBuf::from(engine_tests("fixtures/pssm/q1b_gi.asn"));
        let error = run_collect(&search, Inputs::pssm(&pssm, &subjects), &["6"]).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("a query and a subject record have the Seq-id gi|4504347"),
            "{program}: {error}"
        );
    }
    std::fs::remove_dir_all(&directory).unwrap();
}

/// A query alignment (line LIB item 1.6): `Inputs::msa` and `Inputs::msa_database`, the command
/// line's `-in_msa` of psiblast, against records and a database. The PS6 alignments are those
/// `tests/db_search_manifest/psiblast.tsv` freezes from NCBI 2.17.0; the whole set of them
/// through the API is `tests/db_search_fixtures.py check --library` (rows `db.ps.inmsa.*`).
struct MsaCase {
    /// Under `tests/fixtures/msa`.
    alignment: &'static str,
    /// Under `tests/`: a FASTA subject set (`-subject`), or a database (`-db`).
    subject: &'static str,
    database: bool,
    options: &'static [&'static str],
    /// The frozen NCBI bytes of the outfmt 6 report, a row of the manifest (`db.ps.inmsa.*`), when
    /// the case is that row's.
    frozen: Option<&'static str>,
}

const MSA_CASES: &[MsaCase] = &[
    MsaCase {
        alignment: "ps6_b_db_0.aln",
        subject: "db_search/paja",
        database: true,
        options: &[],
        frozen: Some("db.ps.inmsa.b.db.6"),
    },
    MsaCase {
        alignment: "ps6_b_db_0.aln",
        subject: "fasta/PajaWSV.faa",
        database: false,
        options: &[],
        frozen: None,
    },
    MsaCase {
        alignment: "ps6_b_db_0.aln",
        subject: "db_search/paja",
        database: true,
        options: &["-msa_master_idx", "2"],
        frozen: None,
    },
    MsaCase {
        alignment: "ps6_b_db_0.aln",
        subject: "db_search/paja",
        database: true,
        options: &["-ignore_msa_master"],
        frozen: None,
    },
    // The rounds of an alignment query, with the master chosen.
    MsaCase {
        alignment: "ps6_nc_db_0.aln",
        subject: "db_search/paja",
        database: true,
        options: &["-msa_master_idx", "2", "-num_iterations", "2"],
        frozen: None,
    },
];

impl MsaCase {
    fn file(&self) -> PathBuf {
        PathBuf::from(engine_tests("fixtures/msa")).join(self.alignment)
    }

    /// The command line's report for one `-outfmt`, and its standard error.
    fn cli(&self, extra: &[&str], outfmt: &str) -> (Vec<u8>, Vec<u8>) {
        let output = command(env!("CARGO_BIN_EXE_losat"))
            .args(["psiblast", "-in_msa"])
            .arg(self.file())
            .arg(if self.database { "-db" } else { "-subject" })
            .arg(engine_tests(self.subject))
            .args(self.options)
            .args(extra)
            .arg("-outfmt")
            .arg(outfmt)
            .output()
            .expect("run the command line");
        assert!(
            output.status.success(),
            "{} {:?} -outfmt {outfmt}: {}",
            self.alignment,
            self.options,
            String::from_utf8_lossy(&output.stderr)
        );
        (output.stdout, output.stderr)
    }
}

/// The command line's standard error for a psiblast `-in_msa` run that fails, with `words`.
fn msa_cli_error(alignment: &Path, words: &[&str]) -> Vec<u8> {
    let output = command(env!("CARGO_BIN_EXE_losat"))
        .args(["psiblast", "-in_msa"])
        .arg(alignment)
        .arg("-subject")
        .arg(fixture("PajaWSV.faa"))
        .args(words)
        .args(["-outfmt", "6"])
        .output()
        .unwrap();
    assert!(!output.status.success(), "{alignment:?} {words:?}");
    output.stderr
}

#[test]
fn a_query_alignment_is_the_command_lines() {
    for case in MSA_CASES {
        let context = format!("{} {} {:?}", case.alignment, case.subject, case.options);
        let alignment = case.file();
        let search = validate(Program::Psiblast, case.options).unwrap();
        let subject = engine_tests(case.subject);
        let subjects = match case.database {
            true => Vec::new(),
            false => read_fasta(
                Program::Psiblast,
                Role::Subject,
                &std::fs::read(&subject).unwrap(),
            )
            .unwrap(),
        };
        let inputs = match case.database {
            true => Inputs::msa_database(&alignment, &subject),
            false => Inputs::msa(&alignment, &subjects).named("ignored", &subject),
        };
        let formats = ["6", "0", "7", COLUMNS];
        let outcome = run_collect(&search, inputs, &formats)
            .unwrap_or_else(|error| panic!("{context}: {error}"));
        for (report, outfmt) in outcome.reports.iter().zip(formats) {
            let (stdout, stderr) = case.cli(&[], outfmt);
            assert!(report.bytes == stdout, "{context} -outfmt {outfmt}");
            if outfmt == "6" {
                assert_eq!(
                    outcome.diagnostics, stderr,
                    "{context}: the diagnostics are the command line's"
                );
            }
        }
        // The report NCBI 2.17.0 wrote for the row of the manifest.
        if let Some(row) = case.frozen {
            let frozen = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/db_search")
                .join(format!("{row}.out"));
            assert!(
                outcome.reports[0].bytes == std::fs::read(frozen).unwrap(),
                "{context}: NCBI's bytes"
            );
        }
        // One query, the alignment's master: query 0, with its statistics and rounds.
        assert!(!outcome.hsps.is_empty(), "{context}");
        assert!(outcome.hsps.iter().all(|hsp| hsp.q_idx == 0), "{context}");
        assert_eq!(outcome.query_stats.len(), 1, "{context}");
        assert!(outcome.query_stats[0].is_some(), "{context}");
        let rounds = outcome.hsps.iter().map(|hsp| hsp.round).max().unwrap();
        assert_eq!(outcome.rounds[0].stats.len() as u32, rounds, "{context}");
        if case.options.contains(&"-num_iterations") {
            assert!(rounds > 1, "{context}");
        }
        // One batch with the same records.
        let batches = run_batches(&search, inputs).unwrap();
        assert_eq!(batches.query_count(), 1, "{context}");
        let batches: Vec<_> = batches.map(Result::unwrap).collect();
        assert_eq!(batches.len(), 1, "{context}");
        assert_eq!(batches[0].query_indices, 0..1, "{context}");
        assert_eq!(batches[0].hsps, outcome.hsps, "{context}");
        assert_eq!(batches[0].rounds, outcome.rounds, "{context}");
    }
    // `-ignore_msa_master` and `-msa_master_idx` change the PSSM, so they change the report.
    let report = |options: &[&str]| {
        let search = validate(Program::Psiblast, options).unwrap();
        let alignment = MSA_CASES[0].file();
        let database = engine_tests("db_search/paja");
        run_collect(&search, Inputs::msa_database(&alignment, &database), &["0"])
            .unwrap()
            .reports
            .remove(0)
            .bytes
    };
    let plain = report(&[]);
    assert!(plain != report(&["-msa_master_idx", "2"]));
    assert!(plain != report(&["-ignore_msa_master"]));
}

#[test]
fn a_query_alignment_that_fails_gives_the_command_lines_message() {
    let directory = PathBuf::from(engine_tests("fixtures/msa"));
    let subjects = read_fasta(
        Program::Psiblast,
        Role::Subject,
        &std::fs::read(fixture("PajaWSV.faa")).unwrap(),
    )
    .unwrap();
    let missing = temp_path("none.aln");
    for (alignment, words) in [
        // A file that is none of the formats, a master index past the rows, a missing file.
        (directory.join("ps6_f_unknown.aln"), &[][..]),
        (
            directory.join("ps6_nc_db_0.aln"),
            &["-msa_master_idx", "27"][..],
        ),
        (missing, &[][..]),
    ] {
        let search = validate(Program::Psiblast, words).unwrap();
        let error = run_collect(&search, Inputs::msa(&alignment, &subjects), &["6"]).unwrap_err();
        // The engine's error is the command line's standard error, byte for byte (its last
        // line feed included, as `run_local_blastx`'s errors of the ranges are).
        assert_eq!(
            error.message(),
            msa_cli_error(&alignment, words),
            "{alignment:?} {words:?}: {error}"
        );
    }
    // The alignment is psiblast's: another program refuses it, as it refuses a PSSM it cannot take.
    let alignment = directory.join("ps6_b_db_0.aln");
    for program in [
        Program::Blastp,
        Program::Blastn,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ] {
        let search = validate(program, &[]).unwrap();
        let error = run_collect(&search, Inputs::msa(&alignment, &[]), &["6"]).unwrap_err();
        assert!(
            matches!(error, Error::Argv(_))
                && error
                    .to_string()
                    .contains(&format!("is searched by psiblast, not by {program}")),
            "{program}: {error}"
        );
    }
    // The master options need the alignment, as `-in_msa` is their option on the command line.
    let search = validate(Program::Psiblast, &["-msa_master_idx", "2"]).unwrap();
    let error = run_collect(&search, Inputs::new(&subjects, &subjects), &["6"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error:?}");
    let error = validate(Program::Psiblast, &["-in_msa", "m.aln"]).unwrap_err();
    assert!(error.to_string().contains("(Inputs::msa)"), "{error}");
}

/// Each round's `Gapped` Karlin-Altschul values (lambda, K, H) and effective search space in
/// psiblast's pairwise report of one query, round 1 first.
fn round_footers(report: &str) -> Vec<([f64; 3], i64)> {
    let lines: Vec<&str> = report.lines().collect();
    let mut footers = Vec::new();
    for (at, line) in lines.iter().enumerate() {
        if *line == "Gapped" {
            let values: Vec<f64> = lines[at + 2]
                .split_whitespace()
                .take(3)
                .map(|value| value.parse().unwrap())
                .collect();
            let space = lines[at + 4..]
                .iter()
                .find_map(|line| line.strip_prefix("Effective search space used: "))
                .unwrap()
                .parse()
                .unwrap();
            footers.push(([values[0], values[1], values[2]], space));
        }
    }
    footers
}

/// Whether `value` prints as `printed` with NCBI's three significant digits.
fn prints_as(value: f64, printed: f64) -> bool {
    (value - printed).abs() <= 0.0051 * printed.abs()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1476
// Line RG item Q5: psiblast's rounds in the library: the records carry their round, the rounds of
// a query hold each round's statistics (those its pairwise report prints) and its convergence, and
// round 1 is the search of `-num_iterations 1`.
#[test]
fn psiblast_rounds_carry_their_round_and_its_statistics() {
    let one_round = &CASES[ROUNDS_CASE - 1];
    let rounds_case = &CASES[ROUNDS_CASE];
    let query_bytes = std::fs::read(fixture(rounds_case.query)).unwrap();
    let subject_bytes = std::fs::read(fixture(rounds_case.subject)).unwrap();
    let searches = [
        (
            "subject",
            search_case(one_round, &query_bytes, &subject_bytes, &["0", "7"]),
            search_case(rounds_case, &query_bytes, &subject_bytes, &["0", "7"]),
        ),
        {
            let [one, rounds] = [&DATABASE_CASES[7], &DATABASE_CASES[8]].map(|case| {
                assert_eq!(case.program, Program::Psiblast);
                let (_, records, search) = case.read();
                let db = engine_tests(case.db);
                run_collect(&search, Inputs::database(&records, &db), &["0", "7"]).unwrap()
            });
            ("database", one, rounds)
        },
    ];
    for (context, one, outcome) in searches {
        // One query of two rounds: round 1 is the one-round search.
        assert_eq!(outcome.rounds.len(), 1, "{context}");
        let rounds = &outcome.rounds[0];
        assert_eq!(rounds.stats.len(), 2, "{context}");
        assert!(rounds.converged, "{context}");
        assert!(!one.rounds[0].converged, "{context}");
        assert_eq!(one.rounds[0].stats, [one.query_stats[0]], "{context}");
        assert_eq!(outcome.query_stats, one.query_stats, "{context}");
        assert_eq!(rounds.stats[0], outcome.query_stats[0], "{context}");
        let first: Vec<HspRecord> = outcome
            .hsps
            .iter()
            .filter(|hsp| hsp.round == 1)
            .cloned()
            .collect();
        assert!(first == one.hsps, "{context}: round 1");
        let second = outcome.hsps.len() - first.len();
        assert!(second > 0, "{context}: round 2");
        // Round 2 follows round 1; each round ranks its own HSPs.
        for (index, hsp) in outcome.hsps.iter().enumerate() {
            assert_eq!(hsp.index, index, "{context}");
            let (round, rank) = if index < first.len() {
                (1, index)
            } else {
                (2, index - first.len())
            };
            assert_eq!((hsp.round, hsp.rank), (round, rank), "{context}");
        }
        // The statistics of each round are those that its pairwise report prints.
        let pairwise = String::from_utf8(outcome.reports[0].bytes.clone()).unwrap();
        assert!(pairwise.contains("\nSearch has CONVERGED!\n"), "{context}");
        let footers = round_footers(&pairwise);
        assert_eq!(footers.len(), 2, "{context}");
        for (stats, (printed, space)) in rounds.stats.iter().zip(&footers) {
            let stats = stats.unwrap();
            assert_eq!(stats.eff_space, *space, "{context}");
            for (value, printed) in [stats.lambda, stats.kappa, stats.entropy]
                .iter()
                .zip(printed)
            {
                assert!(prints_as(*value, *printed), "{context}: {value} {printed}");
            }
        }
        assert_ne!(
            rounds.stats[0], rounds.stats[1],
            "{context}: a PSSM's own values"
        );
        // Every row of outfmt 7 lies in its round's part of the report.
        let commented = &outcome.reports[1];
        let text = String::from_utf8(commented.bytes.clone()).unwrap();
        for hsp in &outcome.hsps {
            let row = commented.spans[hsp.index].clone().unwrap();
            let heading = text[..row.start].rfind("# Iteration: ").unwrap();
            let round = text[heading..].lines().next().unwrap();
            assert_eq!(round, format!("# Iteration: {}", hsp.round), "{context}");
        }
    }
}

// Line RG item Q5: the rounds of a query after a record that `-query_loc` skips are the second
// record's, in the whole search and in its batches.
#[test]
fn psiblast_rounds_after_a_skipped_record_are_the_next_records() {
    let case = &DATABASE_CASES[9];
    let (_, records, search) = case.read();
    let db = engine_tests(case.db);
    let inputs = Inputs::database(&records, &db);
    let outcome = run_collect(&search, inputs, &["7"]).unwrap();
    assert_eq!(outcome.rounds.len(), 2);
    assert_eq!(outcome.query_stats[0], None);
    assert_eq!(outcome.rounds[0].stats, [None]);
    assert!(!outcome.rounds[0].converged);
    assert_eq!(outcome.rounds[1].stats.len(), 2);
    assert_eq!(outcome.rounds[1].stats[0], outcome.query_stats[1]);
    assert!(outcome.hsps.iter().all(|hsp| hsp.q_idx == 1));
    assert!(outcome.hsps.iter().any(|hsp| hsp.round == 2));
    let text = String::from_utf8(outcome.reports[0].bytes.clone()).unwrap();
    for hsp in &outcome.hsps {
        let row = outcome.reports[0].spans[hsp.index].clone().unwrap();
        let heading = text[..row.start].rfind("# Iteration: ").unwrap();
        let round = text[heading..].lines().next().unwrap();
        assert_eq!(round, format!("# Iteration: {}", hsp.round));
    }
    let (mut hsps, mut rounds) = (Vec::new(), Vec::new());
    for batch in run_batches(&search, inputs).unwrap() {
        let batch = batch.unwrap();
        hsps.extend(batch.hsps);
        rounds.extend(batch.rounds);
    }
    assert!(hsps == outcome.hsps, "the batches' records");
    assert_eq!(rounds, outcome.rounds, "the batches' rounds");
}

#[test]
fn a_database_record_is_the_database_ordinal_id_in_every_batch() {
    for case in DATABASE_CASES {
        let context = format!("{} {}", case.program, case.db);
        let (_, records, search) = case.read();
        // Enough letters for two query batches of every program but BLASTN (whose batches follow
        // the extensions): each batch names the sequences of its own hit lists.
        let letters: usize = records.iter().map(|record| record.seq().len()).sum();
        let repeated: Vec<FastaRecord> = (0..=40_000 / letters.max(1))
            .flat_map(|copy| {
                records.iter().enumerate().map(move |(index, record)| {
                    FastaRecord::new(
                        format!("Query_{}", copy * 1000 + index + 1),
                        record.title(),
                        record.seq(),
                    )
                })
            })
            .take(if case.program == Program::Blastn {
                8
            } else {
                usize::MAX
            })
            .collect();
        let db = engine_tests(case.db);
        let inputs = Inputs::database(&repeated, &db);
        let whole = run_collect(&search, inputs, &["6"]).unwrap();
        let batches = run_batches(&search, inputs).unwrap();
        assert_eq!(batches.query_count(), repeated.len(), "{context}");
        let (mut hsps, mut stats, mut rounds, mut diagnostics, mut count) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), 0);
        for batch in batches {
            let batch = batch.unwrap_or_else(|error| panic!("{context}: {error}"));
            hsps.extend(batch.hsps);
            stats.extend(batch.query_stats);
            rounds.extend(batch.rounds);
            diagnostics.extend(batch.diagnostics);
            count += 1;
        }
        assert!(
            count >= 2 || case.program == Program::Blastn,
            "{context}: {count} batches"
        );
        assert!(hsps == whole.hsps, "{context}: the records");
        assert_eq!(stats, whole.query_stats, "{context}: the statistics");
        assert_eq!(rounds, whole.rounds, "{context}: the rounds");
        assert_eq!(diagnostics, whole.diagnostics, "{context}: the diagnostics");
    }
}

#[test]
fn a_database_fails_as_on_the_command_line() {
    // The database does not open: NCBI's error, as the command line writes it.
    let case = &DATABASE_CASES[0];
    let (_, records, search) = case.read();
    let missing = engine_tests("seqdb/db/no_such_database");
    let error = run_collect(&search, Inputs::database(&records, &missing), &["6"]).unwrap_err();
    assert!(matches!(error, Error::Search(_)), "{error}");
    let output = command(env!("CARGO_BIN_EXE_losat"))
        .args(["blastn", "-query"])
        .arg(engine_tests(case.query))
        .arg("-db")
        .arg(&missing)
        .output()
        .unwrap();
    assert_eq!(
        error.message(),
        output.stderr,
        "the message of a database that does not open"
    );
    // The batches meet it in their first search.
    let mut batches = run_batches(&search, Inputs::database(&records, &missing)).unwrap();
    assert!(matches!(batches.next(), Some(Err(Error::Search(_)))));
    assert!(batches.next().is_none());
    // Each kind of subject has its options: `validate` holds an argv for either, and the search
    // of the inputs gives the parser's message.
    let only_database = validate(Program::Blastn, &["-taxids", "9606"]).unwrap();
    let error = run_collect(&only_database, Inputs::new(&records, &records), &["6"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    assert!(error.to_string().contains("-taxids"), "{error}");
    let only_subject = validate(Program::Blastn, &["-subject_loc", "1-10"]).unwrap();
    let inputs = Inputs::database(&records, &missing);
    let error = run_collect(&only_subject, inputs, &["6"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    assert!(error.to_string().contains("-subject_loc"), "{error}");
    // Queries without a record are NCBI's warning, after the database opened.
    let db = engine_tests(case.db);
    let outcome = run_collect(&search, Inputs::database(&[], &db), &["6"]).unwrap();
    assert_eq!(outcome.diagnostics, b"Warning: [blastn] Query is Empty!\n");
    assert!(outcome.reports[0].bytes.is_empty());
}

/// `run_makeblastdb` takes the masking options of the command line (`-mask_data`, `-mask_id`,
/// `-mask_desc`): the database it writes has the mask column of the command line's, read back
/// by `blastdbcmd` through the library and through the binary alike.
#[test]
fn a_masked_database_made_through_the_library_equals_the_command_lines() {
    let seqdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb");
    let fasta = seqdb.join("mask_nucl.fna");
    let masks = seqdb.join("mask_writer/dust.asnt");
    let (fasta, masks) = (fasta.to_str().unwrap(), masks.to_str().unwrap());
    // Removed when the test passes, kept when it fails (the databases to compare).
    let directory = TempPath::keep_on_failure(temp_path("masked"));
    std::fs::create_dir_all(&directory).unwrap();
    let library = directory.join("library");
    let binary = directory.join("binary");
    let words = |name: &Path| -> Vec<String> {
        [
            "-in",
            fasta,
            "-dbtype",
            "nucl",
            "-title",
            "masked",
            "-out",
            name.to_str().unwrap(),
            "-mask_data",
            masks,
            "-mask_id",
            "dust1",
            "-mask_desc",
            "the masks of dustmasker",
        ]
        .map(str::to_string)
        .to_vec()
    };
    let argv = words(&library);
    let argv: Vec<&str> = argv.iter().map(String::as_str).collect();
    let (mut log, mut warnings) = (Vec::new(), Vec::new());
    run_makeblastdb(&argv, &mut std::io::empty(), &mut log, &mut warnings).unwrap();
    assert!(
        String::from_utf8_lossy(&log).contains(&format!("Mask file: {masks}\n")),
        "{}",
        String::from_utf8_lossy(&log)
    );
    assert!(warnings.is_empty());
    let cli = command(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .args(words(&binary))
        .output()
        .unwrap();
    assert!(
        cli.status.success(),
        "{}",
        String::from_utf8_lossy(&cli.stderr)
    );
    assert!(cli.stderr.is_empty());
    // The algorithm of `-mask_id` is the first free id from 100 on.
    for name in [&library, &binary] {
        let name = name.to_str().unwrap();
        let words = ["-db", name, "-entry", "all", "-outfmt", "%o %m100"];
        let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
        run_blastdbcmd(&words, &mut std::io::empty(), &mut output, &mut diagnostics).unwrap();
        let text = String::from_utf8_lossy(&output).into_owned();
        assert!(text.starts_with("0 199-281;\n"), "{text}");
        assert!(diagnostics.is_empty());
        let cli = command(env!("CARGO_BIN_EXE_losat"))
            .arg("blastdbcmd")
            .args(words)
            .output()
            .unwrap();
        assert_eq!(output, cli.stdout);
    }
    // The dependencies of the options are the parser's.
    let error = run_makeblastdb(
        &["-dbtype", "nucl", "-mask_id", "x"],
        &mut std::io::empty(),
        &mut std::io::sink(),
        &mut std::io::sink(),
    )
    .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    assert!(error.to_string().contains("mask_data"), "{error}");
}

#[test]
fn a_database_made_through_the_library_is_read_and_searched_through_it() {
    let directory = temp_path("made");
    std::fs::create_dir_all(&directory).unwrap();
    let name = directory.join("made");
    let name = name.to_str().unwrap();
    let fasta = format!(">one a protein\n{PROTEIN}\n>two the same protein\n{PROTEIN}A\n");
    // From FASTA in memory: the command line's standard input.
    let (mut log, mut warnings) = (Vec::new(), Vec::new());
    run_makeblastdb(
        &["-dbtype", "prot", "-title", "made", "-out", name],
        &mut fasta.as_bytes(),
        &mut log,
        &mut warnings,
    )
    .unwrap();
    assert!(
        String::from_utf8_lossy(&log).contains("added 2 sequences"),
        "{}",
        String::from_utf8_lossy(&log)
    );
    assert!(warnings.is_empty());
    // The command line reads what the library wrote, and the two read it alike; a batch file
    // `-` is the caller's reader, the command line's standard input (of a database with
    // accessions: an ordinal id names no sequence of a version 5 database, as in NCBI).
    let accessions = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastdbcmd/db/p5");
    let accessions = accessions.to_str().unwrap();
    for (words, input) in [
        (vec!["-db", name, "-info"], ""),
        (
            vec!["-db", name, "-entry", "all", "-outfmt", "%o %l %t"],
            "",
        ),
        (
            vec!["-db", accessions, "-entry_batch", "-", "-outfmt", "%o %a"],
            "P99999.2\n",
        ),
    ] {
        let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
        run_blastdbcmd(&words, &mut input.as_bytes(), &mut output, &mut diagnostics).unwrap();
        let cli = run_with_stdin(
            command(env!("CARGO_BIN_EXE_losat"))
                .arg("blastdbcmd")
                .args(&words),
            input,
        );
        assert!(cli.status.success(), "{words:?}");
        assert!(!output.is_empty(), "{words:?}");
        assert_eq!(output, cli.stdout, "{words:?}");
        assert_eq!(diagnostics, cli.stderr, "{words:?}");
    }
    // A search of it names its sequences by their ordinal ids.
    let search = validate(Program::Blastp, &["-max_hsps", "1"]).unwrap();
    let query = read_fasta(Program::Blastp, Role::Query, fasta.as_bytes()).unwrap();
    let outcome = run_collect(&search, Inputs::database(&query, name), &["6"]).unwrap();
    let mut pairs: Vec<(usize, usize, usize)> = outcome
        .hsps
        .iter()
        .map(|hsp| (hsp.q_idx, hsp.s_idx, hsp.subject_length))
        .collect();
    pairs.sort_unstable();
    let long = PROTEIN.len();
    assert_eq!(
        pairs,
        [
            (0, 0, long),
            (0, 1, long + 1),
            (1, 0, long),
            (1, 1, long + 1)
        ]
    );
    // A failure's message is what the command line writes to standard error.
    let missing = directory.join("none");
    let missing = missing.to_str().unwrap();
    let words = ["-db", missing, "-info"];
    let mut diagnostics = Vec::new();
    let error = run_blastdbcmd(
        &words,
        &mut std::io::empty(),
        &mut std::io::sink(),
        &mut diagnostics,
    )
    .unwrap_err();
    let cli = command(env!("CARGO_BIN_EXE_losat"))
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    assert!(!cli.status.success());
    assert_eq!(diagnostics, cli.stderr);
    assert!(matches!(error, Error::Search(_)), "{error}");
    assert_eq!(error.message(), cli.stderr.trim_ascii_end(), "{error}");
    // An argv that the parser refuses is the parser's message.
    let error = run_makeblastdb(
        &["-dbtype", "neither"],
        &mut std::io::empty(),
        &mut std::io::sink(),
        &mut std::io::sink(),
    )
    .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    std::fs::remove_dir_all(&directory).unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:254-258
/// A tool's error is the bytes that the command line writes to standard error, also a byte that
/// is not UTF-8: blastdbcmd reads `%é` as the specifier `%` and the first byte of `é` (C3), a
/// `char` that NCBI's stream writes alone into its error.
#[test]
fn a_tool_error_holds_the_bytes_that_the_command_line_writes() {
    let db = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/prot_v5");
    let words = [
        "-db",
        db.to_str().unwrap(),
        "-entry",
        "all",
        "-outfmt",
        "%é",
    ];
    let mut diagnostics = Vec::new();
    let error = run_blastdbcmd(
        &words,
        &mut std::io::empty(),
        &mut std::io::sink(),
        &mut diagnostics,
    )
    .unwrap_err();
    let cli = command(env!("CARGO_BIN_EXE_losat"))
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    assert!(!cli.status.success());
    assert_eq!(diagnostics, cli.stderr);
    let written = cli.stderr.strip_suffix(b"\n").unwrap();
    assert!(
        written.ends_with(b"Unrecognized format specification: '%\xC3'"),
        "{:?}",
        String::from_utf8_lossy(written)
    );
    assert!(matches!(error, Error::Search(_)), "{error:?}");
    assert_eq!(error.message(), written, "{error:?}");
    assert_eq!(error.to_string(), String::from_utf8_lossy(written));
}

#[test]
fn the_api_owns_the_inputs_and_the_outputs() {
    for owned in API_OWNED {
        let flag = format!("-{owned}");
        // `-in_pssm` is psiblast's and tblastn's (`Inputs::pssm`); blastp's parser knows no such
        // option.
        let program = match *owned {
            "in_pssm" => Program::Psiblast,
            _ => Program::Blastp,
        };
        let error = validate(program, &[&flag, "x"]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{flag}: {error}");
        assert!(error.to_string().starts_with(&flag), "{flag}: {error}");
        for program in Program::ALL {
            assert!(
                describe(program)
                    .parameters
                    .iter()
                    .all(|parameter| parameter.flag != flag),
                "{program} describes {flag}"
            );
        }
    }
    // An option the command line rejects is the command line's message; a program the command
    // line has and the API has not is named.
    let error = validate(Program::Tblastx, &["-evalue"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    let error = Program::parse("makeblastdb").unwrap_err();
    assert!(error.to_string().contains("run_makeblastdb"), "{error}");
    let error = Program::parse("dustmasker").unwrap_err();
    assert!(error.to_string().contains("run_dustmasker"), "{error}");
}

/// The output of `command` run with `input` as its standard input.
fn run_with_stdin(command: &mut Command, input: &str) -> std::process::Output {
    use std::io::Write as _;
    let mut child = command
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

// `tests/masker_fixtures.rs` runs every row of the NCBI fixture through the library too
// (`examples/library_command_line`); this is the API's own contract: the caller's reader is the
// standard input, and a failure's message is what the command line writes to standard error.
#[test]
fn dustmasker_masks_sequences_in_memory_as_the_command_line_does() {
    let fasta =
        ">one low complexity\nACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACGTAGCTAGCATCGACT\n";
    for words in [vec![], vec!["-in", "-", "-outfmt", "fasta"]] {
        let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
        run_dustmasker(&words, &mut fasta.as_bytes(), &mut output, &mut diagnostics).unwrap();
        let cli = run_with_stdin(
            command(env!("CARGO_BIN_EXE_losat"))
                .arg("dustmasker")
                .args(&words),
            fasta,
        );
        assert!(cli.status.success(), "{words:?}");
        assert!(!output.is_empty(), "{words:?}");
        assert_eq!(output, cli.stdout, "{words:?}");
        assert_eq!(diagnostics, cli.stderr, "{words:?}");
    }
    // A failure (an `-in` that does not open: NCBI's argument report, exit status 1 there).
    let words = ["-in", "no/such/file.fa"];
    let mut diagnostics = Vec::new();
    let error = run_dustmasker(
        &words,
        &mut std::io::empty(),
        &mut std::io::sink(),
        &mut diagnostics,
    )
    .unwrap_err();
    let cli = run_with_stdin(
        command(env!("CARGO_BIN_EXE_losat"))
            .arg("dustmasker")
            .args(words),
        "",
    );
    assert!(!cli.status.success());
    assert_eq!(diagnostics, cli.stderr);
    assert!(matches!(error, Error::Search(_)), "{error}");
    assert_eq!(error.message(), cli.stderr.trim_ascii_end(), "{error}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:209-214 (the subjects'
// `GetAllSeqs` turns the reader's `eNoDefline` into its own input error, with the NUL of `ends`);
// the queries keep the reader's error. `read_inputs` ends with the command line's message for the
// programs of either reader of the command line (BLASTX reads its files with its own).
#[test]
fn read_inputs_without_a_first_defline_end_with_the_command_line_message() {
    for (program, query, subject) in [
        (Program::Blastn, "ACGTACGTAC", "ACGTACGTAC"),
        (Program::Blastx, "ACGTACGTAC", "MKLVEFMKLV"),
    ] {
        let search = validate(program, &["-parse_deflines"]).unwrap();
        let clean_query = format!(">q1\n{query}\n");
        let clean_subject = format!(">s1\n{subject}\n");
        let bad_subject = format!(";c\n{subject}\n>s2\n{subject}\n");
        let bad_query = format!(";c\n\n{query}\n");
        let message = |error: Error| match error {
            Error::Input(message) => String::from_utf8(message).expect("a UTF-8 message"),
            other => panic!("{program:?}: not an input error: {other}"),
        };
        let error =
            read_inputs(&search, clean_query.as_bytes(), bad_subject.as_bytes()).unwrap_err();
        assert_eq!(
            message(error),
            "LOSAT query/options error: Query input doesn't start with a defline or comment, \
             line 2\0\nPlease refer to the LOSAT user manual.",
            "{program:?}"
        );
        // The subjects are read first: a bad query and a bad subject end with the subjects' error.
        let error = read_inputs(&search, bad_query.as_bytes(), bad_subject.as_bytes()).unwrap_err();
        assert_eq!(
            message(error),
            "LOSAT query/options error: Query input doesn't start with a defline or comment, \
             line 2\0\nPlease refer to the LOSAT user manual.",
            "{program:?}"
        );
        let error =
            read_inputs(&search, bad_query.as_bytes(), clean_subject.as_bytes()).unwrap_err();
        assert_eq!(
            message(error),
            "LOSAT query error: CFastaReader: Input doesn't start with a defline or comment \
             around line 3",
            "{program:?}"
        );
    }
}

// Line LIB E1 (2026-11-21): the message of each producer of `Error` against the command line's
// standard error for the same failure, run in the same test and compared byte for byte by the
// rules of `Error` (NCBI's message with its newlines; LOSAT's refusal without `Error: ` and the
// newline; the parser's text; a reader's or a tool's text without its last newlines).

/// The small records of `tests/fasta/outfmt11` that a program searches in the tests of failures.
fn small_inputs(program: Program) -> (PathBuf, PathBuf) {
    let (query, subject) = match program {
        Program::Blastp | Program::Psiblast => ("arch_bp_q.faa", "arch_bp_s.faa"),
        Program::Blastx => ("arch_bx_q.fna", "arch_bx_s.faa"),
        Program::Tblastn => ("arch_tbn_q.faa", "arch_tbn_s.fna"),
        _ => ("arch_tbx_q.fna", "arch_tbx_s.fna"),
    };
    (
        fixture(&format!("outfmt11/{query}")),
        fixture(&format!("outfmt11/{subject}")),
    )
}

/// `losat <program> <words>` (`command`'s environment).
fn cli_output(program: &str, words: &[&str]) -> std::process::Output {
    command(env!("CARGO_BIN_EXE_losat"))
        .arg(program)
        .args(words)
        .output()
        .expect("run the command line")
}

/// The rule of `Error` by which the command line's standard error `stderr` is the library's
/// `diagnostics` (what the search wrote before the failure) and then the message of `error`:
/// NCBI's message as it is (`"ncbi"`), LOSAT's refusal after `Error: ` and before a newline
/// (`"losat"`), or the reader's error of a file that does not open (`read_input_file`: an
/// `Error::Input` of NCBI's `Command line argument error: `) without its last newlines
/// (`"reader"`); `None` when none gives it.
fn failure_rule(stderr: &[u8], diagnostics: &[u8], error: &Error) -> Option<&'static str> {
    let tail = stderr.strip_prefix(diagnostics)?;
    if tail == error.message() {
        Some("ncbi")
    } else if tail == [b"Error: ", error.message(), b"\n"].concat() {
        Some("losat")
    } else if matches!(error, Error::Input(_))
        && error
            .message()
            .starts_with(b"Command line argument error: ")
        && without_last_newlines(tail) == error.message()
    {
        Some("reader")
    } else {
        None
    }
}

/// `bytes` without its last newlines (the command line's text of a reader's or a tool's error).
fn without_last_newlines(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .rposition(|&byte| byte != b'\n')
        .map_or(0, |at| at + 1);
    &bytes[..end]
}

/// The library's search of the files `query` and `subject` with `options`, as
/// `examples/library_command_line` makes it (`read_input_file`, `read_inputs`, `run`): the result,
/// the report of `outfmt` and the diagnostics.
fn search_of_files(
    program: Program,
    options: &[&str],
    query: &Path,
    subject: &Path,
    outfmt: &str,
) -> (Result<(), Error>, Vec<u8>, Vec<u8>) {
    let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
    let (query_name, subject_name) = (query.to_str().unwrap(), subject.to_str().unwrap());
    let result = validate(program, options).and_then(|search| {
        let query_bytes = read_input_file(Role::Query, query)?;
        let subject_bytes = read_input_file(Role::Subject, subject)?;
        let (queries, subjects) = read_inputs(&search, &query_bytes, &subject_bytes)?;
        let inputs = Inputs::new(&queries, &subjects).named(query_name, subject_name);
        run(
            &search,
            inputs,
            &mut [(outfmt, &mut report)],
            &mut diagnostics,
        )
    });
    (result, report, diagnostics)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
/// `run_to_file` to a file whose writes fail (`/dev/full`, through a link in a directory of the
/// test, where outfmt 13 and 14 write their side files) ends as the command line's `-out` (line
/// LIB E1, fix 1; line RG item T6a, the command line's stream), in every program and every outfmt:
/// the call succeeds where the command line exits 0, with the posts that it writes in the
/// diagnostics; where the command line exits 6, the message is NCBI's text of that end (the
/// failed region's posts, `BLAST failed to write output` and its newline, and for outfmt 11 the
/// error archive's write fault after it); where it aborts, `std::terminate`'s text. The warning
/// of `-max_target_seqs 1`, which the command line writes before the failure, is in the
/// diagnostics and not in the message.
#[cfg(target_os = "linux")]
#[test]
fn a_report_that_cannot_be_written_fails_as_the_command_lines_out() {
    if !Path::new("/dev/full").exists() {
        return;
    }
    let directory = TempPath::new(temp_path("full"));
    std::fs::create_dir_all(&directory).unwrap();
    let full = directory.join("full");
    std::os::unix::fs::symlink("/dev/full", &full).unwrap();
    let full_name = full.to_str().unwrap();
    for program in Program::ALL {
        let (query, subject) = small_inputs(program);
        let (query_name, subject_name) = (query.to_str().unwrap(), subject.to_str().unwrap());
        let search = validate(program, &["-max_target_seqs", "1"]).unwrap();
        let query_bytes = std::fs::read(&query).unwrap();
        let subject_bytes = std::fs::read(&subject).unwrap();
        let (queries, subjects) = read_inputs(&search, &query_bytes, &subject_bytes).unwrap();
        for number in 0..=18 {
            let outfmt = number.to_string();
            let context = format!("{program} -outfmt {outfmt}");
            let cli = cli_output(
                program.name(),
                &[
                    "-query",
                    query_name,
                    "-subject",
                    subject_name,
                    "-max_target_seqs",
                    "1",
                    "-outfmt",
                    outfmt.as_str(),
                    "-out",
                    full_name,
                ],
            );
            let mut diagnostics = Vec::new();
            let inputs = Inputs::new(&queries, &subjects).named(query_name, subject_name);
            let result = run_to_file(&search, inputs, &outfmt, &full, &mut diagnostics);
            let stderr = String::from_utf8_lossy(&cli.stderr);
            let error = match result {
                Ok(()) => {
                    assert_eq!(cli.status.code(), Some(0), "{context}: {stderr}");
                    assert_eq!(diagnostics, cli.stderr, "{context}");
                    continue;
                }
                Err(error) => error,
            };
            assert_ne!(cli.status.code(), Some(0), "{context}: {error:?}");
            assert!(
                failure_rule(&cli.stderr, &diagnostics, &error).is_some(),
                "{context}: the command line's {:?} {stderr:?}, the library's {:?} and {error:?}",
                cli.status.code(),
                String::from_utf8_lossy(&diagnostics),
            );
            if number == 0 {
                assert!(
                    String::from_utf8_lossy(&diagnostics)
                        .contains("Examining 5 or more matches is recommended"),
                    "{context}: {stderr}"
                );
                assert_eq!(cli.status.code(), Some(6), "{context}: {stderr}");
                assert_eq!(
                    error.message(),
                    b"BLAST failed to write output\n",
                    "{context}: {error:?}"
                );
            }
            if cli.status.code() == Some(6) {
                // The posts of the failed region, NCBI's `BLAST failed to write output`, and for
                // outfmt 11 the write fault of the error archive written on the failed stream.
                const FAILED: &[u8] = b"BLAST failed to write output\n";
                let message = error.message();
                let end = message
                    .windows(FAILED.len())
                    .position(|window| window == FAILED)
                    .map(|at| at + FAILED.len());
                assert!(
                    end.is_some_and(|end| end == message.len()
                        || (number == 11 && message[end..].starts_with(b"Error: (802.5) "))),
                    "{context}: {error:?}"
                );
            }
        }
    }
}

/// An organism report (outfmt 18) longer than the `-out` file's buffer, which writes at the
/// insertion that reaches its 8,191st byte (`cli::NcbiReportBuffer`; here the query's title of
/// 9,600 bytes), fails where that write fails, as the command line's `-out`: NCBI's `BLAST failed
/// to write output` and its newline, exit status 6 there (line LIB E1, review N5). A short report,
/// held to the end, passes on both sides
/// (`a_report_that_cannot_be_written_fails_as_the_command_lines_out`).
#[cfg(target_os = "linux")]
#[test]
fn a_long_organism_report_that_cannot_be_written_fails_as_the_command_lines_out() {
    if !Path::new("/dev/full").exists() {
        return;
    }
    let directory = TempPath::new(temp_path("full18"));
    std::fs::create_dir_all(&directory).unwrap();
    let full = directory.join("full");
    std::os::unix::fs::symlink("/dev/full", &full).unwrap();
    let (query, subject) = small_inputs(Program::Blastp);
    let records = std::fs::read_to_string(&query).unwrap();
    let titled = directory.join("titled.faa");
    let title = format!(" {}\n", "abcdefg ".repeat(1200).trim_end());
    std::fs::write(&titled, records.replacen('\n', &title, 1)).unwrap();
    let (query_name, subject_name) = (titled.to_str().unwrap(), subject.to_str().unwrap());
    let cli = cli_output(
        "blastp",
        &[
            "-query",
            query_name,
            "-subject",
            subject_name,
            "-outfmt",
            "18",
            "-out",
            full.to_str().unwrap(),
        ],
    );
    assert_eq!(cli.status.code(), Some(6), "{cli:?}");
    assert!(
        cli.stderr.ends_with(b"BLAST failed to write output\n"),
        "{cli:?}"
    );
    let search = validate(Program::Blastp, &[]).unwrap();
    let query_bytes = std::fs::read(&titled).unwrap();
    let subject_bytes = std::fs::read(&subject).unwrap();
    let (queries, subjects) = read_inputs(&search, &query_bytes, &subject_bytes).unwrap();
    let inputs = Inputs::new(&queries, &subjects).named(query_name, subject_name);
    let mut diagnostics = Vec::new();
    let error = run_to_file(&search, inputs, "18", &full, &mut diagnostics).unwrap_err();
    assert_eq!(
        (
            error.message(),
            failure_rule(&cli.stderr, &diagnostics, &error)
        ),
        (&b"BLAST failed to write output\n"[..], Some("ncbi")),
        "{error:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
/// A directory named as the query or the subject is read as the command line reads it (line LIB
/// E1, fix 2): it opens and its first read fails, which ends the input, so it holds no record
/// (`read_input_file`: blastp's query of a directory is NCBI's `Query is Empty!`); a directory
/// named `.gz`, `.bz2` or `.zst` is refused as not of that format, as the command line refuses
/// it. On Windows a directory does not open, whatever its name (NCBI's `ifstream` there opens
/// none, and the open comes before the format): it is a file that is not accessible. Every
/// program, both sides, against the command line in the same test.
#[test]
fn a_directory_input_is_read_as_the_command_line_reads_it() {
    let directory = TempPath::new(temp_path("directories"));
    let names = ["records", "records.gz", "records.bz2", "records.zst"];
    for name in names {
        std::fs::create_dir_all(directory.join(name)).unwrap();
    }
    for program in Program::ALL {
        let (query, subject) = small_inputs(program);
        for name in names {
            let named = directory.join(name);
            for (query, subject) in [(&named, &subject), (&query, &named)] {
                let (query_name, subject_name) =
                    (query.to_str().unwrap(), subject.to_str().unwrap());
                let context = format!("{program} -query {query_name} -subject {subject_name}");
                let cli = cli_output(
                    program.name(),
                    &[
                        "-query",
                        query_name,
                        "-subject",
                        subject_name,
                        "-outfmt",
                        "6",
                    ],
                );
                let (result, report, diagnostics) =
                    search_of_files(program, &[], query, subject, "6");
                let stderr = String::from_utf8_lossy(&cli.stderr);
                match result {
                    Ok(()) => {
                        assert_eq!(cli.status.code(), Some(0), "{context}: {stderr}");
                        assert_eq!(report, cli.stdout, "{context}");
                        assert_eq!(diagnostics, cli.stderr, "{context}");
                    }
                    Err(error) => {
                        assert_ne!(cli.status.code(), Some(0), "{context}: {error:?}");
                        assert!(cli.stdout.is_empty() && report.is_empty(), "{context}");
                        assert!(
                            failure_rule(&cli.stderr, &diagnostics, &error).is_some(),
                            "{context}: the command line's {stderr:?}, the library's {:?} and \
                             {error:?}",
                            String::from_utf8_lossy(&diagnostics),
                        );
                        if name != "records" {
                            assert!(matches!(error, Error::Input(_)), "{context}: {error:?}");
                        }
                    }
                }
                if cfg!(windows) {
                    let not_accessible = format!(
                        "Command line argument error: Argument \"{}\". File is not \
                         accessible:  `{}'\n",
                        if query == &named { "query" } else { "subject" },
                        named.display()
                    );
                    assert_eq!(
                        (cli.status.code(), stderr.as_ref()),
                        (Some(1), not_accessible.as_str()),
                        "{context}"
                    );
                } else if program == Program::Blastp && name == "records" && query == &named {
                    assert_eq!(
                        (cli.status.code(), cli.stderr.as_slice()),
                        (Some(0), b"Warning: [blastp] Query is Empty!\n".as_slice()),
                        "{context}"
                    );
                }
            }
        }
    }
}

/// `run_makeblastdb`'s argv that its parser refuses is a tool's error of the parser's kind (line
/// LIB E1, fix 4): the command line's standard error goes to the diagnostics, and the message is
/// that text without its last newlines, also where the caller's diagnostics writer fails.
/// `examples/library_command_line` writes that text once, as the command line does (its status is
/// 1 for every failure).
#[test]
fn makeblastdbs_refused_argv_is_a_tools_error() {
    /// A diagnostics writer that takes nothing.
    struct Failing;
    impl std::io::Write for Failing {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("closed"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let example = example_binary("library_command_line");
    for words in [
        &["-nosuch"][..],
        &["-dbtype", "foo", "-in", "x"],
        &["-in", "x"],
    ] {
        let mut diagnostics = Vec::new();
        let error = run_makeblastdb(
            words,
            &mut std::io::empty(),
            &mut std::io::sink(),
            &mut diagnostics,
        )
        .unwrap_err();
        let cli = cli_output("makeblastdb", words);
        assert_eq!(cli.status.code(), Some(1), "{words:?}");
        assert!(matches!(error, Error::Argv(_)), "{words:?}: {error:?}");
        assert_eq!(diagnostics, cli.stderr, "{words:?}");
        assert!(cli.stderr.ends_with(b"\n"), "{words:?}");
        assert_eq!(
            error.message(),
            without_last_newlines(&cli.stderr),
            "{words:?}"
        );
        let failing = run_makeblastdb(
            words,
            &mut std::io::empty(),
            &mut std::io::sink(),
            &mut Failing,
        )
        .unwrap_err();
        assert_eq!(failing, error, "{words:?}: a failing writer");
        let library = command(example.to_str().unwrap())
            .arg("makeblastdb")
            .args(words)
            .output()
            .expect("run the example");
        assert_eq!(
            (library.status.code(), library.stdout, library.stderr),
            (Some(1), cli.stdout, cli.stderr),
            "{words:?}: the example"
        );
    }
}

/// A message of NCBI's is its bytes with its newline (`Error`'s first rule): the engine's error
/// of a `-query_loc` that ends before it starts (exit status 3), and the formatter's of an
/// `-outfmt` out of range (an `Error::Argv`, exit status 255).
#[test]
fn an_ncbi_message_keeps_its_newline() {
    let (query, subject) = small_inputs(Program::Blastp);
    let (query_name, subject_name) = (query.to_str().unwrap(), subject.to_str().unwrap());
    for (options, outfmt, status) in [(&["-query_loc", "5-3"][..], "6", 3), (&[], "99", 255)] {
        let context = format!("{options:?} -outfmt {outfmt}");
        let mut words = vec![
            "-query",
            query_name,
            "-subject",
            subject_name,
            "-outfmt",
            outfmt,
        ];
        words.extend_from_slice(options);
        let cli = cli_output("blastp", &words);
        let (result, _, diagnostics) =
            search_of_files(Program::Blastp, options, &query, &subject, outfmt);
        let error = result.unwrap_err();
        assert_eq!(cli.status.code(), Some(status), "{context}");
        assert!(error.message().ends_with(b"\n"), "{context}: {error:?}");
        assert_eq!(
            [diagnostics.as_slice(), error.message()].concat(),
            cli.stderr,
            "{context}"
        );
        let kind_is_argv = matches!(error, Error::Argv(_));
        assert_eq!(kind_is_argv, outfmt == "99", "{context}: {error:?}");
    }
}

/// LOSAT's own refusal is its text without the `Error: ` and the newline that the command line
/// writes around it (`Error`'s second rule): two query records with one Seq-id under
/// `-parse_deflines`.
#[test]
fn a_losat_refusal_is_its_text_without_error_and_the_newline() {
    let records = TempPath::new(temp_path("same_ids.faa"));
    std::fs::write(&records, ">q1 a\nMKVLAAGIVG\n>q1 b\nMKVLAAGIVG\n").unwrap();
    let (_, subject) = small_inputs(Program::Blastp);
    let words = [
        "-query",
        records.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-parse_deflines",
        "-outfmt",
        "6",
    ];
    let cli = cli_output("blastp", &words);
    let (result, _, diagnostics) = search_of_files(
        Program::Blastp,
        &["-parse_deflines"],
        &records,
        &subject,
        "6",
    );
    let error = result.unwrap_err();
    assert_eq!(cli.status.code(), Some(1));
    assert!(matches!(error, Error::Search(_)), "{error:?}");
    assert!(diagnostics.is_empty());
    assert_eq!(
        cli.stderr,
        [b"Error: ".as_slice(), error.message(), b"\n"].concat(),
        "{error:?}"
    );
}

/// NCBI's reader's report is the command line's standard error without its last newline
/// (`Error`'s third rule; line RG item K4i-2b): a value that the option does not convert and an
/// option that the program does not have, both NCBI's usage and posts with status 1.
#[test]
fn the_readers_report_is_the_command_lines_text() {
    for words in [&["-evalue", "abc"][..], &["-nosuch"]] {
        let error = validate(Program::Blastp, words).unwrap_err();
        let cli = cli_output("blastp", words);
        assert_eq!(cli.status.code(), Some(1), "{words:?}");
        assert!(matches!(error, Error::Argv(_)), "{words:?}: {error:?}");
        assert_eq!(cli.stderr, [error.message(), b"\n"].concat(), "{words:?}");
        assert!(
            cli.stderr
                .ends_with(b"\nError:  Unknown argument: \"nosuch\"\n")
                || cli
                    .stderr
                    .ends_with(b"Argument cannot be converted:  `abc'\n"),
            "{words:?}: {}",
            String::from_utf8_lossy(&cli.stderr)
        );
    }
}

/// A reader's error is the command line's text without its last newline (`Error`'s fourth rule):
/// a NUL in a query record; and LOSAT's refusal of a file named `.gz` that is not gzip data is
/// also without the command line's `Error: `.
#[test]
fn an_input_error_is_the_command_lines_text_without_its_last_newline() {
    let (_, subject) = small_inputs(Program::Blastp);
    let nul = TempPath::new(temp_path("nul.faa"));
    std::fs::write(&nul, b">q\nMK\0V\n").unwrap();
    let junk = TempPath::new(temp_path("junk.faa.gz"));
    std::fs::write(
        &junk,
        b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03garbagegarbagegarbage",
    )
    .unwrap();
    for (query, prefix) in [(&nul, ""), (&junk, "Error: ")] {
        let words = [
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            "6",
        ];
        let cli = cli_output("blastp", &words);
        let (result, _, diagnostics) = search_of_files(Program::Blastp, &[], query, &subject, "6");
        let error = result.unwrap_err();
        assert_eq!(cli.status.code(), Some(1), "{query:?}");
        assert!(matches!(error, Error::Input(_)), "{query:?}: {error:?}");
        assert!(diagnostics.is_empty(), "{query:?}");
        assert_eq!(
            cli.stderr,
            [prefix.as_bytes(), error.message(), b"\n"].concat(),
            "{query:?}: {error:?}"
        );
    }
}

/// NCBI's application start that ends the program (a `-conffile` that does not exist: exit status
/// 255) is the start's text without its last newline.
#[test]
fn a_start_that_ends_the_program_is_its_text_without_the_last_newline() {
    let (query, subject) = small_inputs(Program::Blastp);
    let missing = temp_path("missing.ini");
    let options = ["-conffile", missing.to_str().unwrap()];
    let mut words = vec![
        "-query",
        query.to_str().unwrap(),
        "-subject",
        subject.to_str().unwrap(),
        "-outfmt",
        "6",
    ];
    words.extend(options);
    let cli = cli_output("blastp", &words);
    let (result, _, diagnostics) =
        search_of_files(Program::Blastp, &options, &query, &subject, "6");
    let error = result.unwrap_err();
    assert_eq!(cli.status.code(), Some(255));
    assert!(matches!(error, Error::Search(_)), "{error:?}");
    assert!(diagnostics.is_empty());
    assert_eq!(cli.stderr, [error.message(), b"\n"].concat(), "{error:?}");
}

/// A tool of the library (`losat::api::run_<tool>`).
type ToolCall = fn(
    &[&str],
    &mut dyn std::io::Read,
    &mut dyn std::io::Write,
    &mut dyn std::io::Write,
) -> Result<(), Error>;

// One rule for a help word in the library's tools (line RG item K4i-2b round 2, R-M1): each writes
// the help that the command line prints to its output and returns `Ok`, as the command line exits
// with status 0, makeblastdb included (it returned the help as `Error::Argv` before NCBI's reader
// read its arguments). A search has no output for it: `validate` gives the help as `Error::Argv`.
#[test]
fn a_help_word_is_the_tools_output_and_the_searchs_error() {
    let tools: [(&str, ToolCall); 7] = [
        ("makeblastdb", run_makeblastdb),
        ("blastdbcmd", run_blastdbcmd),
        ("dustmasker", run_dustmasker),
        ("segmasker", run_segmasker),
        ("windowmasker", run_windowmasker),
        ("convert2blastmask", run_convert2blastmask),
        ("makembindex", run_makembindex),
    ];
    for (name, tool) in tools {
        for word in ["-h", "-help"] {
            let cli = command(env!("CARGO_BIN_EXE_losat"))
                .args([name, word])
                .stdin(std::process::Stdio::null())
                .output()
                .unwrap();
            assert_eq!(cli.status.code(), Some(0), "{name} {word}");
            assert!(!cli.stdout.is_empty(), "{name} {word}");
            let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
            tool(
                &[word],
                &mut std::io::empty(),
                &mut output,
                &mut diagnostics,
            )
            .unwrap();
            assert_eq!(
                (output, diagnostics),
                (cli.stdout, cli.stderr),
                "{name} {word}"
            );
        }
    }
    for word in ["-h", "-help"] {
        let cli = command(env!("CARGO_BIN_EXE_losat"))
            .args(["blastp", word])
            .stdin(std::process::Stdio::null())
            .output()
            .unwrap();
        assert_eq!(cli.status.code(), Some(0), "blastp {word}");
        match validate(Program::Blastp, &[word]) {
            Err(Error::Argv(help)) => assert_eq!(help, cli.stdout, "blastp {word}"),
            other => panic!("blastp {word}: {other:?}"),
        }
    }
}

/// A Smith-Waterman alignment of `-use_sw_tback` whose query row starts with 95 gaps (f101, 190
/// letters of nest from its 96th, against nest; line RG item X4b-gap): the records' length, gap
/// openings and gaps are the table's, which count only the columns in the query row's range
/// (190, 0, 0), for `-subject` (BLASTX and TBLASTN map the HSP, BLASTP keeps its end gaps) and
/// for `-db` (where every HSP keeps them), each report the command line's.
#[test]
fn a_query_row_that_starts_with_gaps_counts_as_the_table() {
    let options: &[&str] = &["-use_sw_tback"];
    let check = |context: &str, outcome: &Outcome, cli: Vec<u8>| {
        let table = &outcome.reports[0].bytes;
        assert!(*table == cli, "{context}: the report");
        let rows = rows(table);
        assert_eq!(rows.len(), outcome.hsps.len(), "{context}");
        for (hsp, row) in outcome.hsps.iter().zip(rows) {
            check_columns(context, hsp, &table[row]);
        }
        assert!(
            outcome
                .hsps
                .iter()
                .any(|hsp| (hsp.length, hsp.gap_open, hsp.gaps) == (190, 0, 0)),
            "{context}: the alignment of f101"
        );
    };
    for (program, query, subject) in [
        (Program::Blastp, "swgap/f101.faa", "rps/nest.faa"),
        (Program::Blastx, "swgap/f101.fna", "rps/nest.faa"),
        (Program::Tblastn, "swgap/f101.faa", "rps/nest_p.fna"),
    ] {
        let case = Case {
            program,
            query,
            subject,
            options,
            formats: &[COLUMNS],
        };
        let context = format!("{program} -subject");
        let (query, subject) = (fixture(query), fixture(subject));
        let outcome = api(&case, &query, &subject, &[]);
        check(&context, &outcome, cli(&case, &query, &subject, COLUMNS).0);
    }
    for (program, query, db) in [
        (
            Program::Blastp,
            "fasta/swgap/f101.faa",
            "db_search/swgap/nestaa",
        ),
        (
            Program::Blastx,
            "fasta/swgap/f101.fna",
            "db_search/swgap/nestaa",
        ),
        (
            Program::Tblastn,
            "fasta/swgap/f101.faa",
            "db_search/swgap/nestp",
        ),
    ] {
        let case = DatabaseCase {
            program,
            query,
            db,
            options,
        };
        let context = format!("{program} -db");
        let (_, records, search) = case.read();
        let (db, query_name) = (engine_tests(db), engine_tests(query));
        let inputs = Inputs::database(&records, &db).named(&query_name, "ignored");
        let outcome = run_collect(&search, inputs, &[COLUMNS]).unwrap();
        check(&context, &outcome, case.cli(COLUMNS));
    }
}
