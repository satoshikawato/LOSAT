//! The library API (`LOSAT::api`, docs/losat_library_surface_plan.md §6): a search of records in
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

use LOSAT::api::{
    describe, read_fasta, run, run_collect, validate, Error, FastaRecord, HspRecord, Inputs,
    Outcome, Program, Role, Strand, API_OWNED,
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
];

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

/// A command of this crate's binaries without the environment that changes LOSAT's behavior
/// (`LOSAT_BLASTX_WORKER_LOG`, which the `blastx-worker-probe` feature requires, stays).
fn command(binary: &str) -> Command {
    let mut command = Command::new(binary);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text == "LOSAT_BLASTX_WORKER_LOG" {
            continue;
        }
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLASTDB") {
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

/// The line ranges of a tabular report: one per row.
fn rows(report: &[u8]) -> Vec<std::ops::Range<usize>> {
    let mut rows = Vec::new();
    let mut start = 0;
    for (at, &byte) in report.iter().enumerate() {
        if byte == b'\n' && !report[start..at].starts_with(b"#") {
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
    assert!(shown_methods > 0, "a case with composition-based statistics");
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
     sstart send qframe sframe sstrand score bitscore evalue";

/// One row of `COLUMNS` against the record of its HSP.
fn check_columns(context: &str, hsp: &HspRecord, row: &[u8]) {
    let row = std::str::from_utf8(row).unwrap().trim_end();
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(fields.len(), 18, "{context}: {row}");
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
    // The reports cut these to a few digits (`report::tabular`).
    let bits: f64 = fields[16].parse().unwrap();
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
    let protein_query = matches!(program, Program::Blastp | Program::Tblastn);
    let protein_subject = matches!(program, Program::Blastp | Program::Blastx);
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

/// The command line's outfmt 5 report of `case` (LOSAT Web's BLASTX, which `run_collect` runs,
/// writes no XML).
fn xml_report(case: &Case, query: &Path, subject: &Path) -> String {
    String::from_utf8(cli(case, query, subject, "5").0).unwrap()
}

#[test]
fn the_statistics_are_those_of_the_xml_report() {
    for case in CASES.iter().chain(RANGE_CASES) {
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
    const PROTEIN: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
    const NUCLEOTIDE: &str =
        "ATGAAAACCGCGTATATTGCGAAACAGCGCCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGCCAGCTG\
                              GAAGAACGCCTGGGCCTGATTGAAGTGCAGGCGCCGATTCTGAGCCGCGTGGGCGATGGCACCCAGGAT\
                              AACCTGAGCGGCGCGGAAAAAGCGGTGCAGGTGAAAGTGAAAGCGCTGCCGGATGCGCAG";
    for program in Program::ALL {
        let inputs = |role: Role| -> Vec<Vec<FastaRecord>> {
            let protein = matches!(
                (program, role),
                (Program::Blastp, _)
                    | (Program::Blastx, Role::Subject)
                    | (Program::Tblastn, Role::Query)
            );
            let letters = if protein { PROTEIN } else { NUCLEOTIDE };
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

/// The example binary, built by `cargo test` (every target) or here by `cargo build` when a
/// single test target was asked for.
fn example_binary() -> PathBuf {
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_losat"))
        .parent()
        .unwrap()
        .to_path_buf();
    let example = bin_dir
        .join("examples")
        .join(format!("library_search{}", std::env::consts::EXE_SUFFIX));
    if !example.exists() {
        let status = Command::new(env!("CARGO"))
            .args(["build", "--example", "library_search", "--profile", "test"])
            .arg("--target-dir")
            .arg(bin_dir.parent().unwrap())
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .status()
            .expect("build the example");
        assert!(status.success(), "cargo build --example library_search");
    }
    example
}

#[test]
fn the_example_prints_the_command_lines_report_and_writes_only_its_warnings_to_stderr() {
    let example = example_binary();
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

#[test]
fn the_api_owns_the_inputs_and_the_outputs() {
    for owned in API_OWNED {
        let flag = format!("-{owned}");
        let error = validate(Program::Blastp, &[&flag, "x"]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{flag}: {error}");
        assert!(error.message().starts_with(&flag), "{flag}: {error}");
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
    assert!(error.message().contains("command line only"), "{error}");
}
