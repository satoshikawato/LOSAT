//! The library API (`LOSAT::api`, docs/losat_library_surface_plan.md §6): a search of records in
//! memory writes the bytes of the command line, one report per `-outfmt`, with the command line's
//! warnings; two searches at once give the reports of each alone; the example writes nothing to
//! standard error; the options the API owns are refused.
//!
//! The command line is the oracle here: its reports are the frozen NCBI bytes of the regression
//! fixtures (`fast output regressions`), so a report equal to the command line's is NCBI's.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use LOSAT::api::{
    describe, read_fasta, run_collect, validate, Error, Inputs, Outcome, Program, Role, API_OWNED,
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
        let spans: Vec<_> = outcome
            .hsps
            .iter()
            .map(|hsp| hsp.spans[tabular].clone().expect("a row in outfmt 6"))
            .collect();
        assert_eq!(spans, rows(report), "{context}: the rows of outfmt 6");
        for (index, hsp) in outcome.hsps.iter().enumerate() {
            assert_eq!(hsp.index, index);
            let row = &report[hsp.spans[tabular].clone().unwrap()];
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
                if case.options.contains(&"-num_alignments") {
                    assert!(
                        hsp.spans[pairwise].is_none(),
                        "{context}: another search's report"
                    );
                } else {
                    let bytes = &outcome.reports[pairwise].bytes;
                    match (&hsp.spans[pairwise], &hsp.subject_headings[pairwise]) {
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
                let sections = outcome
                    .hsps
                    .iter()
                    .filter(|hsp| hsp.spans[pairwise].is_some())
                    .count();
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
    // The warnings case has warnings, so the comparison above checked them.
    let case = &CASES[WARNINGS_CASE];
    let outcome = api(case, &fixture(case.query), &fixture(case.subject), &[]);
    assert!(
        String::from_utf8_lossy(&outcome.diagnostics).contains("Ignoring invalid residues"),
        "{}",
        String::from_utf8_lossy(&outcome.diagnostics)
    );
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
