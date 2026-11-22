//! The receivers of the later formats' query warnings (`ReportOutputs::later_warnings`) in
//! process, for each program: what a format after a structured first one gives where it is
//! written, where it failed alone before its reports, and where its report stops.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1452 (cited at
//! `FormatOutput::writes_query_warnings`): a structured report returns before the query's result
//! warnings, and the others post them when the query's report starts.

use super::support::{fixture_sequence, TempFasta};
use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    FormatOutput, FormatRefusal, OutputSink, ReportOutputs,
};
use crate::api::Program;
use crate::cli::{try_parse_from, Cli, Commands};

const PROGRAMS: [&str; 5] = ["blastn", "blastp", "blastx", "tblastn", "tblastx"];

/// One search of a program: its inputs and its options.
struct Search {
    program: &'static str,
    query: TempFasta,
    subject: TempFasta,
    extra: &'static [&'static str],
}

/// What a run of several formats gave: the diagnostics, each format's receiver of the later
/// warnings, and whether each format failed alone.
#[derive(Debug, PartialEq)]
struct Given {
    diagnostics: String,
    later: Vec<String>,
    refused: Vec<bool>,
}

impl Search {
    /// Runs `formats` as LOSAT Web runs the search of the argv's format: the first format fails
    /// the call and the others fail alone. `receives`: the caller takes the later formats'
    /// warnings.
    fn run(&self, formats: &[&str], receives: bool) -> Given {
        self.run_as(formats, receives, FormatRefusal::FailsCall)
    }

    /// Runs `formats`, the first with `first` and the others failing alone (LOSAT Web's other
    /// search of a run of two: every format fails alone).
    fn run_as(&self, formats: &[&str], receives: bool, first: FormatRefusal) -> Given {
        let (query, subject) = (
            self.query.0.display().to_string(),
            self.subject.0.display().to_string(),
        );
        let command = || {
            let mut argv = vec![
                "losat",
                self.program,
                "-query",
                &query,
                "-subject",
                &subject,
            ];
            argv.extend(self.extra);
            let cli: Cli = try_parse_from(argv).expect("argv");
            cli.command
        };
        let program = match self.program {
            "blastn" => Program::Blastn,
            "blastp" => Program::Blastp,
            "blastx" => Program::Blastx,
            "tblastn" => Program::Tblastn,
            _ => Program::Tblastx,
        };
        let crate::api::InputRecords {
            query: queries,
            query_end,
            subject: subjects,
        } = crate::api::input::read_search_inputs(
            program,
            &command(),
            Some(false),
            true,
            &std::fs::read(&self.query.0).expect("query"),
            &std::fs::read(&self.subject.0).expect("subject"),
        )
        .map_err(|error| error.to_string())
        .expect("inputs");
        let mut sinks = vec![Vec::new(); formats.len()];
        let mut diagnostics = Vec::new();
        let mut later = vec![Vec::new(); formats.len()];
        let mut first = Some(first);
        let mut outputs = ReportOutputs {
            other_outfmts: &[],
            formats: formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(outfmt, sink)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(sink),
                    refusal: first.take().unwrap_or(FormatRefusal::FailsAlone),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            later_warnings: receives.then_some(&mut later[..]),
            hits: None,
            observer: None,
        };
        let result = match command() {
            Commands::Blastn(args) => {
                run_local_blastn(args, &queries, &query_end, &subjects, &mut outputs)
            }
            Commands::Blastp(args) => {
                run_local_blastp(args, &queries, &query_end, &subjects, "", "", &mut outputs)
            }
            Commands::Blastx(args) => {
                run_local_blastx(args, &queries, &query_end, &subjects, &mut outputs)
            }
            Commands::Tblastn(args) => {
                run_local_tblastn(args, &queries, &query_end, &subjects, &mut outputs)
            }
            Commands::Tblastx(args) => {
                run_local_tblastx(args, &queries, &query_end, &subjects, &mut outputs)
            }
            _ => unreachable!("a search program"),
        };
        result.expect("run_local");
        let refused = outputs
            .formats
            .iter_mut()
            .map(|format| format.refusal.take_refused().is_some())
            .collect();
        drop(outputs);
        let text = |bytes: Vec<u8>| String::from_utf8(bytes).expect("UTF-8");
        Given {
            diagnostics: text(diagnostics),
            later: later.into_iter().map(text).collect(),
            refused,
        }
    }
}

/// A nucleotide sequence that codes for `protein`, one codon per residue.
fn coding(protein: &[u8]) -> Vec<u8> {
    protein
        .iter()
        .flat_map(|residue| {
            let codon: &[u8; 3] = match residue {
                b'A' => b"GCT",
                b'C' => b"TGT",
                b'D' => b"GAT",
                b'E' => b"GAA",
                b'F' => b"TTT",
                b'G' => b"GGT",
                b'H' => b"CAT",
                b'I' => b"ATT",
                b'K' => b"AAA",
                b'L' => b"CTT",
                b'M' => b"ATG",
                b'N' => b"AAT",
                b'P' => b"CCT",
                b'Q' => b"CAA",
                b'R' => b"CGT",
                b'S' => b"TCT",
                b'T' => b"ACT",
                b'V' => b"GTT",
                b'W' => b"TGG",
                b'Y' => b"TAT",
                other => panic!("no codon for {}", char::from(*other)),
            };
            *codon
        })
        .collect()
}

/// A search of `program` whose first query record has no residues and whose second query hits
/// the subject. `subject_id` names the subject; `masked` puts the first 20 letters of the second
/// query and of the subject in lower case.
fn search(
    program: &'static str,
    subject_id: &str,
    masked: bool,
    extra: &'static [&'static str],
) -> Search {
    search_of(program, subject_id, masked, extra, false)
}

/// `search`, with a third query record, `qlast`, which has no residues either, where
/// `last_empty`.
fn search_of(
    program: &'static str,
    subject_id: &str,
    masked: bool,
    extra: &'static [&'static str],
    last_empty: bool,
) -> Search {
    let protein = fixture_sequence("SicyWSV.faa")[..100].to_vec();
    let nucleotide = fixture_sequence("LC738884.fasta")[..600].to_vec();
    let (mut query, mut subject) = match program {
        "blastn" | "tblastx" => (nucleotide.clone(), nucleotide),
        "blastp" => (protein.clone(), protein),
        "blastx" => (coding(&protein), protein),
        _ => (protein.clone(), coding(&protein)),
    };
    if masked {
        query[..20].make_ascii_lowercase();
        subject[..20].make_ascii_lowercase();
    }
    Search {
        program,
        query: TempFasta::new(
            &format!("later_{program}_query{}.fa", u8::from(last_empty)),
            &[("qempty", &[][..]), ("q1", &query), ("qlast", &[][..])]
                [..2 + usize::from(last_empty)],
        ),
        subject: TempFasta::new(
            &format!("later_{program}_subject.fa"),
            &[(subject_id, &subject)],
        ),
        extra,
    }
}

// A structured first format (outfmt 5) writes no result warning of a query, as the command line's
// run of it. Each later format that writes them (outfmt 6, 7, 0) gives its own to its own
// receiver, as its command line's run writes them; a later structured format gives none, and
// nothing goes anywhere without receivers. A first format that writes them writes them to the
// diagnostics, and the later formats give none (line WQA's sweep, D3 review r3, 2026-11-19).
#[test]
fn each_later_format_gives_its_query_warnings_to_its_own_receiver() {
    for program in PROGRAMS {
        let search = search(program, "s1", false, &[]);
        let warning = format!("Warning: [{program}] Query_1 qempty: Sequence contains no data \n");
        let alone = search.run(&["5"], false);
        assert!(!alone.diagnostics.contains(&warning), "{program}");
        let given = search.run(&["5", "6", "15", "7", "0"], true);
        assert_eq!(given.diagnostics, alone.diagnostics, "{program}");
        assert_eq!(given.refused, [false; 5], "{program}");
        assert_eq!(
            given.later,
            ["", &warning, "", &warning, &warning],
            "{program}"
        );
        let dropped = search.run(&["5", "6", "15", "7", "0"], false);
        assert_eq!(dropped.diagnostics, alone.diagnostics, "{program}");
        let tabular = search.run(&["6"], false);
        assert!(tabular.diagnostics.ends_with(&warning), "{program}");
        let first = search.run(&["6", "5", "0"], true);
        assert_eq!(first.diagnostics, tabular.diagnostics, "{program}");
        assert_eq!(first.later, ["", "", ""], "{program}");
    }
}

// A receiver holds its format's warnings in the queries' order, each once: those of the first and
// of the last query, which have no residues, around the report of the second (review r2 of this
// change, L3, 2026-11-20).
#[test]
fn a_receiver_holds_the_warnings_of_several_queries_in_their_order() {
    for program in PROGRAMS {
        let search = search_of(program, "s1", false, &[], true);
        let warnings = format!(
            "Warning: [{program}] Query_1 qempty: Sequence contains no data \n\
             Warning: [{program}] Query_3 qlast: Sequence contains no data \n"
        );
        let tabular = search.run(&["6"], false);
        assert!(tabular.diagnostics.ends_with(&warnings), "{program}");
        let given = search.run(&["5", "0", "7"], true);
        assert_eq!(given.later, ["", &warnings, &warnings], "{program}");
    }
}

// A later format that failed alone before its reports gives none: the command line's run of it
// stops at the refusal, before the first query's report and its warning. Outfmt 0 is refused for
// the masks of a query that is also a subject (`-parse_deflines -lcase_masking`), or for a
// description id longer than 65 bytes (`report::db_subject::check_pairwise_records`). The formats
// after it give theirs (review r1 of this change, M1, 2026-11-20).
#[test]
fn a_later_format_refused_before_its_reports_gives_no_query_warnings() {
    const LONG_ID: &str =
        "gnl|jjjjjjjjjjjjjjjjjjjj|zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz";
    for program in PROGRAMS {
        let search = match program {
            "blastn" | "blastp" | "tblastx" => {
                search(program, "q1", true, &["-parse_deflines", "-lcase_masking"])
            }
            _ => search(program, LONG_ID, false, &["-parse_deflines"]),
        };
        let given = search.run(&["5", "0"], true);
        assert_eq!(given.refused, [false, true], "{program}");
        assert_eq!(given.later, ["", ""], "{program}");
        let given = search.run(&["5", "0", "6", "7"], true);
        assert_eq!(given.refused, [false, true, false, false], "{program}");
        let warning = &given.later[2];
        assert!(
            warning.starts_with(&format!("Warning: [{program}] "))
                && warning.ends_with(": Sequence contains no data \n")
                && warning.lines().count() == 1,
            "{program}: {warning:?}"
        );
        assert_eq!(given.later, ["", "", warning, warning], "{program}");
    }
}

// A later format whose report stops after a query's warning gave that warning, once, and the
// formats after it give their own: each receiver holds what the command line's run of its format
// wrote before it stopped. BLASTX with HSPs of length 0 (`-ungapped -window_size 0`): NCBI's
// outfmt 6, 7 and 0 throw at the second query's alignments, after the warning of the first, and
// outfmt 5 goes on (review r1 of this change, M2, 2026-11-20).
#[test]
fn a_later_format_whose_report_stops_gave_the_warnings_before_it_once() {
    let search = stopping_blastx("later_stops");
    let warning = "Warning: [blastx] Query_1 qempty: Sequence contains no data \n";
    let given = search.run(&["5", "6", "7", "0"], true);
    assert_eq!(given.refused, [false, true, true, true]);
    assert!(!given.diagnostics.contains(warning));
    assert_eq!(given.later, ["", warning, warning, warning]);
}

// Formats that share the diagnostics (a first format that writes the queries' warnings, or a
// caller without receivers): the first format that is written writes the warnings there, and one
// that stopped after a query's warning wrote it, so the formats after it do not write it again.
// LOSAT Web's other search of a run of two is such a call (outfmt 6 and 7, each failing alone).
// BLASTX's report writer makes the warnings itself and read the batch as not warned off the
// refusal: the format after the one that stopped wrote the warning again (line WQA, 2026-11-20).
#[test]
fn formats_that_share_the_diagnostics_write_a_warning_once_where_one_stops_after_it() {
    let search = stopping_blastx("shared_stops");
    let warning = "Warning: [blastx] Query_1 qempty: Sequence contains no data \n";
    for receives in [false, true] {
        let given = search.run_as(&["6", "7", "0"], receives, FormatRefusal::FailsAlone);
        assert_eq!(given.refused, [true, true, true]);
        assert_eq!(given.diagnostics.matches(warning).count(), 1, "{given:?}");
        assert_eq!(given.later, ["", "", ""]);
    }
}

/// A BLASTX search whose first query record has no residues and whose second query gives an HSP
/// of length 0 (`tests/fasta/outfmt0/x2g_c0054`, line RG's item X2g): outfmt 6, 7 and 0 stop at
/// it, after the warning of the first query, and outfmt 5 goes on.
fn stopping_blastx(name: &str) -> Search {
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/outfmt0");
    let query = super::support::read_records(&fixture.join("x2g_c0054_q.fna")).remove(0);
    let subject = super::support::read_records(&fixture.join("x2g_c0054_s.faa")).remove(0);
    Search {
        program: "blastx",
        query: TempFasta::new(
            &format!("{name}_query.fna"),
            &[("qempty", &[][..]), ("q1", query.seq())],
        ),
        subject: TempFasta::new(&format!("{name}_subject.faa"), &[("s1", subject.seq())]),
        extra: &[
            "-ungapped",
            "-comp_based_stats",
            "0",
            "-window_size",
            "0",
            "-word_size",
            "5",
            "-threshold",
            "1",
            "-evalue",
            "1e5",
            "-seg",
            "no",
        ],
    }
}
