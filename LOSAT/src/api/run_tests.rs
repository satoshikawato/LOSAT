use super::*;
use crate::api::{validate, Program};

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

fn one_protein() -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    (
        vec![FastaRecord::new("Query_1", b"q", PROTEIN)],
        vec![FastaRecord::new("Subject_1", b"s", PROTEIN)],
    )
}

#[test]
fn the_strand_is_ncbis_rule_and_the_ranges_are_half_open_on_the_plus_strand() {
    assert_eq!(strand(true, 9, 1), None);
    assert_eq!(strand(false, 1, 9), Some(Strand::Plus));
    assert_eq!(strand(false, 9, 1), Some(Strand::Minus));
    // A row of one letter: `(s_start - s_end) > 0` is false.
    assert_eq!(strand(false, 5, 5), Some(Strand::Plus));
    let (query, subject) = one_protein();
    let search = validate(Program::Blastp, &[]).unwrap();
    let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
    let mut hsp = outcome.hsps[0].clone();
    (hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end) = (3, 9, 9, 3);
    assert_eq!(hsp.query_range(), 2..9);
    assert_eq!(hsp.subject_range(), 2..9);
    (hsp.s_start, hsp.s_end) = (5, 5);
    assert_eq!(hsp.subject_range(), 4..5);
}

#[test]
fn groups_follow_the_hit_list_sizes() {
    let inputs = Inputs::new(&[], &[]);
    let plain = validate(Program::Blastp, &[])
        .unwrap()
        .command(inputs)
        .unwrap();
    assert_eq!(search_groups(&plain, &[6, 0, 7]), vec![vec![0, 1, 2]]);
    // -num_alignments changes the pairwise formats' hit list size only: two searches, the
    // tabular formats first.
    let split = validate(Program::Blastp, &["-num_alignments", "1"])
        .unwrap()
        .command(inputs)
        .unwrap();
    assert_eq!(search_groups(&split, &[0, 6, 7]), vec![vec![1, 2], vec![0]]);
    // -max_target_seqs sets every format's size.
    let same = validate(Program::Blastp, &["-max_target_seqs", "5"])
        .unwrap()
        .command(inputs)
        .unwrap();
    assert_eq!(search_groups(&same, &[0, 6]), vec![vec![0, 1]]);
}

#[test]
fn a_self_search_gives_one_hsp_with_its_rows_and_sections() {
    let (query, subject) = one_protein();
    let search = validate(Program::Blastp, &[]).unwrap();
    let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6", "0"]).unwrap();
    assert_eq!(outcome.reports.len(), 2);
    assert_eq!(outcome.reports[0].outfmt, "6");
    assert_eq!(outcome.reports[1].outfmt, "0");
    assert!(
        outcome.diagnostics.is_empty(),
        "{:?}",
        String::from_utf8_lossy(&outcome.diagnostics)
    );
    assert_eq!(outcome.hsps.len(), 1);
    let hsp = &outcome.hsps[0];
    assert_eq!((hsp.index, hsp.q_idx, hsp.s_idx, hsp.rank), (0, 0, 0, 0));
    assert_eq!(
        (hsp.q_start, hsp.q_end, hsp.s_start, hsp.s_end),
        (1, 66, 1, 66)
    );
    assert!(hsp.bit_score > 100.0);
    assert!(hsp.e_value < 1e-30);
    assert_eq!(
        hsp.query_aligned.as_deref(),
        Some(std::str::from_utf8(PROTEIN).unwrap())
    );
    // The row in outfmt 6 is the report's one line; the section in outfmt 0 starts at
    // the score line and the heading at the defline.
    let row = outcome.reports[0].spans[0].clone().unwrap();
    let tabular = &outcome.reports[0].bytes;
    assert_eq!(&tabular[row.clone()], tabular.as_slice());
    assert!(tabular.ends_with(b"\n"));
    assert_eq!(tabular.iter().filter(|&&b| b == b'\n').count(), 1);
    let section = outcome.reports[1].spans[0].clone().unwrap();
    let pairwise = &outcome.reports[1].bytes;
    assert!(pairwise[section].starts_with(b" Score ="));
    let heading = outcome.reports[1].subject_headings[0].clone().unwrap();
    assert!(pairwise[heading].starts_with(b">"));
    assert!(outcome.reports[0].subject_headings[0].is_none());
}

#[test]
fn two_searches_give_each_format_its_own_report_and_records_from_the_first() {
    let (query, subject) = one_protein();
    let search = validate(Program::Blastp, &["-num_alignments", "1"]).unwrap();
    let outcome = run_collect(&search, Inputs::new(&query, &subject), &["0", "6"]).unwrap();
    assert_eq!(outcome.reports[0].outfmt, "0");
    assert_eq!(outcome.reports[1].outfmt, "6");
    assert_eq!(outcome.hsps.len(), 1);
    // The records come from the tabular search: a row in report 1, no section in report
    // 0 (another search).
    assert!(outcome.reports[1].spans[0].is_some());
    assert_eq!(outcome.reports[0].spans, [None]);
    // The same bytes as the formats searched alone.
    let alone = run_collect(&search, Inputs::new(&query, &subject), &["0"]).unwrap();
    assert_eq!(alone.reports[0].bytes, outcome.reports[0].bytes);
}

#[test]
fn the_run_writes_to_the_callers_outputs_and_the_errors_are_the_command_lines() {
    let (query, subject) = one_protein();
    let search = validate(Program::Blastp, &[]).unwrap();
    let mut report = Vec::new();
    let mut diagnostics = Vec::new();
    run(
        &search,
        Inputs::new(&query, &subject),
        &mut [("6", &mut report)],
        &mut diagnostics,
    )
    .unwrap();
    let collected = run_collect(&search, Inputs::new(&query, &subject), &["6"]).unwrap();
    assert_eq!(report, collected.reports[0].bytes);
    assert!(diagnostics.is_empty());
    // A format the command line rejects fails before any search, as does a call without one.
    let error = run_collect(&search, Inputs::new(&query, &subject), &["99"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    let error = run_collect(&search, Inputs::new(&query, &subject), &[]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
    // An engine error is the command line's message.
    let empty: Vec<FastaRecord> = Vec::new();
    let error = run_collect(&search, Inputs::new(&query, &empty), &["6"]).unwrap_err();
    assert!(matches!(error, Error::Search(_)), "{error}");
    assert!(
        error.to_string().contains("Empty CBlastQueryVector"),
        "{error}"
    );
}
