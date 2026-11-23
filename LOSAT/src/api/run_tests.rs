use super::*;
use crate::api::{validate, Program};
use crate::blastinput::ncbi_environment::DiagTarget;

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
        .command(inputs, &DiagTarget::Stderr)
        .unwrap();
    assert_eq!(search_groups(&plain, &[6, 0, 7]), vec![vec![0, 1, 2]]);
    // -num_alignments changes the pairwise formats' hit list size only: two searches, the
    // tabular formats first.
    let split = validate(Program::Blastp, &["-num_alignments", "1"])
        .unwrap()
        .command(inputs, &DiagTarget::Stderr)
        .unwrap();
    assert_eq!(search_groups(&split, &[0, 6, 7]), vec![vec![1, 2], vec![0]]);
    // -max_target_seqs sets every format's size.
    let same = validate(Program::Blastp, &["-max_target_seqs", "5"])
        .unwrap()
        .command(inputs, &DiagTarget::Stderr)
        .unwrap();
    assert_eq!(search_groups(&same, &[0, 6]), vec![vec![0, 1]]);
}

// Review of PS15b-1, S1: psiblast's XML cuts each round's results before the round's PSSM, and
// psiblast's archive keeps the uncut lists (`BlastpApp::list_cut`), so these formats take a
// search of their own, the first requested format's first.
#[test]
fn groups_follow_the_cuts_of_the_hit_lists() {
    let inputs = Inputs::new(&[], &[]);
    let command =
        |program, options: &[&str]| validate(program, options).unwrap().command(inputs).unwrap();
    let psiblast = command(Program::Psiblast, &[]);
    assert_eq!(search_groups(&psiblast, &[6, 11]), vec![vec![0], vec![1]]);
    assert_eq!(search_groups(&psiblast, &[11, 6]), vec![vec![0], vec![1]]);
    assert_eq!(search_groups(&psiblast, &[8, 5]), vec![vec![0], vec![1]]);
    assert_eq!(
        search_groups(&psiblast, &[6, 5, 7, 11, 8, 9, 10]),
        vec![vec![0, 2, 4, 5, 6], vec![1], vec![3]]
    );
    // With the hit list sizes: the tabular formats first, then by cut.
    let split = command(Program::Psiblast, &["-num_alignments", "1"]);
    assert_eq!(
        search_groups(&split, &[0, 5, 6, 11]),
        vec![vec![1], vec![2], vec![3], vec![0]]
    );
    // BLASTP's lists are not cut: one search.
    let blastp = command(Program::Blastp, &[]);
    assert_eq!(search_groups(&blastp, &[6, 11, 5]), vec![vec![0, 1, 2]]);
}

// Review of PS15b-1, S1: each report of a call of psiblast formats of other cuts is the command
// line's for its format alone, on the case of the rows `db.ps15b.pk.mts2.{5,8,9}` (NCBI BLAST+
// 2.17.0: outfmt 5 stops after round 2, where 8 has a round 3). The paths are the rows', from the
// package's directory (the tests' working directory).
#[test]
fn psiblast_formats_of_other_cuts_are_each_the_command_lines() {
    let options = [
        "-num_iterations",
        "3",
        "-inclusion_ethresh",
        "5",
        "-max_target_seqs",
        "2",
    ];
    let search = validate(Program::Psiblast, &options).unwrap();
    let query_name = "tests/fasta/psiblast/ps5_pk.faa";
    let bytes = std::fs::read(query_name).unwrap();
    let (records, _) = crate::api::read_inputs(&search, &bytes, b"").unwrap();
    let inputs = Inputs::database(&records, "tests/db_search/paja").named(query_name, "ignored");
    let alone = |outfmt: &str| {
        run_collect(&search, inputs, &[outfmt])
            .unwrap()
            .reports
            .remove(0)
            .bytes
    };
    // The command line's bytes: NCBI's frozen rows, with LOSAT's version (`losat_output_bytes`).
    for outfmt in ["5", "8"] {
        let ncbi = std::fs::read(format!(
            "tests/fixtures/db_search/db.ps15b.pk.mts2.{outfmt}.out"
        ))
        .unwrap();
        assert!(
            alone(outfmt) == crate::wording::losat_output_bytes(&ncbi).as_ref(),
            "-outfmt {outfmt} alone"
        );
    }
    for formats in [["6", "11"], ["11", "6"], ["8", "5"]] {
        let outcome = run_collect(&search, inputs, &formats).unwrap();
        for (report, outfmt) in outcome.reports.iter().zip(formats) {
            assert!(
                report.bytes == alone(outfmt),
                "{formats:?}: -outfmt {outfmt}"
            );
        }
    }
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

/// `PROTEIN` back-translated, one codon per residue: a nucleotide sequence for the programs whose
/// query or subjects are nucleotides.
fn coding_sequence() -> Vec<u8> {
    PROTEIN
        .iter()
        .flat_map(|residue| -> &[u8] {
            match residue {
                b'A' => b"GCT",
                b'R' => b"CGT",
                b'N' => b"AAT",
                b'D' => b"GAT",
                b'Q' => b"CAA",
                b'E' => b"GAA",
                b'G' => b"GGT",
                b'H' => b"CAT",
                b'I' => b"ATT",
                b'L' => b"CTG",
                b'K' => b"AAA",
                b'M' => b"ATG",
                b'F' => b"TTT",
                b'P' => b"CCG",
                b'S' => b"TCT",
                b'T' => b"ACC",
                b'V' => b"GTT",
                b'Y' => b"TAT",
                other => panic!("no codon for {}", char::from(*other)),
            }
            .as_slice()
        })
        .copied()
        .collect()
}

/// A query and a subject of `program`, each the AsnA part (`PROTEIN`) or its coding sequence.
fn self_pair(program: Program) -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    let dna = coding_sequence();
    let (query, subject): (&[u8], &[u8]) = match program {
        Program::Blastn | Program::Tblastx => (&dna, &dna),
        Program::Blastx => (&dna, PROTEIN),
        Program::Tblastn => (PROTEIN, &dna),
        _ => (PROTEIN, PROTEIN),
    };
    (
        vec![FastaRecord::new("Query_1", b"q", query)],
        vec![FastaRecord::new("Subject_1", b"s", subject)],
    )
}

// Line RG item T4a2: only the pairwise report reads `-html`'s links (`blastinput::app::html_links`),
// so a call of outfmt 6 then 0 has them for its second format: that report is the HTML report of
// outfmt 0 alone, in each of the five programs.
#[test]
fn the_html_links_reach_outfmt_0_in_any_place() {
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ] {
        let (query, subject) = self_pair(program);
        let search = validate(program, &["-html"]).unwrap();
        let alone = run_collect(&search, Inputs::new(&query, &subject), &["0"]).unwrap();
        assert!(
            alone.reports[0].bytes.starts_with(b"<HTML>"),
            "{program}: {}",
            String::from_utf8_lossy(&alone.reports[0].bytes)
        );
        let second = run_collect(&search, Inputs::new(&query, &subject), &["6", "0"]).unwrap();
        assert_eq!(
            String::from_utf8_lossy(&second.reports[1].bytes),
            String::from_utf8_lossy(&alone.reports[0].bytes),
            "{program}"
        );
    }
}
