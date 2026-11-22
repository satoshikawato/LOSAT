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

/// The variable that makes `html_posts_child` run (`html_posts_reach_the_diagnostics`), and the
/// formats of its calls (`0`, or `6,0`).
const HTML_POSTS_CHILD: &str = "LOSAT_TEST_HTML_POSTS_CHILD";

// Line RG T4a2 review B3: a `-html` report's reads of `./.ncbirc` post on the call's diagnostics
// (`Outcome::diagnostics`), not on the process's standard error, also where the report is the
// call's second format. A `./.ncbirc` with a repeated key posts at each read: NCBI BLAST+ 2.17.0
// reads it at each reference of the prolog, at the description table and at the alignments of
// each query with hits, and its id links' `InitConfig` reads it once (merged probes
// `T/t4/fix1/merged/{bp_nd0,bn_mb,bx,tn,tx}.ncbi`: blastp 2 references, the others 1). A later
// format's posts follow its program's rule for the lines of a later format's writer
// (`ReportOutputs::later_warnings`): the diagnostics, and in BLASTX the format's receiver, which
// `run_collect` does not give (so none). The file is read in the current directory, which a child
// process of this test binary gets (the tests do not change their own).
#[cfg(not(target_family = "wasm"))]
#[test]
fn html_posts_reach_the_diagnostics() {
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-html-posts-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(
        directory.join(".ncbirc"),
        b"[BLASTFMTUTIL]\nPROTOCOL=https:\nPROTOCOL=https:\n",
    )
    .unwrap();
    for formats in ["0", "6,0"] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "api::run::tests::html_posts_child",
                "--exact",
                "--nocapture",
                "--test-threads",
                "1",
            ])
            .env(HTML_POSTS_CHILD, formats)
            .env_remove("GETURL_DEBUG")
            .env_remove("FMTCFG")
            .env_remove("NCBI")
            .current_dir(&directory)
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            output.status.success(),
            "{stdout}\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(stdout.contains("html posts checked"), "{stdout}");
        // Nothing of the posts went to the process's standard error.
        assert!(
            !String::from_utf8_lossy(&output.stderr).contains("PROTOCOL"),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    std::fs::remove_dir_all(&directory).unwrap();
}

/// The child of `html_posts_reach_the_diagnostics`: without `HTML_POSTS_CHILD` it does nothing.
#[test]
fn html_posts_child() {
    let Some(formats) = std::env::var_os(HTML_POSTS_CHILD) else {
        return;
    };
    let formats = formats.into_string().unwrap();
    let formats: Vec<&str> = formats.split(',').collect();
    for (program, references) in [
        (Program::Blastn, 1),
        (Program::Blastp, 2),
        (Program::Blastx, 1),
        (Program::Tblastn, 1),
        (Program::Tblastx, 1),
    ] {
        let (query, subject) = self_pair(program);
        let search = validate(program, &["-html"]).unwrap();
        let outcome = run_collect(&search, Inputs::new(&query, &subject), &formats).unwrap();
        let post = format!(
            "Warning: (110.8) [{}] Found multiple [BLASTFMTUTIL] PROTOCOL settings; using the one from line 3\n",
            program.name().to_ascii_lowercase()
        );
        // The references, `InitConfig`, the query's description table and its alignments.
        let posts = if program == Program::Blastx && formats.len() > 1 {
            0
        } else {
            references + 3
        };
        assert_eq!(
            String::from_utf8_lossy(&outcome.diagnostics),
            post.repeat(posts),
            "{program} {formats:?}"
        );
    }
    println!("html posts checked");
}
