//! Line RG item T4a2: a `-html` report's reads of `./.ncbirc` post on the
//! call's diagnostics, in the five programs, through `run_collect` and `run` of `losat::api`.
//!
//! The file is read in the current directory, so the checks run in a child process of this test
//! binary that starts in a directory of its own with the file (the tests do not change the
//! directory of their own process, which the other tests share). Integration tests may start
//! processes; `src/` may not (`tests/check_pure_rust_runtime_boundary.py`).

#![cfg(not(target_family = "wasm"))]

use losat::api::{run, run_collect, validate, FastaRecord, Inputs, Program};

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROTEIN: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";

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

/// The variable that makes `html_posts_child` run (`html_posts_reach_the_diagnostics`).
const HTML_POSTS_CHILD: &str = "LOSAT_TEST_HTML_POSTS_CHILD";

// Line RG T4a2 review B3: a `-html` report's reads of `./.ncbirc` post on the call's diagnostics
// (`Outcome::diagnostics`), not on the process's standard error. A `./.ncbirc` with a repeated key
// posts at each read: NCBI BLAST+ 2.17.0 reads it at each reference of the prolog, at the
// description table and at the alignments of each query with hits, and its id links' `InitConfig`
// reads it once, with `GETURL_DEBUG`'s `REGISTRY:` line (merged probes
// `T/t4/fix1/merged/{bp_nd0,bn_mb,bx,tn,tx}.ncbi`: blastp 2 references, the others 1). Fix round 2
// (review c C3 = review d D5): one rule in the five programs, each format's report posts what a
// call of that format alone posts, at any place among the call's formats (BLASTX's later formats
// too), in each report of the same format (formats 0 and 0: `InitConfig`'s first call is each
// report's) and in a call split into two searches (formats 6 and 0 with
// `-num_alignments 5`, whose second search's other texts repeat the first's), through `run_collect`
// and `run`. The file is read in the current directory, which a child process of this test binary
// gets (the tests do not change their own).
#[test]
fn html_posts_reach_the_diagnostics() {
    let directory = std::env::temp_dir().join(format!("losat-html-posts-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    std::fs::write(
        directory.join(".ncbirc"),
        b"[BLASTFMTUTIL]\nPROTOCOL=https:\nPROTOCOL=https:\n",
    )
    .unwrap();
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "library_html_posts::html_posts_child",
            "--exact",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(HTML_POSTS_CHILD, "1")
        .env("GETURL_DEBUG", "1")
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
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("PROTOCOL") && !stderr.contains("REGISTRY"),
        "{stderr}"
    );
    std::fs::remove_dir_all(&directory).unwrap();
}

/// The child of `html_posts_reach_the_diagnostics`: without `HTML_POSTS_CHILD` it does nothing.
#[test]
fn html_posts_child() {
    if std::env::var_os(HTML_POSTS_CHILD).is_none() {
        return;
    }
    for (program, references) in [
        (Program::Blastn, 1),
        (Program::Blastp, 2),
        (Program::Blastx, 1),
        (Program::Tblastn, 1),
        (Program::Tblastx, 1),
    ] {
        let (query, subject) = self_pair(program);
        let inputs = Inputs::new(&query, &subject);
        let call = |options: &[&str], formats: &[&str]| {
            let search = validate(program, options).unwrap();
            run_collect(&search, inputs, formats).unwrap().diagnostics
        };
        let post = format!(
            "Warning: (110.8) [{}] Found multiple [BLASTFMTUTIL] PROTOCOL settings; using the one from line 3\n",
            program.name().to_ascii_lowercase()
        );
        // The references, the query's description table, `InitConfig` at its first id link (the
        // post, then the line), and its alignments.
        let report = format!("{}REGISTRY: .ncbirc\n{post}", post.repeat(references + 2));
        let html = ["-html"];
        let zero = call(&html, &["0"]);
        assert_eq!(String::from_utf8_lossy(&zero), report, "{program} [0]");
        // One search: each format's report posts what a call of its format alone posts. The
        // reports interleave (each format's prolog, then each batch's reports in the call's
        // order), and so do their posts: with one pairwise format, the call's posts are its own;
        // with two, the same lines as the two calls.
        let lines = |text: &[u8]| {
            let mut lines: Vec<String> = String::from_utf8_lossy(text)
                .lines()
                .map(str::to_string)
                .collect();
            lines.sort();
            lines
        };
        for formats in [&["6", "0"][..], &["0", "6"], &["0", "0"], &["18", "0"]] {
            let alone: Vec<u8> = formats
                .iter()
                .flat_map(|&outfmt| call(&html, &[outfmt]))
                .collect();
            let both = call(&html, formats);
            if formats.contains(&"6") {
                assert_eq!(
                    String::from_utf8_lossy(&both),
                    String::from_utf8_lossy(&alone),
                    "{program} {formats:?}"
                );
            }
            assert_eq!(lines(&both), lines(&alone), "{program} {formats:?}");
        }
        // The organism report after outfmt 0 takes `-html` alone (no id links).
        assert!(
            call(&html, &["0", "18"]).starts_with(&zero),
            "{program} [0, 18]"
        );
        // Two searches: the first's diagnostics, then the second's report's posts.
        // (The call splits into two searches: `search_groups` is not public, so the expected
        // value below, the split call's diagnostics = format 6's call + format 0's report, is what
        // shows the split.)
        let split = ["-html", "-num_alignments", "5"];
        let search = validate(program, &split).unwrap();
        let expected = [call(&split, &["6"]), report.clone().into_bytes()].concat();
        assert_eq!(
            String::from_utf8_lossy(&call(&split, &["6", "0"])),
            String::from_utf8_lossy(&expected),
            "{program} [6, 0] split"
        );
        let (mut tabular, mut pairwise, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        run(
            &search,
            inputs,
            &mut [("6", &mut tabular), ("0", &mut pairwise)],
            &mut diagnostics,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8_lossy(&diagnostics),
            String::from_utf8_lossy(&expected),
            "{program} run [6, 0] split"
        );
    }
    println!("html posts checked");
}

/// The variable that makes `later_formats_child` run (`later_formats_post_as_their_runs_alone`).
const LATER_FORMATS_CHILD: &str = "LOSAT_TEST_HTML_LATER_FORMATS_CHILD";

// Line RG T4b review b N3: in a call of several formats, a `-html` report of outfmt 1 to 4 at any
// place gets the posts of NCBI's formatter registries (`report::query_warnings::report_posts`,
// `ReportFormat::writes_html_links`). Here `FMTCFG` names a registry file whose read posts (a
// repeated key) and `GETURL_DEBUG` is set, the run's directory has no `./.ncbirc`: each format's
// call alone posts the command line's standard error of that format alone, each report's bytes
// in the call are those of the call of that format alone, and the call's posts are the same lines
// as those calls' posts (in their order where the call has one pairwise format), in the five
// programs.
// The child process gets the environment (the tests do not change their own).
#[test]
fn later_formats_post_as_their_runs_alone() {
    let directory = std::env::temp_dir().join(format!("losat-html-later-{}", std::process::id()));
    let run_directory = directory.join("run");
    std::fs::create_dir_all(&run_directory).unwrap();
    let registry = directory.join("fmt.rc");
    std::fs::write(
        &registry,
        b"[BLASTFMTUTIL]\nPROTOCOL=https:\nPROTOCOL=https:\n",
    )
    .unwrap();
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "library_html_posts::later_formats_child",
            "--exact",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(LATER_FORMATS_CHILD, "1")
        .env("FMTCFG", &registry)
        .env("GETURL_DEBUG", "1")
        .env_remove("NCBI")
        .current_dir(&run_directory)
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "{stdout}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(stdout.contains("later formats checked"), "{stdout}");
    std::fs::remove_dir_all(&directory).unwrap();
}

/// The child of `later_formats_post_as_their_runs_alone`: without `LATER_FORMATS_CHILD` it does
/// nothing.
#[test]
fn later_formats_child() {
    if std::env::var_os(LATER_FORMATS_CHILD).is_none() {
        return;
    }
    // The records as FASTA (`>` and the title; the local ids are made, `Query_1`, `Subject_1`).
    let fasta = |records: &[FastaRecord]| -> Vec<u8> {
        records
            .iter()
            .flat_map(|record| {
                [
                    b">".as_slice(),
                    &record.title,
                    b"\n",
                    &record.sequence,
                    b"\n",
                ]
                .concat()
            })
            .collect()
    };
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ] {
        let (query, subject) = self_pair(program);
        std::fs::write("query.fa", fasta(&query)).unwrap();
        std::fs::write("subject.fa", fasta(&subject)).unwrap();
        // The command line's standard error of one format alone.
        let command_line = |outfmt: &str| {
            let output = std::process::Command::new(env!("CARGO_BIN_EXE_losat"))
                .args([
                    &program.name().to_ascii_lowercase(),
                    "-query",
                    "query.fa",
                    "-subject",
                    "subject.fa",
                    "-html",
                    "-outfmt",
                    outfmt,
                ])
                .output()
                .unwrap();
            assert!(output.status.success(), "{program} outfmt {outfmt}");
            output.stderr
        };
        let search = validate(program, &["-html"]).unwrap();
        let inputs = Inputs::new(&query, &subject);
        let alone = |outfmt: &str| {
            let mut outcome = run_collect(&search, inputs, &[outfmt]).unwrap();
            (outcome.reports.remove(0).bytes, outcome.diagnostics)
        };
        let (one, one_posts) = alone("1");
        assert!(one.starts_with(b"<HTML>"), "{program}");
        assert!(
            String::from_utf8_lossy(&one_posts).contains("PROTOCOL settings"),
            "{program}: {}",
            String::from_utf8_lossy(&one_posts)
        );
        let lines = |text: &[u8]| {
            let mut lines: Vec<String> = String::from_utf8_lossy(text)
                .lines()
                .map(str::to_string)
                .collect();
            lines.sort();
            lines
        };
        for formats in [
            &["6", "1"][..],
            &["0", "1"],
            &["1", "18"],
            &["3", "0", "4", "2"],
        ] {
            let together = run_collect(&search, inputs, formats).unwrap();
            let mut expected = Vec::new();
            for (at, &format) in formats.iter().enumerate() {
                let (bytes, posts) = alone(format);
                assert_eq!(
                    String::from_utf8_lossy(&posts),
                    String::from_utf8_lossy(&command_line(format)),
                    "{program}: outfmt {format} alone"
                );
                assert_eq!(
                    String::from_utf8_lossy(&together.reports[at].bytes),
                    String::from_utf8_lossy(&bytes),
                    "{program} {formats:?}: outfmt {format}"
                );
                expected.extend(posts);
            }
            // The reports interleave, and so do their posts: with one pairwise format the call's
            // posts are its own; with several, the same lines as their calls alone.
            if formats
                .iter()
                .filter(|&&format| format != "6" && format != "18")
                .count()
                == 1
            {
                assert_eq!(
                    String::from_utf8_lossy(&together.diagnostics),
                    String::from_utf8_lossy(&expected),
                    "{program} {formats:?}"
                );
            }
            assert_eq!(
                lines(&together.diagnostics),
                lines(&expected),
                "{program} {formats:?}"
            );
        }
    }
    println!("later formats checked");
}
