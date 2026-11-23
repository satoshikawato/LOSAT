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
