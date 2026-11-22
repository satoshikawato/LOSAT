//! blast_formatter's command line and `Run` up to the archive (`run_with`). The texts are NCBI
//! BLAST+ 2.17.0's (line RG item F3, probes P2.2xx and P3.122; the frozen rows are in
//! `LOSAT/tests/blast_formatter/`).

use super::*;

fn run_words(line: &[&str], input: &[u8]) -> (i32, String, String) {
    let words: Vec<OsString> = line.iter().map(OsString::from).collect();
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status = run_with(&words, &mut &input[..], &mut out, &mut err);
    (
        status,
        String::from_utf8_lossy(&out).into_owned(),
        String::from_utf8(err).unwrap(),
    )
}

#[test]
fn the_short_usage_is_the_head_of_the_help() {
    let usage = crate::ncbi_app::short_usage(HELP);
    // NCBI's 482 bytes, its version replaced by LOSAT's.
    assert_eq!(usage.len(), 482 + crate::wording::IDENT_GROWTH);
    assert!(usage.starts_with("USAGE\n  blast_formatter [-h] [-help] [-rid BLAST_RID]"));
    assert_eq!(HELP.len(), 6306);
}

// NCBI's multimap of the dependencies (`ncbi_app::args::dependencies`) is the order of NCBI's
// `-xmlhelp` (probe P2.226), which decides which conflict is reported (P2.085, P2.220).
#[test]
fn the_dependencies_are_in_ncbi_s_order() {
    let deps: Vec<(&str, &str)> = crate::ncbi_app::args::dependencies(TOOL.program)
        .into_iter()
        .map(|(first, _, second)| (first, second))
        .collect();
    assert_eq!(
        deps,
        [
            ("archive", "rid"),
            ("max_target_seqs", "num_descriptions"),
            ("max_target_seqs", "num_alignments"),
            ("num_alignments", "max_target_seqs"),
            ("num_descriptions", "max_target_seqs"),
            ("rid", "archive"),
        ]
    );
}

#[test]
fn the_arguments_are_read_into_the_options() {
    let words: Vec<String> = [
        "-archive",
        "a",
        "-outfmt",
        "6 qseqid",
        "-show_gis",
        "-num_alignments",
        "0",
        "-num_descriptions",
        "3",
        "-line_length",
        "40",
        "-html",
        "-sorthits",
        "4",
        "-sorthsps",
        "0",
        "-out",
        "o.txt",
        "-parse_deflines",
    ]
    .iter()
    .map(|word| word.to_string())
    .collect();
    let args = crate::ncbi_app::args::parse(&TOOL, &words).unwrap();
    assert_eq!(
        Options::of(&args),
        Options {
            outfmt: "6 qseqid".to_string(),
            show_gis: true,
            num_descriptions: Some(3),
            num_alignments: Some(0),
            line_length: Some(40),
            html: true,
            sort_hits: Some(4),
            sort_hsps: Some(0),
            max_target_seqs: None,
            out: "o.txt".to_string(),
            parse_deflines: true,
        }
    );
    let args =
        crate::ncbi_app::args::parse(&TOOL, &["-max_target_seqs".to_string(), "7".to_string()])
            .unwrap();
    let options = Options::of(&args);
    assert_eq!((options.outfmt.as_str(), options.out.as_str()), ("0", "-"));
    assert_eq!(options.max_target_seqs, Some(7));
    assert!(!options.show_gis && !options.html && !options.parse_deflines);
}

// NCBI BLAST+ 2.17.0 `blast_formatter` alone (P2.205): `Run` asks for `-rid`'s value.
#[test]
fn without_archive_or_rid_the_rid_has_no_value() {
    assert_eq!(
        run_words(&[], b""),
        (
            1,
            String::new(),
            "Command line argument error: Argument \"rid\". The argument has no value\n"
                .to_string()
        )
    );
}

// P2.221 (`-rid ''`): `CRemoteBlast` refuses the empty RID before it asks NCBI's servers.
#[test]
fn an_empty_rid_is_ncbi_s_engine_error() {
    assert_eq!(
        run_words(&["-rid", ""], b""),
        (
            3,
            String::new(),
            "LOSAT engine error: Empty RID string specified\n".to_string()
        )
    );
}

// Where NCBI asks its servers for the search (P2.219, P2.222), LOSAT refuses `-rid` as the
// search programs refuse `-remote`.
#[test]
fn a_rid_is_refused_as_remote_is() {
    let refusal =
        "error: the NCBI BLAST+ option -rid is not supported by LOSAT's BLAST_FORMATTER\n";
    for line in [
        &["-rid", "ABCDEFGH01R"][..],
        &["-rid", "ABCDEFGH01R", "-outfmt", "6"],
        &["-rid", "X", "-outfmt", "abc", "-out", "/no/such/dir/o"],
    ] {
        assert_eq!(
            run_words(line, b""),
            (2, String::new(), refusal.to_string()),
            "{line:?}"
        );
    }
    // The arguments are checked first (NCBI's bytes: row `rid_archive` of
    // `LOSAT/tests/blast_formatter`), and `-dryrun` ends the run before `Run`.
    let (status, _, err) = run_words(&["-rid", "X", "-archive", "a"], b"");
    assert_eq!(status, 1);
    assert!(
        err.contains("\nError: Argument \"rid\". Incompatible with argument:  `archive'\n"),
        "{err}"
    );
    assert_eq!(
        run_words(&["-rid", "X", "-dryrun"], b""),
        (0, String::new(), String::new())
    );
}

// P3.122 (a missing file) and P3.124 (an unreadable one): `AsInputFile()` throws in `Run`, which
// `CATCH_ALL` reports.
#[test]
fn an_archive_that_does_not_open_is_ncbi_s_argument_error() {
    for path in ["/no/such/dir/archive.asn", ""] {
        let attr = if path.is_empty() {
            String::new()
        } else {
            format!(":  `{path}'")
        };
        assert_eq!(
            run_words(&["-archive", path, "-outfmt", "6"], b""),
            (
                1,
                String::new(),
                format!(
                    "Command line argument error: Argument \"archive\". File is not accessible{attr}\n"
                )
            ),
            "{path:?}"
        );
    }
}

// An input that opens reaches `archive::format_archive`: one that is no archive is NCBI's invalid
// input (probes P3.077-083), in LOSAT's words (`wording`).
#[test]
fn an_archive_that_opens_is_formatted() {
    let invalid = format!(
        "LOSAT query/options error: Invalid input format for BLAST Archive.\n{}\n",
        crate::wording::USER_MANUAL_FOOTER
    );
    assert_eq!(
        run_words(&["-archive", "-", "-outfmt", "8"], b"not an archive\n"),
        (1, String::new(), invalid.clone())
    );
    assert_eq!(
        run_words(&["-archive", "-", "-outfmt", "8"], b""),
        (1, String::new(), invalid)
    );
}

// `CATCH_ALL`'s posts go where the diagnostics go: the `-logfile` file.
#[test]
fn the_run_s_errors_go_to_the_logfile() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_blast_formatter_log_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let log = dir.join("log.txt");
    let log_name = log.to_str().unwrap().to_string();
    let (status, out, err) = run_words(&["-logfile", &log_name], b"");
    assert_eq!((status, out.as_str(), err.as_str()), (1, "", ""));
    assert_eq!(
        std::fs::read_to_string(&log).unwrap(),
        "Command line argument error: Argument \"rid\". The argument has no value\n"
    );
    // LOSAT's refusals stay on the standard error stream.
    let (status, _, err) = run_words(&["-logfile", &log_name, "-rid", "X"], b"");
    assert_eq!(status, 2);
    assert!(
        err.starts_with("error: the NCBI BLAST+ option -rid"),
        "{err}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
