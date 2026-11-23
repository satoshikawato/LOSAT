use super::*;
use crate::cli::{Cli, Commands};

fn parse(extra: &[&str]) -> Result<RpstblastnArgs, String> {
    let mut argv = vec!["losat", "rpstblastn", "-query", "q.fna", "-db", "cdd"];
    argv.extend_from_slice(extra);
    match crate::cli::try_parse_from::<Cli, _, _>(argv) {
        Ok(Cli {
            command: Commands::Rpstblastn(args),
            ..
        }) => Ok(args),
        Ok(_) => Err("not rpstblastn".to_string()),
        Err(error) => Err(crate::cli::render_message(&error)),
    }
}

/// The refusal of `check_options`, once per run.
fn late_refusal(extra: &[&str]) -> String {
    search(extra).check_options().unwrap_err().to_string()
}

fn search(extra: &[&str]) -> TranslatedRpsSearch {
    parse(extra).unwrap().rps_search().unwrap()
}

// ncbiargs.cpp:3829-3834 (cited at `USAGE`): `-help` and `-h` write NCBI's texts, whatever
// follows them; the words before them are read first.
#[test]
fn help_and_usage_are_ncbis_texts() {
    assert_eq!(HELP.len(), 9266);
    assert!(HELP.ends_with("   Compute locally optimal Smith-Waterman alignments?\n\n"));
    let head = &HELP[..HELP.find("OPTIONAL ARGUMENTS").unwrap()];
    assert_eq!(
        USAGE,
        format!("{head}Use '-help' to print detailed descriptions of command line arguments\n")
    );
    assert!(USAGE.contains(&format!(
        "   Translated Reverse Position Specific BLAST {}\n",
        crate::NCBI_BLAST_VERSION
    )));
    for (words, text) in [
        (&["-help"][..], HELP),
        (&["-help", "-no_such_option"], HELP),
        (&["--help"], HELP),
        (&["-h"], USAGE),
        (&["-evalue", "1", "-h", "-xyz", "2"], USAGE),
        // NCBI checks the dependencies of the words after its reader (`x_PostCheck`).
        (
            &["-max_target_seqs", "1", "-num_alignments", "1", "-help"],
            HELP,
        ),
    ] {
        assert_eq!(
            parse(words).unwrap().ncbi_help_text(),
            Some(text),
            "{words:?}"
        );
    }
    assert_eq!(parse(&[]).unwrap().ncbi_help_text(), None);
    // The hidden argument is the reader's: no word of the user names it.
    for words in [&["-ncbi_help", "full"][..], &["--ncbi_help=full"]] {
        assert!(
            parse(words).unwrap_err().contains("unknown option"),
            "{words:?}"
        );
    }
    // A word before the help that NCBI's reader refuses is refused first.
    assert!(parse(&["-no_such_option", "-help"])
        .unwrap_err()
        .contains("unknown option"));
    assert!(parse(&["-evalue", "x", "-help"])
        .unwrap_err()
        .contains("-evalue"));
}

// blast_args.cpp:839-865 (cited at `rps_composition_mode`): RPS-BLAST's letters; every mode
// runs (item PS12b); with -ungapped BLASTX's resolve gives NCBI's error.
#[test]
fn composition_modes_run() {
    for value in ["0", "F", "f", "D", "d"] {
        let args = parse(&["-comp_based_stats", value])
            .unwrap()
            .into_blastx_args()
            .unwrap();
        assert_eq!(args.comp_based_stats, "0", "{value}");
        assert!(args.rpstblastn.is_some());
        assert_eq!(
            search(&["-comp_based_stats", value]).search.requested_cbs,
            0
        );
    }
    for (value, mode) in [("1", 1), ("T", 1), ("t", 1), ("2", 2), ("3", 3)] {
        let args = parse(&["-comp_based_stats", value])
            .unwrap()
            .into_blastx_args()
            .unwrap();
        assert_eq!(args.comp_based_stats, mode.to_string(), "{value}");
        let search = search(&["-comp_based_stats", value]);
        assert_eq!(search.search.requested_cbs, mode, "{value}");
        assert!(search.check_options().is_ok(), "{value}");
    }
    // The default (1).
    let default = search(&[]);
    assert_eq!(default.search.requested_cbs, 1);
    assert!(default.check_options().is_ok());
    // ungapped with the default: NCBI's options error (BLASTX's resolve).
    let error = parse(&["-ungapped"])
        .unwrap()
        .into_blastx_args()
        .unwrap()
        .resolve()
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Composition-adjusted searched are not supported with an ungapped search"),
        "{error}"
    );
    // ungapped without statistics reaches the search, which refuses it.
    assert!(search(&["-ungapped", "-comp_based_stats", "0"]).ungapped);
}

// The refusals of the options that later items run, before the search. Every output format
// that NCBI's rpstblastn runs is BLASTX's report (items PS15a-1, PS15a-2).
#[test]
fn options_of_later_items_are_refused() {
    for outfmt in [
        "1", "4", "5", "8", "9", "11", "12", "13", "14", "15", "16", "18", "20",
    ] {
        parse(&["-outfmt", outfmt, "-out", "x.out"])
            .unwrap()
            .into_blastx_args()
            .unwrap_or_else(|error| panic!("outfmt {outfmt}: {error}"));
    }
    // -use_sw_tback with composition-based statistics (item PS12c): the redo of the frames is
    // NCBI's Smith-Waterman redo.
    assert!(search(&["-use_sw_tback"]).search.use_sw_tback);
    assert!(
        search(&["-use_sw_tback", "-comp_based_stats", "2"])
            .search
            .use_sw_tback
    );
    assert!(!search(&[]).search.use_sw_tback);
    // LOSAT's refusals of options, once per run (`check_options`).
    assert!(late_refusal(&["-sum_stats", "true"]).contains("-sum_stats true"));
    // -searchsp and -dbsize reach BLASTX's options.
    let args = parse(&[
        "-comp_based_stats",
        "0",
        "-searchsp",
        "1000",
        "-dbsize",
        "5000",
    ])
    .unwrap()
    .into_blastx_args()
    .unwrap();
    assert_eq!((args.searchsp, args.dbsize), (Some(1000), Some(5000)));
    // The formats LOSAT writes (1 to 4 through BLASTX's report, line RG item T3c), a custom
    // field list and NCBI's own error of 17 (BLASTX's). 13 and 14 to the standard output go on
    // to NCBI's error too (BLASTX's).
    for outfmt in [
        "0",
        "1",
        "2",
        "3",
        "4",
        "6",
        "7",
        "10",
        "6 qseqid sseqid qframe",
        "17",
        "13",
        "14",
    ] {
        assert!(
            parse(&["-comp_based_stats", "0", "-outfmt", outfmt])
                .unwrap()
                .into_blastx_args()
                .is_ok(),
            "{outfmt}"
        );
    }
    // -use_sw_tback without statistics: accepted, NCBI's search does not read it.
    assert!(parse(&["-use_sw_tback", "-comp_based_stats", "0"])
        .unwrap()
        .into_blastx_args()
        .is_ok());
    for name in ["html", "remote", "entrez_query", "import_search_strategy"] {
        let error = parse(&[&format!("-{name}"), "x"]).unwrap_err();
        assert!(
            error.contains(&format!("-{name} is not supported by LOSAT's RPSTBLASTN")),
            "{name}: {error}"
        );
    }
    // NCBI's options of the other programs: LOSAT's parser (PD exception 1).
    for name in [
        "subject",
        "word_size",
        "matrix",
        "threshold",
        "gapopen",
        "task",
    ] {
        assert!(
            parse(&[&format!("-{name}"), "1"])
                .unwrap_err()
                .contains("unknown option"),
            "{name}"
        );
    }
}

// NCBI's values: the genetic codes, the strands, the threads and modes, the formatting values.
#[test]
fn values_are_ncbis() {
    for extra in [
        &["-query_gencode", "7"][..],
        &["-query_gencode", "0"],
        &["-strand", "foo"],
        &["-num_threads", "-1"],
        &["-mt_mode", "2"],
        &["-line_length", "0"],
        &["-sorthits", "5"],
        &["-max_hsps", "0"],
        &["-max_target_seqs", "0"],
        &["-qcov_hsp_perc", "101"],
        &["-sum_stats", "maybe"],
        &["-xdrop_gap", "nan"],
        &["-window_size", "-1"],
        &["-searchsp", "-1"],
        &["-max_target_seqs", "5", "-num_alignments", "2"],
    ] {
        assert!(parse(extra).is_err(), "{extra:?}");
    }
    for code in [
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 25, 31, 33,
    ] {
        assert_eq!(
            parse(&["-query_gencode", &code.to_string()])
                .unwrap()
                .query_gencode,
            code
        );
    }
}

// What reaches BLASTX's arguments and the RPS search: the query handling, the masks (no SEG
// by default), the counts, the X-drops, `-max_hsps` and `-qcov_hsp_perc` (the RPS
// traceback's), the formatting options, `-num_threads 0`, no sum statistics.
#[test]
fn the_search_is_blastxs_with_the_rps_search() {
    let words = [
        "-comp_based_stats",
        "F",
        "-strand",
        "minus",
        "-query_gencode",
        "4",
        "-query_loc",
        "10-500",
        "-lcase_masking",
        "-max_hsps",
        "2",
        "-qcov_hsp_perc",
        "50",
        "-xdrop_ungap",
        "3",
        "-xdrop_gap_final",
        "0",
        "-line_length",
        "70",
        "-sorthits",
        "1",
        "-num_threads",
        "0",
        "-mt_mode",
        "1",
        "-max_target_seqs",
        "7",
        "-evalue",
        "0.5",
        "-window_size",
        "0",
    ];
    let args = parse(&words).unwrap().into_blastx_args().unwrap();
    assert_eq!(
        (
            args.strand.as_str(),
            args.query_gencode,
            args.query_loc.as_deref()
        ),
        ("minus", 4, Some("10-500"))
    );
    assert!(args.lcase_masking);
    assert_eq!((args.seg.as_str(), args.soft_masking), ("no", false));
    assert_eq!((args.max_hsps, args.qcov_hsp_perc), (None, None));
    assert_eq!((args.max_target_seqs, args.evalue), (Some(7), 0.5));
    assert_eq!((args.num_threads, args.window_size), (1, Some(0)));
    assert_eq!(args.sum_stats, Some(false));
    let rps = search(&words);
    assert_eq!(
        (
            rps.search.max_hsps_per_subject,
            rps.search.query_cov_hsp_perc,
            rps.search.auto_threads
        ),
        (2, 50.0, true)
    );
    assert_eq!(
        (
            rps.search.x_drop_ungapped_bits,
            rps.search.x_drop_gapped_bits,
            rps.search.x_drop_final_bits
        ),
        (Some(3.0), None, Some(0.0))
    );
    // The formatting options are BLASTX's (one owner): the line length reaches the report, the
    // sorts sort nothing (NCBI's rpstblastn hands its formatter the line length only).
    assert_eq!(
        args.formatting,
        FormattingArgs {
            line_length: Some(70),
            sort_hits: Some(1),
            sort_hsps: None,
            hitlist_start: None,
        }
    );
    let resolved = args.resolve().unwrap();
    assert_eq!(resolved.formatting.line_length(), 70);
    assert_eq!(resolved.report_sort(), Default::default());
    // RPS-BLAST prints no word threshold (`blast_rps_options.cpp:60-67`, cited above).
    assert_eq!(resolved.threshold, 0.0);
    assert!(!resolved.sum_stats);
    assert_eq!(resolved.composition, 0);
    assert!(resolved.rpstblastn.is_some());
    // `-sum_stats true` reaches the search, which refuses it.
    let words = [
        "-comp_based_stats",
        "0",
        "-sum_stats",
        "true",
        "-seg",
        "yes",
    ];
    assert!(search(&words).search.sum_stats);
    assert_eq!(
        parse(&words).unwrap().into_blastx_args().unwrap().seg,
        "yes"
    );
}
