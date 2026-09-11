//! Frozen CLI v2 grammar and configuration boundary.
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-207,332-420,838-883,3158-3163
// AddOptionalKey(kArgWordSize, ...); AddOptionalKey(kArgMaxHSPsPerSubject, ...);
// x_TokenizeFilteringArgs(seg_opts, tokens); opt.SetUnifiedP(1);

use LOSAT::algorithm::blastp::args::parse_comp_based_stats;
use LOSAT::blastinput::value_parsers::{
    parse_dust_filtering, parse_seg_filtering, DustSpec, SegSpec,
};
use LOSAT::cli::{try_parse_from, Cli, Commands};

fn parse(program: &str, extra: &[&str]) -> Result<Cli, clap::Error> {
    let mut argv = vec!["losat", program, "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(extra);
    try_parse_from(argv)
}

#[test]
fn canonical_options_for_every_program() {
    for (program, word) in [("blastn", "11"), ("blastp", "3"), ("tblastx", "3")] {
        let cli = parse(
            program,
            &[
                "-outfmt",
                "6",
                "-evalue",
                "0.001",
                "-num_threads",
                "4",
                "-word_size",
                word,
                "-max_target_seqs",
                "9",
            ],
        )
        .unwrap();
        match cli.command {
            Commands::Blastn(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, 0.001);
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Blastp(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, Some(0.001));
                assert_eq!(a.max_target_seqs, 9);
            }
            Commands::Tblastx(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, 0.001);
                assert_eq!(a.max_target_seqs, 9);
            }
        }
    }
}

#[test]
fn legacy_and_unsupported_syntax_is_unknown() {
    for program in ["blastn", "blastp", "tblastx"] {
        for old in [
            "--query",
            "--word-size",
            "--num-threads",
            "--max-target-seqs",
            "--use_sw_tback",
            "--seg=false",
            "--seg",
            "-seg_window",
            "-max-hsps-per-subject",
            "-hitlist_size",
            "-scan_step",
            "-min_diag_separation",
            "-q",
            "-db",
            "-remote",
            "--word_size",
            "--dust-window",
        ] {
            let err = parse(program, &["-outfmt", "6", old]).unwrap_err();
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::UnknownArgument,
                "{program} {old}: {err}"
            );
        }
    }
}

#[test]
fn filtering_has_exactly_one_shared_value_grammar() {
    assert_eq!(parse_seg_filtering("no").unwrap(), SegSpec::No);
    assert_eq!(parse_seg_filtering("yes").unwrap(), SegSpec::Yes);
    assert_eq!(
        parse_seg_filtering("12 2.2 2.5")
            .unwrap()
            .to_ncbi_cli_string(),
        "12 2.2 2.5"
    );
    // NCBI blast_seg.c:2253-2258: clamp negative cutoffs and raise hicut to locut.
    for (spec, expected) in [("12 -1 -2", (0.0, 0.0)), ("12 3 2", (3.0, 3.0))] {
        let params = parse_seg_filtering(spec).unwrap().params().unwrap();
        assert_eq!((params.locut, params.hicut), expected);
    }
    assert_eq!(parse_dust_filtering("no").unwrap(), DustSpec::No);
    assert_eq!(
        parse_dust_filtering("yes").unwrap().params(),
        Some((20, 64, 1))
    );
    assert_eq!(
        parse_dust_filtering("30 32 2").unwrap().params(),
        Some((30, 32, 2))
    );
    for bad in [
        "",
        "true",
        "false",
        "0",
        "1",
        "12 2.2",
        "12 2.2 2.5 4",
        "x 2.2 2.5",
        "0 2.2 2.5",
        "12 NaN 2.5",
        "12 2.2 inf",
    ] {
        assert!(parse_seg_filtering(bad).is_err(), "{bad}");
    }
    for bad in [
        "",
        "true",
        "false",
        "0",
        "20 64",
        "20 64 1 2",
        "x 64 1",
        "20 0 1",
        "NaN 64 1",
        "20 64 -1",
    ] {
        assert!(parse_dust_filtering(bad).is_err(), "{bad}");
    }
    for program in ["blastp", "tblastx"] {
        for value in ["no", "yes", "12 2.2 2.5"] {
            parse(program, &["-outfmt", "6", "-seg", value]).unwrap();
        }
        assert!(parse(program, &["-outfmt", "6", "-seg"]).is_err());
    }
    for value in ["no", "yes", "20 64 1"] {
        parse("blastn", &["-outfmt", "6", "-dust", value]).unwrap();
    }
}

#[test]
fn composition_tokens_round_trip_and_reject_suffix_garbage() {
    for token in [
        "0", "1", "2", "3", "F", "f", "D", "d", "T", "t", "1u", "2u", "3u", "Du", "du", "Tu", "tu",
    ] {
        let parsed = parse_comp_based_stats(token).unwrap();
        assert_eq!(
            parse_comp_based_stats(&parsed.to_ncbi_cli_string()).unwrap(),
            parsed
        );
        parse("blastp", &["-comp_based_stats", token]).unwrap();
    }
    for token in [
        "", "2garbage", "2uu", "2uX", "2U", "0u", "Fu", "4", "22", " 2", "2 ",
    ] {
        assert!(parse_comp_based_stats(token).is_err(), "{token}");
        assert!(
            parse("blastp", &["-comp_based_stats", token]).is_err(),
            "{token}"
        );
    }
}

#[test]
fn defaults_and_task_overrides_remain_distinct() {
    let Commands::Blastp(args) = parse("blastp", &[]).unwrap().command else {
        panic!()
    };
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.evalue, None);
    assert_eq!(args.outfmt, "0");
    assert_eq!(args.max_target_seqs, 500);
    assert_eq!(args.max_hsps_per_subject, None);
    assert_eq!(args.word_size, None);
    assert_eq!(args.resolve().unwrap().max_hsps_per_subject, 0);
    // NCBI api/blast_options_handle.cpp:381-387: Create(eBlastp, locality).
    assert_eq!(args.task, "blastp");
    assert_eq!(args.resolve().unwrap().evalue, 10.0);
    let Commands::Blastp(a) = parse("blastp", &["-task", "blastp"]).unwrap().command else {
        panic!()
    };
    assert_eq!(a.task, "blastp");
    assert_eq!(a.word_size, None);
    assert_eq!(a.resolve().unwrap().word_size, 3);
    assert_eq!(a.evalue, None);
    assert_eq!(a.resolve().unwrap().evalue, 10.0);
    let Commands::Blastp(a) = parse("blastp", &["-task", "blastp", "-evalue", "42"])
        .unwrap()
        .command
    else {
        panic!()
    };
    assert_eq!(a.resolve().unwrap().evalue, 42.0);
    for task in ["blastp-short", "blastp-fast"] {
        let error = parse("blastp", &["-task", task]).unwrap_err();
        assert_eq!(error.kind(), clap::error::ErrorKind::InvalidValue);
        assert!(error.to_string().contains(task));
    }
    for program in ["blastn", "tblastx"] {
        assert!(parse(program, &[]).unwrap_err().to_string().contains("0"));
        parse(program, &["-outfmt", "6"]).unwrap();
    }
    assert!(parse("blastn", &["-outfmt", "6", "-task", "dc-megablast"]).is_err());
}

#[test]
fn numeric_values_are_validated_before_io() {
    for program in ["blastn", "blastp", "tblastx"] {
        for (key, values) in [
            ("-num_threads", vec!["0", "-1"]),
            ("-word_size", vec!["0", "-1"]),
            ("-max_target_seqs", vec!["0", "-1"]),
            ("-evalue", vec!["NaN", "inf", "-inf", "-1"]),
        ] {
            for value in values {
                assert!(
                    parse(program, &["-outfmt", "6", key, value]).is_err(),
                    "{program} {key} {value}"
                );
            }
        }
    }
    for program in ["blastn", "blastp"] {
        assert!(parse(program, &["-outfmt", "6", "-max_hsps", "0"]).is_err());
        parse(program, &["-outfmt", "6", "-max_hsps", "1"]).unwrap();
    }
    for value in ["2", "4"] {
        assert!(parse("tblastx", &["-outfmt", "6", "-word_size", value]).is_err());
    }
    // NCBI core/blast_options.c:1322-1334: word_size is bounded by 4 and 100.
    for value in ["3", "101"] {
        assert!(parse("blastn", &["-outfmt", "6", "-word_size", value]).is_err());
    }
    for value in ["4", "100"] {
        parse("blastn", &["-outfmt", "6", "-word_size", value]).unwrap();
    }
    for program in ["blastp", "tblastx"] {
        assert!(parse(program, &["-outfmt", "6", "-threshold", "0"]).is_err());
        assert!(parse(program, &["-outfmt", "6", "-window_size", "2147483648"]).is_err());
        parse(program, &["-outfmt", "6", "-window_size", "0"]).unwrap();
    }
    for value in ["0", "7", "8", "32", "255"] {
        assert!(parse("tblastx", &["-outfmt", "6", "-db_gencode", value]).is_err());
    }
}

#[test]
fn flags_and_value_tokens_are_not_reinterpreted() {
    // NCBI api/blast_advprot_options.cpp:58: SetSmithWatermanMode(false).
    // v0.1.0 does not expose the unported Smith-Waterman path.
    let error = parse("blastp", &["-use_sw_tback"]).unwrap_err();
    assert_eq!(error.kind(), clap::error::ErrorKind::UnknownArgument);
    for value in [
        "-use_sw_tback=true",
        "-use_sw_tback=false",
        "-use_sw_tback=yes",
    ] {
        assert!(parse("blastp", &[value]).is_err());
    }
    assert!(parse("blastp", &["-use_sw_tback", "true"]).is_err());
    let cli: Cli = try_parse_from([
        "losat",
        "blastp",
        "-query",
        "-num_threads",
        "-subject",
        "-query",
    ])
    .unwrap();
    let Commands::Blastp(a) = cli.command else {
        panic!()
    };
    assert_eq!(a.query.to_str(), Some("-num_threads"));
    assert_eq!(a.subject.to_str(), Some("-query"));
    assert!(
        try_parse_from::<Cli, _, _>(["losat", "blastp", "-query", "-", "-subject", "s"]).is_err()
    );
}

#[test]
fn output_capabilities_fail_explicitly_and_help_is_canonical() {
    parse("blastp", &["-outfmt", "6 qseqid sseqid pident length"]).unwrap();
    for (program, specs) in [
        ("blastp", vec!["5", "0 qseqid", "6 unknown"]),
        ("blastn", vec!["0", "6 qseqid"]),
        ("tblastx", vec!["0", "7", "6 qseqid"]),
    ] {
        for spec in specs {
            assert!(parse(program, &["-outfmt", spec]).is_err());
        }
    }
    for program in ["blastn", "blastp", "tblastx"] {
        for help in ["-help", "--help"] {
            let err = try_parse_from::<Cli, _, _>(["losat", program, help]).unwrap_err();
            assert_eq!(err.kind(), clap::error::ErrorKind::DisplayHelp);
            let text = LOSAT::cli::render_message(&err);
            assert!(text.contains("-query <PATH>"), "{text}");
            assert!(text.contains("-num_threads"));
            assert!(
                !text.contains("--query") && !text.contains("--word") && !text.contains("--num")
            );
            assert!(!text.contains("seg_window") && !text.contains("dust_window"));
            if program == "blastp" {
                for hidden in ["blastp-short", "blastp-fast", "use_sw_tback"] {
                    assert!(!text.contains(hidden), "{text}");
                }
            }
        }
    }
}
