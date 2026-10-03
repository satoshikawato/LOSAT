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
    for (program, word) in [
        ("blastn", "11"),
        ("blastp", "3"),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:54-55
        // ```c++
        //     static const char kDefaultTask[] = "blastx";
        //     SetTask(kDefaultTask);
        // ```
        ("blastx", "3"),
        ("tblastx", "3"),
        ("tblastn", "3"),
    ] {
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:54-55
            // ```c++
            //     static const char kDefaultTask[] = "blastx";
            //     SetTask(kDefaultTask);
            // ```
            Commands::Blastx(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, 0.001);
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Blastn(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, Some(0.001));
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
                assert_eq!(a.max_target_seqs, Some(9));
            }
            // NCBI tblastn_args.cpp:55-62: SetTask(kDefaultTask);
            // blast_args.cpp:2726-2731: AddOptionalKey(kArgMaxTargetSequences, ...);
            Commands::Tblastn(a) => {
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
            // BLASTN names the NCBI options that it does not implement (AGENTS.md rule 2).
            if program == "blastn" && matches!(old, "-db" | "-remote") {
                assert!(
                    err.to_string().contains("not supported by LOSAT's BLASTN"),
                    "{err}"
                );
                continue;
            }
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
    // NCBI blast_filter.c:1147-1154: a value that is not above 0 keeps the default
    // (12, 2.2, 2.5); then blast_seg.c:2253-2258 raises hicut to locut.
    for (spec, expected) in [
        ("12 -1 -2", (12, 2.2, 2.5)),
        ("12 0 2.5", (12, 2.2, 2.5)),
        ("0 2.2 2.5", (12, 2.2, 2.5)),
        ("-5 2.2 0", (12, 2.2, 2.5)),
        ("12 3 2", (12, 3.0, 3.0)),
        ("15 2.5 3.0", (15, 2.5, 3.0)),
    ] {
        let params = parse_seg_filtering(spec).unwrap().params().unwrap();
        assert_eq!(
            (params.window, params.locut, params.hicut),
            expected,
            "{spec}"
        );
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
        "12  2.2 2.5",
        " 12 2.2 2.5",
        "12 2.2 2.5 ",
        "12\t2.2 2.5",
        "2147483648 2.2 2.5",
        "12 NaN 2.5",
        "12 2.2 inf",
    ] {
        assert!(parse_seg_filtering(bad).is_err(), "{bad}");
    }
    // NCBI blast_args.cpp:375-384,410-426: single spaces separate three signed numbers;
    // the masker replaces values out of range with its defaults (utils/dust.rs).
    for bad in [
        "",
        "true",
        "false",
        "Yes",
        "0",
        "20 64",
        "20 64 1 2",
        "20  64 1",
        " 20 64 1",
        "20 64 1 ",
        "20\u{a0}64 1",
        "x 64 1",
        "0x14 64 1",
        "NaN 64 1",
    ] {
        assert!(parse_dust_filtering(bad).is_err(), "{bad}");
    }
    for (spec, params) in [
        ("-1 64 1", (u32::MAX, 64, 1)),
        ("20 0 1", (20, 0, 1)),
        ("20 64 -1", (20, 64, u32::MAX as usize)),
    ] {
        assert_eq!(parse_dust_filtering(spec).unwrap().params(), Some(params));
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
    // NCBI blast_args.cpp:2800-2803: the default -outfmt is 0, which BLASTN and TBLASTX
    // implement.
    parse("blastn", &[]).unwrap();
    parse("tblastx", &[]).unwrap();
    for program in ["blastn", "tblastx"] {
        parse(program, &["-outfmt", "6"]).unwrap();
    }
    // NCBI blast_options_handle.cpp:344-380: dc-megablast and blastn-short are tasks of
    // blastn (Session SD); rmblastn stays an explicit rejection.
    parse("blastn", &["-outfmt", "6", "-task", "dc-megablast"]).unwrap();
    parse("blastn", &["-outfmt", "6", "-task", "blastn-short"]).unwrap();
    assert!(parse("blastn", &["-outfmt", "6", "-task", "rmblastn"]).is_err());
    // NCBI blast_args.cpp:708-730: -template_type and -template_length require each other
    // and take coding/optimal/coding_and_optimal and 16/18/21.
    parse(
        "blastn",
        &[
            "-task",
            "dc-megablast",
            "-template_type",
            "optimal",
            "-template_length",
            "21",
        ],
    )
    .unwrap();
    for words in [
        &["-template_type", "coding"][..],
        &["-template_length", "18"],
        &["-template_type", "Coding", "-template_length", "18"],
        &["-template_type", "coding", "-template_length", "17"],
        // blast_input_aux.hpp:214-222,239: the set is checked with a base-10 conversion.
        &["-template_type", "coding", "-template_length", "0x12"],
    ] {
        assert!(parse("blastn", words).is_err(), "{words:?}");
    }
    for length in ["018", "+18", "000000000000000018"] {
        parse(
            "blastn",
            &["-template_type", "coding", "-template_length", length],
        )
        .unwrap();
    }
    let words = ["-template_type", "coding", "-template_length", " 18"];
    assert!(parse("blastn", &words).is_err(), "{words:?}");
}

#[test]
fn numeric_values_are_validated_before_io() {
    for program in ["blastn", "blastp", "tblastx"] {
        // NCBI rejects an e-value that does not start with a digit, a point or a sign
        // (`inf`, `nan`, ` 1`) and one with trailing text; BLASTN's other forms are below.
        let evalues = if program == "blastn" {
            vec!["NaN", "inf", "nan", " 1", "1 ", "1,5", ""]
        } else {
            vec!["NaN", "inf", "-inf", "-1"]
        };
        for (key, values) in [
            ("-num_threads", vec!["0", "-1"]),
            ("-word_size", vec!["0", "-1"]),
            ("-max_target_seqs", vec!["0", "-1"]),
            ("-evalue", evalues),
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
    // NCBI blast_args.cpp:168-170: the argument is 4 or more; blast_options.c:1326-1333:
    // an option check rejects more than 100, with NCBI's message.
    assert!(parse("blastn", &["-outfmt", "6", "-word_size", "3"]).is_err());
    // Forms that NCBI reads (strtod) and LOSAT does not: rejected explicitly.
    for evalue in ["0x10", "0x1p-3", "1e"] {
        let error = parse("blastn", &["-outfmt", "6", "-evalue", evalue])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("not supported by LOSAT's BLASTN"),
            "-evalue {evalue}: {error}"
        );
    }
    // NCBI's option check rejects 0 or less and accepts infinity and NaN (E2g: NCBI
    // searches with them, as LOSAT does).
    for evalue in [
        "-1", "-inf", "0", "+inf", "+nan", "-nan", ".5e-400", "1e400",
    ] {
        let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-evalue", evalue])
            .unwrap()
            .command
        else {
            unreachable!("blastn")
        };
        let checked = LOSAT::algorithm::blastn::scoring::check_scoring_options(&args);
        let limited = LOSAT::algorithm::blastn::scoring::check_losat_limits(&args);
        if matches!(evalue, "+inf" | "+nan" | "-nan" | "1e400") {
            assert!(checked.is_ok(), "-evalue {evalue}");
            assert!(limited.is_ok(), "-evalue {evalue}");
        } else {
            let error = checked.unwrap_err().to_string();
            assert!(
                error.contains("expect value or cutoff score must be greater than zero"),
                "-evalue {evalue}: {error}"
            );
        }
    }
    // NCBI's integer arguments (ncbiargs.cpp:118-129) read hexadecimal after `0x`.
    let Commands::Blastn(args) = parse(
        "blastn",
        &[
            "-outfmt",
            "6",
            "-task",
            "blastn",
            "-reward",
            "0x2",
            "-penalty",
            "-3",
            "-gapopen",
            "0X5",
            "-gapextend",
            "0x2",
            "-word_size",
            "0xB",
            "-max_target_seqs",
            "0x10",
        ],
    )
    .unwrap()
    .command
    else {
        unreachable!("blastn")
    };
    assert_eq!(
        (args.reward, args.gap_open, args.gap_extend, args.word_size),
        (Some(2), Some(5), Some(2), Some(11))
    );
    assert_eq!(args.max_target_seqs, Some(16));
    // A constrained argument is read again with NStr::StringToDouble, which does not read
    // a bare 0x (blast_input_aux.hpp:110-113); -gapopen has no constraint.
    for (key, value) in [("-reward", "0x"), ("-penalty", "0X"), ("-word_size", "0x")] {
        assert!(
            parse("blastn", &["-outfmt", "6", key, value]).is_err(),
            "{key} {value}"
        );
    }
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-gapopen", "0x"])
        .unwrap()
        .command
    else {
        unreachable!("blastn")
    };
    assert_eq!(args.gap_open, Some(0));
    // NCBI's other blastn tasks and options are rejected explicitly (AGENTS.md rule 2).
    for extra in [
        &["-task", "rmblastn"][..],
        &["-strand", "plus"],
        &["-ungapped"],
        &["-h"],
        &["-version"],
    ] {
        let mut words = vec!["-outfmt", "6"];
        words.extend_from_slice(extra);
        let error = LOSAT::cli::render_message(&parse("blastn", &words).unwrap_err());
        assert!(
            error.contains("not supported by LOSAT's BLASTN"),
            "{extra:?}: {error}"
        );
    }
    assert!(parse("blastn", &["-outfmt", "6", "-task", "BLASTN"]).is_err());
    // NCBI blast_input_aux.hpp:79-87 (CDirEntry::GetName, ncbifile.cpp:298-312,358-363,
    // 465-472): the name after the trailing separators are removed.
    let dotted = format!("{}/.", "o".repeat(256));
    parse("blastn", &["-outfmt", "6", "-out", &dotted]).unwrap();
    assert!(parse(
        "blastn",
        &["-outfmt", "6", "-out", &format!("{}/", "o".repeat(256))]
    )
    .is_err());
    // NCBI blast_input_aux.hpp:79-87: an -out file name is shorter than 256 bytes.
    let long = format!("dir/{}", "o".repeat(256));
    assert!(parse("blastn", &["-outfmt", "6", "-out", &long]).is_err());
    parse("blastn", &["-outfmt", "6", "-out", &long[..long.len() - 1]]).unwrap();
    for (key, value) in [
        ("-max_target_seqs", "0x"),
        ("-max_target_seqs", "+0x1"),
        ("-penalty", "-0x3"),
        ("-word_size", "0x3"),
        ("-reward", "0x80000000"),
        ("-gapopen", "0x5g"),
    ] {
        assert!(
            parse("blastn", &["-outfmt", "6", key, value]).is_err(),
            "{key} {value}"
        );
    }
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-penalty", "0x0"])
        .unwrap()
        .command
    else {
        unreachable!("blastn")
    };
    assert!(
        LOSAT::algorithm::blastn::scoring::check_scoring_options(&args)
            .unwrap_err()
            .to_string()
            .contains("BLASTN penalty must be negative")
    );
    for (value, accepted) in [("4", true), ("100", true), ("101", false)] {
        let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-word_size", value])
            .unwrap()
            .command
        else {
            unreachable!("blastn")
        };
        let checked = LOSAT::algorithm::blastn::scoring::check_scoring_options(&args);
        assert_eq!(checked.is_ok(), accepted, "-word_size {value}");
        if !accepted {
            assert!(checked
                .unwrap_err()
                .to_string()
                .contains("Word-size must be less than or equal to 100"));
        }
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
        ("tblastx", vec!["5", "6 qseqid", "7 std"]),
    ] {
        for spec in specs {
            assert!(parse(program, &["-outfmt", spec]).is_err());
        }
    }
    // TBLASTX implements the pairwise report and both tabular formats (session S08).
    for spec in ["0", "6", "7"] {
        parse("tblastx", &["-outfmt", spec]).unwrap();
    }
    // BLASTN parses -outfmt when it sets the options, as NCBI (blast_args.cpp:2801-2851).
    for spec in ["5", "6 qseqid", "0x6"] {
        parse("blastn", &["-outfmt", spec]).unwrap();
        assert!(LOSAT::algorithm::blastn::hsp::parse_blastn_output_format(spec).is_err());
    }
    // NCBI blast_args.cpp:2845-2851: a custom specification counts only for the tabular
    // formats, so BLASTN's outfmt 0 ignores it.
    parse("blastn", &["-outfmt", "0 qseqid"]).unwrap();
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
