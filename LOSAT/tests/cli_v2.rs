//! Frozen CLI v2 grammar and configuration boundary.
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-207,332-420,838-883,3158-3163
// AddOptionalKey(kArgWordSize, ...); AddOptionalKey(kArgMaxHSPsPerSubject, ...);
// x_TokenizeFilteringArgs(seg_opts, tokens); opt.SetUnifiedP(1);

use LOSAT::algorithm::blastp::args::parse_comp_based_stats;
use LOSAT::blastinput::value_parsers::{
    parse_dust_filtering, parse_seg_filtering, DustSpec, SegSpec,
};
use LOSAT::cli::{try_parse_from, Cli, Commands};

fn parse_argv(argv: Vec<&str>) -> Result<Cli, clap::Error> {
    try_parse_from(argv)
}

fn parse(program: &str, extra: &[&str]) -> Result<Cli, clap::Error> {
    let mut argv = vec!["losat", program, "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(extra);
    try_parse_from(argv)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2263-2299 (the
// taxid lists exclude each other); blast_args.cpp:2370-2381 (`-subject` excludes them).
#[test]
fn taxid_options_go_with_a_database_only() {
    let lists: [&[&str]; 5] = [
        &["-taxids", "9606"],
        &["-taxidlist", "t.txt"],
        &["-negative_taxids", "9606"],
        &["-negative_taxidlist", "t.txt"],
        &["-no_taxid_expansion"],
    ];
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        for words in lists {
            let mut argv = vec!["losat", program, "-query", "q.fa", "-db", "db"];
            argv.extend(words);
            assert!(parse_argv(argv).is_ok(), "{program} -db {words:?}");
            let error = LOSAT::cli::render_message(&parse(program, words).unwrap_err());
            assert!(
                error.contains("cannot be used with"),
                "{program} {words:?}: {error}"
            );
        }
        let argv = vec![
            "losat",
            program,
            "-query",
            "q.fa",
            "-db",
            "db",
            "-taxids",
            "1",
            "-negative_taxids",
            "2",
        ];
        assert!(parse_argv(argv).is_err(), "{program}: two lists");
        let argv = vec![
            "losat",
            program,
            "-query",
            "q.fa",
            "-db",
            "db",
            "-taxids",
            "1",
            "-no_taxid_expansion",
        ];
        assert!(
            parse_argv(argv).is_ok(),
            "{program}: a list without expansion"
        );
    }
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:58-116
        // psiblast's options are blastp's search options without -task and -ungapped.
        ("psiblast", "3"),
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
                assert_eq!(a.max_target_seqs, Some(9));
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
                assert_eq!(a.evalue, Some(0.001));
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Psiblast(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, Some(0.001));
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Makedb(_) => unreachable!("{program} is a search program"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:225-246
// ```c++
//     arg_desc->AddDefaultKey(kInput, "input_file",
//                             "Input file/database name",
//                             CArgDescriptions::eInputFile, "-");
//     arg_desc->AddDefaultKey("input_type", "type",
//                              "Type of the data specified in input_file",
//                              CArgDescriptions::eString, "fasta");
//     arg_desc->SetConstraint("input_type", &(*new CArgAllow_Strings,
//                                             "fasta", "blastdb",
//                                             "asn1_bin",
//                                             "asn1_txt"));
//
//     arg_desc->AddKey(kArgDbType, "molecule_type",
//                      "Molecule type of target db", CArgDescriptions::eString);
//     arg_desc->SetConstraint(kArgDbType, &(*new CArgAllow_Strings,
//                                         "nucl", "prot"));
//
//     arg_desc->SetCurrentGroup("Configuration options");
//     arg_desc->AddOptionalKey(kArgDbTitle, "database_title",
//                              "Title for BLAST database\n" + dflt,
//                              CArgDescriptions::eString);
//
//     arg_desc->AddFlag("parse_seqids",
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:293-294
// ```c++
//     arg_desc->AddOptionalKey(kOutput, "database_name",
//                              "Name of BLAST database to be created\n" + dflt +
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:304-305
// ```c++
//     arg_desc->AddDefaultKey("max_file_sz", "number_of_bytes",
//                             "Maximum file size for BLAST database files",
// ```
#[test]
fn makedb_takes_makeblastdb_s_option_names() {
    let makedb = |extra: &[&str]| {
        let mut argv = vec!["losat", "makeblastdb"];
        argv.extend(extra);
        try_parse_from(argv)
    };
    let args = |cli: Cli| match cli.command {
        Commands::Makedb(args) => args,
        _ => panic!("not makedb"),
    };

    let given = args(
        makedb(&[
            "-in",
            "x.faa",
            "-dbtype",
            "prot",
            "-out",
            "db",
            "-title",
            "a title",
            "-max_file_sz",
            "1GB",
        ])
        .unwrap(),
    );
    assert_eq!(given.input, "x.faa");
    assert_eq!(given.dbtype, "prot");
    assert_eq!(given.out.as_deref(), Some("db"));
    assert_eq!(given.title.as_deref(), Some("a title"));
    assert_eq!(given.max_file_sz, "1GB");

    // makeblastdb's defaults.
    let defaults = args(makedb(&["-dbtype", "nucl"]).unwrap());
    assert_eq!(defaults.input, "-");
    assert_eq!(defaults.input_type, "fasta");
    assert_eq!(defaults.out, None);
    assert_eq!(defaults.title, None);
    assert_eq!(defaults.max_file_sz, "3GB");
    assert_eq!(defaults.blastdb_version, 5);
    assert_eq!(defaults.logfile, None);

    // `-dbtype` is mandatory and takes makeblastdb's two values.
    assert!(makedb(&["-in", "x.faa"]).is_err());
    assert!(makedb(&["-in", "x.faa", "-dbtype", "protein"]).is_err());
    assert!(makedb(&["-in", "x.faa", "-dbtype", "prot", "-blastdb_version", "3"]).is_err());

    assert!(!defaults.hash_index);
    assert_eq!(defaults.taxid, None);

    // `-hash_index` and `-taxid` (makeblastdb.cpp:249-251, 315-318: a flag, and an integer
    // of at least 0 that excludes `-taxid_map`).
    let given = args(makedb(&["-dbtype", "prot", "-hash_index", "-taxid", "9606"]).unwrap());
    assert!(given.hash_index);
    assert_eq!(given.taxid, Some(9606));
    assert!(makedb(&["-dbtype", "prot", "-taxid", "-1"]).is_err());
    assert!(makedb(&["-dbtype", "prot", "-taxid", "abc"]).is_err());

    // `-parse_seqids` is a flag; `-taxid_map` takes a file and requires `-parse_seqids`, and
    // `-taxid` excludes it (makeblastdb.cpp:246-247, 322-329). The parser's errors end with
    // exit status 2 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES; makeblastdb prints its usage, 1).
    assert!(!defaults.parse_seqids);
    assert_eq!(defaults.taxid_map, None);
    let given = args(
        makedb(&[
            "-dbtype",
            "prot",
            "-parse_seqids",
            "-taxid_map",
            "map.txt",
            "-hash_index",
        ])
        .unwrap(),
    );
    assert!(given.parse_seqids && given.hash_index);
    assert_eq!(given.taxid_map.as_deref(), Some("map.txt"));
    let given = args(makedb(&["-dbtype", "prot", "-parse_seqids", "-taxid", "9606"]).unwrap());
    assert!(given.parse_seqids);
    assert_eq!(given.taxid, Some(9606));
    for words in [
        &["-dbtype", "prot", "-taxid_map", "map.txt"][..],
        &[
            "-dbtype",
            "prot",
            "-parse_seqids",
            "-taxid_map",
            "map.txt",
            "-taxid",
            "9606",
        ],
        &[
            "-dbtype",
            "prot",
            "-taxid",
            "1",
            "-taxid_map",
            "map.txt",
            "-parse_seqids",
        ],
        &["-dbtype", "prot", "-parse_seqids", "-taxid_map"],
    ] {
        let err = makedb(words).unwrap_err();
        assert_eq!(err.exit_code(), 2, "{words:?}: {err}");
        assert!(
            !err.to_string().contains("not supported"),
            "{words:?}: {err}"
        );
    }

    // The options of makeblastdb that LOSAT does not implement are named (`-h` and
    // `-version` are served: `makeblastdb_is_a_name_of_makedb`).
    for option in ["-mask_data", "-gi_mask", "-randomize", "-oid_masks"] {
        let err = makedb(&["-in", "x.faa", "-dbtype", "prot", option]).unwrap_err();
        assert!(
            err.to_string().contains(&format!(
                "the NCBI BLAST+ option {option} is not supported by LOSAT's makeblastdb"
            )),
            "{option}: {err}"
        );
    }

    // `-help` prints LOSAT's usage of makeblastdb, under makeblastdb's option names.
    for help in ["-help", "--help"] {
        let err = makedb(&[help]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::DisplayHelp);
        let text = LOSAT::cli::render_message(&err);
        assert!(text.contains("-dbtype <molecule_type>"), "{text}");
        assert!(text.contains("-max_file_sz"), "{text}");
        assert!(
            !text.contains("--in") && !text.contains("--dbtype"),
            "{text}"
        );
    }

    // The options of the search programs are unknown to makeblastdb, and makeblastdb's
    // options are unknown to the search programs.
    assert!(makedb(&["-dbtype", "prot", "-evalue", "1"]).is_err());
    assert!(makedb(&["-dbtype", "prot", "-query", "q.fa"]).is_err());
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        for option in [
            &["-dbtype", "prot"][..],
            &["-parse_seqids"],
            &["-in", "x.faa"],
        ] {
            let err = parse(program, option).unwrap_err();
            assert!(!err.to_string().contains("makeblastdb"), "{program}: {err}");
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
            // The programs name the NCBI options that they do not implement (AGENTS.md
            // rule 2). BLASTP, BLASTN, TBLASTN, BLASTX and TBLASTX have `-db` (it takes a value,
            // which is missing here).
            if matches!(
                program,
                "blastp" | "blastn" | "tblastn" | "blastx" | "tblastx"
            ) && old == "-db"
            {
                assert!(err.to_string().contains("-db requires one value"), "{err}");
                continue;
            } else if matches!(old, "-db" | "-remote") {
                assert!(
                    err.to_string().contains(&format!(
                        "not supported by LOSAT's {}",
                        program.to_ascii_uppercase()
                    )),
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
fn composition_tokens_follow_ncbis_first_character_rule() {
    // Web ABI v1 keeps its own strict grammar (frozen).
    for token in [
        "0", "1", "2", "3", "F", "f", "D", "d", "T", "t", "1u", "2u", "3u", "Du", "du", "Tu", "tu",
    ] {
        let parsed = parse_comp_based_stats(token).unwrap();
        assert_eq!(
            parse_comp_based_stats(&parsed.to_ncbi_cli_string()).unwrap(),
            parsed
        );
    }
    for token in [
        "", "2garbage", "2uu", "2uX", "2U", "0u", "Fu", "4", "22", " 2", "2 ",
    ] {
        assert!(parse_comp_based_stats(token).is_err(), "{token}");
    }
    // NCBI blast_args.cpp:834-892: the CLI reads the first character (any other is mode 0)
    // and, for BLASTP, a `u` or `U` second character of an adjusting mode.
    use LOSAT::blastinput::app::{parse_comp_based_stats as ncbi_cbs, CompositionMode};
    for (token, mode, unified) in [
        ("2garbage", CompositionMode::CompositionMatrixAdjust, false),
        ("t", CompositionMode::CompositionMatrixAdjust, false),
        ("Du", CompositionMode::CompositionMatrixAdjust, true),
        ("2U", CompositionMode::CompositionMatrixAdjust, true),
        ("3x", CompositionMode::CompoForceFullMatrixAdjust, false),
        ("1", CompositionMode::CompositionBasedStats, false),
        ("0u", CompositionMode::NoCompositionBasedStats, false),
        ("4", CompositionMode::NoCompositionBasedStats, false),
        ("", CompositionMode::NoCompositionBasedStats, false),
        (" 2", CompositionMode::NoCompositionBasedStats, false),
    ] {
        assert_eq!(
            ncbi_cbs(token, true, false).unwrap(),
            (mode, unified),
            "{token}"
        );
        assert_eq!(
            ncbi_cbs(token, false, false).unwrap(),
            (mode, false),
            "{token}"
        );
        let Commands::Blastp(args) = parse("blastp", &["-comp_based_stats", token])
            .unwrap()
            .command
        else {
            unreachable!("blastp")
        };
        // LOSAT's BLASTP runs mode 2 without unified P-values and mode 0 (the ordinary
        // traceback), and rejects the others.
        let checked = LOSAT::algorithm::blastp::blast_engine::check_options(&args);
        if (mode == CompositionMode::CompositionMatrixAdjust && !unified)
            || mode == CompositionMode::NoCompositionBasedStats
        {
            checked.unwrap();
        } else {
            assert!(
                checked
                    .unwrap_err()
                    .to_string()
                    .contains("not supported by LOSAT's BLASTP"),
                "{token}"
            );
        }
    }
    assert!(ncbi_cbs("2", true, true)
        .unwrap_err()
        .to_string()
        .contains("Composition-adjusted searched are not supported with an ungapped search"));
    ncbi_cbs("F", true, true).unwrap();
}

#[test]
fn defaults_and_task_overrides_remain_distinct() {
    let Commands::Blastp(args) = parse("blastp", &[]).unwrap().command else {
        panic!()
    };
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.evalue, None);
    assert_eq!(args.outfmt, "0");
    assert_eq!(args.max_target_seqs, None);
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
    // NCBI api/blast_options_handle.cpp:381-402: blastp-fast (word size 5, threshold 20)
    // runs; blastp-short (word size 2, PAM30) is rejected when the arguments are parsed.
    let Commands::Blastp(a) = parse("blastp", &["-task", "blastp-fast"]).unwrap().command else {
        panic!()
    };
    let resolved = a.resolve().unwrap();
    assert_eq!((resolved.word_size, resolved.threshold), (5, 20.0));
    LOSAT::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    let error = parse("blastp", &["-task", "blastp-short"]).unwrap_err();
    assert_eq!(error.kind(), clap::error::ErrorKind::ValueValidation);
    assert_eq!(error.exit_code(), 2);
    let text = LOSAT::cli::render_message(&error);
    assert!(
        text.contains("invalid value 'blastp-short' for '-task <TASK>'")
            && text.contains("-word_size 2 (the lookup table) is not supported by LOSAT's BLASTP"),
        "{text}"
    );
    // The search entrance still refuses the task (a defensive check) when it is set by hand.
    let Commands::Blastp(mut a) = parse("blastp", &[]).unwrap().command else {
        panic!()
    };
    a.task = "blastp-short".into();
    assert_eq!(a.resolve().unwrap().evalue, 20000.0);
    assert!(LOSAT::algorithm::blastp::blast_engine::check_options(&a)
        .unwrap_err()
        .to_string()
        .contains("not supported by LOSAT's BLASTP"));
    // TBLASTN takes both of its tasks (tblastn-fast: the compressed lookup table).
    parse("tblastn", &["-task", "tblastn"]).unwrap();
    parse("tblastn", &["-task", "tblastn-fast"]).unwrap();
    // A task that is not an NCBI task keeps clap's list of the accepted names.
    for (program, other, listed) in [
        ("blastp", "blastp-x", "blastp-fast"),
        ("tblastn", "tblastn-x", "tblastn-fast"),
    ] {
        let text = LOSAT::cli::render_message(&parse(program, &["-task", other]).unwrap_err());
        assert!(
            text.contains("possible values") && text.contains(listed),
            "{text}"
        );
        assert!(!text.contains("-short"), "{text}");
    }
    // blast_args.cpp:3158-3163 (CArgAllowStringSet): the task names are case-sensitive.
    assert!(parse("blastp", &["-task", "BLASTP"]).is_err());
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

// NCBI blast_args.cpp:2665-2685: -num_descriptions and -num_alignments are integers of 0 or
// more (`CArgAllowValuesGreaterThanOrEqual(0)`), which NCBI reads with the sign and
// hexadecimal forms; blast_args.cpp:2732-2738: -max_target_seqs excludes both.
#[test]
fn num_descriptions_and_num_alignments_are_integers_of_zero_or_more() {
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        for key in ["-num_descriptions", "-num_alignments"] {
            for value in ["0", "5", "+5", "0x10", "2147483647"] {
                parse(program, &["-outfmt", "6", key, value])
                    .unwrap_or_else(|error| panic!("{program} {key} {value}: {error}"));
            }
            for value in ["-1", "abc", "1e1", "2147483648", ""] {
                assert!(
                    parse(program, &["-outfmt", "6", key, value]).is_err(),
                    "{program} {key} {value}"
                );
            }
            // The exclusion is NCBI's `SetDependency(eExcludes)`.
            let error = parse(program, &["-max_target_seqs", "5", key, "5"])
                .err()
                .unwrap_or_else(|| panic!("{program} -max_target_seqs with {key}"));
            assert_eq!(error.exit_code(), 2, "{program} {key}");
        }
        parse(
            program,
            &[
                "-num_descriptions",
                "3",
                "-num_alignments",
                "2",
                "-outfmt",
                "0",
            ],
        )
        .unwrap();
        parse(program, &["-max_target_seqs", "5"]).unwrap();
    }
    let Commands::Blastp(args) = parse(
        "blastp",
        &["-num_descriptions", "0x10", "-num_alignments", "+5"],
    )
    .unwrap()
    .command
    else {
        unreachable!("blastp")
    };
    assert_eq!(
        (args.num_descriptions, args.num_alignments),
        (Some(16), Some(5))
    );
}

#[test]
fn numeric_values_are_validated_before_io() {
    for program in ["blastn", "blastp", "tblastx"] {
        // NCBI rejects an e-value that does not start with a digit, a point or a sign
        // (`inf`, `nan`, ` 1`) and one with trailing text; BLASTN's other forms are below.
        let evalues = vec!["NaN", "inf", "nan", " 1", "1 ", "1,5", ""];
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
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        assert!(parse(program, &["-outfmt", "6", "-max_hsps", "0"]).is_err());
        parse(program, &["-outfmt", "6", "-max_hsps", "1"]).unwrap();
        // NCBI blast_args.cpp:196-201: -qcov_hsp_perc is a real number in 0..100.
        for value in ["0", "50.5", "1e1", "100"] {
            parse(program, &["-outfmt", "6", "-qcov_hsp_perc", value]).unwrap();
        }
        for value in ["-1", "100.5", "inf", "nan", "abc"] {
            assert!(
                parse(program, &["-outfmt", "6", "-qcov_hsp_perc", value]).is_err(),
                "{program} -qcov_hsp_perc {value}"
            );
        }
    }
    // NCBI tblastx runs word sizes 2 to 4 (blast_options.c:1342-1364), which LOSAT's
    // TBLASTX rejects where the search starts; 5 or more is NCBI's error.
    for (value, message) in [
        ("2", "not supported by LOSAT's TBLASTX"),
        ("4", "not supported by LOSAT's TBLASTX"),
        ("5", "Word-size must be less than 6 for protein comparison"),
    ] {
        let Commands::Tblastx(args) = parse("tblastx", &["-outfmt", "6", "-word_size", value])
            .unwrap()
            .command
        else {
            unreachable!("tblastx")
        };
        let error = LOSAT::algorithm::tblastx::blast_engine::check_options(&args)
            .unwrap_err()
            .to_string();
        assert!(error.contains(message), "{value}: {error}");
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
    for extra in [&["-task", "rmblastn"][..], &["-version"]] {
        let mut words = vec!["-outfmt", "6"];
        words.extend_from_slice(extra);
        let error = LOSAT::cli::render_message(&parse("blastn", &words).unwrap_err());
        assert!(
            error.contains("not supported by LOSAT's BLASTN"),
            "{extra:?}: {error}"
        );
    }
    // `-ungapped` is NCBI's flag (blast_args.cpp:917) and LOSAT's BLASTN searches with it (O4a).
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-ungapped"])
        .unwrap()
        .command
    else {
        unreachable!("blastn")
    };
    assert!(args.ungapped);
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6"]).unwrap().command else {
        unreachable!("blastn")
    };
    assert!(!args.ungapped);
    assert!(parse("blastn", &["-outfmt", "6", "-task", "BLASTN"]).is_err());
    // `-strand` is NCBI's (blast_args.cpp:1953-1955): plus, minus and both are accepted, any
    // other value is a clap error (LOSAT's rule for argument errors, rc 2).
    for (value, want) in [("plus", "plus"), ("minus", "minus"), ("both", "both")] {
        let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-strand", value])
            .unwrap()
            .command
        else {
            unreachable!("blastn")
        };
        assert_eq!(args.strand, want);
    }
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6"]).unwrap().command else {
        unreachable!("blastn")
    };
    assert_eq!(args.strand, "both");
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-strand=minus"])
        .unwrap()
        .command
    else {
        unreachable!("blastn")
    };
    assert_eq!(args.strand, "minus");
    for value in ["foo", "PLUS", "Both", "", "1"] {
        assert!(
            parse("blastn", &["-outfmt", "6", "-strand", value]).is_err(),
            "{value:?}"
        );
    }
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
    for program in ["blastn", "blastp", "tblastx"] {
        assert!(parse(program, &["-outfmt", "6", "-window_size", "2147483648"]).is_err());
        parse(program, &["-outfmt", "6", "-window_size", "0"]).unwrap();
    }
    // TBLASTX checks its options after the files are opened, with NCBI's messages.
    for (words, message) in [
        (&["-threshold", "0"][..], "Non-zero threshold required"),
        (
            &["-evalue", "-inf"],
            "expect value or cutoff score must be greater than zero",
        ),
        (
            &["-evalue", "0"],
            "expect value or cutoff score must be greater than zero",
        ),
        (
            &["-seg", "12 2.2"],
            "Invalid number of arguments to filtering option",
        ),
    ] {
        let Commands::Tblastx(args) = parse("tblastx", words).unwrap().command else {
            unreachable!("tblastx")
        };
        let error = LOSAT::algorithm::tblastx::blast_engine::check_options(&args)
            .unwrap_err()
            .to_string();
        assert!(error.contains(message), "{words:?}: {error}");
    }
    // BLASTP checks its options after the files are opened, with NCBI's messages
    // (blast_options.c BLAST_ValidateOptions).
    for (words, message) in [
        (&["-threshold", "0"][..], "Non-zero threshold required"),
        (
            &["-evalue", "0"],
            "expect value or cutoff score must be greater than zero",
        ),
        (
            &["-evalue", "-inf"],
            "expect value or cutoff score must be greater than zero",
        ),
        (
            &["-word_size", "8"],
            "Word-size must be less than 8 for a tblastn, blastp or blastx search",
        ),
        (
            &["-matrix", "FOO"],
            "FOO is not a supported matrix, supported matrices are:",
        ),
        (
            &["-matrix", "PAM30", "-gapopen", "1"],
            "Gap existence and extension values of 1 and 1 not supported for PAM30",
        ),
        (
            &["-seg", "12 2.2"],
            "Invalid number of arguments to filtering option",
        ),
    ] {
        let Commands::Blastp(args) = parse("blastp", words).unwrap().command else {
            unreachable!("blastp")
        };
        let error = args.resolve().map(|_| ()).unwrap_err().to_string();
        assert!(error.contains(message), "{words:?}: {error}");
    }
    for evalue in ["0x10", "0x1p-3", "1e"] {
        let error = parse("blastp", &["-evalue", evalue])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("not supported by LOSAT's BLASTP"),
            "{evalue}: {error}"
        );
    }
    for value in ["0", "7", "8", "32", "255"] {
        assert!(parse("tblastx", &["-outfmt", "6", "-db_gencode", value]).is_err());
    }
}

#[test]
fn flags_and_value_tokens_are_not_reinterpreted() {
    // NCBI blast_args.cpp:803-806: -use_sw_tback is a flag; LOSAT's BLASTP does not
    // implement the Smith-Waterman traceback and rejects it where the search starts.
    let Commands::Blastp(a) = parse("blastp", &["-use_sw_tback"]).unwrap().command else {
        panic!()
    };
    assert!(LOSAT::algorithm::blastp::blast_engine::check_options(&a)
        .unwrap_err()
        .to_string()
        .contains("not supported by LOSAT's BLASTP"));
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
    assert_eq!(
        a.subject.as_deref().and_then(|path| path.to_str()),
        Some("-query")
    );
    // NCBI's `-` is standard input (blast_args.cpp:3425-3427, ncbiargs.cpp).
    let cli =
        try_parse_from::<Cli, _, _>(["losat", "blastp", "-query", "-", "-subject", "s"]).unwrap();
    let Commands::Blastp(a) = cli.command else {
        panic!()
    };
    assert_eq!(a.query.to_str(), Some("-"));
}

/// What a program does with a `-outfmt` value (`OUTFMT_ROWS`).
#[derive(Debug, Clone, Copy)]
enum Outfmt {
    /// LOSAT writes the format.
    Writes(LOSAT::blastinput::app::ReportFormat),
    /// NCBI's error when it reads the value, before the option handlers: exit status and
    /// standard error (blast_args.cpp:2801-2851).
    Ncbi(i32, &'static str),
    /// NCBI's error after the input and output files are opened (its formatting handler or the
    /// XInclude check): the message of "BLAST query/options error: ...".
    NcbiLater(&'static str),
    /// A format LOSAT writes whose custom specification NCBI's formatting handler then rejects
    /// (`s_ValidateCustomDelim`): the message of "BLAST query/options error: ...".
    WritesThenNcbi(LOSAT::blastinput::app::ReportFormat, &'static str),
    /// LOSAT's rejection before the inputs are read; `{P}` is the program's name.
    Losat(&'static str),
}

const ALL: &[&str] = &["blastn", "blastp", "blastx", "tblastn", "tblastx"];
const NOT_BLASTN: &[&str] = &["blastp", "blastx", "tblastn", "tblastx"];
/// The programs that write every field, a delimiter and outfmt 10 and 20 with the shared tabular
/// writer: all of them (line CLI C2: BLASTN in PR 1, TBLASTN and TBLASTX in PR 2, BLASTP and
/// BLASTX in PR 3).
const SHARED: &[&str] = ALL;
const UNSUPPORTED: Outfmt = Outfmt::Losat("output format {N} is not supported by LOSAT's {P}");
const OUT_OF_RANGE: Outfmt = Outfmt::Ncbi(255, "Error: Formatting choice is out of range\n");
const INVALID_DELIMITER: Outfmt = Outfmt::Ncbi(
    1,
    "BLAST query/options error: Delimiter format is invalid. Valid format is delim=<delimiter value>\nPlease refer to the BLAST+ user manual.\n",
);

/// The programs, `-outfmt` value, `-out` file and result of each case. A change that makes a
/// program write a format edits or deletes the rows of that format only; the rows of a line OF
/// track are grouped. `{N}` in a LOSAT message is the format number.
#[rustfmt::skip]
const OUTFMT_ROWS: &[(&[&str], &str, Option<&str>, Outfmt)] = {
    use LOSAT::blastinput::app::ReportFormat::*;
    use Outfmt::*;
    &[
        // Written today by every program; NCBI drops the specification of a format that is
        // not tabular (blast_args.cpp:2844-2851) and reads the number as `NStr::StringToInt`.
        (ALL, "0", None, Writes(Pairwise)),
        (ALL, "6", None, Writes(Tabular)),
        (ALL, "7", None, Writes(TabularWithComments)),
        (ALL, " 6\t", None, Writes(Tabular)),
        (ALL, "+6", None, Writes(Tabular)),
        (ALL, "07", None, Writes(TabularWithComments)),
        (ALL, "0 qseqid", None, Writes(Pairwise)),
        (ALL, "0 std", None, Writes(Pairwise)),
        (ALL, "0 delim=,", None, Writes(Pairwise)),
        (ALL, "6 delim=", None, Writes(Tabular)),
        (ALL, "7 delim= ", None, Writes(TabularWithComments)),
        // Custom fields (every program writes every field and a delimiter with the shared
        // tabular writer, line CLI C2) and delimiters.
        (SHARED, "6 std", None, Writes(Tabular)),
        (SHARED, "07 std", None, Writes(TabularWithComments)),
        (SHARED, "6 qseqid sseqid pident length", None, Writes(Tabular)),
        (SHARED, "7  qaccver", None, Writes(TabularWithComments)),
        (SHARED, "6 delim=,", None, Writes(Tabular)),
        (SHARED, "6 delimiter=;", None, Writes(Tabular)),
        (SHARED, "6 delim=, std", None, Writes(Tabular)),
        (SHARED, "7 delim=| qseqid sseqid", None, Writes(TabularWithComments)),
        (SHARED, "6 delim=; qseqid", None, Writes(Tabular)),
        // NCBI's check of a custom delimiter against the fields that separate their entries with
        // it, searched as text (blast_args.cpp:2753-2792).
        (SHARED, "6 delim=; qseqid staxids", None, WritesThenNcbi(Tabular, "Your custom record separator (;) is also used by the format specifier (staxids) to separate multiple entries. Please use a different record separator (delim keyword).")),
        (SHARED, "6 delim=; qseqid xstaxidsx", None, WritesThenNcbi(Tabular, "Your custom record separator (;) is also used by the format specifier (staxids) to separate multiple entries. Please use a different record separator (delim keyword).")),
        (SHARED, "7 delim=/ frames", None, WritesThenNcbi(TabularWithComments, "Your custom record separator (/) is also used by the format specifier (frames) to separate multiple entries. Please use a different record separator (delim keyword).")),
        (SHARED, "10 delim=<> salltitles", None, WritesThenNcbi(CommaSeparatedValues, "Your custom record separator (<>) is also used by the format specifier (salltitles) to separate multiple entries. Please use a different record separator (delim keyword).")),
        // NCBI's errors when it reads the value (the delimiter is read before the number).
        (ALL, "abc", None, Ncbi(1, "BLAST query/options error: 'abc' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")),
        (ALL, "0x6", None, Ncbi(1, "BLAST query/options error: '0x6' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")),
        (ALL, "6\u{a0}", None, Ncbi(1, "BLAST query/options error: '6\u{a0}' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")),
        (ALL, "6\tqseqid", None, Ncbi(1, "BLAST query/options error: '6\tqseqid' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")),
        (ALL, " ", None, Ncbi(1, "BLAST query/options error: '' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")),
        (ALL, "-1", None, OUT_OF_RANGE),
        (ALL, "22", None, OUT_OF_RANGE),
        (ALL, "99", None, OUT_OF_RANGE),
        (ALL, "6 delim", None, INVALID_DELIMITER),
        (ALL, "0 delim", None, INVALID_DELIMITER),
        (ALL, "abc delim", None, INVALID_DELIMITER),
        (ALL, "99 delimiter qaccver", None, INVALID_DELIMITER),
        (ALL, "17 delim SQ", None, INVALID_DELIMITER),
        // NCBI's errors after the inputs are opened: the formats of other programs, and the
        // XInclude formats written to standard output.
        (NOT_BLASTN, "17", None, NcbiLater("SAM format is only applicable to blastn")),
        (NOT_BLASTN, "17 qseqid", None, NcbiLater("SAM format is only applicable to blastn")),
        (ALL, "19", None, NcbiLater("AIRR rearrangement format is only applicable to igblastn")),
        (ALL, "21", None, NcbiLater("FASTA output format is only applicable to magicblast")),
        (ALL, "13", None, NcbiLater("Please provide a file name for outfmt 13.")),
        (ALL, "13", Some("-"), NcbiLater("Please provide a file name for outfmt 13.")),
        (ALL, "14", None, NcbiLater("Please provide a file name for outfmt 14.")),
        // Line OF track T (text formats): not written yet.
        (ALL, "1", None, UNSUPPORTED),
        (ALL, "2", None, UNSUPPORTED),
        (ALL, "3", None, UNSUPPORTED),
        (ALL, "4", None, UNSUPPORTED),
        (ALL, "18", None, UNSUPPORTED),
        // Line CLI C2: BLASTN (PR 1), TBLASTN and TBLASTX (PR 2), BLASTP and BLASTX (PR 3) write 10
        // and 20 with the shared tabular writer.
        (SHARED, "10", None, Writes(CommaSeparatedValues)),
        (SHARED, "10 delim=; qseqid sseqid", None, Writes(CommaSeparatedValues)),
        (SHARED, "20", None, Writes(CommaSeparatedValuesWithHeader)),
        (SHARED, " 20 qseqid sseqid qcovus ", None, Writes(CommaSeparatedValuesWithHeader)),
        // Line OF track S (structured formats): 8 and 12 in every program (OF3b), 9 (OF3c), 5 in
        // every program (BLASTN and BLASTP: OF4a; BLASTX, TBLASTN and TBLASTX: OF4b), 15 and 16 in every program (BLASTP and BLASTN: OF5a; BLASTX, TBLASTN and TBLASTX: OF5b), NCBI dropping a custom specification
        // (blast_args.cpp:2844-2851); the others not yet.
        (ALL, "8", None, Writes(AsnText)),
        (ALL, "12", None, Writes(JsonSeqalign)),
        (ALL, " 8 std", None, Writes(AsnText)),
        (ALL, "12 delim=,", None, Writes(JsonSeqalign)),
        (ALL, "5", None, Writes(Xml)),
        (ALL, " 5 std", None, Writes(Xml)),
        (ALL, "9", None, Writes(AsnBinary)),
        (ALL, " 9 std", None, Writes(AsnBinary)),
        // 11 (the BLAST archive, OF8a and OF8b): every program; NCBI ignores a specification
        // after it.
        (ALL, "11", None, Writes(ArchiveFormat)),
        (ALL, " 11 qseqid", None, Writes(ArchiveFormat)),
        // 13 and 14 (BLAST JSON and XML2 as side files and a master file at `-out`): every
        // program (BLASTP and BLASTN: OF5c part 1; BLASTX, TBLASTN and TBLASTX: part 2).
        (ALL, "13", Some("o.json"), Writes(Json)),
        (ALL, "14", Some("o.xml"), Writes(Xml2)),
        (ALL, " 14 std", Some("o"), Writes(Xml2)),
        (ALL, "15", None, Writes(JsonS)),
        (ALL, "16", None, Writes(Xml2S)),
        (ALL, " 16 std", None, Writes(Xml2S)),
        // SAM (OF6): BLASTN alone writes it, reading the words SQ and SR and ignoring the rest and a
        // delimiter (`CBlast_SAM_Formatter::x_ProcessCustomSpec`).
        (&["blastn"], "17", None, Writes(Sam)),
        (&["blastn"], "17 qseqid", None, Writes(Sam)),
        (&["blastn"], " 17 SQ SR ", None, Writes(Sam)),
        (&["blastn"], "17 delim=, SQ", None, Writes(Sam)),
    ]
};

fn program_of(name: &str) -> LOSAT::blastinput::app::Program {
    use LOSAT::blastinput::app::Program;
    match name {
        "blastn" => Program::Blastn,
        "blastp" => Program::Blastp,
        "blastx" => Program::Blastx,
        "tblastn" => Program::Tblastn,
        "tblastx" => Program::Tblastx,
        _ => unreachable!("{name}"),
    }
}

// Every program reads -outfmt as NCBI when it sets the options (blast_args.cpp:2801-2851,
// `blastinput/app.rs`): the argument parser accepts any value, and the format's report, NCBI's
// errors and LOSAT's rejections follow `OUTFMT_ROWS`.
#[test]
fn output_formats_follow_the_format_table() {
    use LOSAT::blastinput::app::{
        formatting_handler_check, parse_formatting_string, report_format, xinclude_check,
    };
    let native = |error: anyhow::Error| match error.downcast::<LOSAT::cli::NativeError>() {
        Ok(native) => (native.exit, native.message),
        Err(other) => panic!("not NCBI's error: {other}"),
    };
    for &(programs, value, out, expected) in OUTFMT_ROWS {
        for &name in programs {
            let case = format!("{name} -outfmt {value:?} -out {out:?}");
            let mut words = vec!["-outfmt", value];
            if let Some(out) = out {
                words.extend(["-out", out]);
            }
            parse(name, &words).unwrap_or_else(|error| panic!("{case}: {error}"));
            let program = program_of(name);
            let choice = match (parse_formatting_string(value), expected) {
                (Err(error), Outfmt::Ncbi(exit, message)) => {
                    assert_eq!(native(error), (exit, message.to_string()), "{case}");
                    continue;
                }
                (Ok(choice), Outfmt::Ncbi(..)) => panic!("{case}: read as {choice:?}"),
                (Ok(choice), _) => choice,
                (Err(error), _) => panic!("{case}: {error:#}"),
            };
            let report = report_format(&choice, program, out.map(std::path::Path::new));
            match expected {
                Outfmt::Writes(format) => assert_eq!(report.unwrap(), Some(format), "{case}"),
                Outfmt::Losat(message) => assert_eq!(
                    report.unwrap_err().to_string(),
                    message
                        .replace("{P}", program.name())
                        .replace("{N}", &choice.number.to_string()),
                    "{case}"
                ),
                Outfmt::NcbiLater(message) => {
                    assert_eq!(report.unwrap(), None, "{case}");
                    let error = formatting_handler_check(&choice, program)
                        .and_then(|()| xinclude_check(&choice))
                        .unwrap_err();
                    assert_eq!(
                        native(error),
                        (
                            1,
                            format!("BLAST query/options error: {message}\nPlease refer to the BLAST+ user manual.\n")
                        ),
                        "{case}"
                    );
                }
                Outfmt::WritesThenNcbi(format, message) => {
                    assert_eq!(report.unwrap(), Some(format), "{case}");
                    let error = formatting_handler_check(&choice, program).unwrap_err();
                    assert_eq!(
                        native(error),
                        (
                            1,
                            format!("BLAST query/options error: {message}\nPlease refer to the BLAST+ user manual.\n")
                        ),
                        "{case}"
                    );
                }
                Outfmt::Ncbi(..) => unreachable!(),
            }
        }
    }
}

#[test]
fn output_capabilities_fail_explicitly_and_help_is_canonical() {
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
            // NCBI's BLASTP tasks and flags are options (blast_args.cpp:803-806,3158-3163).
            if program == "blastp" {
                for option in ["blastp-fast", "-use_sw_tback"] {
                    assert!(text.contains(option), "{text}");
                }
                assert!(!text.contains("blastp-short"), "{text}");
            }
        }
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3292-3333
// ```c++
//     arg_desc.AddDefaultKey(kArgCullingLimit, "int_value",
//                            "If the query range of a hit is enveloped by that of at least this many "
//                            "higher-scoring hits, delete the hit",
//                            CArgDescriptions::eInteger,
//                            NStr::IntToString(kDfltArgCullingLimit));
//     arg_desc.SetConstraint(kArgCullingLimit, new CArgAllowValuesGreaterThanOrEqual(0));
// ```
// BLASTN, BLASTP and TBLASTN read `-culling_limit` as a non-negative integer (a hexadecimal
// value too, as NCBI's argument parser does).
#[test]
fn culling_limit_is_read_as_ncbi_does() {
    for program in ["blastn", "blastp", "tblastn"] {
        let limit = |words: &[&str]| -> Result<i32, String> {
            match parse(program, words) {
                Ok(cli) => match cli.command {
                    Commands::Blastn(args) => Ok(args.culling_limit_value()),
                    Commands::Blastp(args) => Ok(args.culling_limit_value()),
                    Commands::Tblastn(args) => Ok(args.culling_limit_value()),
                    _ => unreachable!("{program}"),
                },
                Err(error) => Err(LOSAT::cli::render_message(&error)),
            }
        };
        assert_eq!(limit(&[]), Ok(0), "{program}");
        assert_eq!(limit(&["-culling_limit", "0"]), Ok(0), "{program}");
        assert_eq!(limit(&["-culling_limit", "1"]), Ok(1), "{program}");
        assert_eq!(limit(&["-culling_limit", "0x2"]), Ok(2), "{program}");
        assert_eq!(
            limit(&["-culling_limit", "2147483647"]),
            Ok(2147483647),
            "{program}"
        );
        for bad in ["-1", "1.5", "abc", "2147483648"] {
            let error = limit(&["-culling_limit", bad]).unwrap_err();
            assert!(
                !error.contains("is not supported by"),
                "{program} {bad}: {error}"
            );
        }
    }
}

// NCBI blast_args.cpp:2663: -show_gis is a flag of every search program.
#[test]
fn show_gis_is_a_flag_of_every_search_program() {
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        assert!(parse(program, &["-show_gis"]).is_ok(), "{program}");
        assert!(parse(program, &["-show_gis", "true"]).is_err(), "{program}");
    }
    let Commands::Blastn(args) = parse("blastn", &["-show_gis"]).unwrap().command else {
        panic!()
    };
    assert!(args.show_gis);
    let Commands::Blastn(args) = parse("blastn", &[]).unwrap().command else {
        panic!()
    };
    assert!(!args.show_gis);
}

// NCBI blastp 2.17.0+ -help options that LOSAT's BLASTP does not implement, and the NCBI C++
// Toolkit's standard options (ncbiargs.cpp:78-87), are rejected explicitly.
#[test]
fn unported_ncbi_options_are_rejected_explicitly() {
    for (program, words) in [
        ("blastn", &["-db_soft_mask", "11"][..]),
        ("blastn", &["-db_hard_mask", "11"]),
        ("tblastn", &["-db_soft_mask", "11"]),
        ("tblastx", &["-db_soft_mask", "11"]),
        ("blastp", &["-gilist", "ids.gil"]),
        ("blastp", &["-seqidlist", "ids.txt"]),
        ("blastp", &["-negative_seqidlist", "ids.txt"]),
        // Line DI-8 ported the taxid options to BLASTP's database search, not to PSIBLAST's.
        ("psiblast", &["-taxids", "9606"]),
        ("psiblast", &["-no_taxid_expansion"]),
        ("blastp", &["-db_soft_mask", "21"]),
        ("blastp", &["-db_hard_mask", "21"]),
        ("blastp", &["-lcase_masking"]),
        ("blastp", &["-version"]),
        ("blastp", &["-xdrop_gap_final", "30"]),
        ("blastn", &["-logfile", "log"]),
        ("blastp", &["-xmlhelp"]),
        ("tblastn", &["-version-full"]),
        ("tblastx", &["-help-full"]),
        // NCBI ncbiargs.cpp:88: -dryrun, which NCBI accepts without searching.
        ("blastp", &["-dryrun"]),
        ("tblastn", &["-dryrun"]),
        ("tblastx", &["-dryrun"]),
        ("blastn", &["-dryrun"]),
    ] {
        let error = LOSAT::cli::render_message(&parse(program, words).unwrap_err());
        assert!(
            error.contains(&format!(
                "is not supported by LOSAT's {}",
                program.to_ascii_uppercase()
            )),
            "{program} {words:?}: {error}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
// ```c++
//     arg_desc.SetCurrentGroup("General search options");
//     // database filename
//     if (m_IsIgBlast){
//         arg_desc.AddOptionalKey(kArgDb, "database_name", "Optional additional database name",
//                                 CArgDescriptions::eString);
//     } else {
//         arg_desc.AddOptionalKey(kArgDb, "database_name", "BLAST database name",
//                                 CArgDescriptions::eString);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
// ```c++
//     // DB size
//     if (!m_IsMapper) {
//         arg_desc.SetCurrentGroup("Statistical options");
//         arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
//                                 "Effective length of the database ",
//                                 CArgDescriptions::eInt8);
//     }
// ```
// BLASTP's `-db NAME` and `-dbsize N`: `-db` excludes `-subject` and `-subject_loc` (a
// conflict of LOSAT's parser, where NCBI prints its usage), `-dbsize` is an Int8 read with
// either of them.
#[test]
fn blastp_takes_a_database_and_its_size() {
    let blastp = |words: &[&str]| {
        let mut argv = vec!["losat", "blastp", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv).map(|cli| match cli.command {
            Commands::Blastp(args) => args,
            _ => unreachable!("blastp"),
        })
    };
    let args = blastp(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    // NCBI keeps the name without the white space at its ends (blast_args.cpp:2473).
    let resolved = args.resolve().unwrap();
    assert_eq!(resolved.db.as_deref(), Some("dir/name"));
    assert_eq!(resolved.dbsize, 3_000_000_000);
    // An empty name is a name (NCBI's "Database name is required." comes from the open).
    assert_eq!(blastp(&["-db", ""]).unwrap().db.as_deref(), Some(""));
    assert_eq!(
        blastp(&["-db", "  "])
            .unwrap()
            .resolve()
            .unwrap()
            .db
            .as_deref(),
        Some("")
    );
    assert_eq!(blastp(&["-db=a b"]).unwrap().db.as_deref(), Some("a b"));
    // `-dbsize` is read with `-subject` too, and takes NCBI's Int8 forms.
    for (word, value) in [
        ("0", 0),
        ("-5", -5),
        ("+7", 7),
        ("0x10", 16),
        ("9223372036854775807", i64::MAX),
        ("-9223372036854775808", i64::MIN),
    ] {
        let args = blastp(&["-subject", "s.fa", "-dbsize", word]).unwrap();
        assert_eq!(args.dbsize, Some(value), "{word}");
        assert_eq!(args.resolve().unwrap().dbsize, value, "{word}");
    }
    // Without `-dbsize` the options keep a database length of 0.
    let resolved = blastp(&["-subject", "s.fa"]).unwrap().resolve().unwrap();
    assert_eq!((resolved.db, resolved.dbsize), (None, 0));
    for word in ["abc", "1.5", "", "9223372036854775808", "1e3"] {
        let error = blastp(&["-db", "d", "-dbsize", word]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ValueValidation,
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = blastp(words).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{words:?}: {error}"
        );
        let text = LOSAT::cli::render_message(&error);
        assert!(text.contains("'-db <NAME>'"), "{text}");
    }
    // LOSAT Web has no database files: its validation (`check_options`) rejects `-db`.
    let error =
        LOSAT::algorithm::blastp::blast_engine::check_options(&blastp(&["-db", "d"]).unwrap())
            .unwrap_err();
    assert!(
        format!("{error:#}").contains("-db (a database search) is not supported by LOSAT Web yet"),
        "{error:#}"
    );
}

// BLASTN's `-db NAME` and `-dbsize N` (blast_args.cpp:2191-2199 and 2231-2237, quoted
// above): as BLASTP's, read by the handler of the database arguments of every program
// (blast_args.cpp:2472-2476 and 2570-2572).
#[test]
fn blastn_takes_a_database_and_its_size() {
    let blastn = |words: &[&str]| {
        let mut argv = vec!["losat", "blastn", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv).map(|cli| match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn"),
        })
    };
    let args = blastn(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    // NCBI keeps the name without the white space at its ends (blast_args.cpp:2473).
    assert_eq!(args.database_name(), Some("dir/name"));
    assert_eq!(args.db_length(), 3_000_000_000);
    assert_eq!(blastn(&["-db", "  "]).unwrap().database_name(), Some(""));
    // `-dbsize` is read with `-subject` too, and takes NCBI's Int8 forms.
    for (word, value) in [
        ("0", 0),
        ("-5", -5),
        ("0x10", 16),
        ("9223372036854775807", i64::MAX),
    ] {
        let args = blastn(&["-subject", "s.fa", "-dbsize", word]).unwrap();
        assert_eq!(args.db_length(), value, "{word}");
    }
    // Without `-dbsize` the options keep a database length of 0.
    let args = blastn(&["-subject", "s.fa"]).unwrap();
    assert_eq!((args.database_name(), args.db_length()), (None, 0));
    for word in ["abc", "1.5", "", "9223372036854775808"] {
        let error = blastn(&["-db", "d", "-dbsize", word]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ValueValidation,
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = blastn(words).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{words:?}: {error}"
        );
    }
    // LOSAT Web has no database files: its validation (`check_scoring`) rejects `-db`.
    let error = LOSAT::algorithm::blastn::scoring::check_scoring(&blastn(&["-db", "d"]).unwrap())
        .unwrap_err();
    assert!(
        format!("{error:#}").contains("-db (a database search) is not supported by LOSAT Web yet"),
        "{error:#}"
    );
}

// TBLASTN's `-db NAME` and `-dbsize N` (blast_args.cpp:2191-2199 and 2231-2237, quoted
// above), as BLASTN's.
#[test]
fn tblastn_takes_a_database_and_its_size() {
    let tblastn = |words: &[&str]| {
        let mut argv = vec!["losat", "tblastn", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv).map(|cli| match cli.command {
            Commands::Tblastn(args) => args,
            _ => unreachable!("tblastn"),
        })
    };
    let args = tblastn(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    assert_eq!(args.database_name(), Some("dir/name"));
    for (word, value) in [("-5", -5), ("0x10", 16), ("9223372036854775807", i64::MAX)] {
        let args = tblastn(&["-subject", "s.fa", "-dbsize", word]).unwrap();
        assert_eq!(args.dbsize, Some(value), "{word}");
    }
    for word in ["abc", "1.5", "", "9223372036854775808"] {
        let error = tblastn(&["-db", "d", "-dbsize", word]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ValueValidation,
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = tblastn(words).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{words:?}: {error}"
        );
    }
    // LOSAT Web has no database files: its validation (`check_options`) rejects `-db`.
    let error =
        LOSAT::algorithm::tblastn::check_options(&tblastn(&["-db", "d"]).unwrap()).unwrap_err();
    assert!(
        format!("{error:#}").contains("-db (a database search) is not supported by LOSAT Web yet"),
        "{error:#}"
    );
}

// TBLASTX's `-db NAME` and `-dbsize N`, as TBLASTN's: `-db` excludes `-subject` and
// `-subject_loc`, `-dbsize` is an Int8 read with either of them.
#[test]
fn tblastx_takes_a_database_and_its_size() {
    let tblastx = |words: &[&str]| {
        let mut argv = vec!["losat", "tblastx", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv).map(|cli| match cli.command {
            Commands::Tblastx(args) => args,
            _ => unreachable!("tblastx"),
        })
    };
    let args = tblastx(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    assert_eq!(args.database_name(), Some("dir/name"));
    for (word, value) in [("-5", -5), ("0x10", 16), ("9223372036854775807", i64::MAX)] {
        let args = tblastx(&["-subject", "s.fa", "-dbsize", word]).unwrap();
        assert_eq!(args.dbsize, Some(value), "{word}");
    }
    for word in ["abc", "1.5", "", "9223372036854775808"] {
        let error = tblastx(&["-db", "d", "-dbsize", word]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ValueValidation,
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = tblastx(words).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{words:?}: {error}"
        );
    }
    // LOSAT Web has no database files: its validation (`check_options`) rejects `-db`.
    let error =
        LOSAT::algorithm::tblastx::blast_engine::check_options(&tblastx(&["-db", "d"]).unwrap())
            .unwrap_err();
    assert!(
        format!("{error:#}").contains("-db (a database search) is not supported by LOSAT Web yet"),
        "{error:#}"
    );
}

#[test]
fn values_on_which_ncbi_crashes_are_rejected_with_the_options() {
    // Decision D12 (AUTHORITY.md): an -evalue of DBL_MAX or more crashes NCBI's blastp and
    // tblastn on some inputs; blast_hits.c:43-70: a -max_target_seqs whose preliminary size wraps to a
    // negative int crashes NCBI.
    let blastp = |words: &[&str]| {
        let Commands::Blastp(args) = parse("blastp", words).unwrap().command else {
            unreachable!("blastp")
        };
        LOSAT::algorithm::blastp::blast_engine::check_options(&args)
    };
    let tblastn = |words: &[&str]| {
        let Commands::Tblastn(args) = parse("tblastn", words).unwrap().command else {
            unreachable!("tblastn")
        };
        LOSAT::algorithm::tblastn::check_options(&args)
    };
    for words in [
        &["-evalue", "+inf"][..],
        &["-evalue", "1e999"],
        &["-evalue", "1.7976931348623157e308"],
        &["-max_target_seqs", "1073741799"],
    ] {
        let error = format!("{:#}", blastp(words).unwrap_err());
        assert!(
            error.contains("not supported by LOSAT's BLASTP"),
            "{words:?}: {error}"
        );
    }
    for words in [
        &["-evalue", "1e308"][..],
        &["-evalue", "1.7976931348623156e308"],
        &["-max_target_seqs", "2147483647"],
        &["-max_target_seqs", "1073741798"],
    ] {
        blastp(words).unwrap();
    }
    let error = format!("{:#}", tblastn(&["-evalue", "+inf"]).unwrap_err());
    assert!(
        error.contains("not supported by LOSAT's TBLASTN"),
        "{error}"
    );
    let error = format!(
        "{:#}",
        tblastn(&["-evalue", "1.7976931348623157e308"]).unwrap_err()
    );
    assert!(
        error.contains("not supported by LOSAT's TBLASTN"),
        "{error}"
    );
    tblastn(&["-evalue", "1e308"]).unwrap();
}

#[test]
fn toolkit_words_in_an_option_value_are_rejected_as_options() {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001: the application reads
    // -version, -version-full*, -dryrun, -logfile and -conffile anywhere before `--`, also
    // in an option's value (`-out -version` prints the version). LOSAT rejects them (D9)
    // in that place too, instead of writing a file named "-version". The parser refuses
    // `-version` (the Web engine's); the command line prints it before parsing
    // (`ncbi_version_request_reads_the_word_before_the_arguments`).
    for (program, upper) in [
        ("blastn", "BLASTN"),
        ("blastp", "BLASTP"),
        ("tblastn", "TBLASTN"),
        ("tblastx", "TBLASTX"),
    ] {
        for words in [
            &["-out", "-version"][..],
            &["-out", "-dryrun"],
            &["-evalue", "-version-full"],
            &["-out", "-logfile"],
            &["-out", "-conffile=x.ini"],
        ] {
            let error = parse(program, words).unwrap_err().to_string();
            assert!(
                error.contains(&format!("is not supported by LOSAT's {upper}")),
                "{program} {words:?}: {error}"
            );
        }
        assert!(parse(program, &["-out", "o.txt"]).is_ok(), "{program}");
    }
}

#[test]
fn a_last_double_dash_changes_nothing_and_a_word_after_it_is_an_error() {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2866-2872: `--` starts the
    // positional arguments, which the BLAST programs do not have.
    for program in ["blastn", "blastp", "tblastn", "tblastx"] {
        assert!(parse(program, &["-evalue", "1", "--"]).is_ok(), "{program}");
        assert!(
            parse(program, &["--", "-outfmt", "6"]).is_err(),
            "{program}"
        );
        assert!(parse(program, &["--", "--"]).is_err(), "{program}");
        // As the value of an option, `--` is that value.
        assert!(parse(program, &["-out", "--"]).is_ok(), "{program}");
    }
    // makeblastdb reads its arguments with the same parser: it builds the database
    // when `--` is the last word.
    let makedb = |words: &[&str]| {
        let mut argv = vec!["losat", "makeblastdb", "-dbtype", "prot"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv)
    };
    assert!(makedb(&["-in", "x.faa", "--"]).is_ok());
    assert!(makedb(&["--", "-in", "x.faa"]).is_err());
    assert!(makedb(&["--", "--"]).is_err());
    assert!(makedb(&["-out", "--"]).is_ok());
}

// BLASTX's `-db NAME` and `-dbsize N` (stage DB-5): `-db` excludes `-subject` (a conflict of
// LOSAT's parser, where NCBI prints its usage), `-dbsize` is an Int8 read with `-db` and with
// `-subject` and becomes the options' database length (`SetDbLength`).
#[test]
fn blastx_takes_a_database_and_its_size() {
    let blastx = |words: &[&str]| {
        let mut argv = vec!["losat", "blastx", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv).map(|cli| match cli.command {
            Commands::Blastx(args) => args,
            _ => unreachable!("blastx"),
        })
    };
    let args = blastx(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    assert_eq!(args.database_name(), Some("dir/name"));
    assert_eq!(args.resolve().unwrap().db_length, 3_000_000_000);
    for (word, value) in [("-5", -5), ("0x10", 16), ("9223372036854775807", i64::MAX)] {
        let args = blastx(&["-subject", "s.fa", "-dbsize", word]).unwrap();
        assert_eq!(args.dbsize, Some(value), "{word}");
        assert_eq!(args.resolve().unwrap().db_length, value, "{word}");
    }
    assert_eq!(
        blastx(&["-subject", "s.fa"])
            .unwrap()
            .resolve()
            .unwrap()
            .db_length,
        0
    );
    for word in ["abc", "1.5", "", "9223372036854775808"] {
        let error = blastx(&["-db", "d", "-dbsize", word]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ValueValidation,
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
    ] {
        let error = blastx(words).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::ArgumentConflict,
            "{words:?}: {error}"
        );
    }
    // `-subject_loc` stays outside BLASTX's scope, with `-db` too.
    let error = blastx(&["-db", "d", "-subject_loc", "1-10"]).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("unsupported BLASTX option '-subject_loc'"),
        "{error}"
    );
}

#[test]
fn makedb_is_not_a_command_and_makeblastdb_checks_apply() {
    let makedb = |words: &[&str]| {
        let mut argv = vec!["losat", "makeblastdb", "-dbtype", "prot"];
        argv.extend(words);
        try_parse_from::<Cli, _, _>(argv)
    };
    // The old name is gone: an unknown command, with every spelling of its arguments.
    for words in [
        &["makedb", "-in", "x.faa", "-dbtype", "prot"][..],
        &["makedb", "-h"],
        &["makedb", "-version"],
    ] {
        let mut argv = vec!["losat"];
        argv.extend(words);
        let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::UnknownArgument,
            "{words:?}"
        );
        assert_ne!(error.exit_code(), 0);
        assert!(error.to_string().contains("makedb"), "{words:?}: {error}");
    }
    let Commands::Makedb(args) = makedb(&["-in", "x.faa", "-out", "db"]).unwrap().command else {
        unreachable!("makeblastdb")
    };
    assert_eq!(
        (args.input.as_str(), args.out.as_deref()),
        ("x.faa", Some("db"))
    );
    // The checks that name the command apply.
    let error = makedb(&["-mask_data"]).unwrap_err().to_string();
    assert!(
        error.contains("the NCBI BLAST+ option -mask_data is not supported by LOSAT's makeblastdb"),
        "{error}"
    );
    assert!(makedb(&["-parse_seqids"]).is_ok());
    assert!(makedb(&["-in", "x.faa", "--"]).is_ok());
    assert!(makedb(&["--", "-in", "x.faa"]).is_err());
    assert!(makedb(&["-evalue", "1"]).is_err());
    let error = makedb(&["-h"]).unwrap_err();
    assert_eq!(error.kind(), clap::error::ErrorKind::DisplayHelp);
    assert_eq!(error.exit_code(), 0);
}

#[test]
fn short_help_h_is_display_help_for_every_program() {
    for program in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "makeblastdb",
    ] {
        let error = try_parse_from::<Cli, _, _>(["losat", program, "-h"]).unwrap_err();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::DisplayHelp,
            "{program}"
        );
        assert_eq!(error.exit_code(), 0);
        // As an option's value, `-h` is the value.
        assert!(
            parse(program, &["-out", "-h"]).is_ok() || program == "makeblastdb",
            "{program}"
        );
    }
}

#[test]
fn ncbi_version_request_reads_the_word_before_the_arguments() {
    use LOSAT::cli::ncbi_version_request;
    let version = |words: &[&str]| ncbi_version_request(words.iter().copied());
    let text = version(&["losat", "blastn", "-version"]).unwrap();
    assert!(
        text.starts_with("blastn: 2.17.0+\n Package: LOSAT "),
        "{text}"
    );
    assert!(
        text.ends_with(", output compatible with blast 2.17.0\n"),
        "{text}"
    );
    assert!(version(&["losat", "makeblastdb", "-x", "-version"])
        .unwrap()
        .starts_with("makeblastdb: 2.17.0+\n"));
    // Where another option expects its value, and next to invalid arguments.
    assert!(version(&["losat", "blastx", "-out", "-version"]).is_some());
    assert!(version(&["losat", "blastp", "-evalue", "x", "-version", "-bogus"]).is_some());
    // Not a request: after `--`, a longer word, no program, the top level.
    assert!(version(&["losat", "blastn", "--", "-version"]).is_none());
    assert!(version(&["losat", "blastn", "-version-full"]).is_none());
    assert!(version(&["losat", "blastn", "-version=1"]).is_none());
    assert!(version(&["losat", "-version"]).is_none());
    assert!(version(&["losat", "--version"]).is_none());
    assert!(version(&["losat", "nosuch", "-version"]).is_none());
    assert!(version(&["losat"]).is_none());
}

#[test]
fn mt_mode_takes_zero_to_two_and_needs_num_threads() {
    for program in ["blastn", "blastp", "blastx", "tblastn"] {
        for mode in ["0", "1", "2"] {
            assert!(
                parse(program, &["-mt_mode", mode, "-num_threads", "2"]).is_ok(),
                "{program} {mode}"
            );
            // clap's `requires` is not met by the default of `-num_threads`.
            let error = parse(program, &["-mt_mode", mode]).unwrap_err().to_string();
            assert!(error.contains("num_threads"), "{program} {mode}: {error}");
        }
        for mode in ["3", "-1", "x", "1.5", ""] {
            assert!(
                parse(program, &["-mt_mode", mode, "-num_threads", "2"]).is_err(),
                "{program} {mode:?}"
            );
        }
    }
    // NCBI's tblastx and psiblast have no -mt_mode.
    for program in ["tblastx", "psiblast"] {
        let error = parse(program, &["-mt_mode", "1", "-num_threads", "2"])
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("unknown option or argument"),
            "{program}: {error}"
        );
    }
}

#[test]
fn use_index_false_parses_and_true_is_an_unsupported_error() {
    for value in ["false", "False", "f", "no", "0"] {
        assert!(parse("blastn", &["-use_index", value]).is_ok(), "{value}");
    }
    for value in ["true", "T", "yes", "1"] {
        let error =
            LOSAT::cli::render_message(&parse("blastn", &["-use_index", value]).unwrap_err());
        assert!(
            error.contains(
                "-use_index true (the megablast index) is not supported by LOSAT's BLASTN"
            ),
            "{value}: {error}"
        );
    }
    assert!(parse("blastn", &["-use_index", "maybe"]).is_err());
    // NCBI's other programs have no -use_index; `-index_name` stays unported.
    assert!(parse("blastp", &["-use_index", "false"]).is_err());
    let error = parse("blastn", &["-index_name", "x"])
        .unwrap_err()
        .to_string();
    assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
}
