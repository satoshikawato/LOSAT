//! Frozen CLI v2 grammar and configuration boundary.
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-207,332-420,838-883,3158-3163

use crate::algorithm::blastp::args::parse_comp_based_stats;
use crate::blastinput::ncbi_environment::DiagTarget;
use crate::blastinput::value_parsers::{
    parse_dust_filtering, parse_seg_filtering, DustSpec, SegSpec,
};
use crate::cli::{try_parse_from, ArgvError, Cli, Commands};

fn parse_argv(argv: Vec<&str>) -> Result<Cli, ArgvError> {
    try_parse_from(argv)
}

fn parse(program: &str, extra: &[&str]) -> Result<Cli, ArgvError> {
    let mut argv = vec!["losat", program, "-query", "q.fa", "-subject", "s.fa"];
    argv.extend(extra);
    try_parse_from(argv)
}

/// The exit status of the command line whose reading ended with `error`: 0 after a printed help,
/// 1 after an error of NCBI's argument reader or a refusal of LOSAT.
fn status(error: &ArgvError) -> i32 {
    match error {
        ArgvError::Printed(_) => 0,
        ArgvError::Report(_) | ArgvError::Refused(_) => 1,
    }
}

/// The error that NCBI's argument reader writes for `error` (its usage, the line of `=` and its
/// posts), which must be one.
fn report(error: &ArgvError) -> String {
    match error {
        ArgvError::Report(bytes) => String::from_utf8_lossy(bytes).into_owned(),
        other => panic!("not an error of NCBI's argument reader: {other:?}"),
    }
}

/// The text that `-help` (or `-h`, `-help-full`) printed, which `error` must be.
fn printed(error: &ArgvError) -> &str {
    match error {
        ArgvError::Printed(text) => text,
        other => panic!("not a printed help: {other:?}"),
    }
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
    for program in [
        "blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast",
    ] {
        for words in lists {
            let mut argv = vec!["losat", program, "-query", "q.fa", "-db", "db"];
            argv.extend(words);
            assert!(parse_argv(argv).is_ok(), "{program} -db {words:?}");
            let error = (parse(program, words).unwrap_err()).to_string();
            assert!(
                error.contains("Incompatible with argument:  `"),
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2211-2221,2287,2289,2302-2319,2364-2370
// Line DI-8b: `-seqidlist` and `-negative_seqidlist` go with `-db` only, exclude each other,
// every taxid list and `-no_taxid_expansion` (NCBI's argument reader's exclusions, exit status 1,
// as the taxid lists' are).
#[test]
fn seqid_list_options_go_with_a_database_only() {
    let lists: [&[&str]; 2] = [
        &["-seqidlist", "ids.txt"],
        &["-negative_seqidlist", "ids.txt"],
    ];
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        for words in lists {
            let mut argv = vec!["losat", program, "-query", "q.fa", "-db", "db"];
            argv.extend(words);
            assert!(parse_argv(argv).is_ok(), "{program} -db {words:?}");
            let error = (parse(program, words).unwrap_err()).to_string();
            assert!(
                error.contains("Incompatible with argument:  `"),
                "{program} {words:?}: {error}"
            );
            for other in [
                &["-negative_seqidlist", "o.txt"][..],
                &["-seqidlist", "o.txt"],
                &["-taxids", "9606"],
                &["-taxidlist", "t.txt"],
                &["-negative_taxids", "9606"],
                &["-negative_taxidlist", "t.txt"],
                &["-no_taxid_expansion"],
            ] {
                if other[0] == words[0] {
                    continue;
                }
                let mut argv = vec!["losat", program, "-query", "q.fa", "-db", "db"];
                argv.extend(words);
                argv.extend(other);
                assert!(parse_argv(argv).is_err(), "{program} {words:?} {other:?}");
            }
        }
    }
}

#[test]
fn canonical_options_for_every_program() {
    for (program, word) in [
        ("blastn", "11"),
        ("blastp", "3"),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:54-55
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:55-62
            // blast_args.cpp:2726-2731 registers the optional `kArgMaxTargetSequences` key.
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
            Commands::Rpsblast(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, Some(0.001));
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Rpstblastn(a) => {
                assert_eq!(a.num_threads, 4);
                assert_eq!(a.evalue, Some(0.001));
                assert_eq!(a.max_target_seqs, Some(9));
            }
            Commands::Makedb(_) | Commands::Rpsbproc(_) => {
                unreachable!("{program} is a search program")
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:225-246
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:293-294
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:304-305
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
    assert!(!defaults.randomize);

    // The hidden flag `-randomize` is read (its check against the input type is the build's).
    assert!(args(makedb(&["-dbtype", "nucl", "-randomize"]).unwrap()).randomize);

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
    // `-taxid` excludes it (makeblastdb.cpp:246-247, 322-329). NCBI's argument reader reports
    // them with makeblastdb's usage, exit status 1.
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
        assert_eq!(status(&err), 1, "{words:?}: {err}");
        assert!(
            !err.to_string().contains("not supported"),
            "{words:?}: {err}"
        );
    }
    // Without the mandatory `-dbtype`, NCBI's reader reports an exclusion first (the first walk
    // of `x_PostCheck`), and a missing required argument only after the mandatory one, whose
    // error is the parser's.
    let err = makedb(&[
        "-in",
        "x.fa",
        "-taxid",
        "1",
        "-taxid_map",
        "m",
        "-parse_seqids",
    ])
    .unwrap_err();
    assert_eq!(status(&err), 1, "{err}");
    assert!(
        err.to_string()
            .contains("Argument \"taxid_map\". Incompatible with argument:  `taxid'"),
        "{err}"
    );
    let err = makedb(&["-in", "x.fa", "-taxid_map", "m"]).unwrap_err();
    assert!(
        report(&err).contains("Argument \"dbtype\". Mandatory value is missing"),
        "{err}"
    );

    // `-mask_data` takes a list of files, `-mask_id` needs `-mask_data` and `-mask_desc` needs
    // `-mask_id` (makeblastdb.cpp:259-277); the reader reports them as those of `-taxid_map`
    // above.
    assert_eq!(defaults.mask_data, None);
    assert_eq!((defaults.mask_id, defaults.mask_desc), (None, None));
    let given = args(
        makedb(&[
            "-dbtype",
            "prot",
            "-mask_data",
            "a.asnt,b.asnb",
            "-mask_id",
            "x,y",
            "-mask_desc",
            "d1,d2",
        ])
        .unwrap(),
    );
    assert_eq!(given.mask_data.as_deref(), Some("a.asnt,b.asnb"));
    assert_eq!(given.mask_id.as_deref(), Some("x,y"));
    assert_eq!(given.mask_desc.as_deref(), Some("d1,d2"));
    for words in [
        &["-dbtype", "prot", "-mask_id", "x"][..],
        &["-dbtype", "prot", "-mask_data", "a", "-mask_desc", "d"],
        &["-dbtype", "prot", "-mask_desc", "d", "-mask_id", "x"],
        &["-dbtype", "prot", "-mask_data"],
    ] {
        let err = makedb(words).unwrap_err();
        assert_eq!(status(&err), 1, "{words:?}: {err}");
        assert!(
            !err.to_string().contains("not supported"),
            "{words:?}: {err}"
        );
    }

    // `-gi_mask` needs `-parse_seqids` and excludes `-mask_id`; `-gi_mask_name` needs
    // `-mask_data` and `-gi_mask` (makeblastdb.cpp:279-288): the reader's errors, exit status 1.
    assert!(!defaults.gi_mask);
    assert_eq!(defaults.gi_mask_name, None);
    let given = args(
        makedb(&[
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-mask_data",
            "a.asnb,b.asnb",
            "-gi_mask",
            "-gi_mask_name",
            "x,y",
        ])
        .unwrap(),
    );
    assert!(given.gi_mask);
    assert_eq!(given.gi_mask_name.as_deref(), Some("x,y"));
    for words in [
        &["-dbtype", "nucl", "-gi_mask"][..],
        &[
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-gi_mask",
            "-mask_data",
            "a",
            "-mask_id",
            "x",
        ],
        &["-dbtype", "nucl", "-gi_mask_name", "x"],
        &[
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-gi_mask",
            "-gi_mask_name",
            "x",
        ],
        &["-dbtype", "nucl", "-mask_data", "a", "-gi_mask_name", "x"],
        &[
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-mask_data",
            "a",
            "-gi_mask",
            "-gi_mask_name",
        ],
    ] {
        let err = makedb(words).unwrap_err();
        assert_eq!(status(&err), 1, "{words:?}: {err}");
        assert!(
            !err.to_string().contains("not supported"),
            "{words:?}: {err}"
        );
    }
    // The reader names the first dependency that NCBI's `x_PostCheck` meets, also where an
    // argument that one requires is excluded by another (NCBI BLAST+ 2.17.0:
    // `tests/seqdb/gi_mask/expected/e03.stderr`, `e04.stderr`).
    for (words, sentence) in [
        (
            &[
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-mask_data",
                "a",
                "-gi_mask",
                "-mask_desc",
                "d",
            ][..],
            "Argument \"mask_id\". Must be specified, as it is required by argument:  `mask_desc'",
        ),
        (
            &[
                "-dbtype",
                "nucl",
                "-mask_data",
                "a",
                "-mask_id",
                "y",
                "-gi_mask_name",
                "x",
            ],
            "Argument \"gi_mask\". Must be specified, as it is required by argument:  \
             `gi_mask_name'",
        ),
        (
            &[
                "-dbtype",
                "prot",
                "-taxid",
                "1",
                "-taxid_map",
                "m",
                "-parse_seqids",
            ],
            "Argument \"taxid_map\". Incompatible with argument:  `taxid'",
        ),
    ] {
        let err = makedb(words).unwrap_err();
        assert_eq!(status(&err), 1, "{words:?}: {err}");
        assert!(err.to_string().contains(sentence), "{words:?}: {err}");
    }

    // `-help` prints makeblastdb's help, NCBI's text; `--help` is not one of its words.
    let err = makedb(&["-help"]).unwrap_err();
    let text = printed(&err);
    assert!(text.contains(" -dbtype <String, `nucl', `prot'>"), "{text}");
    assert!(text.contains(" -max_file_sz <String>"), "{text}");
    assert!(
        !text.contains("--in") && !text.contains("--dbtype"),
        "{text}"
    );
    assert_eq!(status(&makedb(&["--help"]).unwrap_err()), 1);

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
            // rule 2). `-db` takes a value, which is missing here.
            if old == "-db" {
                assert!(
                    report(&err).contains("Argument \"-db\". Value is missing"),
                    "{err}"
                );
                continue;
            } else if old == "-remote" {
                assert_eq!(
                    err,
                    ArgvError::Refused(format!(
                        "the NCBI BLAST+ option -remote is not supported by LOSAT's {}",
                        program.to_ascii_uppercase()
                    ))
                );
                continue;
            }
            let text = report(&err);
            assert!(
                text.contains("Unknown argument") || text.contains("Invalid argument"),
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1147-1154: a value that is not above 0 keeps the default
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:375-384,410-426: single spaces separate three signed numbers;
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:834-892: the CLI reads the first character (any other is mode 0)
    // and, for BLASTP, a `u` or `U` second character of an adjusting mode.
    use crate::blastinput::app::{parse_comp_based_stats as ncbi_cbs, CompositionMode};
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
        // LOSAT's BLASTP runs every mode (modes 1 and 3 since series SC, SC5), with unified
        // P-values too (item X5 of line RG).
        crate::algorithm::blastp::blast_engine::check_options(&args).unwrap();
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:381-387
    // Creates the default `eBlastp` options handle for the chosen locality.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:381-402: blastp-fast (word size 5, threshold 20) and
    // blastp-short (word size 2, PAM30 9/1 with its threshold 16 and window 15, e-value 20000)
    // run, and so does -word_size 2 of each task (item SC3).
    let Commands::Blastp(a) = parse("blastp", &["-task", "blastp-fast"]).unwrap().command else {
        panic!()
    };
    let resolved = a.resolve().unwrap();
    assert_eq!((resolved.word_size, resolved.threshold), (5, 20.0));
    crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    let Commands::Blastp(a) = parse("blastp", &["-task", "blastp-short"]).unwrap().command else {
        panic!()
    };
    let resolved = a.resolve().unwrap();
    assert_eq!(
        (
            resolved.word_size,
            resolved.threshold,
            resolved.window_size,
            resolved.evalue,
            (resolved.scoring.gap_open, resolved.scoring.gap_extend)
        ),
        (2, 16.0, 15, 20000.0, (9, 1))
    );
    crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    for (task, threshold) in [("blastp", 11.0), ("blastp-fast", 20.0)] {
        let Commands::Blastp(a) = parse("blastp", &["-task", task, "-word_size", "2"])
            .unwrap()
            .command
        else {
            panic!()
        };
        let resolved = a.resolve().unwrap();
        assert_eq!(
            (resolved.word_size, resolved.threshold),
            (2, threshold),
            "{task}"
        );
        crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    }
    // Word size 4 keeps the lookup table with each task's threshold, and runs (item SC8b).
    for (task, threshold) in [
        ("blastp", 11.0),
        ("blastp-fast", 20.0),
        ("blastp-short", 16.0),
    ] {
        let Commands::Blastp(a) = parse("blastp", &["-task", task, "-word_size", "4"])
            .unwrap()
            .command
        else {
            panic!()
        };
        let resolved = a.resolve().unwrap();
        assert_eq!(
            (resolved.word_size, resolved.threshold),
            (4, threshold),
            "{task}"
        );
        crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    }
    // TBLASTN takes both of its tasks (tblastn-fast: the compressed lookup table).
    parse("tblastn", &["-task", "tblastn"]).unwrap();
    parse("tblastn", &["-task", "tblastn-fast"]).unwrap();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:429-433: blastx-fast is word size 5 (the compressed
    // lookup table, lookup type 4) with threshold 20; -word_size 5 sets its own threshold over
    // the task's, -word_size 3 keeps the task's on the standard table
    // (blast_args.cpp:291-297, blast_options_local_priv.hpp:612-619).
    for (extra, word_size, threshold, lookup_type) in [
        (&["-task", "blastx"][..], 3, 12.0, 3),
        (&["-task", "blastx-fast"][..], 5, 20.0, 4),
        (&["-word_size", "5"][..], 5, 19.3, 4),
        (&["-task", "blastx-fast", "-word_size", "5"][..], 5, 19.3, 4),
        (&["-task", "blastx-fast", "-word_size", "3"][..], 3, 20.0, 3),
        (
            &["-task", "blastx-fast", "-threshold", "15"][..],
            5,
            15.0,
            4,
        ),
        (
            &[
                "-task",
                "blastx-fast",
                "-word_size",
                "3",
                "-threshold",
                "14",
            ][..],
            3,
            14.0,
            3,
        ),
    ] {
        let Commands::Blastx(a) = parse("blastx", extra).unwrap().command else {
            panic!()
        };
        let resolved = a.resolve().unwrap();
        assert_eq!(
            (resolved.word_size, resolved.threshold, resolved.lookup_type),
            (word_size, threshold, lookup_type),
            "{extra:?}"
        );
    }
    // A task that is not an NCBI task gets NCBI's list of the permissible names.
    for (program, other, listed) in [
        ("blastp", "blastp-x", &["blastp-fast", "blastp-short"][..]),
        ("blastx", "blastx-x", &["blastx-fast"][..]),
        ("tblastn", "tblastn-x", &["tblastn-fast"][..]),
    ] {
        let text = report(&parse(program, &["-task", other]).unwrap_err());
        assert!(
            text.contains("Permissible values") && listed.iter().all(|task| text.contains(task)),
            "{text}"
        );
        if program != "blastp" {
            assert!(!text.contains("-short"), "{text}");
        }
    }
    // blast_args.cpp:105 (CArgAllowStringSet): the task names are case-sensitive.
    assert!(parse("blastp", &["-task", "BLASTP"]).is_err());
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803: the default -outfmt is 0, which BLASTN and TBLASTX
    // implement.
    parse("blastn", &[]).unwrap();
    parse("tblastx", &[]).unwrap();
    for program in ["blastn", "tblastx"] {
        parse(program, &["-outfmt", "6"]).unwrap();
    }
    // NCBI blast_options_handle.cpp:344-370: dc-megablast, blastn-short (Session SD) and rmblastn
    // (blastn with reward 0 and penalty 0, line RG X6) are tasks of blastn.
    parse("blastn", &["-outfmt", "6", "-task", "dc-megablast"]).unwrap();
    parse("blastn", &["-outfmt", "6", "-task", "blastn-short"]).unwrap();
    parse("blastn", &["-outfmt", "6", "-task", "rmblastn"]).unwrap();
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685: -num_descriptions and -num_alignments are integers of 0 or
// more (`CArgAllowValuesGreaterThanOrEqual` with 0), which NCBI reads with the sign and
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
            // The exclusion is NCBI's `SetDependency` with `eExcludes`.
            let error = parse(program, &["-max_target_seqs", "5", key, "5"])
                .err()
                .unwrap_or_else(|| panic!("{program} -max_target_seqs with {key}"));
            assert_eq!(status(&error), 1, "{program} {key}");
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201: -qcov_hsp_perc is a real number in 0..100.
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
    // NCBI tblastx runs word sizes 2 to 4 (c++/src/algo/blast/core/blast_options.c:1335-1362), as LOSAT's TBLASTX
    // does (series SC, SC10); 5 or more is NCBI's error.
    for (value, message) in [
        ("2", None),
        ("3", None),
        ("4", None),
        (
            "5",
            Some("Word-size must be less than 6 for protein comparison"),
        ),
    ] {
        let Commands::Tblastx(args) = parse("tblastx", &["-outfmt", "6", "-word_size", value])
            .unwrap()
            .command
        else {
            unreachable!("tblastx")
        };
        let checked = crate::algorithm::tblastx::blast_engine::check_options(&args);
        match message {
            None => assert!(checked.is_ok(), "{value}: {checked:?}"),
            Some(message) => {
                let error = checked.unwrap_err().to_string();
                assert!(error.contains(message), "{value}: {error}");
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:168-170: the argument is 4 or more; c++/src/algo/blast/core/blast_options.c:1326-1333:
    // an option check rejects more than 100, with NCBI's message.
    assert!(parse("blastn", &["-outfmt", "6", "-word_size", "3"]).is_err());
    // NCBI's `CArg_Double` reads glibc's hexadecimal numbers and a final `e`
    // (`StringToDoublePosix`); NCBI 2.17.0 searches with each (line SW, SW8).
    for (evalue, read) in [("0x10", 16.0), ("0x1p-3", 0.125), ("1e", 1.0)] {
        let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-evalue", evalue])
            .unwrap()
            .command
        else {
            unreachable!("blastn")
        };
        assert_eq!(args.evalue, Some(read), "-evalue {evalue}");
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
        let checked = crate::algorithm::blastn::scoring::check_scoring_options(&args);
        let limited = crate::algorithm::blastn::scoring::check_losat_limits(&args);
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
    // NCBI's integer arguments (c++/src/corelib/ncbiargs.cpp:118-129) read hexadecimal after `0x`.
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
    // NCBI's rmblastn task is blastn with reward 0 and penalty 0 (matrix-only scoring,
    // blast_options_handle.cpp:365-370); every task of NCBI's blastn parses.
    let Commands::Blastn(args) = parse("blastn", &["-outfmt", "6", "-task", "rmblastn"])
        .unwrap()
        .command
    else {
        unreachable!("blastn")
    };
    assert_eq!(args.task, "rmblastn");
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
    // other value is NCBI's reader's error (exit status 1).
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87 (CDirEntry::GetName, c++/src/corelib/ncbifile.cpp:298-312,358-363,
    // 465-472): the name after the trailing separators are removed.
    let dotted = format!("{}/.", "o".repeat(256));
    parse("blastn", &["-outfmt", "6", "-out", &dotted]).unwrap();
    assert!(parse(
        "blastn",
        &["-outfmt", "6", "-out", &format!("{}/", "o".repeat(256))]
    )
    .is_err());
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87: an -out file name is shorter than 256 bytes.
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
        crate::algorithm::blastn::scoring::check_scoring_options(&args)
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
        let checked = crate::algorithm::blastn::scoring::check_scoring_options(&args);
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
    // `-min_raw_gapped_score` is BLASTN's only (blast_args.cpp:701-704, `eInteger`, no constraint).
    for value in ["15", "0", "-5", "2147483647", "-2147483648"] {
        parse("blastn", &["-min_raw_gapped_score", value]).unwrap();
    }
    for value in ["2147483648", "1.5", "abc"] {
        assert!(parse("blastn", &["-min_raw_gapped_score", value]).is_err());
    }
    for program in ["blastp", "blastx", "tblastn", "tblastx"] {
        assert!(parse(program, &["-min_raw_gapped_score", "15"]).is_err());
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
        let error = crate::algorithm::tblastx::blast_engine::check_options(&args)
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
    for (evalue, read) in [("0x10", 16.0), ("0x1p-3", 0.125), ("1e", 1.0)] {
        let Commands::Blastp(args) = parse("blastp", &["-evalue", evalue]).unwrap().command else {
            unreachable!("blastp")
        };
        assert_eq!(args.evalue, Some(read), "{evalue}");
    }
    for value in ["0", "7", "8", "32", "255"] {
        assert!(parse("tblastx", &["-outfmt", "6", "-db_gencode", value]).is_err());
    }
}

#[test]
fn flags_and_value_tokens_are_not_reinterpreted() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806: -use_sw_tback is a flag; LOSAT's BLASTP runs the
    // Smith-Waterman redo with composition-based statistics (item X4) and without them
    // (`-comp_based_stats 0`, `-matrix IDENTITY`, item X4-c0).
    let Commands::Blastp(a) = parse("blastp", &["-use_sw_tback"]).unwrap().command else {
        panic!()
    };
    assert!(a.use_sw_tback);
    crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    let Commands::Blastp(a) = parse("blastp", &["-use_sw_tback", "-comp_based_stats", "0"])
        .unwrap()
        .command
    else {
        panic!()
    };
    crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    // IDENTITY's reset of the default mode 2 (`identity_cbs_reset`); `-task blastp-short` keeps
    // the default mode 2 (it changes the matrix, the gap costs, the threshold, the word size and the
    // filter, `blast_options_handle.cpp:389-395`) and runs.
    let Commands::Blastp(a) = parse("blastp", &["-use_sw_tback", "-matrix", "IDENTITY"])
        .unwrap()
        .command
    else {
        panic!()
    };
    crate::algorithm::blastp::blast_engine::check_options(&a).unwrap();
    let Commands::Blastp(a) = parse("blastp", &["-use_sw_tback", "-task", "blastp-short"])
        .unwrap()
        .command
    else {
        panic!()
    };
    assert!(crate::algorithm::blastp::blast_engine::check_options(&a).is_ok());
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
    // NCBI's `-` is standard input (c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427, c++/src/corelib/ncbiargs.cpp).
    let cli = try_parse_from(["losat", "blastp", "-query", "-", "-subject", "s"]).unwrap();
    let Commands::Blastp(a) = cli.command else {
        panic!()
    };
    assert_eq!(a.query.to_str(), Some("-"));
}

/// What a program does with a `-outfmt` value (`OUTFMT_ROWS`).
#[derive(Debug, Clone, Copy)]
enum Outfmt {
    /// LOSAT writes the format.
    Writes(crate::blastinput::app::ReportFormat),
    /// NCBI's error when it reads the value, before the option handlers: exit status and
    /// standard error (blast_args.cpp:2801-2851).
    Ncbi(i32, &'static str),
    /// NCBI's error after the input and output files are opened (its formatting handler or the
    /// XInclude check): the message of "BLAST query/options error: ...".
    NcbiLater(&'static str),
    /// A format LOSAT writes whose custom specification NCBI's formatting handler then rejects
    /// (`s_ValidateCustomDelim`): the message of "BLAST query/options error: ...".
    WritesThenNcbi(crate::blastinput::app::ReportFormat, &'static str),
}

const ALL: &[&str] = &["blastn", "blastp", "blastx", "tblastn", "tblastx"];
const NOT_BLASTN: &[&str] = &["blastp", "blastx", "tblastn", "tblastx"];
/// The programs that write every field, a delimiter and outfmt 10 and 20 with the shared tabular
/// writer: all of them (line CLI C2: BLASTN in PR 1, TBLASTN and TBLASTX in PR 2, BLASTP and
/// BLASTX in PR 3).
const SHARED: &[&str] = ALL;
const OUT_OF_RANGE: Outfmt = Outfmt::Ncbi(255, "Error: Formatting choice is out of range\n");
const INVALID_DELIMITER: Outfmt = Outfmt::Ncbi(
    1,
    "LOSAT query/options error: Delimiter format is invalid. Valid format is delim=<delimiter value>\nPlease refer to the LOSAT user manual.\n",
);

/// The programs, `-outfmt` value, `-out` file and result of each case. A change that makes a
/// program write a format edits or deletes the rows of that format only; the rows of a line OF
/// track are grouped. `{N}` in a LOSAT message is the format number.
#[rustfmt::skip]
const OUTFMT_ROWS: &[(&[&str], &str, Option<&str>, Outfmt)] = {
    use crate::blastinput::app::ReportFormat::*;
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
        (ALL, "abc", None, Ncbi(1, "LOSAT query/options error: 'abc' is not a valid output format\nPlease refer to the LOSAT user manual.\n")),
        (ALL, "0x6", None, Ncbi(1, "LOSAT query/options error: '0x6' is not a valid output format\nPlease refer to the LOSAT user manual.\n")),
        (ALL, "6\u{a0}", None, Ncbi(1, "LOSAT query/options error: '6\u{a0}' is not a valid output format\nPlease refer to the LOSAT user manual.\n")),
        (ALL, "6\tqseqid", None, Ncbi(1, "LOSAT query/options error: '6\tqseqid' is not a valid output format\nPlease refer to the LOSAT user manual.\n")),
        (ALL, " ", None, Ncbi(1, "LOSAT query/options error: '' is not a valid output format\nPlease refer to the LOSAT user manual.\n")),
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
        // Line RG T3 (the query-anchored formats): BLASTP (T3a), BLASTN and TBLASTN (T3b), BLASTX
        // and TBLASTX (T3c) write them.
        (ALL, "1", None, Writes(QueryAnchoredIdentities)),
        (ALL, "2", None, Writes(QueryAnchoredNoIdentities)),
        (ALL, "3", None, Writes(FlatQueryAnchoredIdentities)),
        (ALL, "4", None, Writes(FlatQueryAnchoredNoIdentities)),
        (ALL, "1 qseqid", None, Writes(QueryAnchoredIdentities)),
        // Line RG T2: the organism report in every program.
        (SHARED, "18", None, Writes(TaxFormat)),
        (SHARED, "18 qseqid", None, Writes(TaxFormat)),
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

fn program_of(name: &str) -> crate::blastinput::app::Program {
    use crate::blastinput::app::Program;
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
    use crate::blastinput::app::{
        formatting_handler_check, parse_formatting_string, report_format, xinclude_check,
    };
    let native = |error: anyhow::Error| match error.downcast::<crate::cli::NativeError>() {
        Ok(native) => (native.exit, String::from_utf8(native.message).unwrap()),
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
                Outfmt::NcbiLater(message) => {
                    assert_eq!(report.unwrap(), None, "{case}");
                    let error = formatting_handler_check(&choice, program)
                        .and_then(|()| xinclude_check(&choice))
                        .unwrap_err();
                    assert_eq!(
                        native(error),
                        (
                            1,
                            format!("LOSAT query/options error: {message}\nPlease refer to the LOSAT user manual.\n")
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
                            format!("LOSAT query/options error: {message}\nPlease refer to the LOSAT user manual.\n")
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
        let err = try_parse_from(["losat", program, "-help"]).unwrap_err();
        let text = printed(&err);
        assert!(text.contains("\n -query <File_In>\n"), "{text}");
        assert!(text.contains("\n -num_threads <Integer, >=1>\n"), "{text}");
        assert!(!text.contains("--query") && !text.contains("--word") && !text.contains("--num"));
        assert!(!text.contains("seg_window") && !text.contains("dust_window"));
        // NCBI's BLASTP tasks and flags are options (c++/src/algo/blast/blastinput/blast_args.cpp:804-806,3158-3163).
        if program == "blastp" {
            for option in ["blastp-fast", "blastp-short", "-use_sw_tback"] {
                assert!(text.contains(option), "{text}");
            }
        }
        // `--help` is not one of NCBI's words.
        assert_eq!(
            status(&try_parse_from(["losat", program, "--help"]).unwrap_err()),
            1
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3302,3338-3340
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
                Err(error) => Err(error.to_string()),
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1961: -parse_deflines is a flag of every search program.
#[test]
fn parse_deflines_is_a_flag_of_every_search_program() {
    for program in [
        "blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast",
    ] {
        assert!(parse(program, &["-parse_deflines"]).is_ok(), "{program}");
        assert!(
            parse(program, &["-parse_deflines", "true"]).is_err(),
            "{program}"
        );
    }
    let Commands::Blastn(args) = parse("blastn", &["-parse_deflines"]).unwrap().command else {
        panic!("blastn");
    };
    assert!(args.parse_deflines);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663: -show_gis is a flag of every search program.
#[test]
fn show_gis_is_a_flag_of_every_search_program() {
    for program in [
        "blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast",
    ] {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2692-2693
// `-html` is a flag of the five search programs (item T4a2 of line RG); psiblast keeps its refusal.
#[test]
fn html_is_a_flag_of_the_five_search_programs_and_psiblast_refuses_it() {
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        assert!(parse(program, &["-html"]).is_ok(), "{program}");
        assert!(parse(program, &["-html", "true"]).is_err(), "{program}");
        let given = match parse(program, &["-html"]).unwrap().command {
            Commands::Blastn(args) => args.html,
            Commands::Blastp(args) => args.html,
            Commands::Blastx(args) => args.html,
            Commands::Tblastn(args) => args.html,
            Commands::Tblastx(args) => args.html,
            _ => panic!("{program}"),
        };
        let absent = match parse(program, &[]).unwrap().command {
            Commands::Blastn(args) => args.html,
            Commands::Blastp(args) => args.html,
            Commands::Blastx(args) => args.html,
            Commands::Tblastn(args) => args.html,
            Commands::Tblastx(args) => args.html,
            _ => panic!("{program}"),
        };
        assert!(given && !absent, "{program}");
    }
    let error = parse("psiblast", &["-html"]).unwrap_err().to_string();
    assert!(
        error.contains("the NCBI BLAST+ option -html is not supported by LOSAT's PSIBLAST"),
        "{error}"
    );
}

// NCBI blastp 2.17.0+ -help options that LOSAT's BLASTP does not implement, and the NCBI C++
// Toolkit's standard options (ncbiargs.cpp:78-87), are rejected explicitly.
#[test]
fn unported_ncbi_options_are_rejected_explicitly() {
    for (program, words) in [
        ("tblastn", &["-import_search_strategy", "s.asn"][..]),
        ("tblastx", &["-import_search_strategy", "s.asn"]),
        ("blastp", &["-import_search_strategy", "s.asn"]),
    ] {
        // With a database (`-import_search_strategy` remains refused).
        let mut argv = vec!["losat", program, "-query", "q.fa", "-db", "db"];
        argv.extend(words);
        let error = (parse_argv(argv).unwrap_err()).to_string();
        assert!(
            error.contains(&format!(
                "is not supported by LOSAT's {}",
                program.to_ascii_uppercase()
            )),
            "{program} {words:?}: {error}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2342-2356 (cited at
// `blastinput::db_mask::DbMaskArgs`)
// BLASTN's (series NP, M4a-2), BLASTP's, BLASTX's (series NP, M4b), TBLASTN's and TBLASTX's (M4b
// PR 2) `-db_soft_mask` and `-db_hard_mask`: strings, each excluding the other, `-subject` and
// `-subject_loc` (NCBI's reader's exclusions; probes o04, o05, o06 of the M4a-2 survey); neither
// is refused (`ncbi_app::args::refusal`, series NP, M4a-2b).
#[test]
fn the_nucleotide_and_protein_database_programs_take_the_database_masks() {
    let masks = |program: &str, words: &[&str]| {
        let mut argv = vec!["losat", program, "-query", "q.fa"];
        argv.extend(words);
        parse_argv(argv).map(|cli| match cli.command {
            Commands::Blastn(args) => args.db_mask,
            Commands::Blastp(args) => args.db_mask,
            Commands::Blastx(args) => args.db_mask,
            Commands::Tblastn(args) => args.db_mask,
            Commands::Tblastx(args) => args.db_mask,
            _ => panic!("{program}"),
        })
    };
    for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
        let found = masks(program, &["-db", "d", "-db_soft_mask", "+11"]).unwrap();
        assert_eq!(found.db_soft_mask.as_deref(), Some("+11"), "{program}");
        assert_eq!(found.db_hard_mask, None, "{program}");
        let found = masks(program, &["-db", "d", "-db_hard_mask", "seg"]).unwrap();
        assert_eq!(found.db_hard_mask.as_deref(), Some("seg"), "{program}");
        for words in [
            &["-db", "d", "-db_soft_mask", "11", "-db_hard_mask", "11"][..],
            &["-subject", "s.fa", "-db_soft_mask", "11"],
            &["-subject", "s.fa", "-db_hard_mask", "11"],
            &["-subject_loc", "1-5", "-db_soft_mask", "11"],
            &["-subject_loc", "1-5", "-db_hard_mask", "11"],
        ] {
            let error = masks(program, words).unwrap_err();
            assert!(
                report(&error).contains("Incompatible with argument"),
                "{program} {words:?}"
            );
        }
        for name in ["db_hard_mask", "db_soft_mask"] {
            let refusal = crate::ncbi_app::args::refusal(program, name);
            assert_eq!(refusal, None, "{program} -{name}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:61-62,
// blastx_args.cpp:62-63, psiblast_args.cpp:58, rpsblast_args.cpp:86, rpstblastn_args.cpp:61
// Only blastp's and blastx's database arguments turn the database masking on; NCBI's psiblast,
// rpsblast and rpstblastn leave it off, so `-db_soft_mask` and `-db_hard_mask` are NCBI's unknown arguments there (exit
// 1, `Unknown argument`): LOSAT's parser's unknown option, not a refusal by name (series NP,
// M4b).
#[test]
fn psiblast_rpsblast_and_rpstblastn_have_no_database_masks() {
    for program in ["psiblast", "rpsblast", "rpstblastn"] {
        for flag in ["-db_soft_mask", "-db_hard_mask"] {
            let argv = vec!["losat", program, "-query", "q.fa", "-db", "d", flag, "21"];
            let error = parse_argv(argv).unwrap_err();
            assert!(
                report(&error).contains("Unknown argument"),
                "{program} {flag}"
            );
            let text = error.to_string();
            assert!(
                !text.contains("is not supported by LOSAT"),
                "{program} {flag}: {text}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
// BLASTP's `-db NAME` and `-dbsize N`: `-db` excludes `-subject` and `-subject_loc` (a
// conflict of LOSAT's parser, where NCBI prints its usage), `-dbsize` is an Int8 read with
// either of them.
#[test]
fn blastp_takes_a_database_and_its_size() {
    let blastp = |words: &[&str]| {
        let mut argv = vec!["losat", "blastp", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from(argv).map(|cli| match cli.command {
            Commands::Blastp(args) => args,
            _ => unreachable!("blastp"),
        })
    };
    let args = blastp(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    // NCBI keeps the name without the white space at its ends (c++/src/algo/blast/blastinput/blast_args.cpp:2473).
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
        assert!(
            report(&error).contains("Argument cannot be converted"),
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = blastp(words).unwrap_err();
        assert!(
            report(&error).contains("Incompatible with argument"),
            "{words:?}: {error}"
        );
        let text = error.to_string();
        assert!(
            text.contains("Argument \"db\".") || text.contains("`db'"),
            "{text}"
        );
    }
    // LOSAT Web's validation (`check_options`) checks `-db` as the command line does: the run opens
    // the database from the files that the host gives (line WQA D12a).
    crate::algorithm::blastp::blast_engine::check_options(&blastp(&["-db", "d"]).unwrap()).unwrap();
}

// BLASTN's `-db NAME` and `-dbsize N` (blast_args.cpp:2191-2199 and 2231-2237, quoted
// above): as BLASTP's, read by the handler of the database arguments of every program
// (blast_args.cpp:2472-2476 and 2570-2572).
#[test]
fn blastn_takes_a_database_and_its_size() {
    let blastn = |words: &[&str]| {
        let mut argv = vec!["losat", "blastn", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from(argv).map(|cli| match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn"),
        })
    };
    let args = blastn(&["-db", " dir/name\t", "-dbsize", "3000000000"]).unwrap();
    assert_eq!(args.db.as_deref(), Some(" dir/name\t"));
    assert_eq!(args.subject, None);
    assert_eq!(args.dbsize, Some(3_000_000_000));
    // NCBI keeps the name without the white space at its ends (c++/src/algo/blast/blastinput/blast_args.cpp:2473).
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
        assert!(
            report(&error).contains("Argument cannot be converted"),
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = blastn(words).unwrap_err();
        assert!(
            report(&error).contains("Incompatible with argument"),
            "{words:?}: {error}"
        );
    }
}

// TBLASTN's `-db NAME` and `-dbsize N` (blast_args.cpp:2191-2199 and 2231-2237, quoted
// above), as BLASTN's.
#[test]
fn tblastn_takes_a_database_and_its_size() {
    let tblastn = |words: &[&str]| {
        let mut argv = vec!["losat", "tblastn", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from(argv).map(|cli| match cli.command {
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
        assert!(
            report(&error).contains("Argument cannot be converted"),
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = tblastn(words).unwrap_err();
        assert!(
            report(&error).contains("Incompatible with argument"),
            "{words:?}: {error}"
        );
    }
    // LOSAT Web's validation (`check_options`) checks `-db` as the command line does: the run opens
    // the database from the files that the host gives (line WQA D12a).
    crate::algorithm::tblastn::check_options(&tblastn(&["-db", "d"]).unwrap()).unwrap();
}

// TBLASTX's `-db NAME` and `-dbsize N`, as TBLASTN's: `-db` excludes `-subject` and
// `-subject_loc`, `-dbsize` is an Int8 read with either of them.
#[test]
fn tblastx_takes_a_database_and_its_size() {
    let tblastx = |words: &[&str]| {
        let mut argv = vec!["losat", "tblastx", "-query", "q.fa"];
        argv.extend(words);
        try_parse_from(argv).map(|cli| match cli.command {
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
        assert!(
            report(&error).contains("Argument cannot be converted"),
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
        &["-db", "d", "-subject_loc", "1-10"],
    ] {
        let error = tblastx(words).unwrap_err();
        assert!(
            report(&error).contains("Incompatible with argument"),
            "{words:?}: {error}"
        );
    }
    // LOSAT Web's validation (`check_options`) checks `-db` as the command line does: the run opens
    // the database from the files that the host gives (line WQA D12a).
    crate::algorithm::tblastx::blast_engine::check_options(&tblastx(&["-db", "d"]).unwrap())
        .unwrap();
}

#[test]
fn evalues_on_which_ncbi_crashes_are_rejected_and_large_hit_list_sizes_run() {
    // Decision D12 (AUTHORITY.md): an infinite -evalue, and DBL_MAX with the composition-based
    // redo, which crashes NCBI's blastp and tblastn on some inputs
    // (`protein_options_check::check_evalue_limit`); without the redo NCBI runs DBL_MAX (line
    // O4c's review D1). blast_hits.c:43-70: a -max_target_seqs whose preliminary size passes
    // NCBI's int wraps or crashes NCBI; LOSAT computes the size without the overflow and runs
    // (PD-LOSAT-NCBI-DEFECTS exception 19).
    let blastp = |words: &[&str]| {
        let Commands::Blastp(args) = parse("blastp", words).unwrap().command else {
            unreachable!("blastp")
        };
        crate::algorithm::blastp::blast_engine::check_options(&args)
    };
    let tblastn = |words: &[&str]| {
        let Commands::Tblastn(args) = parse("tblastn", words).unwrap().command else {
            unreachable!("tblastn")
        };
        crate::algorithm::tblastn::check_options(&args)
    };
    for words in [
        &["-evalue", "+inf"][..],
        &["-evalue", "1e999"],
        &["-evalue", "1e999", "-comp_based_stats", "0"],
        &["-evalue", "1.7976931348623157e308"],
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
        &[
            "-evalue",
            "1.7976931348623157e308",
            "-comp_based_stats",
            "0",
        ],
        &["-max_target_seqs", "2147483647"],
        &["-max_target_seqs", "1073741798"],
        // Past NCBI's `int` overflow of the preliminary hit list size (LOSAT computes it without).
        &["-max_target_seqs", "1073741799"],
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
    for words in [
        &[
            "-evalue",
            "1.7976931348623157e308",
            "-comp_based_stats",
            "0",
        ][..],
        &[
            "-evalue",
            "1.7976931348623157e308",
            "-ungapped",
            "-comp_based_stats",
            "0",
        ],
    ] {
        tblastn(words).unwrap();
    }
    for words in [
        &["-evalue", "+inf", "-comp_based_stats", "0"][..],
        &["-evalue", "1e999", "-ungapped", "-comp_based_stats", "0"],
    ] {
        let error = format!("{:#}", tblastn(words).unwrap_err());
        assert!(error.contains("an infinite -evalue"), "{words:?}: {error}");
    }
}

#[test]
fn version_words_in_an_option_value_are_read_as_ncbi_reads_them() {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001: the application reads
    // -version, -version-full*, -dryrun, -logfile and -conffile anywhere before `--`, also
    // in an option's value (`-out -version` prints the version). The command line prints the
    // versions before parsing (`ncbi_version_request_reads_the_word_before_the_arguments`), and
    // so does the Web engine's ABI v2 (line WQA's D14, 2026-11-20); the parser's callers without
    // an application start (the library API) refuse them by name in that place too, instead of
    // writing a file named "-version" (line RG item E3). `-version=1` is an unknown option: NCBI's
    // reader's error (`Unknown argument: "version"`, exit status 1).
    for (program, upper) in [
        ("blastn", "BLASTN"),
        ("blastp", "BLASTP"),
        ("tblastn", "TBLASTN"),
        ("tblastx", "TBLASTX"),
    ] {
        for words in [
            &["-out", "-version"][..],
            &["-evalue", "-version-full"],
            &["-version-full-xml"],
            &["-out", "o.txt", "-version-full-json"],
        ] {
            let error = parse(program, words).unwrap_err().to_string();
            assert!(
                error.contains(&format!(
                    "the option {} runs on the command line and in LOSAT Web of LOSAT's {upper} only",
                    words.last().unwrap()
                )),
                "{program} {words:?}: {error}"
            );
        }
        for word in ["-version=1", "-version-full=1"] {
            let error = parse(program, &[word]).unwrap_err().to_string();
            let name = word.trim_start_matches('-').trim_end_matches("=1");
            assert!(
                error.contains(&format!("Unknown argument: \"{name}\"")),
                "{program} {word}: {error}"
            );
        }
        assert!(parse(program, &["-out", "o.txt"]).is_ok(), "{program}");
    }
}

// No NCBI counterpart: the command line's entry makes the application start of the program that
// `argv[1]` names when it is one of `ENTRY_NAMES`, and runs a parsed command with that start
// (`cli::main_entry`, which holds the start of every command it parses).
#[test]
fn every_program_is_an_entry_name_on_the_reader() {
    for (name, _) in crate::cli::PROGRAMS {
        assert!(
            crate::cli::ENTRY_NAMES.contains(&name),
            "the program {name} is not one of ENTRY_NAMES"
        );
        assert!(
            crate::cli::reader::tool(name).is_some(),
            "the program {name} is not on the reader"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (cited at
// `ncbi_environment::preparse`): `-logfile FILE` and `-conffile FILE` go to the argument reader
// with their values, a last `-logfile` or `-conffile` without one is dropped. The command line's
// application start reads them (`cli::main_entry`), and so does the start
// of each library call (`api::start`, `try_parse_with_own_start`, line RG item E1w-2) and of each
// run of LOSAT Web's ABI v2 (`web_adapter::try_parse_cli`, item E1w-3a); the parser of ABI v1,
// which runs none, refuses them (`toolkit_file_word_refusal`).
#[test]
fn toolkit_file_options_need_an_application_start() {
    let library = |program: &str, extra: &[&str]| {
        let mut argv = vec!["losat", program, "-query", "q.fa", "-subject", "s.fa"];
        argv.extend(extra);
        crate::cli::try_parse_with_own_start(argv, &DiagTarget::Stderr)
    };
    for (program, upper) in [
        ("blastn", "BLASTN"),
        ("blastp", "BLASTP"),
        ("blastx", "BLASTX"),
        ("tblastn", "TBLASTN"),
        ("tblastx", "TBLASTX"),
    ] {
        for (words, name) in [
            (&["-logfile", "log.txt"][..], "-logfile"),
            (&["-logfile=log.txt"], "-logfile"),
            (&["-conffile", "x.ini"], "-conffile"),
            (&["-conffile=x.ini"], "-conffile"),
            (&["-out", "-conffile=x.ini"], "-conffile"),
            // `-dryrun` is the value of `-logfile`.
            (&["-logfile", "-dryrun"], "-logfile"),
            // `-dryrun` (line RG item E1b): the command line's entry ends the run after the
            // arguments; in an option's value too, as NCBI reads it there.
            (&["-dryrun"], "-dryrun"),
            (&["-evalue", "1", "-dryrun"], "-dryrun"),
            (&["-out", "-dryrun"], "-dryrun"),
        ] {
            let error = parse(program, words).unwrap_err().to_string();
            assert!(
                error.contains(&format!(
                    "the option {name} runs on the command line, in the library API and in the Web engine's ABI v2 of LOSAT's {upper} only: the Web engine's ABI v1 runs no NCBI application start, which reads it"
                )),
                "{program} {words:?}: {error}"
            );
        }
        // The library's parser takes them, as the command line's does after its start (`-out`
        // without a value aside: the library refuses `-out` itself, `api::search`).
        for words in [
            &["-logfile", "log.txt"][..],
            &["-logfile=log.txt"],
            &["-conffile", "x.ini"],
            &["-conffile=x.ini"],
            &["-logfile", "-dryrun"],
            &["-dryrun"],
            &["-evalue", "1", "-dryrun"],
            &["-dryrun", "-logfile", "log.txt", "-conffile", "x.ini"],
        ] {
            let parsed = library(program, words);
            assert!(parsed.is_ok(), "{program} {words:?}: {:?}", parsed.err());
        }
        // A last word without its value is dropped, as NCBI drops it.
        assert!(
            parse(program, &["-evalue", "1", "-logfile"]).is_ok(),
            "{program}"
        );
        assert!(
            parse(program, &["-evalue", "1", "-conffile"]).is_ok(),
            "{program}"
        );
        // `-out` then has no value: the parser's error.
        assert!(parse(program, &["-out", "-logfile"]).is_err(), "{program}");
        // The versions, which the command line prints before parsing (line RG item E3), are
        // refused by the library's parser and the Web engine's ABI v1 alike (line RG item E1w-2); a
        // run of the Web engine's ABI v2 prints them before its parse, as the command line does
        // (`api::start::CallerStart::version`; line WQA's D14, 2026-11-20).
        for word in ["-version", "-version-full-json"] {
            let text = format!(
                "the option {word} runs on the command line and in LOSAT Web of LOSAT's {upper} only: the library API does not write NCBI's version text"
            );
            for error in [
                parse(program, &["-out", word]).unwrap_err(),
                library(program, &["-out", word]).unwrap_err(),
            ] {
                assert!(
                    error.to_string().contains(&text),
                    "{program} {word}: {error}"
                );
            }
        }
    }
    // makeblastdb: the Web engine's parser refuses `-conffile` and `-dryrun` (makeblastdb takes
    // `-logfile`, its own argument), the library's takes them (`api::run_makeblastdb` runs the
    // start); a last word without value is dropped.
    let argv = |words: &[&'static str]| {
        let mut argv = vec!["losat", "makeblastdb", "-in", "x.faa", "-dbtype", "prot"];
        argv.extend_from_slice(words);
        argv
    };
    let makedb = |words: &[&'static str]| try_parse_from(argv(words));
    let library_makedb = |words: &[&'static str]| {
        crate::cli::try_parse_with_own_start(argv(words), &DiagTarget::Stderr)
    };
    for (words, name) in [
        (&["-conffile", "c.ini"][..], "-conffile"),
        (&["-dryrun"], "-dryrun"),
        (&["-logfile", "log.txt", "-dryrun"], "-dryrun"),
    ] {
        let error = makedb(words).unwrap_err().to_string();
        assert!(
            error.contains(&format!(
                "the option {name} runs on the command line and in the library API of LOSAT's makeblastdb only: the Web engine's ABI v1 runs no NCBI application start, which reads it"
            )),
            "{words:?}: {error}"
        );
        let parsed = library_makedb(words);
        assert!(parsed.is_ok(), "{words:?}: {:?}", parsed.err());
    }
    assert!(makedb(&["-logfile", "log.txt"]).is_ok());
    assert!(makedb(&["-conffile"]).is_ok());
    assert!(makedb(&["-logfile"]).is_ok());
    assert!(library_makedb(&["-conffile"]).is_ok());
    // So do the versions, which the command line prints before parsing (line RG item E3): the
    // library's makeblastdb prints them before it parses (`api::run_makeblastdb`).
    for word in [
        "-version",
        "-version-full",
        "-version-full-xml",
        "-version-full-json",
    ] {
        for error in [
            makedb(&["-title", word]).unwrap_err(),
            library_makedb(&["-title", word]).unwrap_err(),
        ] {
            assert!(
                error.to_string().contains(&format!(
                    "the option {word} runs on the command line and in the library API of LOSAT's makeblastdb only: the Web engine does not write NCBI's version text"
                )),
                "{error}"
            );
        }
    }
    // The programs that the library API does not run: the refusals name the command line alone.
    for (program, upper) in [("rpsblast", "RPSBLAST"), ("rpstblastn", "RPSTBLASTN")] {
        let refused = |word: &'static str| {
            try_parse_from(["losat", program, "-query", "q.fa", "-db", "d", word])
                .unwrap_err()
                .to_string()
        };
        let error = refused("-dryrun");
        assert!(
            error.contains(&format!(
                "the option -dryrun runs on the command line of LOSAT's {upper} only: the Web engine's ABI v1 runs no NCBI application start, which reads it"
            )),
            "{program}: {error}"
        );
        let error = refused("-version");
        assert!(
            error.contains(&format!(
                "the option -version runs on the command line of LOSAT's {upper} only: the Web engine does not write NCBI's version text"
            )),
            "{program}: {error}"
        );
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
        try_parse_from(argv)
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
        try_parse_from(argv).map(|cli| match cli.command {
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
        assert!(
            report(&error).contains("Argument cannot be converted"),
            "{word}"
        );
    }
    for words in [
        &["-db", "d", "-subject", "s.fa"][..],
        &["-subject", "s.fa", "-db", "d"],
    ] {
        let error = blastx(words).unwrap_err();
        assert!(
            report(&error).contains("Incompatible with argument"),
            "{words:?}: {error}"
        );
    }
    // `-db` excludes `-subject_loc` (blast_args.cpp:2377-2381), as for `-subject`.
    let error = blastx(&["-db", "d", "-subject_loc", "1-10"]).unwrap_err();
    assert!(
        report(&error).contains("Incompatible with argument"),
        "{error}"
    );
}

#[test]
fn makedb_is_not_a_command_and_makeblastdb_checks_apply() {
    let makedb = |words: &[&str]| {
        let mut argv = vec!["losat", "makeblastdb", "-dbtype", "prot"];
        argv.extend(words);
        try_parse_from(argv)
    };
    // The old name is gone: an unknown command, with every spelling of its arguments.
    for words in [
        &["makedb", "-in", "x.faa", "-dbtype", "prot"][..],
        &["makedb", "-h"],
        &["makedb", "-version"],
    ] {
        let mut argv = vec!["losat"];
        argv.extend(words);
        let error = try_parse_from(argv).unwrap_err();
        assert!(
            report(&error).contains("unknown option or argument 'makedb'"),
            "{words:?}: {error}"
        );
    }
    let Commands::Makedb(args) = makedb(&["-in", "x.faa", "-out", "db"]).unwrap().command else {
        unreachable!("makeblastdb")
    };
    assert_eq!(
        (args.input.as_str(), args.out.as_deref()),
        ("x.faa", Some("db"))
    );
    // The checks that name the command apply: makeblastdb's dependencies.
    let error = makedb(&["-gi_mask"]).unwrap_err().to_string();
    assert!(
        error.contains(
            "Argument \"parse_seqids\". Must be specified, as it is required by argument:  \
             `gi_mask'"
        ),
        "{error}"
    );
    assert!(makedb(&["-parse_seqids"]).is_ok());
    assert!(makedb(&["-in", "x.faa", "--"]).is_ok());
    assert!(makedb(&["--", "-in", "x.faa"]).is_err());
    assert!(makedb(&["-evalue", "1"]).is_err());
    let error = makedb(&["-h"]).unwrap_err();
    assert!(
        printed(&error).starts_with("USAGE\n  makeblastdb "),
        "{error}"
    );
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
        let error = try_parse_from(["losat", program, "-h"]).unwrap_err();
        assert!(
            printed(&error).starts_with(&format!("USAGE\n  {program} ")),
            "{program}: {error}"
        );
        // As an option's value, `-h` is the value.
        assert!(
            parse(program, &["-out", "-h"]).is_ok() || program == "makeblastdb",
            "{program}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3020-3023, 2840-2852 (cited at
// `cli::try_parse_from`)
// NCBI's reader stops at its help flags, the words before them read; `-help-full` is `-help` with
// the hidden arguments, which LOSAT's help does not list (NCBI BLAST+ 2.17.0: the same bytes in
// every place but makeblastdb's hidden `-randomize`; line RG item E2).
#[test]
fn help_full_is_help_and_the_reader_stops_at_a_help_flag() {
    // The help that a command line prints (NCBI's text, `cli::reader`).
    let help = |program: &str, words: &[&str]| -> String {
        let mut argv = vec!["losat", program];
        argv.extend(words);
        let error = try_parse_from(argv)
            .err()
            .unwrap_or_else(|| panic!("{program} {words:?}: no help"));
        printed(&error).to_string()
    };
    for (program, value_option) in [
        ("blastn", "-out"),
        ("blastp", "-out"),
        ("blastx", "-out"),
        ("tblastn", "-out"),
        ("tblastx", "-out"),
        ("psiblast", "-out"),
        ("rpsblast", "-out"),
        ("rpstblastn", "-out"),
        ("makeblastdb", "-title"),
        ("rpsbproc", "-o"),
    ] {
        let full = help(program, &["-help"]);
        let short = help(program, &["-h"]);
        let hidden = help(program, &["-help-full"]);
        if program == "makeblastdb" {
            assert!(hidden.contains("\n -randomize\n"), "{hidden}");
            assert!(!full.contains("-randomize"), "{full}");
        } else {
            assert_eq!(hidden, full, "{program}");
        }
        for words in [
            &["-help-full", "-h"][..],
            &["-help-full", "-help-full"],
            &["-help-full", "-xmlhelp"],
            &["-help-full", "-bogus"],
            &["-help-full", value_option],
            &["-help-full", "-evalue", "abc"],
        ] {
            assert_eq!(help(program, words), hidden, "{program} {words:?}");
        }
        for words in [
            &["-help", "-h"][..],
            &["-help", "-bogus"],
            &["-help", value_option],
        ] {
            assert_eq!(help(program, words), full, "{program} {words:?}");
        }
        for words in [
            &["-h", "-help-full"][..],
            &["-h", "-help"],
            &["-h", "-xmlhelp"],
            &["-h", "-bogus"],
            &["-h", value_option],
        ] {
            assert_eq!(help(program, words), short, "{program} {words:?}");
        }
        // A word before the flag is read first.
        let error = try_parse_from(["losat", program, "-bogus", "-help-full"])
            .err()
            .unwrap_or_else(|| panic!("{program}: -bogus read"));
        assert_ne!(status(&error), 0, "{program}");
    }
    // As an option's value, `-help-full` is the value; outside a program it is no flag.
    assert!(parse("blastn", &["-out", "-help-full"]).is_ok());
    assert!(try_parse_from(["losat", "-help-full"]).is_err());
    assert_ne!(
        try_parse_from(["losat", "-help-full"])
            .err()
            .map(|error| status(&error)),
        Some(0)
    );
    // `-xmlhelp` prints NCBI's XML on the command line (line RG item E2x, `ncbi_xmlhelp_tests`)
    // and, for the programs it runs, in the Web engine's ABI v2 (line WQA's D14, 2026-11-20); the
    // library API, which these tests are (no application start), refuses it by name.
    for (program, place) in [
        (
            "blastn",
            "the command line and in LOSAT Web of LOSAT's BLASTN",
        ),
        ("makeblastdb", "the command line of LOSAT's"),
        ("rpsbproc", "the command line of LOSAT's"),
    ] {
        let error = try_parse_from(["losat", program, "-xmlhelp"]).unwrap_err();
        let text = error.to_string();
        assert!(
            text.contains(&format!("the option -xmlhelp runs on {place}")),
            "{program}: {text}"
        );
    }
}

#[test]
fn ncbi_version_request_reads_the_word_before_the_arguments() {
    use crate::cli::ncbi_version_request;
    let version = |words: &[&str]| ncbi_version_request(words.iter().copied());
    // LOSAT's version where NCBI's programs print their own (`crate::wording::LOSAT_VERSION`).
    let losat = format!("LOSAT {}", env!("CARGO_PKG_VERSION"));
    let text = version(&["losat", "blastn", "-version"]).unwrap();
    assert!(
        text.starts_with(&format!("blastn: {losat}\n Package: LOSAT ")),
        "{text}"
    );
    assert!(
        text.ends_with(", output compatible with blast 2.17.0\n"),
        "{text}"
    );
    assert!(version(&["losat", "makeblastdb", "-x", "-version"])
        .unwrap()
        .starts_with(&format!("makeblastdb: {losat}\n")));
    // Where another option expects its value, and next to invalid arguments.
    assert!(version(&["losat", "blastx", "-out", "-version"]).is_some());
    assert!(version(&["losat", "blastp", "-evalue", "x", "-version", "-bogus"]).is_some());
    // Not a request: after `--`, a longer word, no program, the top level.
    assert!(version(&["losat", "blastn", "--", "-version"]).is_none());
    // The value of `-logfile` or `-conffile` (NCBI logs to a file named `-version`).
    assert!(version(&["losat", "blastn", "-logfile", "-version"]).is_none());
    assert!(version(&["losat", "makeblastdb", "-conffile", "-version"]).is_none());
    assert!(version(&["losat", "blastn", "-logfile=x", "-version"]).is_some());
    assert!(version(&["losat", "blastn", "-version=1"]).is_none());
    assert!(version(&["losat", "blastn", "-version-full=1"]).is_none());
    // The full versions (line RG item E3, `ncbi_version::print`): the first version word is
    // printed.
    let full = version(&[
        "losat",
        "blastn",
        "-evalue",
        "x",
        "-version-full",
        "-version",
    ])
    .unwrap();
    assert!(
        full.starts_with(&format!("blastn: {losat}\n Package: LOSAT "))
            && full.ends_with(
                ", output compatible with blast 2.17.0\n Package-Config:  \n GI-64bit:  TRUE\n"
            ),
        "{full}"
    );
    assert!(version(&["losat", "makeblastdb", "-title", "-version-full-xml"])
        .unwrap()
        .contains(&format!(
            "<appname>makeblastdb</appname>\n<version_info major=\"{}\" minor=\"{}\" patch_level=\"{}\"",
            env!("CARGO_PKG_VERSION_MAJOR"),
            env!("CARGO_PKG_VERSION_MINOR"),
            env!("CARGO_PKG_VERSION_PATCH")
        )));
    assert!(version(&["losat", "rpsbproc", "-version-full-json", "-version"])
        .unwrap()
        .contains(&format!(
            "\"appname\": \"rpsbproc\",\n    \"version_info\": {{\"major\": {}, \"minor\": {}, \"patch_level\": {}",
            env!("CARGO_PKG_VERSION_MAJOR"),
            env!("CARGO_PKG_VERSION_MINOR"),
            env!("CARGO_PKG_VERSION_PATCH")
        )));
    assert!(version(&["losat", "rpstblastn", "-version-full"])
        .unwrap()
        .starts_with(&format!("rpstblastn: {losat}\n")));
    assert!(version(&["losat", "blastn", "--", "-version-full"]).is_none());
    assert!(version(&["losat", "blastn", "-logfile", "-version-full"]).is_none());
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
            // NCBI's reader: `-mt_mode` requires `-num_threads` given (its default does not count).
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
            error.contains("Unknown argument: \"mt_mode\""),
            "{program}: {error}"
        );
    }
}

#[test]
fn use_index_parses_as_an_ncbi_boolean() {
    for value in [
        "false", "False", "f", "no", "0", "true", "T", "yes", "1", "on",
    ] {
        assert!(parse("blastn", &["-use_index", value]).is_ok(), "{value}");
    }
    assert!(parse("blastn", &["-use_index", "maybe"]).is_err());
    // NCBI's other programs have no -use_index nor -index_name.
    assert!(parse("blastp", &["-use_index", "false"]).is_err());
    assert!(parse("blastp", &["-index_name", "x"]).is_err());
    assert!(parse("blastn", &["-index_name", "x"]).is_ok());
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
// NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:946
#[test]
fn native_error_exception_truncates_what_at_first_nul() {
    let site = crate::ncbi_exception::ThrowSite {
        file: "src/objects/seq/../seqloc/Seq_id.cpp",
        line: 2501,
        kind: "CSeqIdException::eFormat",
        function: Some("ncbi::objects::CSeq_id::Set"),
    };
    let plain = crate::ncbi_exception::NcbiException::new(site, b"Malformatted ID foo".to_vec());
    let plain_err = crate::cli::NativeError::exception(&plain);
    assert_eq!(plain_err.exit, 255);
    assert_eq!(plain_err.message, b"Error: Malformatted ID foo\n\n");

    let with_nul = crate::ncbi_exception::NcbiException::new(
        site,
        b"Malformatted ID \x00=\x00\x00\x00".to_vec(),
    );
    let nul_err = crate::cli::NativeError::exception(&with_nul);
    assert_eq!(nul_err.exit, 255);
    assert_eq!(nul_err.message, b"Error: Malformatted ID \n");
}
