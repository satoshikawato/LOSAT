use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:336-342 and
// blast_aux.cpp:662-673 (cited at `TaskParser::check`): a task read in lower case, LOSAT's
// refusal for another program's NCBI task, NCBI's error for a name NCBI does not know.
#[test]
fn task_check_reads_the_name_as_ncbis_options_factory() {
    let mut task = String::from("BLASTP-Fast");
    BLASTP_TASK.check(&mut task).unwrap();
    assert_eq!(task, "blastp-fast");
    let message = |name: &str| {
        let mut task = name.to_string();
        BLASTP_TASK.check(&mut task).unwrap_err().to_string()
    };
    assert_eq!(
        message("blastx"),
        "the task blastx is not supported by LOSAT's BLASTP (use blastp, blastp-fast or blastp-short)"
    );
    assert_eq!(message("BlastQ"), "'blastq' is not a supported task");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137
// A search strategy can carry a negative `WindowSize` (the command line takes 0 and more) and NCBI
// runs it: the resolved window keeps the sign (the report's footer prints it, and the search is
// one-hit), where it was 0.
#[test]
fn a_negative_window_of_a_strategy_resolves_signed() {
    let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
    let TestCommand::Blastp(mut args) = cli.command;
    for window in [-1, -40, i32::MIN] {
        args.window_size = Some(window);
        assert_eq!(args.resolve().unwrap().window_size, window);
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3152-3171 (`s_TruncateSpaces`
// with `eTrunc_Both`: `isspace` at both ends; `seqdb::db_argument_name`).
#[test]
fn the_database_name_loses_the_white_space_at_its_ends() {
    let resolved_name = |given: &str| {
        let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-db", given]);
        let TestCommand::Blastp(args) = cli.command;
        args.resolve().unwrap().db
    };
    for (given, kept) in [
        ("db", "db"),
        (" db ", "db"),
        ("\tdir/db\n", "dir/db"),
        ("\r\x0b\x0c db", "db"),
        ("a b", "a b"),
        (" a  b ", "a  b"),
        ("", ""),
        ("   ", ""),
        ("\u{a0}db", "\u{a0}db"),
    ] {
        assert_eq!(resolved_name(given).as_deref(), Some(kept), "{given:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-349: parse the same single string as the public CLI.
fn parse_cli<'a>(args: impl IntoIterator<Item = &'a str>) -> TestCli {
    match crate::cli::try_parse_from(args).unwrap().command {
        crate::cli::Commands::Blastp(args) => TestCli {
            command: TestCommand::Blastp(args),
        },
        command => panic!("not blastp: {command:?}"),
    }
}

#[derive(Debug)]
struct TestCli {
    command: TestCommand,
}

#[derive(Debug)]
enum TestCommand {
    Blastp(BlastpArgs),
}

#[test]
fn test_default_blastp_options_resolve_to_ncbi_defaults() {
    let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
    let TestCommand::Blastp(args) = cli.command;
    assert_eq!(args.evalue, None);
    assert_eq!(args.threshold, None);
    assert_eq!(args.word_size, None);
    assert_eq!(args.window_size, None);
    assert_eq!(args.matrix, None);
    assert_eq!(args.gap_open, None);
    assert_eq!(args.gap_extend, None);
    assert_eq!(args.comp_based_stats, "2");
    assert_eq!(args.seg, None);
    assert!(!args.ungapped);
    assert!(!args.use_sw_tback);

    let resolved = args.resolve().expect("resolved blastp args");
    assert_eq!(resolved.task, "blastp");
    assert_eq!(resolved.evalue, 10.0);
    assert_eq!(resolved.word_size, 3);
    assert_eq!(
        resolved.lookup_table_type,
        BlastpLookupTableType::AaLookupTable
    );
    assert_eq!(resolved.threshold, 11.0);
    assert_eq!(resolved.window_size, 40);
    assert_eq!(resolved.scoring.matrix, ScoringMatrix::Blosum62);
    assert_eq!(resolved.scoring.gap_open, 11);
    assert_eq!(resolved.scoring.gap_extend, 1);
    assert_eq!(
        resolved.comp_based_stats,
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: false,
        }
    );
    assert_eq!(resolved.seg, BlastpSegSpec::No);
    assert_eq!(resolved.max_target_seqs, 500);
    assert_eq!(resolved.max_hsps_per_subject, 0);
    assert!(!resolved.ungapped);
    assert!(!resolved.use_sw_tback);
    assert!(!resolved.chaining);
    assert!(!resolved.soft_masking);
    assert_eq!(resolved.outfmt, "0");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-155 (IDENTITY's
// search runs mode 0 whatever mode was asked; NCBI 2.17.0's footers and reports of
// `-matrix IDENTITY` with modes 1, 2, 3, T and 2u equal those of mode 0, item SC4b)
#[test]
fn identity_resets_the_composition_mode_to_0() {
    let no_cbs = BlastpCompBasedStats {
        mode: BlastpCompositionMode::NoCompositionBasedStats,
        unified_p: false,
    };
    for (matrix, cbs, reset) in [
        ("IDENTITY", "2", true),
        ("identity", "1", true),
        ("Identity", "3", true),
        ("IDENTITY", "T", true),
        ("IDENTITY", "2u", true),
        ("IDENTITY", "0", false),
        ("IDENTITY", "F", false),
        ("IDENTITY", "0u", false),
    ] {
        let cli = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-matrix",
            matrix,
            "-comp_based_stats",
            cbs,
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.comp_based_stats, no_cbs, "{matrix} {cbs}");
        assert_eq!(resolved.identity_cbs_reset, reset, "{matrix} {cbs}");
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Identity);
        assert_eq!(
            (resolved.scoring.gap_open, resolved.scoring.gap_extend),
            (15, 2)
        );
        assert_eq!((resolved.threshold, resolved.window_size), (27.0, 40));
    }
    // The other matrices keep their mode.
    let cli = parse_cli([
        "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-matrix", "PAM30",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp args");
    assert!(!resolved.identity_cbs_reset);
    assert_eq!(
        resolved.comp_based_stats.mode,
        BlastpCompositionMode::CompositionMatrixAdjust
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:498-501 (an ungapped search
// takes any `-matrix` name; the score block of the first query batch fails for an unknown
// one, `unknown_matrix`). NCBI BLAST+ 2.17.0 exits 139 on `IDENTITYX` (a name that starts
// with IDENTITY, as BLASTX's `-matrix` refusal says), which LOSAT refuses. A gapped search
// stops at the matrix error of the options.
#[test]
fn ungapped_search_takes_an_unknown_matrix_name() {
    let resolve = |extra: &[&str]| {
        let mut words = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
        words.extend_from_slice(extra);
        let cli = parse_cli(words);
        let TestCommand::Blastp(args) = cli.command;
        args.resolve()
    };
    let unknown = resolve(&["-ungapped", "-comp_based_stats", "0", "-matrix", "FOO"])
        .expect("an ungapped search with an unknown matrix");
    assert!(unknown.unknown_matrix);
    assert!(
        !unknown.has_gumbel(),
        "no Gumbel block in an ungapped search"
    );
    assert!(!unknown.old_fsc, "`old_fsc` is -searchsp or OLD_FSC only");
    let known = resolve(&["-ungapped", "-comp_based_stats", "0"]).expect("resolved");
    assert!(!known.unknown_matrix);
    assert!(!known.has_gumbel());
    assert!(known.matrix_file.is_none());
    let searchsp =
        resolve(&["-ungapped", "-comp_based_stats", "0", "-searchsp", "1000"]).expect("resolved");
    assert!(searchsp.old_fsc && !searchsp.has_gumbel());
    let refused = resolve(&[
        "-ungapped",
        "-comp_based_stats",
        "0",
        "-matrix",
        "IDENTITYX",
    ])
    .expect_err("IDENTITYX");
    assert_eq!(
        refused.to_string(),
        "-matrix IDENTITYX is read as IDENTITY without IDENTITY's statistics table (NCBI BLAST+ 2.17.0 aborts); this is not supported by LOSAT's BLASTP"
    );
    // A gapped search keeps its own default: the Gumbel block.
    let gapped = resolve(&["-comp_based_stats", "0"]).expect("resolved");
    assert!(!gapped.old_fsc && gapped.has_gumbel());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
// `-soft_masking` takes NCBI's boolean spellings (`NStr::StringToBool`); a task does not
// change it (`blastp-short`'s `ClearFilterOptions` comes before the handler).
#[test]
fn soft_masking_sets_mask_at_hash() {
    let resolve = |extra: &[&str]| {
        let mut argv = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
        argv.extend_from_slice(extra);
        let TestCommand::Blastp(args) = parse_cli(argv).command;
        args.resolve().expect("resolved blastp args").soft_masking
    };
    for value in ["true", "T", "1", "YES", "on"] {
        assert!(resolve(&["-soft_masking", value]), "{value}");
    }
    for value in ["false", "F", "0", "no", "off"] {
        assert!(!resolve(&["-soft_masking", value]), "{value}");
    }
    assert!(resolve(&["-task", "blastp-short", "-soft_masking", "true"]));
    assert!(resolve(&[
        "-task",
        "blastp-fast",
        "-soft_masking",
        "true",
        "-seg",
        "yes"
    ]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973 (`blastinput::report_counts`): the hit list size and the
// counts of outfmt 0 from -num_descriptions / -num_alignments; "No hits are being saved".
#[test]
fn num_descriptions_and_num_alignments_set_the_counts() {
    let resolve = |extra: &[&str]| {
        let mut argv = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
        argv.extend_from_slice(extra);
        let cli = parse_cli(argv);
        let TestCommand::Blastp(args) = cli.command;
        args.resolve()
    };
    let counts = |extra: &[&str]| {
        let resolved = resolve(extra).unwrap();
        (
            resolved.max_target_seqs,
            resolved.num_descriptions,
            resolved.num_alignments,
        )
    };
    assert_eq!(counts(&[]), (500, 500, 250));
    assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
    assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
    assert_eq!(
        counts(&["-num_descriptions", "2", "-num_alignments", "4"]),
        (4, 2, 4)
    );
    assert_eq!(counts(&["-max_target_seqs", "7"]), (7, 7, 7));
    assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
    assert_eq!(
        counts(&["-outfmt", "6", "-num_descriptions", "2"]),
        (500, 500, 500)
    );
    let error = resolve(&["-num_descriptions", "0", "-num_alignments", "0"]).unwrap_err();
    assert_eq!(
        error.to_string(),
        "LOSAT query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n"
    );
    assert!(resolve(&["-outfmt", "7", "-num_alignments", "0"]).is_err());
    assert!(resolve(&["-num_descriptions", "0"]).is_ok());
}

#[test]
fn test_blastp_short_task_resolves_ncbi_defaults() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:395-399: blastp-fast sets word size 5, the
    // compressed lookup and BLAST_WORD_THRESHOLD_BLASTP_FAST (20), which is kept
    // because (int)20 is not BLASTP's default threshold (blast_args.cpp:586-600).
    let cli = parse_cli([
        "losat",
        "blastp",
        "-query",
        "q.faa",
        "-subject",
        "s.faa",
        "-task",
        "blastp-short",
    ]);
    let TestCommand::Blastp(mut args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp-short args");
    assert_eq!(resolved.task, "blastp-short");
    assert_eq!(resolved.evalue, 20000.0);
    assert_eq!(resolved.scoring.matrix, ScoringMatrix::Pam30);
    assert_eq!(resolved.scoring.gap_open, 9);
    assert_eq!(resolved.scoring.gap_extend, 1);
    assert_eq!(resolved.word_size, 2);
    assert_eq!(resolved.threshold, 16.0);
    assert_eq!(resolved.window_size, 15);
    assert_eq!(resolved.seg, BlastpSegSpec::No);
    assert!(!resolved.chaining);
    args.evalue = Some(42.0);
    assert_eq!(args.resolve().unwrap().evalue, 42.0);
}

#[test]
fn test_blastp_fast_task_resolves_ncbi_defaults() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:395-399: blastp-fast sets word size 5, the
    // compressed lookup and BLAST_WORD_THRESHOLD_BLASTP_FAST (20), which is kept
    // because (int)20 is not BLASTP's default threshold (blast_args.cpp:586-600).
    let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
    let TestCommand::Blastp(mut args) = cli.command;
    args.task = "blastp-fast".into();
    let resolved = args.resolve().expect("resolved blastp-fast args");
    assert_eq!(resolved.task, "blastp-fast");
    assert_eq!(resolved.evalue, 10.0);
    assert_eq!(resolved.word_size, 5);
    assert_eq!(resolved.threshold, 20.0);
    assert_eq!(
        resolved.lookup_table_type,
        BlastpLookupTableType::CompressedAaLookupTable
    );
    assert!(resolved.chaining);
    args.evalue = Some(42.0);
    assert_eq!(args.resolve().unwrap().evalue, 42.0);
}

#[test]
fn test_seg_yes_no_and_custom_parser() {
    let cli = parse_cli([
        "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-seg", "no",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    assert_eq!(args.resolve().unwrap().seg, BlastpSegSpec::No);

    let cli = parse_cli([
        "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-seg", "yes",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    assert_eq!(args.resolve().unwrap().seg, BlastpSegSpec::Yes);

    let cli = parse_cli([
        "losat",
        "blastp",
        "-query",
        "q.faa",
        "-subject",
        "s.faa",
        "-seg",
        "15 2.5 3.0",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    assert_eq!(
        args.resolve().unwrap().seg,
        BlastpSegSpec::WindowLocutHicut {
            window: 15,
            locut: 2.5,
            hicut: 3.0,
        }
    );
}

#[test]
fn test_comp_based_stats_parser_supports_ncbi_modes() {
    assert_eq!(
        parse_comp_based_stats("0").unwrap(),
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::NoCompositionBasedStats,
            unified_p: false,
        }
    );
    assert_eq!(
        parse_comp_based_stats("1").unwrap(),
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionBasedStats,
            unified_p: false,
        }
    );
    assert_eq!(
        parse_comp_based_stats("2u").unwrap(),
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: true,
        }
    );
    assert_eq!(
        parse_comp_based_stats("3").unwrap(),
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::ForceFullMatrixAdjust,
            unified_p: false,
        }
    );
    assert_eq!(
        parse_comp_based_stats("T").unwrap(),
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: false,
        }
    );
}

#[test]
fn test_matrix_resolution_updates_gap_threshold_and_window_defaults() {
    let cli = parse_cli([
        "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-matrix", "BLOSUM45",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp args");
    assert_eq!(resolved.scoring.matrix, ScoringMatrix::Blosum45);
    assert_eq!(resolved.scoring.gap_open, 14);
    assert_eq!(resolved.scoring.gap_extend, 2);
    assert_eq!(resolved.threshold, 14.0);
    assert_eq!(resolved.window_size, 60);
}

#[test]
fn test_explicit_gap_values_override_matrix_defaults_independently() {
    let cli = parse_cli([
        "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-matrix", "PAM70", "-gapopen",
        "11",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp args");
    assert_eq!(resolved.scoring.matrix, ScoringMatrix::Pam70);
    assert_eq!(resolved.scoring.gap_open, 11);
    assert_eq!(resolved.scoring.gap_extend, 1);
}

#[test]
fn test_word_size_above_four_selects_compressed_lookup_thresholds() {
    let cli = parse_cli([
        "losat",
        "blastp",
        "-query",
        "q.faa",
        "-subject",
        "s.faa",
        "-word_size",
        "5",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp args");
    assert_eq!(
        resolved.lookup_table_type,
        BlastpLookupTableType::CompressedAaLookupTable
    );
    assert_eq!(resolved.threshold, 19.3);

    let cli = parse_cli([
        "losat",
        "blastp",
        "-query",
        "q.faa",
        "-subject",
        "s.faa",
        "-word_size",
        "7",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let resolved = args.resolve().expect("resolved blastp args");
    assert_eq!(resolved.threshold, 20.25);
}

#[test]
fn test_ungapped_with_composition_based_stats_is_rejected() {
    let cli = parse_cli([
        "losat",
        "blastp",
        "-query",
        "q.faa",
        "-subject",
        "s.faa",
        "-ungapped",
        "-comp_based_stats",
        "2",
    ]);
    let TestCommand::Blastp(args) = cli.command;
    let err = args
        .resolve()
        .expect_err("ungapped comp-based stats should fail");
    assert!(err
        .to_string()
        .contains("Composition-adjusted searched are not supported with an ungapped search"));
}

#[test]
fn test_comp_based_stats_cli_string_matches_ncbi_modes() {
    assert_eq!(
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::NoCompositionBasedStats,
            unified_p: false,
        }
        .to_ncbi_cli_string(),
        "0"
    );
    assert_eq!(
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionBasedStats,
            unified_p: false,
        }
        .to_ncbi_cli_string(),
        "1"
    );
    assert_eq!(
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: true,
        }
        .to_ncbi_cli_string(),
        "2u"
    );
    assert_eq!(
        BlastpCompBasedStats {
            mode: BlastpCompositionMode::ForceFullMatrixAdjust,
            unified_p: false,
        }
        .to_ncbi_cli_string(),
        "3"
    );
}

#[test]
fn test_seg_cli_string_matches_ncbi_forms() {
    assert_eq!(BlastpSegSpec::No.to_ncbi_cli_string(), "no");
    assert_eq!(BlastpSegSpec::Yes.to_ncbi_cli_string(), "yes");
    assert_eq!(
        BlastpSegSpec::WindowLocutHicut {
            window: 12,
            locut: 2.2,
            hicut: 2.5,
        }
        .to_ncbi_cli_string(),
        "12 2.2 2.5"
    );
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:276-286,
// api/blast_prot_options.cpp:100-111 and core/blast_options.c:608-618 (cited at the fields
// and at `validate_protein_options`): the X-drops reach the resolved arguments (7, 15 and 25
// bits at every task), and an ungapped X-drop of 0 or below is NCBI's options error.
#[test]
fn the_x_drops_reach_the_resolved_arguments() {
    let resolved = |words: &[&str]| {
        let mut argv = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
        argv.extend_from_slice(words);
        let TestCommand::Blastp(args) = parse_cli(argv).command;
        args.resolve()
    };
    let bits = |args: &ResolvedBlastpArgs| {
        (
            args.x_drop_ungapped_bits,
            args.x_drop_gapped_bits,
            args.x_drop_final_bits,
        )
    };
    for task in ["blastp", "blastp-fast", "blastp-short"] {
        let args = resolved(&["-task", task]).unwrap();
        assert_eq!(bits(&args), (7.0, 15.0, 25.0), "{task}");
    }
    let args = resolved(&[
        "-xdrop_ungap",
        "3",
        "-xdrop_gap",
        "0",
        "-xdrop_gap_final",
        "-1e300",
    ])
    .unwrap();
    assert_eq!(bits(&args), (3.0, 0.0, -1e300));
    for value in ["0", "-1", "-0.0", "-inf"] {
        let error = resolved(&["-xdrop_ungap", value]).unwrap_err();
        assert!(
            format!("{error:#}").contains("x_dropoff must be greater than zero"),
            "{value}: {error:#}"
        );
    }
    // `<= 0.0` is false for NaN, as in C: NCBI searches with it.
    assert!(resolved(&["-xdrop_ungap", "+nan"]).is_ok());
}

// Review of PS15b-1, B1/S2 and S1: psiblast's results carry no query masks in any round
// (`report::record::carries_query_masks`), BLASTP's and RPS-BLAST's do; the formatters' cuts of the
// hit lists (`BlastpApp::list_cut`): the archive keeps them, psiblast's XML cuts them before the
// round's PSSM, the other formats of psiblast and rpsblast after it, BLASTP's are never cut.
#[test]
fn the_apps_masks_and_cuts() {
    use crate::report::record::carries_query_masks;
    assert!(!carries_query_masks(BlastpApp::Psiblast.program_type()));
    assert!(carries_query_masks(BlastpApp::Blastp.program_type()));
    assert!(carries_query_masks(BlastpApp::Rpsblast.program_type()));
    // Formats 0, 5, 8, 11, 12, 13, 16 and 18: the organism report (18) cuts nothing in either
    // app; psiblast's XML2 and JSON documents (13 to 16) cut before the PSSM (line PS, PS15b-2).
    for (app, cuts) in [
        (
            BlastpApp::Psiblast,
            [
                ListCut::AfterPssm,
                ListCut::BeforePssm,
                ListCut::AfterPssm,
                ListCut::Uncut,
                ListCut::AfterPssm,
                ListCut::BeforePssm,
                ListCut::BeforePssm,
                ListCut::Uncut,
            ],
        ),
        (
            BlastpApp::Rpsblast,
            [
                ListCut::AfterPssm,
                ListCut::AfterPssm,
                ListCut::AfterPssm,
                ListCut::Uncut,
                ListCut::AfterPssm,
                ListCut::AfterPssm,
                ListCut::AfterPssm,
                ListCut::Uncut,
            ],
        ),
        (BlastpApp::Blastp, [ListCut::Uncut; 8]),
    ] {
        for (outfmt, cut) in [0, 5, 8, 11, 12, 13, 16, 18].into_iter().zip(cuts) {
            assert_eq!(app.list_cut(outfmt), cut, "{app:?} -outfmt {outfmt}");
        }
    }
}
