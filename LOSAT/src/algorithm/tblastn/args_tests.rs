use super::*;

use crate::cli::{try_parse_from, Cli, Commands};

// TBLASTN's `-db_gencode`: the genetic codes of NCBI's reader and LOSAT's code 32
// (PD-TLOSAN-LOCAL-GENCODE-32); a value that is no code is NCBI's argument error.
#[test]
fn tblastn_genetic_codes_are_read_as_ncbi_reads_them() {
    for (value, taken) in [
        ("1", true),
        ("+5", true),
        ("005", true),
        ("11", true),
        ("32", true),
        ("33", true),
        ("0", false),
        ("7", false),
        ("8", false),
        ("255", false),
        ("abc", false),
        ("", false),
    ] {
        assert_eq!(parse(&["-db_gencode", value]).is_ok(), taken, "{value:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942
#[test]
fn lcase_masking_flag_is_explicit() {
    assert!(!parse(&[]).unwrap().lcase_masking);
    assert!(parse(&["-lcase_masking"]).unwrap().lcase_masking);
}

fn parse(extra: &[&str]) -> std::result::Result<TblastnArgs, crate::cli::ArgvError> {
    let mut argv = vec!["losat", "tblastn", "-query", "q.faa", "-subject", "s.fna"];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv)?;
    let Commands::Tblastn(args) = cli.command else {
        unreachable!()
    };
    Ok(args)
}

fn error(extra: &[&str]) -> String {
    parse(extra).unwrap().resolve().unwrap_err().to_string()
}

fn limit(extra: &[&str]) -> String {
    check_options(&parse(extra).unwrap())
        .unwrap_err()
        .to_string()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973 (`blastinput::report_counts`).
#[test]
fn num_descriptions_and_num_alignments_set_the_counts() {
    let counts = |extra: &[&str]| {
        let resolved = parse(extra).unwrap().resolve().unwrap();
        (
            resolved.hitlist_size,
            resolved.num_descriptions,
            resolved.num_alignments,
        )
    };
    assert_eq!(counts(&[]), (500, 500, 250));
    assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
    assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
    assert_eq!(counts(&["-max_target_seqs", "7"]), (7, 7, 7));
    assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
    assert_eq!(
        error(&["-num_descriptions", "0", "-num_alignments", "0"]),
        "LOSAT query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastn_options.cpp:54-85 and c++/src/algo/blast/api/blast_prot_options.cpp:85-141:
#[test]
fn ordinary_defaults_and_code32_are_accepted() {
    let args = parse(&[]).unwrap();
    assert_eq!(args.outfmt, "0");
    assert_eq!(args.comp_based_stats, "2");
    let resolved = args.resolve().unwrap();
    assert_eq!(resolved.task, "tblastn");
    assert_eq!(resolved.db_gencode, 1);
    assert_eq!(resolved.word_size, 3);
    assert_eq!(resolved.threshold, 13.0);
    assert_eq!(resolved.window_size, 40);
    assert_eq!((resolved.gap_open, resolved.gap_extend), (11, 1));
    assert_eq!(resolved.evalue, 10.0);
    assert_eq!(resolved.matrix_name, "BLOSUM62");
    assert!(resolved.sum_stats && !resolved.soft_masking);
    assert_eq!(resolved.max_intron_length, 0);
    assert_eq!(resolved.hitlist_size, 500);
    check_options(&args).unwrap();
    let resolved = parse(&["-soft_masking", "true", "-sum_stats", "off"])
        .unwrap()
        .resolve()
        .unwrap();
    assert!(resolved.soft_masking && !resolved.sum_stats);
    assert!(parse(&["-sum_stats", "maybe"]).is_err());
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:586-623: the matrix's suggestion (+2 for a translated
    // subject) replaces only the default 13.
    let resolved = parse(&["-matrix", "BLOSUM45"]).unwrap().resolve().unwrap();
    assert_eq!(
        (
            resolved.threshold,
            resolved.window_size,
            resolved.gap_open,
            resolved.gap_extend
        ),
        (16.0, 60, 14, 2)
    );
    for id in [
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        31, 32, 33,
    ] {
        let id_string = id.to_string();
        assert_eq!(parse(&["-db_gencode", &id_string]).unwrap().db_gencode, id);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1056 has 26 CLI values; gc.prt:340-347
// defines ID 32 for the TBLASTN-only product decision.
#[test]
fn invalid_codes_are_explicit_errors() {
    for id in ["0", "7", "8", "17", "20", "34", "255", "256", "abc", "-1"] {
        let error = parse(&["-db_gencode", id]).unwrap_err().to_string();
        assert!(error.contains("genetic code"), "{id}: {error}");
    }
}

// NCBI's checks of the options (c++/src/algo/blast/core/blast_options.c BLAST_ValidateOptions, the option
// handlers of blast_args.cpp) and LOSAT's limits, where NCBI starts the search.
#[test]
fn ncbi_checks_and_losat_limits_are_explicit() {
    assert!(error(&["-matrix", "BAD"]).contains("BAD is not a supported matrix"));
    assert!(error(&["-gapopen", "1", "-gapextend", "1"])
        .contains("Gap existence and extension values of 1 and 1 not supported for BLOSUM62"));
    assert!(error(&["-matrix", "IDENTITY", "-word_size", "6"])
        .contains("Word size larger than 5 is not supported for the identity scoring matrix"));
    assert!(
        error(&["-evalue", "0"]).contains("expect value or cutoff score must be greater than zero")
    );
    assert!(error(&["-word_size", "8"]).contains("Word-size must be less than 8"));
    assert!(error(&["-threshold", "0"]).contains("Non-zero threshold required"));
    assert!(error(&["-ungapped"]).contains("Composition-adjusted searched"));
    assert!(error(&["-seg", "12 2.2"]).contains("Invalid number of arguments"));
    assert!(parse(&["-word_size", "1"]).is_err());
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:255-275 and c++/src/algo/blast/core/blast_options.c:1174-1235: a matrix without gap costs
    // takes its best pair, and its suggested threshold (+2) and window.
    let resolved = parse(&["-matrix", "PAM30"]).unwrap().resolve().unwrap();
    assert_eq!(
        (
            resolved.gap_open,
            resolved.gap_extend,
            resolved.threshold,
            resolved.window_size
        ),
        (9, 1, 18.0, 15)
    );
    // `-ungapped` with composition-based statistics F (0): NCBI's ungapped search, which
    // LOSAT runs (line O4c).
    assert!(check_options(&parse(&["-ungapped", "-comp_based_stats", "F"]).unwrap()).is_ok());
    // NCBI's tblastn crashes with IDENTITY's compressed lookup table.
    for extra in [
        &["-task", "tblastn-fast", "-matrix", "IDENTITY"][..],
        &["-matrix", "identity", "-word_size", "5"],
    ] {
        assert!(
            limit(extra).contains("with which NCBI BLAST+'s tblastn crashes"),
            "{extra:?}"
        );
    }
    // Every matrix of NCBI's tables with every gap cost pair of its table, word sizes 2 to 4,
    // and composition modes 0 and 2 (series SC, SC1) and 1 and 3 (SC5).
    for extra in [
        &["-matrix", "PAM30"][..],
        &["-comp_based_stats", "1"],
        &["-comp_based_stats", "3"],
        &["-matrix", "PAM30", "-comp_based_stats", "1"],
        &[
            "-task",
            "tblastn-fast",
            "-matrix",
            "BLOSUM45",
            "-comp_based_stats",
            "3",
        ],
        &["-matrix", "pam250", "-comp_based_stats", "0"],
        &["-matrix", "BLOSUM45", "-gapopen", "13", "-gapextend", "3"],
        &["-matrix", "BLOSUM45", "-word_size", "3"],
        &["-matrix", "BLOSUM80", "-gapopen", "25", "-gapextend", "2"],
        &["-gapopen", "13", "-gapextend", "1"],
        &["-word_size", "2"],
        &["-word_size", "4", "-matrix", "PAM70"],
        &[
            "-matrix",
            "BLOSUM90",
            "-threshold",
            "11",
            "-window_size",
            "60",
        ],
        // The compressed lookup table with every matrix but IDENTITY (series SC, SC11).
        &["-task", "tblastn-fast", "-matrix", "BLOSUM80"],
        &[
            "-matrix",
            "pam30",
            "-word_size",
            "5",
            "-comp_based_stats",
            "0",
        ],
        &["-matrix", "BLOSUM45", "-word_size", "6", "-threshold", "15"],
        &[
            "-task",
            "tblastn-fast",
            "-matrix",
            "PAM250",
            "-word_size",
            "7",
        ],
        // The one-hit word finder (line RG, X2; aa_ungapped.c:214-224).
        &["-window_size", "0"],
        &["-word_size", "5", "-window_size", "0"],
        &["-task", "tblastn-fast", "-window_size", "0"],
    ] {
        check_options(&parse(extra).unwrap()).unwrap();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:606-620: the matrix's suggested threshold replaces tblastn's 13
    // only, not the compressed table's 19.3, 21 or 20.25 (nor tblastn-fast's 20); the
    // window is the matrix's (BLOSUM45's 60).
    let settings = search_settings(
        &parse(&["-matrix", "BLOSUM45", "-word_size", "6"])
            .unwrap()
            .resolve()
            .unwrap(),
    )
    .unwrap();
    assert_eq!(
        (
            settings.scoring.word_size,
            settings.scoring.threshold,
            settings.scoring.window
        ),
        (6, 21.0, 60)
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:836-879: the typed mode is the search's (`compo_adjust_mode`).
    for (value, mode) in [
        ("0", BlastCompoAdjustMode::NoCompositionBasedStats),
        ("F", BlastCompoAdjustMode::NoCompositionBasedStats),
        ("1", BlastCompoAdjustMode::CompositionBasedStats),
        ("2", BlastCompoAdjustMode::CompositionMatrixAdjust),
        ("T", BlastCompoAdjustMode::CompositionMatrixAdjust),
        ("3", BlastCompoAdjustMode::ForceFullMatrixAdjust),
    ] {
        let settings = search_settings(
            &parse(&["-comp_based_stats", value])
                .unwrap()
                .resolve()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(settings.compo_adjust_mode, mode, "{value}");
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-155: IDENTITY is searched in mode 0 with a warning (for the
    // first query batch), whatever mode is asked.
    for (extra, reset) in [
        (&["-matrix", "IDENTITY"][..], true),
        (&["-matrix", "identity", "-comp_based_stats", "1"], true),
        (&["-matrix", "IDENTITY", "-comp_based_stats", "3"], true),
        (&["-matrix", "IDENTITY", "-comp_based_stats", "0"], false),
        (&["-matrix", "BLOSUM62"], false),
    ] {
        let settings = search_settings(&parse(extra).unwrap().resolve().unwrap()).unwrap();
        assert_eq!(
            (
                settings.identity_cbs_reset,
                settings.compo_adjust_mode.is_enabled()
            ),
            (reset, !reset && extra.len() == 2),
            "{extra:?}"
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:436-447: tblastn-fast sets word size 5, the
    // compressed lookup table and BLAST_WORD_THRESHOLD_BLASTP_FAST (20); with BLOSUM62 11/1
    // the table of word size 3, 5, 6 or 7 takes any threshold and window size.
    let fast = parse(&["-task", "tblastn-fast"]).unwrap();
    let resolved = fast.resolve().unwrap();
    assert_eq!(
        (
            resolved.word_size,
            resolved.compressed_lookup,
            resolved.threshold
        ),
        (5, true, 20.0)
    );
    check_options(&fast).unwrap();
    for extra in [
        &[
            "-task",
            "tblastn-fast",
            "-word_size",
            "5",
            "-window_size",
            "15",
            "-threshold",
            "100",
        ][..],
        &["-word_size", "6"],
        &["-word_size", "7", "-threshold", "15"],
        &["-task", "tblastn-fast", "-word_size", "3"],
        &["-word_size", "3", "-window_size", "15", "-threshold", "100"],
        &["-threshold", "13.5", "-window_size", "1"],
    ] {
        check_options(&parse(extra).unwrap()).unwrap();
    }
    let settings =
        search_settings(&parse(&["-word_size", "7"]).unwrap().resolve().unwrap()).unwrap();
    assert_eq!(
        (settings.scoring.word_size, settings.scoring.threshold),
        (7, 20.25)
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:834-892: the first character chooses the mode (any other is
    // mode 0); a `u` asks for unified P-values only for blastp.
    for value in ["bogus", "0", "F", "2u", "D", "t"] {
        check_options(&parse(&["-comp_based_stats", value]).unwrap()).unwrap();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70: a hit list size whose
    // preliminary size overflows NCBI's `int` (2^30 without composition-based statistics, 2^30 - 25
    // with them, up to 2^31 - 1) runs; LOSAT computes the size without overflow.
    for extra in [
        &["-comp_based_stats", "0", "-max_target_seqs", "1073741824"][..],
        &["-comp_based_stats", "0", "-max_target_seqs", "2147483598"],
        &["-max_target_seqs", "1073741799"],
        &["-max_target_seqs", "2147483647"],
    ] {
        check_options(&parse(extra).unwrap()).unwrap();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:313-319: -qcov_hsp_perc and -max_hsps.
    let resolved = parse(&["-max_hsps", "2", "-qcov_hsp_perc", "50.5"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!((resolved.max_hsps, resolved.query_cov_hsp_perc), (2, 50.5));
    let resolved = parse(&[]).unwrap().resolve().unwrap();
    assert_eq!((resolved.max_hsps, resolved.query_cov_hsp_perc), (0, 0.0));
    for extra in [
        &["-max_hsps", "0"][..],
        &["-qcov_hsp_perc", "100.5"],
        &["-qcov_hsp_perc", "nan"],
    ] {
        assert!(parse(extra).is_err(), "{extra:?}");
    }
    for extra in [&["-gilist", "ids.gil"][..], &["-remote"]] {
        let error = (parse(extra).unwrap_err()).to_string();
        assert!(
            error.contains("is not supported by LOSAT's TBLASTN"),
            "{extra:?}: {error}"
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806,885
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
    // `-use_sw_tback` is a flag that sets the Smith-Waterman traceback, with every
    // `-comp_based_stats` (the redo runs with mode 0 too).
    for cbs in ["0", "1", "2", "3"] {
        let args = parse(&["-use_sw_tback", "-comp_based_stats", cbs]).unwrap();
        check_options(&args).unwrap();
        let resolved = args.resolve().unwrap();
        assert!(resolved.use_sw_tback, "{cbs}");
    }
    assert!(!parse(&[]).unwrap().resolve().unwrap().use_sw_tback);
    assert!(parse(&["-use_sw_tback", "true"]).is_err());
}

/// The tblastn arguments `-in_pssm <psiblast's frozen PSSM> -subject s.fna <extra>`.
fn parse_pssm(extra: &[&str]) -> std::result::Result<TblastnArgs, crate::cli::ArgvError> {
    let pssm = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
    );
    let mut argv = vec!["losat", "tblastn", "-in_pssm", pssm, "-subject", "s.fna"];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv)?;
    let Commands::Tblastn(args) = cli.command else {
        unreachable!()
    };
    Ok(args)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:138-149, c++/src/algo/blast/api/psiblast_options.cpp:69-75 and c++/src/algo/blast/blastinput/blast_args.cpp:866-868,
// 1050-1055: PSI-TBLASTN is psiblast's handle with tblastn's threshold and genetic code 1, no
// sum statistics, composition-based statistics 0 whatever `-comp_based_stats` says, and
// `-task tblastn-fast` without effect (line PS, DECISIONS 35).
#[test]
fn in_pssm_takes_psi_tblastns_options() {
    let resolved = parse_pssm(&[]).unwrap().resolve().unwrap();
    assert!(resolved.pssm.is_some());
    assert!(!resolved.sum_stats && !resolved.compressed_lookup);
    assert_eq!((resolved.word_size, resolved.threshold), (3, 13.0));
    assert_eq!(resolved.db_gencode, 1);
    assert_eq!(
        resolved.composition_mode,
        crate::blastinput::app::CompositionMode::NoCompositionBasedStats
    );
    for value in ["0", "1", "2", "3", "D", "t", "x", "", "2u"] {
        let resolved = parse_pssm(&["-comp_based_stats", value])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!(
            resolved.composition_mode,
            crate::blastinput::app::CompositionMode::NoCompositionBasedStats,
            "{value:?}"
        );
    }
    let fast = parse_pssm(&["-task", "tblastn-fast"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!(
        (fast.word_size, fast.threshold, fast.compressed_lookup),
        (3, 13.0, false)
    );
    let resolved = parse_pssm(&["-db_gencode", "4", "-sum_stats", "true"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!((resolved.db_gencode, resolved.sum_stats), (1, true));
    // The matrix's suggestion replaces the default threshold 13 (+2 for a translated subject).
    let resolved = parse_pssm(&["-matrix", "BLOSUM45"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!((resolved.threshold, resolved.window_size), (16.0, 60));
    // A matrix query keeps tblastn's handle.
    assert!(parse(&[]).unwrap().resolve().unwrap().pssm.is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:910-936 (`std_matrix_only` for a program that is neither blastp nor
// tblastn) and 1358-1363 (the word size of PSI-TBLASTN); the code 32 is not NCBI's
// (blast_args.cpp:1011-1015); `-query` and `-query_loc` exclude `-in_pssm`
// (blast_args.cpp:1241-1247).
#[test]
fn in_pssm_errors_are_ncbis_and_losats() {
    let error = |extra: &[&str]| {
        parse_pssm(extra)
            .unwrap()
            .resolve()
            .unwrap_err()
            .to_string()
    };
    for size in ["5", "6", "7", "8"] {
        assert_eq!(
            error(&["-word_size", size]),
            "LOSAT query/options error: Word-size must be less than 6 for protein comparison\nPlease refer to the LOSAT user manual.\n",
            "{size}"
        );
    }
    assert!(error(&["-task", "tblastn-fast", "-word_size", "6"])
        .contains("Word-size must be less than 6 for protein comparison"));
    let identity = error(&["-matrix", "IDENTITY"]);
    assert!(
        identity.starts_with(
            "LOSAT query/options error: IDENTITY is not a supported matrix, supported matrices are:\nBLOSUM80 \n"
        ),
        "{identity}"
    );
    assert!(!identity[40..].contains("IDENTITY"), "{identity}");
    assert!(error(&["-db_gencode", "32"]).contains("-db_gencode 32 with -in_pssm"));
    assert!(error(&["-threshold", "0"]).contains("Non-zero threshold required"));
    let pssm = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
    );
    for extra in [&["-query", "q.faa"][..], &["-query_loc", "1-10"]] {
        let mut argv = vec!["losat", "tblastn", "-in_pssm", pssm, "-subject", "s.fna"];
        argv.extend_from_slice(extra);
        assert!(try_parse_from(argv).is_err(), "{extra:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:781-815: the largest intron is read by the linking of the HSPs
// only (`core::link_hsps::translated_longest_intron`), so every value NCBI's command line
// takes is searched, with sum statistics and without them (the default of `-in_pssm`,
// `-sum_stats false`), and a negative one is the command line's error.
#[test]
fn the_largest_intron_is_accepted_with_and_without_sum_statistics() {
    let settings = |args: TblastnArgs| check_search_limits(&args.resolve().unwrap());
    for value in ["0", "1", "4", "5", "100", "2147483647"] {
        for args in [
            parse(&["-max_intron_length", value]).unwrap(),
            parse(&["-max_intron_length", value, "-sum_stats", "false"]).unwrap(),
            parse_pssm(&["-max_intron_length", value]).unwrap(),
            parse_pssm(&["-max_intron_length", value, "-sum_stats", "true"]).unwrap(),
        ] {
            assert!(settings(args).is_ok(), "{value}");
        }
    }
    let resolved = parse(&["-max_intron_length", "2147483647"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!(resolved.max_intron_length, i32::MAX);
    for value in ["-1", "2147483648", "1.5", "abc"] {
        assert!(parse(&["-max_intron_length", value]).is_err(), "{value}");
    }
}

// blast_args.cpp:276-278 and blast_options.c:608-618: `-xdrop_ungap` reaches the search's
// options (7 bits by default), one of 0 or below is NCBI's options error, and a gapped X-drop
// whose raw score NCBI's kernels cannot size is refused (`core::x_drop`).
#[test]
fn the_x_drops_reach_the_search_or_ncbis_errors() {
    let resolved = parse(&[]).unwrap().resolve().unwrap();
    assert_eq!(
        (
            resolved.xdrop_ungap,
            resolved.xdrop_gap,
            resolved.xdrop_gap_final
        ),
        (7.0, 15.0, 25.0)
    );
    let resolved = parse(&["-xdrop_ungap", "3.5", "-xdrop_gap", "1e300"])
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!((resolved.xdrop_ungap, resolved.xdrop_gap), (3.5, 1e300));
    assert_eq!(
        search_settings(&resolved).unwrap().scoring.word_xdrop_bits,
        3.5
    );
    for value in ["0", "-1", "-0.0", "-inf"] {
        let error = parse(&["-xdrop_ungap", value])
            .unwrap()
            .resolve()
            .unwrap_err();
        assert!(
            format!("{error:#}").contains("x_dropoff must be greater than zero"),
            "{value}: {error:#}"
        );
    }
    // 4.2e8 bits are 1.09e9 raw at BLOSUM62 11/1's Lambda, past the kernels' 2^30 - 103: the
    // preliminary X-drop at every mode, the final one only without the redo (mode 0, or
    // IDENTITY's reset). LOSAT refuses them where the search starts (`check_search_limits`),
    // which a run that NCBI stops before its first batch never reaches.
    for extra in [
        &["-xdrop_gap", "4.2e8"][..],
        &["-xdrop_gap_final", "4.2e8", "-comp_based_stats", "0"],
    ] {
        let resolved = parse(extra).unwrap().resolve().unwrap();
        assert!(search_settings(&resolved).is_ok());
        let error = check_search_limits(&resolved).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(&format!("an {} of 420000000.0 bits", extra[0])),
            "{extra:?}: {error:#}"
        );
    }
    for cbs in ["1", "2", "3"] {
        let resolved = parse(&["-xdrop_gap_final", "4.2e8", "-comp_based_stats", cbs])
            .unwrap()
            .resolve()
            .unwrap();
        assert!(check_search_limits(&resolved).is_ok());
    }
}

// NCBI BLAST+ 2.17.0 with `-ungapped -comp_based_stats 0` (line O4c): a `-matrix` name that
// NCBI's tables do not have reaches the setup of the first batch, which fails with NCBI's
// engine error (`SetupEngineError::unknown_matrix`); one that the matrix loader reads as
// IDENTITY without a table aborts there, which LOSAT refuses; a gapped search refuses both
// with NCBI's options error first.
#[test]
fn an_ungapped_search_keeps_a_matrix_of_no_table_for_its_setup() {
    use crate::algorithm::tblastx::args::TblastxMatrix;
    let ungapped = |matrix: &str| {
        parse(&["-ungapped", "-comp_based_stats", "0", "-matrix", matrix])
            .unwrap()
            .resolve()
            .unwrap()
    };
    let resolved = ungapped("FOO");
    assert_eq!(resolved.unknown_matrix, Some(TblastxMatrix::Unknown));
    assert!(check_search_limits(&resolved).is_ok());
    let scoring = search_settings(&resolved).unwrap().scoring;
    assert!(scoring.unknown_matrix);
    assert_eq!(scoring.matrix, ScoringMatrix::Blosum62);
    for matrix in ["BLOSUM62", "pam30", "IDENTITY"] {
        let resolved = ungapped(matrix);
        assert_eq!(resolved.unknown_matrix, None, "{matrix}");
        assert!(!search_settings(&resolved).unwrap().scoring.unknown_matrix);
    }
    let resolved = ungapped("identityx");
    assert_eq!(resolved.unknown_matrix, Some(TblastxMatrix::Aborts));
    let error = check_search_limits(&resolved).unwrap_err().to_string();
    assert!(
        error.contains("-matrix identityx is read as IDENTITY") && error.contains("TBLASTN"),
        "{error}"
    );
    for matrix in ["FOO", "IDENTITYX"] {
        assert!(
            parse(&["-matrix", matrix]).unwrap().resolve().is_err(),
            "{matrix}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328,3341-3346
// `-best_hit_overhang` / `-best_hit_score_edge`: above 0 and below 0.5, the one not given takes
// 0.1, neither given is no filter, and `-culling_limit` excludes both (NCBI refuses the command
// line; LOSAT's parser too).
#[test]
fn best_hit_options_are_read_as_ncbi_reads_them() {
    use crate::core::hspfilter_besthit::BestHitParams;
    let best_hit = |extra: &[&str]| parse(extra).unwrap().resolve().unwrap().best_hit;
    assert_eq!(best_hit(&[]), None);
    assert_eq!(
        best_hit(&["-best_hit_overhang", "0.25"]),
        Some(BestHitParams {
            overhang: 0.25,
            score_edge: 0.1
        })
    );
    assert_eq!(
        best_hit(&["-best_hit_score_edge", "0.05", "-best_hit_overhang", "0.3"]),
        Some(BestHitParams {
            overhang: 0.3,
            score_edge: 0.05
        })
    );
    for value in ["0", "0.5", "-0.1", "nan", "x"] {
        assert!(parse(&["-best_hit_overhang", value]).is_err(), "{value}");
        assert!(parse(&["-best_hit_score_edge", value]).is_err(), "{value}");
    }
    assert!(parse(&["-culling_limit", "1", "-best_hit_overhang", "0.1"]).is_err());
    assert!(parse(&["-culling_limit", "0", "-best_hit_score_edge", "0.1"]).is_err());
}
