use super::super::pairwise::query_ungapped_karlin;
use super::*;

fn args(words: &[&str]) -> BlastnArgs {
    let argv = [
        "losat", "blastn", "-query", "q", "-subject", "s", "-outfmt", "6",
    ];
    let cli: crate::cli::Cli =
        crate::cli::try_parse_from(argv.iter().chain(words)).expect("valid arguments");
    let crate::cli::Commands::Blastn(args) = cli.command else {
        panic!("blastn")
    };
    args
}

fn message(result: anyhow::Result<()>) -> String {
    String::from_utf8(
        result
            .expect_err("rejected")
            .downcast::<NativeError>()
            .expect("an NCBI error")
            .message,
    )
    .unwrap()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1247-1261,1399-1413
#[test]
fn discontiguous_templates_are_checked_after_the_word_size_limit() {
    let word = "Invalid discontiguous template parameters: word size must be either 11 or 12";
    let table = "Invalid lookup table type for discontiguous Mega BLAST";
    assert!(check_scoring_options(&args(&["-task", "dc-megablast"])).is_ok());
    assert!(check_scoring_options(&args(&["-task", "dc-megablast", "-word_size", "12"])).is_ok());
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-word_size",
        "13"
    ])))
    .contains(word));
    // The template applies to any task: megablast keeps its MB table, blastn and
    // blastn-short have the standard table.
    let template = ["-template_type", "coding", "-template_length", "18"];
    let with = |words: &[&str]| args(&[words, &template[..]].concat());
    assert!(check_scoring_options(&with(&["-word_size", "11"])).is_ok());
    assert!(message(check_scoring_options(&with(&[]))).contains(word));
    assert!(message(check_scoring_options(&with(&["-task", "blastn"]))).contains(table));
    assert!(message(check_scoring_options(&with(&["-task", "blastn-short"]))).contains(word));
    assert!(message(check_scoring_options(&with(&[
        "-task",
        "blastn-short",
        "-word_size",
        "12"
    ])))
    .contains(table));
    // The word size limit (100) comes first, the zero gap costs last.
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-word_size",
        "101"
    ])))
    .contains("Word-size must be less than or equal to 100"));
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-word_size",
        "13",
        "-gapopen",
        "0",
        "-gapextend",
        "0"
    ])))
    .contains(word));
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-gapopen",
        "0",
        "-gapextend",
        "0"
    ])))
    .contains("Greedy extension must be used"));
}

// DI-6: `-no_greedy` makes megablast's gap costs 0/0 invalid (blast_options.c:1699-1711),
// and an off-diagonal range needs a window (blast_options.c:1583-1596), checked after the
// lookup table's options.
#[test]
fn no_greedy_and_off_diagonal_range_are_checked() {
    assert!(message(check_scoring_options(&args(&["-no_greedy"])))
        .contains("Greedy extension must be used"));
    assert!(
        check_scoring_options(&args(&["-no_greedy", "-gapopen", "5", "-gapextend", "2"])).is_ok()
    );
    assert!(check_scoring_options(&args(&["-task", "blastn", "-no_greedy"])).is_ok());
    assert_eq!(
        message(check_scoring_options(&args(&["-off_diagonal_range", "2"]))),
        "LOSAT query/options error: off_diagonal_range is only useful in 2-hit algorithm\nPlease refer to the LOSAT user manual.\n"
    );
    assert!(check_scoring_options(&args(&["-off_diagonal_range", "0"])).is_ok());
    assert!(
        check_scoring_options(&args(&["-off_diagonal_range", "2", "-window_size", "40"])).is_ok()
    );
    assert!(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-off_diagonal_range",
        "2"
    ]))
    .is_ok());
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "dc-megablast",
        "-window_size",
        "0",
        "-off_diagonal_range",
        "2"
    ])))
    .contains("off_diagonal_range is only useful"));
    // The lookup table's check comes first.
    assert!(message(check_scoring_options(&args(&[
        "-off_diagonal_range",
        "2",
        "-word_size",
        "101"
    ])))
    .contains("Word-size must be less than"));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518: the e-value check needs a cutoff score of 0 or below too.
#[test]
fn min_raw_gapped_score_above_zero_allows_an_evalue_of_zero() {
    const REFUSED: &str = "LOSAT query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the LOSAT user manual.\n";
    for words in [
        &["-evalue", "0"][..],
        &["-evalue", "0", "-min_raw_gapped_score", "0"],
        &["-evalue", "-1", "-min_raw_gapped_score", "-5"],
    ] {
        assert_eq!(message(check_scoring_options(&args(words))), REFUSED);
    }
    for words in [
        &["-evalue", "0", "-min_raw_gapped_score", "1"][..],
        &["-evalue", "-1", "-min_raw_gapped_score", "40"],
        &["-min_raw_gapped_score", "0"],
        &["-min_raw_gapped_score", "-5"],
    ] {
        assert!(check_scoring_options(&args(words)).is_ok(), "{words:?}");
    }
}

#[test]
fn options_are_checked_in_ncbi_order() {
    assert!(check_scoring_options(&args(&[])).is_ok());
    assert!(check_scoring_options(&args(&["-task", "blastn"])).is_ok());
    assert_eq!(
        message(check_scoring_options(&args(&[
            "-task", "blastn", "-penalty", "0", "-gapopen", "0", "-gapextend", "0"
        ]))),
        "LOSAT query/options error: BLASTN penalty must be negative\nPlease refer to the LOSAT user manual.\n"
    );
    assert!(message(check_scoring_options(&args(&["-gapopen", "3"])))
        .contains("BLASTN gap extension penalty cannot be 0"));
    assert!(message(check_scoring_options(&args(&[
        "-task",
        "blastn",
        "-gapopen",
        "0",
        "-gapextend",
        "0"
    ])))
    .contains("Greedy extension must be used"));
    // LOSAT's limits are separate: NCBI's checks accept these options.
    for words in [
        &["-reward", "0"][..],
        &["-reward", "32767", "-penalty", "-1"],
        &["-reward", "16384", "-penalty", "-32768"],
        // 98303 is 32767 as NCBI's 16-bit reward.
        &["-reward", "98303", "-penalty", "-2"],
    ] {
        let parsed = args(words);
        assert!(check_scoring_options(&parsed).is_ok(), "{words:?}");
        let error = check_losat_limits(&parsed).unwrap_err().to_string();
        assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
    }
    // Scores below the 16-bit limits are accepted (E2g V1), and so are the infinite
    // and NaN e-values that NCBI searches with (E2g), and reward 0 with penalty 0
    // (matrix-only scoring, line RG X6).
    for words in [
        &["-reward", "5000", "-penalty", "-1"][..],
        &["-reward", "32766", "-penalty", "-32767"],
        &["-reward", "0", "-penalty", "0"],
        &["-task", "rmblastn"],
        &["-task", "blastn", "-reward", "0", "-penalty", "0"],
        &["-evalue", "+inf"],
        &["-evalue", "+nan"],
        &["-evalue", "-nan(7)"],
    ] {
        assert!(check_losat_limits(&args(words)).is_ok(), "{words:?}");
    }
}

/// The checks of the options that a run makes (`process_options`): the counts of the formats
/// `[0]` from `report_counts_of`, then `check_scoring_options` and `check_losat_limits`.
fn run_checks(words: &[&str]) -> anyhow::Result<()> {
    let mut args = args(words);
    args.report_counts = report_counts_of(&args, &[0]).unwrap().0;
    check_scoring_options(&args)?;
    check_losat_limits(&args)
}

/// `run_checks` for a run of the formats `outfmts`, gapped or not (`ungapped`).
fn run_checks_with(words: &[&str], outfmts: &[i32], ungapped: bool) -> anyhow::Result<()> {
    let args = args(words);
    let (counts, _) = report_counts_of(&args, outfmts).unwrap();
    check_scoring_options_with(&args, &counts, ungapped)?;
    check_losat_limits(&args)
}

#[test]
fn the_counts_a_run_derives_give_the_refusals_of_the_options() {
    let error = message(run_checks(&[
        "-num_descriptions",
        "0",
        "-num_alignments",
        "0",
    ]));
    assert_eq!(
        error,
        "LOSAT query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n"
    );
    // The hit list of format 0 is the larger of the two numbers.
    assert!(run_checks(&["-num_alignments", "0"]).is_ok());
    assert!(run_checks(&["-num_descriptions", "0"]).is_ok());
    assert!(run_checks(&["-num_descriptions", "1", "-num_alignments", "0"]).is_ok());
    // A hit list size whose preliminary size overflows NCBI's `int` runs (LOSAT computes the
    // size without overflow; PD-LOSAT-NCBI-DEFECTS, hit list sizes past `Int4`).
    for words in [
        &["-num_alignments", "2147483598"][..],
        &["-num_descriptions", "2147483647"],
        &["-max_target_seqs", "2147483647"],
    ] {
        assert!(run_checks(words).is_ok(), "{words:?}");
    }
    // Without the counts of a run (`args.report_counts` at its default), no hit list is empty.
    assert!(
        check_scoring_options(&args(&["-num_descriptions", "0", "-num_alignments", "0"])).is_ok()
    );
}

#[test]
fn table_errors_repeat_for_each_query_of_the_batch() {
    let error = karlin_error("Substitution scores 1 and -6 are not supported", [false; 2]);
    assert_eq!(
        error.downcast::<NativeError>().unwrap().message,
        b"LOSAT engine error: Error: Substitution scores 1 and -6 are not supported Error: Substitution scores 1 and -6 are not supported \n"
    );
    // A query without letters keeps its warning after the error (NCBI 2.17.0).
    let error = karlin_error(
        "Substitution scores 1 and -6 are not supported",
        [true, false],
    );
    assert_eq!(
        error.downcast::<NativeError>().unwrap().message,
        b"LOSAT engine error: Error: Substitution scores 1 and -6 are not supported Warning: Sequence contains no data Error: Substitution scores 1 and -6 are not supported \n"
    );
}

#[test]
fn no_valid_query_needs_no_gapped_block() {
    let spec = NuclScoringSpec {
        reward: 1,
        penalty: -6,
        gap_open: 0,
        gap_extend: 0,
    };
    assert!(context_blocks(&[None, None], &spec).is_ok());
    let ungapped = query_ungapped_karlin(b"ACGTACGT", 1, -6);
    assert!(context_blocks(&[ungapped, ungapped], &spec).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3955-3963, 3989-3992, 4024-4029
// (oracle `-ungapped`: 2/-3 beta -2, 1/-1 beta -2, 2/-2 and 4/-6 beta 0, 1/-6 and 3/-1 no table
// row: alpha = beta = 0 and a length adjustment of 0, effective search space 377 x 373)
#[test]
fn ungapped_blocks_take_alpha_beta_from_the_ungapped_block() {
    let spec = |reward, penalty| NuclScoringSpec {
        reward,
        penalty,
        gap_open: 0,
        gap_extend: 0,
    };
    let block = KarlinParams {
        lambda: 0.634,
        k: 0.408,
        h: 0.912,
        alpha: 7.0,
        beta: 7.0,
    };
    for (reward, penalty, beta) in [
        (2, -3, -2.0),
        (1, -1, -2.0),
        (2, -2, 0.0),
        (4, -6, 0.0),
        (1, -2, 0.0),
    ] {
        let (blocks, round_down) =
            context_blocks_ungapped(&[Some(block), None], &spec(reward, penalty));
        assert!(!round_down);
        assert!(blocks[1].is_none());
        let blocks = blocks[0].unwrap();
        assert_eq!(blocks.ungapped.lambda, block.lambda);
        assert_eq!(blocks.gapped.lambda, block.lambda);
        assert_eq!(blocks.gapped.k, block.k);
        assert_eq!(blocks.gapped.h, block.h);
        assert_eq!(
            blocks.gapped.alpha,
            block.lambda / block.h,
            "{reward}/{penalty}"
        );
        assert_eq!(blocks.gapped.beta, beta, "{reward}/{penalty}");
    }
    // A pair with no table row (1/-6, 3/-1) is no error: alpha = beta = 0.
    for (reward, penalty) in [(1, -6), (3, -1)] {
        let (blocks, _) = context_blocks_ungapped(&[Some(block)], &spec(reward, penalty));
        let blocks = blocks[0].unwrap();
        assert_eq!((blocks.gapped.alpha, blocks.gapped.beta), (0.0, 0.0));
        assert!(context_blocks(&[Some(block)], &spec(reward, penalty)).is_err());
    }
    // No valid context: no blocks, as `context_blocks`.
    assert!(context_blocks_ungapped(&[None, None], &spec(1, -6))
        .0
        .iter()
        .all(Option::is_none));
}

/// Oracle (`t2`, 377 x 373, 1/-6 `-ungapped`): the effective search space is 377 x 373 =
/// 140,621, so the length adjustment of alpha = beta = 0 is 0.
#[test]
fn ungapped_blocks_without_a_table_row_have_no_length_adjustment() {
    let spec = NuclScoringSpec {
        reward: 1,
        penalty: -6,
        gap_open: 0,
        gap_extend: 0,
    };
    let block = KarlinParams {
        lambda: 1.39,
        k: 0.749,
        h: 1.38,
        alpha: 0.0,
        beta: 0.0,
    };
    let (blocks, _) = context_blocks_ungapped(&[Some(block)], &spec);
    let adjustment = crate::stats::length_adjustment::compute_length_adjustment_ncbi(
        377,
        373,
        1,
        &blocks[0].unwrap().gapped,
    );
    assert_eq!(adjustment.length_adjustment, 0);
}

#[test]
fn ungapped_option_checks_follow_gapped_calculation() {
    let zero_extension = args(&["-gapopen", "3"]);
    let zero_costs = args(&["-task", "blastn", "-gapopen", "0", "-gapextend", "0"]);
    let counts = ReportCounts::default();
    // "gap extension penalty cannot be 0" only in a gapped search.
    assert!(check_scoring_options_with(&zero_extension, &counts, false).is_err());
    assert!(check_scoring_options_with(&zero_extension, &counts, true).is_ok());
    // "Greedy extension must be used ..." is raised in both.
    assert!(
        message(check_scoring_options_with(&zero_costs, &counts, true))
            .contains("Greedy extension must be used")
    );
    // Gap costs without a table row (gapped: exit 3) and an unsupported pair are accepted.
    for words in [
        &["-reward", "1", "-penalty", "-6"][..],
        &["-reward", "3", "-penalty", "-1"],
        &["-task", "blastn", "-gapopen", "3", "-gapextend", "7"],
        &["-gapopen", "5", "-gapextend", "0"],
        &["-gapopen", "40000", "-gapextend", "40000"],
    ] {
        assert!(run_checks_with(words, &[6], true).is_ok(), "{words:?}");
    }
    // (A gapped search of 1/-6 stops at its table: `context_blocks`,
    // `no_valid_query_needs_no_gapped_block`.)
    // A reward of 0 with a negative penalty (megablast's -2) stays refused; with penalty 0 it
    // is matrix-only scoring.
    assert!(run_checks_with(&["-reward", "0"], &[6], true).is_err());
    assert!(run_checks_with(&["-reward", "0", "-penalty", "0"], &[6], true).is_ok());
    // The preliminary hit list size of 2147483598 runs, ungapped (the size itself, as NCBI) and
    // gapped (without NCBI's overflow).
    assert!(run_checks_with(&["-num_alignments", "2147483598"], &[0], true).is_ok());
    assert!(run_checks_with(&["-num_alignments", "2147483598"], &[0], false).is_ok());
}

// NCBI 2.17.0's greedy offsets (`greedy_xdrop_offset`, x86-64 `Int4`) and the refusals of
// the X-drops that NCBI cannot run (line RG X1b probe `b2.*`).
#[test]
fn gapped_x_drops_that_ncbis_kernels_cannot_run_are_refused() {
    use super::super::alignment::greedy::greedy_xdrop_offset;
    // megablast 1/-2, linear (0/0): the X-drop doubles, (2x + 1) / 6 + 1.
    assert_eq!(greedy_xdrop_offset(13, 1, -2, 0, 0), 5);
    assert_eq!(greedy_xdrop_offset(-1, 1, -2, 0, 0), 1);
    assert_eq!(greedy_xdrop_offset(-3, 1, -2, 0, 0), 1);
    // -8 bits = raw -4: the offset 0 reads `max_score` before NCBI writes it.
    assert_eq!(greedy_xdrop_offset(-4, 1, -2, 0, 0), 0);
    // INT_MIN (past Int4) doubles to 0; 2^30 doubles to INT_MIN; 2^31 - 1 doubles to -2.
    assert_eq!(greedy_xdrop_offset(i32::MIN, 1, -2, 0, 0), 1);
    assert!(greedy_xdrop_offset(1 << 30, 1, -2, 0, 0) <= 0);
    assert_eq!(greedy_xdrop_offset(i32::MAX, 1, -2, 0, 0), 1);
    // An even reward is not doubled: INT_MIN and INT_MAX (whose sum wraps) leave no offset.
    assert_eq!(greedy_xdrop_offset(-1, 2, -3, 0, 0), 1);
    assert!(greedy_xdrop_offset(i32::MIN, 2, -3, 0, 0) <= 0);
    assert!(greedy_xdrop_offset(i32::MAX, 2, -3, 0, 0) <= 0);
    // Affine 1/-2 2/2: doubled 2/4 4/4, common factor gdb3(6, 4, 5) = 1.
    assert_eq!(greedy_xdrop_offset(13, 1, -2, 2, 2), 28);
    assert_eq!(greedy_xdrop_offset(i32::MIN, 1, -2, 2, 2), 2);
    // Large megablast gap costs above 32767 up to i32::MAX (widened to 64 bits, no overflow).
    assert_eq!(greedy_xdrop_offset(13, 1, -2, 32_768, 32_768), 28);
    assert_eq!(greedy_xdrop_offset(13, 1, -2, 40_000, 40_000), 28);
    assert_eq!(greedy_xdrop_offset(13, 1, -2, 1_000_000, 1_000_000), 28);
    assert_eq!(
        greedy_xdrop_offset(13, 1, -2, 2_000_000_000, 2_000_000_000),
        28
    );
    assert_eq!(greedy_xdrop_offset(13, 1, -2, i32::MAX, i32::MAX), 28);
    assert_eq!(greedy_xdrop_offset(13, 1, -2, i32::MAX, 2), 28);
    // Even reward with high GCD (1000/-2000 40000/2000 -> gcd 500): (13000 + 500) / 500 + 1 = 28.
    assert_eq!(greedy_xdrop_offset(13_000, 1000, -2000, 40_000, 2000), 28);
    assert_eq!(greedy_xdrop_offset(20, 2, -3, i32::MAX, i32::MAX), 22);

    let check = |use_dp, (gapped, finale), (reward, penalty), (gap_open, gap_extend)| {
        check_gapped_x_dropoffs(
            use_dp,
            false,
            (25.0, 100.0),
            (gapped, finale),
            (reward, penalty),
            (gap_open, gap_extend),
            "BLASTN",
        )
    };
    let refusal = |result: anyhow::Result<()>| result.expect_err("refused").to_string();
    // Greedy: the defaults, INT_MIN at either stage, large gap costs, and the top of Int4 run.
    assert!(check(false, (13, 54), (1, -2), (0, 0)).is_ok());
    assert!(check(false, (13, 54), (1, -2), (40_000, 40_000)).is_ok());
    assert!(check(false, (13, 54), (1, -2), (i32::MAX, i32::MAX)).is_ok());
    assert!(check(false, (13, 54), (1, -2), (i32::MAX, 2)).is_ok());
    assert!(check(false, (20, 80), (2, -3), (i32::MAX, i32::MAX)).is_ok());
    assert!(check(false, (i32::MIN, 54), (1, -2), (0, 0)).is_ok());
    assert!(check(false, (13, i32::MIN), (1, -2), (0, 0)).is_ok());
    assert!(
        refusal(check(false, (-4, 54), (1, -2), (0, 0))).starts_with("an -xdrop_gap of 25.0 bits")
    );
    assert!(refusal(check(false, (13, 1 << 30), (1, -2), (0, 0)))
        .starts_with("an -xdrop_gap_final of 100.0 bits"));
    // The final X-drop that the preliminary one set is named -xdrop_gap.
    assert!(refusal(check(false, (1 << 30, 1 << 30), (1, -2), (0, 0)))
        .starts_with("an -xdrop_gap of 25.0 bits"));
    assert!(check(false, (i32::MIN, 54), (2, -3), (0, 0)).is_err());
    // -xdrop_gap -3e9 (raw -1624563704) doubles to 1045839888: its offset 174306649 passes
    // the 19 cells that the final X-drop's 54 gives NCBI's array (NCBI: SIGSEGV). With a
    // final X-drop as large in NCBI's arithmetic (-3e9 or 1e9 bits), NCBI runs it.
    assert!(refusal(check(false, (-1_624_563_704, 54), (1, -2), (0, 0)))
        .starts_with("an -xdrop_gap of 25.0 bits"));
    assert!(check(false, (-1_624_563_704, -1_624_563_704), (1, -2), (0, 0)).is_ok());
    assert!(check(false, (-1_624_563_704, 541_521_234), (1, -2), (0, 0)).is_ok());
    // Dynamic programming: only the final X-drop, against the Int4 sizes of `ALIGN_EX`; the
    // preliminary one is never past the HSP's ungapped score, INT_MIN rises to the gap costs.
    assert!(check(true, (i32::MIN, i32::MIN), (2, -3), (5, 2)).is_ok());
    assert!(check(true, (-100, 2_147_483_441), (2, -3), (5, 2)).is_ok());
    assert!(check(true, (-100, 2_147_483_442), (2, -3), (5, 2)).is_err());
    assert!(check(true, (13, (1 << 30) - 104), (1, -2), (2, 1)).is_ok());
    assert!(refusal(check(true, (13, (1 << 30) - 103), (1, -2), (2, 1)))
        .contains("-xdrop_gap_final of 100.0 bits"));
    // Matrix-only scoring (option reward 0 and penalty 0, engine pair 1/-3): the options are
    // the raw X-drops (blast_parameters.c:472-482), checked as dynamic programming's final
    // one, which keeps its own option's name when it equals the preliminary one.
    let matrix_only = |(gap_bits, final_bits): (f64, f64), (gap_open, gap_extend)| {
        check_gapped_x_dropoffs(
            false,
            true,
            (gap_bits, final_bits),
            (gap_bits as i32, final_bits as i32),
            (1, -3),
            (gap_open, gap_extend),
            "BLASTN",
        )
    };
    assert!(matrix_only((30.0, 100.0), (5, 2)).is_ok());
    assert!(matrix_only((30.0, 0.0), (5, 2)).is_ok());
    assert!(matrix_only((-100.0, 2_147_483_441.0), (5, 2)).is_ok());
    assert_eq!(
        refusal(matrix_only((-100.0, 2_147_483_442.0), (5, 2))),
        "an -xdrop_gap_final of 2147483442.0 (a raw score with reward 0 and penalty 0), whose raw X-drop makes the sizes of NCBI BLAST+'s gapped alignment kernels leave 32-bit integers (undefined behaviour), is not supported by LOSAT's BLASTN"
    );
    assert!(
        refusal(matrix_only((2_147_483_442.0, 2_147_483_442.0), (5, 2)))
            .starts_with("an -xdrop_gap_final of 2147483442.0 (a raw score")
    );
}

// NCBI 2.17.0 (line RG X6 probes, ws:~/losat-bench/work/rgx6): megablast's greedy extension
// with reward 0 and penalty 0 divides by zero with gap costs 0 and 0 (SIGFPE, exit 136, also
// with -ungapped) and does not end with other gap costs; -no_greedy and -ungapped with other
// gap costs run, as every other task does.
#[test]
fn matrix_only_scoring_with_the_greedy_extension_is_refused() {
    let refused = |words: &[&str], ungapped| {
        let error = check_greedy_matrix_only(&args(words), ungapped)
            .expect_err("refused")
            .to_string();
        assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
    };
    let runs = |words: &[&str], ungapped| {
        assert!(
            check_greedy_matrix_only(&args(words), ungapped).is_ok(),
            "{words:?} {ungapped}"
        );
    };
    let r0p0 = ["-reward", "0", "-penalty", "0"];
    let affine = [&r0p0[..], &["-gapopen", "5", "-gapextend", "2"]].concat();
    refused(&r0p0, false);
    refused(&r0p0, true);
    refused(&affine, false);
    runs(&affine, true);
    runs(&[&affine[..], &["-no_greedy"]].concat(), false);
    for task in ["blastn", "dc-megablast", "blastn-short"] {
        runs(&[&r0p0[..], &["-task", task]].concat(), false);
    }
    runs(&["-task", "rmblastn"], false);
    runs(&["-task", "rmblastn"], true);
    // megablast's own scores extend greedily.
    runs(&[], false);
}
