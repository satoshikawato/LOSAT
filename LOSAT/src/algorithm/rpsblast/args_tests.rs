use super::*;
use crate::cli::{Cli, Commands};

fn parse(extra: &[&str]) -> Result<RpsblastArgs, String> {
    let mut argv = vec!["losat", "rpsblast", "-query", "q.faa", "-db", "cdd"];
    argv.extend_from_slice(extra);
    match crate::cli::try_parse_from::<Cli, _, _>(argv) {
        Ok(Cli {
            command: Commands::Rpsblast(args),
            ..
        }) => Ok(args),
        Ok(_) => Err("not rpsblast".to_string()),
        Err(error) => Err(error.to_string()),
    }
}

#[test]
fn the_composition_modes_read_rps_blasts_letters() {
    for (value, mode) in [
        ("0", 0),
        ("F", 0),
        ("f", 0),
        ("D", 0),
        ("d", 0),
        ("1", 1),
        ("T", 1),
        ("t", 1),
        ("2", 2),
        ("3", 3),
        ("x", 0),
        ("", 0),
        ("1u", 1),
    ] {
        assert_eq!(rps_composition_mode(value), mode, "{value:?}");
    }
}

#[test]
fn the_search_is_blastps_with_rpsblasts_program() {
    let args = parse(&[
        "-comp_based_stats",
        "F",
        "-num_threads",
        "0",
        "-mt_mode",
        "1",
    ])
    .unwrap()
    .into_blastp_args()
    .unwrap();
    assert_eq!(args.app, BlastpApp::Rpsblast);
    assert_eq!(args.comp_based_stats, "0");
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.rps.as_ref().unwrap().requested_cbs, 0);
    let args = parse(&["-max_hsps", "2", "-qcov_hsp_perc", "50"])
        .unwrap()
        .into_blastp_args()
        .unwrap();
    // RPS-BLAST's defaults (blastinput_unit_test.cpp:3296-3303 `RpsBlastCBS`): composition-based
    // statistics 1 and no SEG (`blast_rps_options.cpp:68-73`).
    assert_eq!(args.seg, None);
    assert_eq!(
        args.resolve().unwrap().seg,
        crate::algorithm::blastp::args::BlastpSegSpec::No
    );
    let rps = args.rps.unwrap();
    assert_eq!(rps.requested_cbs, 1);
    assert_eq!(
        (rps.max_hsps_per_subject, rps.query_cov_hsp_perc),
        (2, 50.0)
    );
    assert_eq!(args.max_hsps_per_subject, None);
}

#[test]
fn options_of_later_items_are_refused() {
    for extra in [&["-outfmt", "8"][..], &["-outfmt", "1"]] {
        let error = parse(extra).unwrap().into_blastp_args().unwrap_err();
        assert!(
            error
                .to_string()
                .contains("not supported by LOSAT's RPSBLAST"),
            "{extra:?}: {error}"
        );
    }
    let error = parse(&["-html", "1"]).unwrap_err();
    assert!(
        error.contains("-html is not supported by LOSAT's RPSBLAST"),
        "{error}"
    );
    assert!(parse(&["-subject", "s.faa"]).is_err());
    assert!(parse(&["-mt_mode", "2"]).is_err());
}

// blast_args.cpp:234-238 and 303-307 (cited at the field), ncbiargs.cpp:118-129: `-searchsp`
// takes NCBI's `Int8` forms at 0 or above and, at any value (0 too), the search runs without
// the Gumbel block; the same parser as psiblast's and blastp's.
#[test]
fn searchsp_takes_ncbis_int8_forms_and_drops_the_gumbel_block() {
    for (value, space) in [
        ("+1000000000", 1_000_000_000),
        ("0x3B9ACA00", 1_000_000_000),
        ("9223372036854775807", i64::MAX),
        ("0", 0),
    ] {
        let args = parse(&["-searchsp", value]).unwrap();
        assert_eq!(args.searchsp, Some(space), "{value}");
        let resolved = args.into_blastp_args().unwrap().resolve().unwrap();
        assert_eq!(resolved.searchsp, Some(space), "{value}");
        assert!(resolved.old_fsc, "{value}");
    }
    assert_eq!(parse(&[]).unwrap().searchsp, None);
    for value in [" 1", "1e6", "1.5", "", "-1", "-5", "9223372036854775808"] {
        assert!(
            parse(&["-searchsp", value]).is_err(),
            "-searchsp {value:?} is refused"
        );
    }
}

// blast_args.cpp:209-244, 332-342, 1939-1942, 2683-2721 (cited at the fields): the masks, the
// X-drops, the formatting options and `-sum_stats false` reach blastp's arguments.
#[test]
fn the_options_of_ps11c_1_reach_blastps_arguments() {
    let args = parse(&[
        "-seg",
        "10 1.8 2.1",
        "-soft_masking",
        "true",
        "-lcase_masking",
        "-xdrop_ungap",
        "3",
        "-xdrop_gap",
        "-1",
        "-xdrop_gap_final",
        "0",
        "-line_length",
        "1",
        "-sorthits",
        "4",
        "-sorthsps",
        "0",
        "-sum_stats",
        "F",
    ])
    .unwrap()
    .into_blastp_args()
    .unwrap();
    assert_eq!(args.seg.as_deref(), Some("10 1.8 2.1"));
    assert_eq!(args.soft_masking, Some(true));
    assert!(args.lcase_masking);
    assert_eq!(
        args.formatting,
        FormattingArgs {
            line_length: Some(1),
            sort_hits: Some(4),
            sort_hsps: Some(0),
        }
    );
    let rps = args.rps.as_ref().unwrap();
    assert_eq!(
        (
            rps.x_drop_ungapped_bits,
            rps.x_drop_gapped_bits,
            rps.x_drop_final_bits
        ),
        (Some(3.0), Some(-1.0), Some(0.0))
    );
    assert!(!rps.sum_stats);
    let resolved = args.resolve().unwrap();
    assert_eq!(resolved.line_length, 1);
    assert!(resolved.lcase_masking && resolved.soft_masking);
    // Nothing given: NCBI's defaults.
    let args = parse(&[]).unwrap().into_blastp_args().unwrap();
    assert_eq!(args.formatting, FormattingArgs::default());
    assert_eq!(args.resolve().unwrap().line_length, 60);
    let rps = args.rps.unwrap();
    assert_eq!(
        (
            rps.x_drop_ungapped_bits,
            rps.x_drop_gapped_bits,
            rps.x_drop_final_bits
        ),
        (None, None, None)
    );
}

// blast_args.cpp:2683-2721: `-line_length` is 1 or more, `-sorthits` and `-sorthsps` 0 to 4
// (LOSAT's parser text and rc 2, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 1); the X-drops are
// NCBI's doubles (`nan` and `inf` are its conversion errors).
#[test]
fn the_formatting_and_x_drop_values_are_ncbis() {
    for extra in [
        &["-line_length", "0"][..],
        &["-line_length", "-3"],
        &["-sorthits", "5"],
        &["-sorthsps", "-1"],
        &["-sorthsps", "x"],
        &["-xdrop_gap", "nan"],
        &["-xdrop_ungap", "inf"],
        &["-sum_stats", "maybe"],
    ] {
        assert!(parse(extra).is_err(), "{extra:?}");
    }
    assert_eq!(
        parse(&["-line_length", "2147483647"])
            .unwrap()
            .formatting
            .line_length,
        Some(2_147_483_647)
    );
    // `-sum_stats true` reaches the search, which refuses it (`rps_hooks::rps_search`).
    assert!(
        parse(&["-sum_stats", "true"])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .rps
            .unwrap()
            .sum_stats
    );
    assert_eq!(
        parse(&["-xdrop_gap", "1e300"]).unwrap().xdrop_gap,
        Some(1e300)
    );
}

// blast_options.c:608-618 (cited at `validate_protein_options`): an ungapped X-drop of 0 or
// below is NCBI's options error, made when the options are validated (after the SEG
// option's own error, which its handler raises).
#[test]
fn an_ungapped_x_drop_of_zero_or_below_is_ncbis_options_error() {
    for value in ["0", "-1", "-0.0"] {
        let error = parse(&["-xdrop_ungap", value])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "LOSAT query/options error: x_dropoff must be greater than zero\nPlease refer to the LOSAT user manual.\n",
            "{value}"
        );
    }
    let error = parse(&["-xdrop_ungap", "0", "-seg", "maybe"])
        .unwrap()
        .into_blastp_args()
        .unwrap()
        .resolve()
        .unwrap_err();
    assert!(
        format!("{error:#}").contains("Invalid number of arguments to filtering option"),
        "{error:#}"
    );
    // The hit saving options are checked after the X-drop.
    let error = parse(&[
        "-xdrop_ungap",
        "0",
        "-num_descriptions",
        "0",
        "-num_alignments",
        "0",
    ])
    .unwrap()
    .into_blastp_args()
    .unwrap()
    .resolve()
    .unwrap_err();
    assert!(
        format!("{error:#}").contains("x_dropoff must be greater than zero"),
        "{error:#}"
    );
    assert!(parse(&["-xdrop_ungap", "1e-300"])
        .unwrap()
        .into_blastp_args()
        .unwrap()
        .resolve()
        .is_ok());
}
