use super::*;
use crate::cli::{Cli, Commands};

fn parse(extra: &[&str]) -> Result<PsiblastArgs, String> {
    let mut argv = vec!["losat", "psiblast", "-query", "q.faa", "-subject", "s.faa"];
    argv.extend_from_slice(extra);
    match crate::cli::try_parse_from(argv) {
        Ok(Cli {
            command: Commands::Psiblast(args),
        }) => Ok(args),
        Ok(other) => panic!("not psiblast: {other:?}"),
        Err(error) => Err(error.to_string()),
    }
}

// psiblast_args.cpp:44-118 and psiblast_options.cpp:46-60: blastp's task options, the
// program psiblast, and the PSSM engine's defaults (inclusion threshold 0.002, pseudo-count 0,
// one round).
#[test]
fn the_first_round_is_blastps_search_with_psiblasts_program() {
    let args = parse(&[]).unwrap();
    assert_eq!(args.num_iterations, 1);
    assert_eq!(args.pseudocount, 0);
    assert_eq!(args.inclusion_ethresh, 0.002);
    let blastp = args.into_blastp_args().unwrap();
    assert_eq!(blastp.app, BlastpApp::Psiblast);
    assert_eq!(blastp.task, "blastp");
    assert!(!blastp.ungapped);
    assert_eq!(blastp.inclusion_ethresh, Some(0.002));
    let resolved = blastp.resolve().unwrap();
    assert_eq!(resolved.app, BlastpApp::Psiblast);
    assert_eq!(resolved.inclusion_ethresh, 0.002);
    assert_eq!(resolved.query_batch_size().unwrap(), 1);
    // traceback_stage.cpp:216-225: the preliminary hit list size of CBS (1050 for 500).
    assert_eq!(resolved.max_target_seqs, 500);
    assert_eq!(resolved.traceback_hitlist_size(), 1050);
    let resolved = parse(&["-max_target_seqs", "600", "-inclusion_ethresh", "1e-300"])
        .unwrap()
        .into_blastp_args()
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!(resolved.traceback_hitlist_size(), 1250);
    assert_eq!(resolved.inclusion_ethresh, 1e-300);
    // blastp keeps its hit list size and its threshold in the traceback.
    let mut blastp = parse(&[]).unwrap().into_blastp_args().unwrap();
    blastp.app = BlastpApp::Blastp;
    blastp.inclusion_ethresh = None;
    let resolved = blastp.resolve().unwrap();
    assert_eq!(resolved.traceback_hitlist_size(), 500);
    assert_eq!(resolved.inclusion_ethresh, 0.002);
}

// blast_args.cpp:234-238 and 303-307, ncbiargs.cpp:118-129: `-searchsp` takes NCBI's `Int8`
// forms (a sign, hex, leading zeros) at 0 or above, and at any value (0 too) the search runs
// without the Gumbel block (NCBI 2.17.0 probes `x.sp.{plus,hex,lead0}`, `err.sp.*` and
// `x.sp.{space,trail}` of line PS, item PS5c).
#[test]
fn searchsp_takes_ncbis_int8_forms_and_drops_the_gumbel_block() {
    for (value, space) in [
        ("+1000000000", 1_000_000_000),
        ("0x3B9ACA00", 1_000_000_000),
        ("0001000000000", 1_000_000_000),
        ("9223372036854775807", i64::MAX),
        ("0", 0),
    ] {
        let args = parse(&["-searchsp", value]).unwrap();
        assert_eq!(args.searchsp, Some(space), "{value}");
        let resolved = args.into_blastp_args().unwrap().resolve().unwrap();
        assert_eq!(resolved.searchsp, Some(space), "{value}");
        assert!(resolved.old_fsc, "{value}");
    }
    let blastp = parse(&[]).unwrap().into_blastp_args().unwrap();
    assert_eq!(blastp.searchsp, None);
    for value in [" 1", "1 ", "1e6", "1.5", "", "-1", "9223372036854775808"] {
        assert!(
            parse(&["-searchsp", value]).is_err(),
            "-searchsp {value:?} is refused"
        );
    }
}

// psiblast_app.cpp:246-321 (`rounds.rs`): the rounds of a FASTA query run (0: until
// convergence), and those of `-in_pssm` (PS5b, with one round too when a PSSM file is
// written), with `-query_loc` too (Q4a).
#[test]
fn rounds_run_and_the_later_output_formats_are_refused() {
    for (rounds, plan) in [("0", Some(0)), ("2", Some(2)), ("1", None)] {
        let blastp = parse(&["-num_iterations", rounds])
            .unwrap()
            .into_blastp_args()
            .unwrap();
        assert_eq!(blastp.psiblast_pssm.as_ref().unwrap().rounds, plan);
        let resolved = blastp.resolve().unwrap();
        assert_eq!(
            resolved.psi_rounds.as_ref().map(|plan| plan.num_iterations),
            plan
        );
    }
    let in_pssm = |extra: &[&str]| {
        let mut argv = vec![
            "losat", "psiblast", "-in_pssm", "x.asn", "-subject", "s.faa",
        ];
        argv.extend_from_slice(extra);
        let Ok(Cli {
            command: Commands::Psiblast(args),
        }) = crate::cli::try_parse_from(argv)
        else {
            panic!("psiblast -in_pssm {extra:?}");
        };
        args.into_blastp_args()
            .unwrap()
            .psiblast_pssm
            .unwrap()
            .rounds
    };
    assert_eq!(in_pssm(&["-num_iterations", "3"]), Some(3));
    assert_eq!(in_pssm(&[]), None);
    assert_eq!(in_pssm(&["-out_pssm", "f"]), Some(1));
    assert_eq!(in_pssm(&["-out_ascii_pssm", "f"]), Some(1));
    let ranged = parse(&["-num_iterations", "2", "-query_loc", "1-20"])
        .unwrap()
        .into_blastp_args()
        .unwrap()
        .resolve()
        .unwrap();
    assert_eq!(
        ranged.psi_rounds.as_ref().map(|plan| plan.num_iterations),
        Some(2),
        "the rounds of a range run"
    );
    assert!(ranged.query_range.is_some(), "the range reaches the rounds");
    for masks in [
        &["-seg", "yes"][..],
        &["-soft_masking", "true"],
        &["-seg", "yes", "-soft_masking", "true"],
    ] {
        let masked = parse(&[&["-num_iterations", "3"][..], masks].concat())
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!(
            masked.psi_rounds.as_ref().map(|plan| plan.num_iterations),
            Some(3),
            "the rounds of a masked query run ({masks:?})"
        );
    }
    // Items PS15b-1 and PS15b-2: BLAST XML, the Seq-annots, the archive, the query-anchored
    // reports, the JSON Seq-aligns and the organism report run; the XML2 and JSON documents (the
    // second part of PS15b-2) are refused.
    for format in ["1", "2", "3", "4", "5", "8", "9", "11", "12", "18"] {
        assert!(
            parse(&["-outfmt", format])
                .unwrap()
                .into_blastp_args()
                .is_ok(),
            "outfmt {format}"
        );
    }
    for format in ["15", "16"] {
        let error = parse(&["-outfmt", format])
            .unwrap()
            .into_blastp_args()
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            format!("output format {format} is not supported by LOSAT's PSIBLAST yet")
        );
    }
    // Outfmt 13 and 14 to a file are refused here; to standard output they go on to NCBI's
    // own error (`Please provide a file name for outfmt 13.`, blast_app_util.cpp:888-899),
    // which comes later (rows db.dlm.ps.o13, .o14).
    for format in ["13", "14"] {
        let error = parse(&["-outfmt", format, "-out", "r.out"])
            .unwrap()
            .into_blastp_args()
            .unwrap_err()
            .to_string();
        assert_eq!(
            error,
            format!("output format {format} is not supported by LOSAT's PSIBLAST yet")
        );
        assert!(parse(&["-outfmt", format])
            .unwrap()
            .into_blastp_args()
            .is_ok());
    }
    for format in [
        "0",
        "6",
        "7",
        "6 qseqid sseqid evalue",
        "6 delim=| std",
        "10",
        "20 qseqid sseqid",
        "17",
    ] {
        assert!(parse(&["-outfmt", format])
            .unwrap()
            .into_blastp_args()
            .is_ok());
    }
    // With `-phi_pattern` (items PS8a, PS8b): PHI-BLAST's reports run (0 to 7, 10, 11, 18, 20),
    // and formats 8, 9 and 12 go on to NCBI's engine error.
    let with_pattern = |format: &str| {
        parse(&["-phi_pattern", "p.pat", "-outfmt", format])
            .unwrap()
            .into_blastp_args()
            .map(|_| ())
            .map_err(|error| error.to_string())
    };
    for format in [
        "0",
        "6",
        "7",
        "10",
        "20 qseqid sseqid",
        "8",
        "9",
        "12",
        "5",
        "11",
        "1",
        "2",
        "3",
        "4",
        "18",
    ] {
        assert_eq!(with_pattern(format), Ok(()), "outfmt {format}");
    }
    for format in ["15", "16"] {
        assert_eq!(
            with_pattern(format),
            Err(format!(
                "output format {format} is not supported by LOSAT's PSIBLAST yet"
            )),
            "phi outfmt {format}"
        );
    }
    // NCBI's own error of the format string comes first (NCBI 2.17.0 `psiblast -outfmt 99`:
    // "Error: Formatting choice is out of range", exit 255).
    assert!(parse(&["-outfmt", "99"])
        .unwrap()
        .into_blastp_args()
        .unwrap_err()
        .to_string()
        .contains("Formatting choice is out of range"));
}

#[test]
fn blastps_database_filters_and_id_options_reach_the_search() {
    let db = |extra: &[&str]| {
        let mut argv = vec!["losat", "psiblast", "-query", "q.faa", "-db", "d"];
        argv.extend_from_slice(extra);
        match crate::cli::try_parse_from(argv) {
            Ok(Cli {
                command: Commands::Psiblast(args),
            }) => Ok(args.into_blastp_args().unwrap()),
            Ok(other) => panic!("not psiblast: {other:?}"),
            Err(error) => Err(error.to_string()),
        }
    };
    let blastp = db(&[
        "-taxids",
        "9606",
        "-no_taxid_expansion",
        "-show_gis",
        "-parse_deflines",
        "-lcase_masking",
        "-num_iterations",
        "3",
    ])
    .unwrap();
    assert_eq!(
        blastp.db_filter.taxids.as_deref(),
        Some(std::ffi::OsStr::new("9606"))
    );
    assert!(blastp.db_filter.no_taxid_expansion);
    assert!(blastp.show_gis && blastp.parse_deflines && blastp.lcase_masking);
    let blastp = db(&["-negative_seqidlist", "l.txt"]).unwrap();
    assert_eq!(
        blastp.db_filter.negative_seqidlist.as_deref(),
        Some(std::ffi::OsStr::new("l.txt"))
    );
    let blastp = db(&[]).unwrap();
    assert!(!blastp.show_gis && !blastp.parse_deflines && !blastp.lcase_masking);
    assert!(blastp.db_filter.taxids.is_none() && blastp.db_filter.seqidlist.is_none());
    // The filters exclude `-subject` and `-subject_loc`, as blastp's do.
    for extra in [
        &["-subject", "s.faa", "-taxids", "9606"][..],
        &["-subject_loc", "1-10", "-seqidlist", "l.txt"][..],
    ] {
        let mut argv = vec!["losat", "psiblast", "-query", "q.faa"];
        argv.extend_from_slice(extra);
        let error = crate::cli::try_parse_from(argv).unwrap_err();
        assert!(
            error.to_string().contains("Incompatible with argument:  `"),
            "{extra:?}"
        );
    }
}

#[test]
fn the_x_drops_gap_trigger_culling_and_sum_stats_reach_blastps_search() {
    let blastp = |extra: &[&str]| parse(extra).unwrap().into_blastp_args();
    let args = blastp(&[
        "-culling_limit",
        "2",
        "-xdrop_ungap",
        "3.5",
        "-xdrop_gap",
        "1e9",
        "-xdrop_gap_final",
        "-1",
        "-gap_trigger",
        "-5",
        "-sum_stats",
        "false",
    ])
    .unwrap();
    assert_eq!(args.culling_limit, Some(2));
    assert_eq!(
        (args.xdrop_ungap, args.xdrop_gap, args.xdrop_gap_final),
        (Some(3.5), Some(1e9), Some(-1.0))
    );
    assert_eq!(args.gap_trigger, Some(-5.0));
    // NCBI's default gap trigger of a protein search (BLAST_GAP_TRIGGER_PROT).
    let args = blastp(&[]).unwrap();
    assert_eq!(args.gap_trigger, Some(22.0));
    assert_eq!(args.culling_limit, None);
    assert!(!args.sum_stats);
    // Sum statistics link the HSPs, in every spelling of NCBI's boolean.
    for spelling in ["true", "T", "t", "1", "yes"] {
        assert!(
            blastp(&["-sum_stats", spelling]).unwrap().sum_stats,
            "{spelling}"
        );
    }
    assert!(!blastp(&["-sum_stats", "false"]).unwrap().sum_stats);
    // NCBI's constraint of -culling_limit (0 or more).
    assert!(parse(&["-culling_limit", "-1"]).is_err());
}

#[test]
fn psiblasts_arguments_that_losat_does_not_run_are_named() {
    // `-html` (line RG item T4c2) is accepted and forwarded to blastp's search.
    assert!(!parse(&[]).unwrap().into_blastp_args().unwrap().html);
    assert!(parse(&["-html"]).unwrap().into_blastp_args().unwrap().html);
    // `-in_pssm` (PS4a) is an argument: with `-query` it is NCBI's reader's exclusion
    // (blast_args.cpp:1240-1247).
    let error = parse(&["-in_pssm", "x"]).unwrap_err();
    assert!(error.contains("Incompatible with argument:  `"), "{error}");
    // `-phi_pattern` (PS8a) is an argument; `-in_msa` excludes it (blast_args.cpp:1173-1192).
    let args = parse(&["-phi_pattern", "p.pat"]).unwrap();
    assert_eq!(
        args.phi_pattern.as_deref(),
        Some(std::path::Path::new("p.pat"))
    );
    let error = parse(&["-phi_pattern", "p.pat", "-in_msa", "a.aln"]).unwrap_err();
    assert!(error.contains("Incompatible with argument:  `"), "{error}");
    // Arguments that psiblast does not have are unknown, as in NCBI's reader.
    for option in ["-task", "-ungapped", "-db_soft_mask", "-mt_mode"] {
        let error = parse(&[option, "1"]).unwrap_err();
        assert!(
            error.contains(&format!("Unknown argument: \"{}\"", &option[1..])),
            "{option}: {error}"
        );
    }
}
