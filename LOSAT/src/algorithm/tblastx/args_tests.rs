//! Unit tests for tblastx/args.rs

use crate::algorithm::tblastx::TblastxArgs;
use std::path::PathBuf;

fn parse_args(args: &[&str]) -> TblastxArgs {
    let mut all_args = vec!["losat".to_string(), "tblastx".to_string()];
    all_args.extend(args.iter().map(|s| s.to_string()));

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660: select the implemented tabular formatter explicitly.
    all_args.extend(["-outfmt".into(), "6".into()]);
    let cli: crate::cli::Cli = crate::cli::try_parse_from(all_args).unwrap();
    let crate::cli::Commands::Tblastx(args) = cli.command else {
        panic!("wrong program")
    };
    args
}

#[test]
fn test_default_values() {
    let args = parse_args(&["-query", "query.fasta", "-subject", "subject.fasta"]);

    assert_eq!(args.evalue, 10.0);
    assert_eq!(args.threshold, None);
    assert_eq!(args.word_threshold(), 13.0);
    assert_eq!(args.matrix, None);
    assert_eq!(args.matrix_name(), "BLOSUM62");
    assert_eq!(args.window(), 40);
    assert_eq!(args.word_size, None);
    assert_eq!(args.word_size(), 3);
    assert_eq!(args.culling_limit, None);
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.query_gencode, 1);
    assert_eq!(args.db_gencode, 1);
    // An omitted -max_target_seqs keeps NCBI's default hit list size of 500 (and the
    // pairwise report's 250 alignments).
    assert_eq!(args.max_target_seqs, None);
    let params = args.seg_spec().unwrap().params().unwrap();
    assert_eq!((params.window, params.locut, params.hicut), (12, 2.2, 2.5));
}

#[test]
fn test_custom_evalue() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-evalue",
        "1e-5",
    ]);
    assert_eq!(args.evalue, 1e-5);
}

#[test]
fn test_custom_threshold() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-threshold",
        "20",
    ]);
    assert_eq!(args.threshold, Some(20.0));
    assert_eq!(args.word_threshold(), 20.0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1174-1235 (BLAST_GetSuggestedThreshold and
// BLAST_GetSuggestedWindowSize) and c++/src/algo/blast/blastinput/blast_args.cpp:586-623
// (TBLASTX takes the suggestion, plus 2, unless -threshold is given): PAM30 reads 18 and 15;
// the name is kept as typed.
#[test]
fn test_matrix_suggests_threshold_and_window() {
    use crate::algorithm::tblastx::args::TblastxMatrix;
    use crate::config::ScoringMatrix;
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-matrix",
        "pam30",
    ]);
    assert_eq!(args.matrix.as_deref(), Some("pam30"));
    assert_eq!(args.matrix_name(), "pam30");
    assert!(matches!(
        args.scoring_matrix(),
        TblastxMatrix::Table(ScoringMatrix::Pam30)
    ));
    assert_eq!((args.word_threshold(), args.window()), (18.0, 15));
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-matrix",
        "PAM30",
        "-threshold",
        "13",
        "-window_size",
        "40",
    ]);
    assert_eq!((args.word_threshold(), args.window()), (13.0, 40));
}

#[test]
fn test_custom_word_size() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-word_size",
        "3",
    ]);
    assert_eq!(args.word_size, Some(3));
}

#[test]
fn test_custom_num_threads() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-num_threads",
        "4",
    ]);
    assert_eq!(args.num_threads, 4);
}

#[test]
fn test_custom_gencode() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-query_gencode",
        "11",
        "-db_gencode",
        "11",
    ]);
    assert_eq!(args.query_gencode, 11);
    assert_eq!(args.db_gencode, 11);
}

#[test]
fn test_custom_max_target_seqs() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-max_target_seqs",
        "1000",
    ]);
    assert_eq!(args.max_target_seqs, Some(1000));
}

#[test]
fn test_seg_options() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-seg",
        "15 1.9 2.3",
    ]);
    let params = args.seg_spec().unwrap().params().unwrap();
    assert_eq!((params.window, params.locut, params.hicut), (15, 1.9, 2.3));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
// `-soft_masking` takes NCBI's boolean spellings (`NStr::StringToBool`); its default is false
// (`kDfltArgLookupTableMaskingOnlyProt`).
#[test]
fn test_soft_masking_sets_mask_at_hash() {
    let soft = |extra: &[&str]| {
        let mut argv = vec!["-query", "query.fasta", "-subject", "subject.fasta"];
        argv.extend_from_slice(extra);
        parse_args(&argv).soft_masking()
    };
    assert!(!soft(&[]));
    for value in ["true", "T", "1", "YES", "on"] {
        assert!(soft(&["-soft_masking", value]), "{value}");
    }
    for value in ["false", "F", "0", "no", "off"] {
        assert!(!soft(&["-soft_masking", value]), "{value}");
    }
    assert!(soft(&["-seg", "no", "-soft_masking", "true"]));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-244 and 321-323
// (cited at `TblastxArgs::sum_stats`); c++/src/algo/blast/api/tblastx_options.cpp:76-81 (the
// default, true)
// `-sum_stats` takes NCBI's boolean spellings (`NStr::StringToBool`); omitted, the search links
// with sum statistics and the archive records nothing (`None`).
#[test]
fn test_sum_stats() {
    let parse = |extra: &[&str]| {
        let mut argv = vec![
            "losat", "tblastx", "-query", "q.fa", "-subject", "s.fa", "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        crate::cli::try_parse_from(argv).map(|cli: crate::cli::Cli| match cli.command {
            crate::cli::Commands::Tblastx(args) => args,
            _ => panic!("wrong program"),
        })
    };
    let omitted = parse(&[]).unwrap();
    assert_eq!(omitted.sum_stats, None);
    assert!(omitted.sum_stats());
    for value in ["true", "TRUE", "t", "T", "1", "yes", "on"] {
        let args = parse(&["-sum_stats", value]).unwrap();
        assert_eq!(args.sum_stats, Some(true), "{value}");
        assert!(args.sum_stats(), "{value}");
    }
    for value in ["false", "F", "f", "0", "no", "off"] {
        let args = parse(&["-sum_stats", value]).unwrap();
        assert_eq!(args.sum_stats, Some(false), "{value}");
        assert!(!args.sum_stats(), "{value}");
    }
    for value in ["maybe", "2", ""] {
        assert!(parse(&["-sum_stats", value]).is_err(), "{value}");
    }
    assert!(parse(&["-sum_stats"]).is_err());
}

#[test]
fn test_output_path() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-out",
        "output.txt",
    ]);
    assert_eq!(args.out, Some(PathBuf::from("output.txt")));
}

#[test]
fn test_query_and_subject_paths() {
    let args = parse_args(&["-query", "query.fasta", "-subject", "subject.fasta"]);
    assert_eq!(args.query, PathBuf::from("query.fasta"));
    assert_eq!(args.subject, Some(PathBuf::from("subject.fasta")));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-944
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1525-1534 (cited at
// `check_ncbi_options`)
// `-max_intron_length` is an integer of 0 or more (NCBI's argument errors for the others are the
// parser's); 0 runs as without it, any other value is NCBI's options error for TBLASTX.
#[test]
fn test_max_intron_length() {
    let parse = |extra: &[&str]| {
        let mut argv = vec![
            "losat", "tblastx", "-query", "q.fa", "-subject", "s.fa", "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        crate::cli::try_parse_from(argv).map(|cli: crate::cli::Cli| match cli.command {
            crate::cli::Commands::Tblastx(args) => args,
            _ => panic!("wrong program"),
        })
    };
    let check = |args: &TblastxArgs| crate::algorithm::tblastx::blast_engine::check_options(args);
    let omitted = parse(&[]).unwrap();
    assert_eq!(omitted.max_intron_length, None);
    assert!(check(&omitted).is_ok());
    let zero = parse(&["-max_intron_length", "0"]).unwrap();
    assert_eq!(zero.max_intron_length, Some(0));
    assert!(check(&zero).is_ok());
    let uneven = crate::blastinput::app::options_error(
        "Uneven gap linking of HSPs is allowed for blastx, tblastn, and psitblastn only",
    )
    .to_string();
    for value in ["1", "4", "100", "2147483647"] {
        let args = parse(&["-max_intron_length", value]).unwrap();
        assert_eq!(check(&args).unwrap_err().to_string(), uneven, "{value}");
    }
    // The e-value's check comes first (blast_options.c:1518-1523).
    let both = parse(&["-max_intron_length", "100", "-evalue", "0"]).unwrap();
    assert!(check(&both)
        .unwrap_err()
        .to_string()
        .contains("expect value or cutoff score must be greater than zero"));
    for value in ["-1", "1.5", "abc", "2147483648", ""] {
        assert!(parse(&["-max_intron_length", value]).is_err(), "{value}");
    }
}
