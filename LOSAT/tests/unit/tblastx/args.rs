//! Unit tests for tblastx/args.rs

use clap::{Args, Command, FromArgMatches};
use std::path::PathBuf;
use LOSAT::algorithm::tblastx::TblastxArgs;

fn parse_args(args: &[&str]) -> TblastxArgs {
    let mut all_args = vec!["losat".to_string(), "tblastx".to_string()];
    all_args.extend(args.iter().map(|s| s.to_string()));

    // NCBI blast_args.cpp:2657-2660: select the implemented tabular formatter explicitly.
    all_args.extend(["-outfmt".into(), "6".into()]);
    let cli: LOSAT::cli::Cli = LOSAT::cli::try_parse_from(all_args).unwrap();
    let LOSAT::cli::Commands::Tblastx(args) = cli.command else {
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
    assert_eq!(args.word_size, 3);
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

// NCBI c++/src/algo/blast/core/blast_options.c:1174-1235 (BLAST_GetSuggestedThreshold and
// BLAST_GetSuggestedWindowSize) and c++/src/algo/blast/blastinput/blast_args.cpp:586-623
// (TBLASTX takes the suggestion, plus 2, unless -threshold is given): PAM30 reads 18 and 15;
// the name is kept as typed.
#[test]
fn test_matrix_suggests_threshold_and_window() {
    use LOSAT::algorithm::tblastx::args::TblastxMatrix;
    use LOSAT::config::ScoringMatrix;
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
    assert_eq!(args.word_size, 3);
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
// ```c++
// CFilteringArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
// {
//     if (args[kArgLookupTableMaskingOnly]) {
//         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
//     }
// ```
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
