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
    assert_eq!(args.threshold, 13);
    assert_eq!(args.word_size, 3);
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.query_gencode, 1);
    assert_eq!(args.db_gencode, 1);
    assert_eq!(args.max_target_seqs, 500);
    let params = args.seg.params().unwrap();
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
    assert_eq!(args.threshold, 20);
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
    assert_eq!(args.max_target_seqs, 1000);
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
    let params = args.seg.params().unwrap();
    assert_eq!((params.window, params.locut, params.hicut), (15, 1.9, 2.3));
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
    assert_eq!(args.subject, PathBuf::from("subject.fasta"));
}
