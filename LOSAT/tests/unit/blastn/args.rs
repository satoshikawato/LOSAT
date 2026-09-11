//! Unit tests for blastn/args.rs

use clap::{Args, Command, FromArgMatches};
use std::path::PathBuf;
use LOSAT::algorithm::blastn::BlastnArgs;

fn parse_args(args: &[&str]) -> BlastnArgs {
    let mut all_args = vec!["losat".to_string(), "blastn".to_string()];
    all_args.extend(args.iter().map(|s| s.to_string()));

    // NCBI blast_args.cpp:2657-2660: select the implemented tabular formatter explicitly.
    all_args.extend(["-outfmt".into(), "6".into()]);
    let cli: LOSAT::cli::Cli = LOSAT::cli::try_parse_from(all_args).unwrap();
    let LOSAT::cli::Commands::Blastn(args) = cli.command else {
        panic!("wrong program")
    };
    args
}

#[test]
fn test_default_values() {
    let args = parse_args(&["-query", "query.fasta", "-subject", "subject.fasta"]);

    assert_eq!(args.task, "megablast");
    assert_eq!(args.word_size, 28);
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.evalue, 10.0);
    assert_eq!(args.max_target_seqs, Some(500));
    assert_eq!(args.reward, 1);
    assert_eq!(args.penalty, -2);
    assert_eq!(args.gap_open, 0);
    assert_eq!(args.gap_extend, 0);
    assert_eq!(args.dust.params(), Some((20, 64, 1)));
    assert_eq!(args.verbose, false);
    assert_eq!(args.scan_step, 0);
}

#[test]
fn test_custom_task() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-task",
        "blastn",
    ]);
    assert_eq!(args.task, "blastn");
}

#[test]
fn test_custom_word_size() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-word_size",
        "11",
    ]);
    assert_eq!(args.word_size, 11);
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
fn test_custom_scoring_parameters() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-reward",
        "2",
        "-penalty=-3",
        "-gapopen",
        "5",
        "-gapextend",
        "2",
    ]);
    assert_eq!(args.reward, 2);
    assert_eq!(args.penalty, -3);
    assert_eq!(args.gap_open, 5);
    assert_eq!(args.gap_extend, 2);
}

#[test]
fn test_dust_options() {
    // Note: -dust is a bool flag, so we can't set it to false directly
    // We'll test the other dust options instead
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-dust",
        "30 32 2",
    ]);
    // dust defaults to true
    assert_eq!(args.dust.params(), Some((30, 32, 2)));
}

#[test]
fn test_verbose_flag() {
    let args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-verbose",
    ]);
    assert_eq!(args.verbose, true);
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
