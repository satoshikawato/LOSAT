//! Unit tests for blastn/args.rs

use LOSAT::algorithm::blastn::BlastnArgs;
use clap::{Command, FromArgMatches, Args};
use std::path::PathBuf;

fn parse_args(args: &[&str]) -> BlastnArgs {
    let mut all_args = vec!["losat".to_string(), "blastn".to_string()];
    all_args.extend(args.iter().map(|s| s.to_string()));
    
    // Create a command and add BlastnArgs as arguments
    // Use the same approach as main.rs: create a subcommand
    let cmd = Command::new("losat")
        .subcommand(BlastnArgs::augment_args(Command::new("blastn")));
    
    let matches = cmd.get_matches_from(all_args);
    let sub_matches = matches.subcommand_matches("blastn").unwrap();
    
    BlastnArgs::from_arg_matches(sub_matches).unwrap()
}

#[test]
fn test_default_values() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta"]);
    
    assert_eq!(args.task, "megablast");
    assert_eq!(args.word_size, 28);
    assert_eq!(args.num_threads, 0);
    assert_eq!(args.evalue, 10.0);
    assert_eq!(args.max_target_seqs, 500);
    assert_eq!(args.reward, 1);
    assert_eq!(args.penalty, -2);
    assert_eq!(args.gap_open, 0);
    assert_eq!(args.gap_extend, 0);
    assert_eq!(args.dust, true);
    assert_eq!(args.dust_level, 20);
    assert_eq!(args.dust_window, 64);
    assert_eq!(args.dust_linker, 1);
    assert_eq!(args.verbose, false);
    assert_eq!(args.chain, false);
    assert_eq!(args.scan_step, 0);
}

#[test]
fn test_custom_task() {
    let args = parse_args(&["-q", "query.fasta", "-s", "subject.fasta", "--task", "blastn"]);
    assert_eq!(args.task, "blastn");
}

#[test]
fn test_custom_word_size() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--word-size", "11"]);
    assert_eq!(args.word_size, 11);
}

#[test]
fn test_custom_num_threads() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-n", "4"]);
    assert_eq!(args.num_threads, 4);
}

#[test]
fn test_custom_evalue() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--evalue", "1e-5"]);
    assert_eq!(args.evalue, 1e-5);
}

#[test]
fn test_custom_max_target_seqs() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--max-target-seqs", "1000"]);
    assert_eq!(args.max_target_seqs, 1000);
}

#[test]
fn test_custom_scoring_parameters() {
    let args = parse_args(&[
        "blastn", "-q", "query.fasta", "-s", "subject.fasta",
        "--reward", "2",
        "--penalty", "-3",
        "--gap-open", "5",
        "--gap-extend", "2",
    ]);
    assert_eq!(args.reward, 2);
    assert_eq!(args.penalty, -3);
    assert_eq!(args.gap_open, 5);
    assert_eq!(args.gap_extend, 2);
}

#[test]
fn test_dust_options() {
    let args = parse_args(&[
        "blastn", "-q", "query.fasta", "-s", "subject.fasta",
        "--dust", "false",
        "--dust-level", "30",
        "--dust-window", "32",
        "--dust-linker", "2",
    ]);
    assert_eq!(args.dust, false);
    assert_eq!(args.dust_level, 30);
    assert_eq!(args.dust_window, 32);
    assert_eq!(args.dust_linker, 2);
}

#[test]
fn test_verbose_flag() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-v"]);
    assert_eq!(args.verbose, true);
}

#[test]
fn test_chain_flag() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--chain"]);
    assert_eq!(args.chain, true);
}

#[test]
fn test_scan_step() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--scan-step", "8"]);
    assert_eq!(args.scan_step, 8);
}

#[test]
fn test_output_path() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-o", "output.txt"]);
    assert_eq!(args.out, Some(PathBuf::from("output.txt")));
}

#[test]
fn test_query_and_subject_paths() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta"]);
    assert_eq!(args.query, PathBuf::from("query.fasta"));
    assert_eq!(args.subject, PathBuf::from("subject.fasta"));
}

