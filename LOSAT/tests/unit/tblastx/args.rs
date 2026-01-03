//! Unit tests for tblastx/args.rs

use LOSAT::algorithm::tblastx::TblastxArgs;
use clap::{Command, FromArgMatches, Args};
use std::path::PathBuf;

fn parse_args(args: &[&str]) -> TblastxArgs {
    let mut all_args = vec!["losat".to_string(), "tblastx".to_string()];
    all_args.extend(args.iter().map(|s| s.to_string()));
    
    // Create a command and add TblastxArgs as arguments
    // Use the same approach as main.rs: create a subcommand
    let cmd = Command::new("losat")
        .subcommand(TblastxArgs::augment_args(Command::new("tblastx")));
    
    let matches = cmd.get_matches_from(all_args);
    let sub_matches = matches.subcommand_matches("tblastx").unwrap();
    
    TblastxArgs::from_arg_matches(sub_matches).unwrap()
}

#[test]
fn test_default_values() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta"]);
    
    assert_eq!(args.evalue, 10.0);
    assert_eq!(args.threshold, 13);
    assert_eq!(args.word_size, 3);
    assert_eq!(args.num_threads, 0);
    assert_eq!(args.query_gencode, 1);
    assert_eq!(args.db_gencode, 1);
    assert_eq!(args.max_target_seqs, 500);
    assert_eq!(args.seg, true);
    assert_eq!(args.seg_window, 12);
    assert_eq!(args.seg_locut, 2.2);
    assert_eq!(args.seg_hicut, 2.5);
}

#[test]
fn test_custom_evalue() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-e", "1e-5"]);
    assert_eq!(args.evalue, 1e-5);
}

#[test]
fn test_custom_threshold() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-t", "20"]);
    assert_eq!(args.threshold, 20);
}

#[test]
fn test_custom_word_size() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-w", "2"]);
    assert_eq!(args.word_size, 2);
}

#[test]
fn test_custom_num_threads() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "-n", "4"]);
    assert_eq!(args.num_threads, 4);
}

#[test]
fn test_custom_gencode() {
    let args = parse_args(&[
        "-q", "query.fasta", "-s", "subject.fasta",
        "--query-gencode", "11",
        "--db-gencode", "11",
    ]);
    assert_eq!(args.query_gencode, 11);
    assert_eq!(args.db_gencode, 11);
}

#[test]
fn test_custom_max_target_seqs() {
    let args = parse_args(&[ "-q", "query.fasta", "-s", "subject.fasta", "--max-target-seqs", "1000"]);
    assert_eq!(args.max_target_seqs, 1000);
}

#[test]
fn test_seg_options() {
    let args = parse_args(&[
        "-q", "query.fasta", "-s", "subject.fasta",
        "--seg=false",
        "--seg-window", "15",
        "--seg-locut", "1.9",
        "--seg-hicut", "2.3",
    ]);
    assert_eq!(args.seg, false);
    assert_eq!(args.seg_window, 15);
    assert_eq!(args.seg_locut, 1.9);
    assert_eq!(args.seg_hicut, 2.3);
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

