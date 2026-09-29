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
    // NCBI blast_args.cpp:166-170,648-660: the word size and scores are optional keys;
    // an omitted one keeps the task's default (coordination.rs).
    assert_eq!(args.word_size, None);
    assert_eq!(args.num_threads, 1);
    assert_eq!(args.evalue, 10.0);
    // NCBI blast_args.cpp:2913-2927: an omitted -max_target_seqs keeps the defaults of
    // the hit list (500) and of the pairwise alignments (250), so it stays unset here.
    assert_eq!(args.max_target_seqs, None);
    assert_eq!(args.reward, None);
    assert_eq!(args.penalty, None);
    assert_eq!(args.gap_open, None);
    assert_eq!(args.gap_extend, None);
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
    assert_eq!(args.word_size, Some(11));
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
    assert_eq!(args.reward, Some(2));
    assert_eq!(args.penalty, Some(-3));
    assert_eq!(args.gap_open, Some(5));
    assert_eq!(args.gap_extend, Some(2));
}

#[test]
fn given_scores_replace_the_task_defaults_even_when_equal_to_megablast_ones() {
    use LOSAT::algorithm::blastn::coordination::{
        determine_effective_word_size, determine_scoring_params,
    };
    let base = [
        "-query", "q.fasta", "-subject", "s.fasta", "-task", "blastn",
    ];
    let args = parse_args(&base);
    assert_eq!(determine_scoring_params(&args), (2, -3, 5, 2));
    assert_eq!(determine_effective_word_size(&args), 11);
    let mut given = base.to_vec();
    given.extend([
        "-reward",
        "1",
        "-penalty",
        "-2",
        "-gapopen",
        "0",
        "-gapextend",
        "0",
        "-word_size",
        "28",
    ]);
    let args = parse_args(&given);
    assert_eq!(determine_scoring_params(&args), (1, -2, 0, 0));
    assert_eq!(determine_effective_word_size(&args), 28);
    let args = parse_args(&["-query", "q.fasta", "-subject", "s.fasta", "-gapopen=-1"]);
    assert_eq!(determine_scoring_params(&args), (1, -2, -1, 0));
}

#[test]
fn test_dust_options() {
    // Note: -dust is a bool flag, so we can't set it to false directly
    // We'll test the other dust options instead
    let mut args = parse_args(&[
        "-query",
        "query.fasta",
        "-subject",
        "subject.fasta",
        "-dust",
        "30 32 2",
    ]);
    // NCBI reads -dust in its filtering handler (`resolve_dust`); dust defaults to true.
    args.resolve_dust().unwrap();
    assert_eq!(args.dust.params(), Some((30, 32, 2)));
}

#[test]
fn options_without_an_ncbi_blastn_equivalent_are_rejected() {
    // AGENTS.md rule 5: NCBI blastn has none of these (-limit_lookup and
    // -max_db_word_count are magicblast's; NCBI has no -verbose or -min_hit_length).
    for extra in [
        &["-verbose"][..],
        &["-limit_lookup"],
        &["-max_db_word_count", "30"],
        &["-min_hit_length", "10"],
    ] {
        let mut argv = vec!["losat", "blastn", "-query", "q", "-subject", "s"];
        argv.extend(extra);
        assert!(
            LOSAT::cli::try_parse_from::<LOSAT::cli::Cli, _, _>(argv).is_err(),
            "{extra:?}"
        );
    }
}

#[test]
fn query_defaults_to_standard_input() {
    // NCBI cmdline_flags.cpp:47: const string kDfltArgQuery("-");
    let args = parse_args(&["-subject", "subject.fasta"]);
    assert_eq!(args.query, PathBuf::from("-"));
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
