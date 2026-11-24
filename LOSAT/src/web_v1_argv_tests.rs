//! v1's frozen command (`frozen_command`, a builder) is the command that clap's derive of the
//! command line made when v1 was frozen, and v1 reads every request as it read it with the derive
//! (line RG item K4i-3c, its proof before the derive leaves; DECISIONS K4i-3-1 Q4): the same
//! arguments (ids, long names, value names, actions, defaults, value types) and usage, the same
//! list of long names for its messages (`LONGS`), and, on the frozen cases and a sweep of every
//! argument, the same options or the same message.

use std::ffi::OsString;

use clap::{CommandFactory, FromArgMatches};

use super::*;
use crate::cli::Cli;

/// What the arguments of a program's command are, without their help texts.
fn shape(command: &clap::Command, program: &str) -> (Vec<String>, String) {
    let mut command = command.clone();
    command.build();
    let scope = command.find_subcommand(program).unwrap();
    let args = scope
        .get_arguments()
        .map(|arg| {
            format!(
                "{} {:?} {:?} {:?} {:?} {:?} {:?} {} {} {:?}",
                arg.get_id(),
                arg.get_long(),
                arg.get_short(),
                arg.get_value_names(),
                arg.get_action(),
                arg.get_num_args(),
                arg.get_default_values(),
                arg.is_required_set(),
                arg.is_hide_set(),
                arg.get_value_parser().type_id(),
            )
        })
        .collect();
    let mut scope = scope.clone();
    (args, scope.render_usage().to_string())
}

#[test]
fn the_frozen_command_is_the_derive_s() {
    for program in ["blastx", "tblastn"] {
        assert_eq!(
            shape(&frozen_command(program), program),
            shape(&crate::cli::command_of(program), program),
            "{program}"
        );
    }
}

#[test]
fn the_long_names_are_the_derive_s() {
    let command = Cli::command();
    let derive: Vec<&str> = std::iter::once(&command)
        .chain(command.get_subcommands())
        .flat_map(|scope| scope.get_arguments().filter_map(clap::Arg::get_long))
        .collect();
    assert_eq!(LONGS, derive.as_slice());
}

/// v1's reading of `argv` with the derive (before line RG item K4i-3c): the derive's command with
/// v1's parsers (and `mutate`), the derive's conversion of the words it took, and its messages.
fn with_derive(
    program: &str,
    mutate: fn(clap::Command) -> clap::Command,
    argv: Vec<OsString>,
) -> Result<String, String> {
    let command = mutate(v1_command(crate::cli::command_of(program), program));
    let read = || -> Result<Commands, clap::Error> {
        let (program, given) = checked_words(&command, argv)?;
        let mut translated: Vec<OsString> = vec!["losat".into(), program.into()];
        for (name, value) in given {
            match value {
                Some(value) => {
                    let mut word = OsString::from(format!("--{name}="));
                    word.push(value);
                    translated.push(word);
                }
                None => translated.push(format!("--{name}").into()),
            }
        }
        let mut matches = command.clone().try_get_matches_from(&translated)?;
        Cli::from_arg_matches_mut(&mut matches)
            .map(|cli| cli.command)
            .map_err(|error| error.format(&mut command.clone()))
    };
    read()
        .map(|command| format!("{command:?}"))
        .map_err(|error| render(&error))
}

/// v1's reading of `argv` now (`parse`).
fn with_builder(
    program: &str,
    mutate: fn(clap::Command) -> clap::Command,
    argv: Vec<OsString>,
) -> Result<String, String> {
    parse(mutate(command(program)), argv).map(|command| format!("{command:?}"))
}

fn blastx_outfmt(command: clap::Command) -> clap::Command {
    command.mut_subcommand("blastx", |blastx| {
        blastx.mut_arg("outfmt", |outfmt| {
            outfmt.value_parser(crate::algorithm::blastx::web::v1_outfmt)
        })
    })
}

fn unchanged(command: clap::Command) -> clap::Command {
    command
}

/// Values for the sweep of an argument that reads its value as `reads`.
fn values(reads: Reads, long: &str) -> Vec<&'static str> {
    match reads {
        Flag => vec![],
        Text => vec![
            "x",
            "",
            "0",
            "6",
            "7",
            "7 qseqid sseqid",
            "11",
            "both",
            "plus",
            "minus",
            "blastx",
            "blastx-fast",
            "tblastn",
            "tblastn-fast",
            "BLOSUM62",
            "1-10",
            "yes",
            "no",
            "20 64 1",
            "0x2",
        ],
        Path => vec!["x", "-", "in/q.fa", "--evalue", "---word_size"],
        Integer | Count | Small => vec![
            "-1",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "11",
            "12",
            "33",
            "255",
            "256",
            "x",
            "0x10",
            "0x",
            "2147483648",
            "+3",
            " 3",
            "1e1",
        ],
        Int8 => vec!["-1", "0", "1", "x", "0x10", "9223372036854775808", "1.5"],
        Double => vec![
            "-1", "0", "0.25", "0.5", "1", "50", "100", "101", "x", "1e-5", "inf", "nan", "0x1p2",
            "2,5",
        ],
        Boolean => vec!["true", "false", "T", "f", "1", "0", "x", ""],
    }
    .into_iter()
    .chain((long == "query_loc").then_some("---evalue"))
    .collect()
}

fn argv(program: &str, extra: &[&str]) -> Vec<OsString> {
    ["losat", program, "-query", "query", "-subject", "subject"]
        .iter()
        .chain(extra)
        .map(OsString::from)
        .collect()
}

#[test]
fn v1_reads_every_request_as_with_the_derive() {
    let mut cases: Vec<(&str, Vec<String>)> = Vec::new();
    for (program, args) in [("blastx", BLASTX), ("tblastn", TBLASTN)] {
        cases.push((program, vec![]));
        for &(_, long, _, reads, _) in args {
            if matches!(long, "query" | "subject") {
                continue;
            }
            let flag = format!("-{long}");
            match reads {
                Flag => {
                    cases.push((program, vec![flag.clone()]));
                    cases.push((program, vec![format!("{flag}=true")]));
                    cases.push((program, vec![flag.clone(), flag.clone()]));
                }
                _ => {
                    for value in values(reads, long) {
                        cases.push((program, vec![flag.clone(), value.to_string()]));
                        cases.push((program, vec![format!("{flag}={value}")]));
                    }
                    cases.push((program, vec![flag.clone()]));
                    cases.push((
                        program,
                        vec![flag.clone(), "1".into(), flag.clone(), "2".into()],
                    ));
                }
            }
        }
        for (name, excluded) in V1_EXCLUSIONS {
            for other in *excluded {
                let words = |name: &str| -> Vec<String> {
                    let reads = args.iter().find(|arg| arg.1 == name).map(|arg| arg.3);
                    match reads {
                        Some(Flag) => vec![format!("-{name}")],
                        _ => vec![format!("-{name}"), "x".into()],
                    }
                };
                let mut both = words(name);
                both.extend(words(other));
                cases.push((program, both));
            }
        }
        for extra in [
            &["-mt_mode", "1"][..],
            &["-mt_mode", "1", "-num_threads", "2"],
            &["-in_pssm", "x"],
            &["-query_loc", "1-5"],
            &["-nosuch"],
            &["--evalue", "1"],
            &["-evalue"],
            &["--"],
            &["-evalue", "1", "--"],
            &["-logfile", "x"],
            &["-dryrun"],
            &["-version"],
            &["-db", "x", "-dbsize", "5"],
            &["-outfmt", "7 qseqid sseqid evalue", "-evalue", "1e-5"],
            &["-outfmt", "6 nosuch"],
            &["-word_size", "1"],
            &["-matrix", "FOO"],
        ] {
            cases.push((program, extra.iter().map(|word| word.to_string()).collect()));
        }
    }
    let mut same = 0;
    let mut differing = Vec::new();
    for (program, extra) in &cases {
        let extra: Vec<&str> = extra.iter().map(String::as_str).collect();
        let mutate = match *program {
            "blastx" => blastx_outfmt,
            _ => unchanged,
        };
        let derive = with_derive(program, mutate, argv(program, &extra));
        let builder = with_builder(program, mutate, argv(program, &extra));
        match derive == builder {
            true => same += 1,
            false => differing.push(format!(
                "{program} {extra:?}:\n  derive  {derive:?}\n  builder {builder:?}"
            )),
        }
    }
    eprintln!("v1 cases: {same} the same, {} differing", differing.len());
    assert!(
        differing.is_empty(),
        "{}",
        differing[..differing.len().min(20)].join("\n")
    );
}
