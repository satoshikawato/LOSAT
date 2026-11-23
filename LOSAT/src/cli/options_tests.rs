//! The options of the ten programs on NCBI's reader (`options`, each program's `from_args`) are
//! those that clap converted from the same arguments (`convert`), field by field (line RG item
//! K4i-3c, its proof commit): on every accepted case of the frozen syntax cases
//! (`tests/unit/ncbi_syntax_cases.tsv.gz`) and of NCBI's verdicts (`tests/unit/ncbi_args_verdicts.tsv`),
//! and on a sweep of every argument of every program that LOSAT does not refuse by name, given
//! (values of its kind, at the bounds of its constraint) and not given.

use std::collections::BTreeSet;
use std::ffi::OsString;
use std::io::Read;

use super::{convert, options, refusal_after_read, tool, TOOLS};
use crate::blastinput::ncbi_environment::{preparse, Preparsed};
use crate::ncbi_app::args::{described, parse, Args, Dep, Kind};
use crate::ncbi_app::Tool;

/// The arguments that NCBI's reader reads from `words` (after the toolkit's words that the
/// application start takes), and whether the words are a dry run; `None` where the reader refuses
/// them or prints a text.
fn read(tool: &Tool, words: &[OsString]) -> Option<(Args, bool)> {
    match preparse(words) {
        Preparsed::Words { words, dryrun, .. } => {
            parse(tool, &words).ok().map(|args| (args, dryrun))
        }
        Preparsed::Version(_) => None,
    }
}

/// How the two conversions of `words` compare: `None` where the reader or LOSAT's refusals end the
/// words before a conversion, else whether they give the same options (or the same error).
fn compare(tool: &Tool, words: &[OsString]) -> Option<Result<(), String>> {
    let (args, dry_run) = read(tool, words)?;
    if !dry_run && refusal_after_read(tool, &args).is_some() {
        return None;
    }
    let clap = convert(tool, &args, dry_run).map(|command| format!("{command:?}"));
    let ours = options(tool, &args, dry_run).map(|command| format!("{command:?}"));
    Some(match clap == ours {
        true => Ok(()),
        false => Err(format!(
            "{} {words:?}:\n  clap {clap:?}\n  ours {ours:?}",
            tool.program
        )),
    })
}

/// The counts of a set of cases: compared and equal, ended before a conversion, and the
/// differences.
#[derive(Default)]
struct Counts {
    same: usize,
    ended: usize,
    differing: Vec<String>,
}

impl Counts {
    fn add(&mut self, compared: Option<Result<(), String>>) {
        match compared {
            Some(Ok(())) => self.same += 1,
            Some(Err(difference)) => self.differing.push(difference),
            None => self.ended += 1,
        }
    }

    fn assert_same(&self, what: &str) {
        eprintln!(
            "{what}: {} the same, {} ended before a conversion, {} differing",
            self.same,
            self.ended,
            self.differing.len()
        );
        assert!(
            self.differing.is_empty(),
            "{what}: {} differ:\n{}",
            self.differing.len(),
            self.differing[..self.differing.len().min(20)].join("\n")
        );
    }
}

fn words(words: &[&str]) -> Vec<OsString> {
    words.iter().map(OsString::from).collect()
}

#[test]
fn the_options_are_clap_s_on_every_accepted_syntax_case() {
    let gz = include_bytes!("../../tests/unit/ncbi_syntax_cases.tsv.gz");
    let mut text = String::new();
    flate2::read::GzDecoder::new(&gz[..])
        .read_to_string(&mut text)
        .unwrap();
    let mut counts = Counts::default();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (program, kind) = (fields[1], fields[3]);
        if kind != "accepted" {
            continue;
        }
        let tool = tool(program).expect("a program on the reader");
        counts.add(compare(tool, &words(&fields[5..])));
    }
    assert!(counts.same > 0);
    counts.assert_same("accepted syntax cases");
}

#[test]
fn the_options_are_clap_s_on_every_accepted_verdict() {
    let table = include_str!("../../tests/unit/ncbi_args_verdicts.tsv");
    let mut counts = Counts::default();
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (program, code) = (fields[1], fields[2]);
        let Some(tool) = tool(program) else {
            continue;
        };
        if code != "ok" {
            continue;
        }
        let count: usize = fields[3].parse().unwrap();
        counts.add(compare(tool, &words(&fields[4 + count..])));
    }
    assert!(counts.same > 0);
    counts.assert_same("accepted verdicts");
}

/// Values of each kind for the sweep: NCBI's reader keeps those that its constraint allows.
fn values(kind: Kind, choices: &[String]) -> Vec<OsString> {
    let mut values: Vec<OsString> = choices.iter().map(OsString::from).collect();
    let more: &[&str] = match kind {
        Kind::Flag => &[],
        Kind::Boolean => &["true", "false", "t", "F", "1", "0", "yes", "No"],
        Kind::Integer => &[
            "-2147483648",
            "-1",
            "0",
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "11",
            "16",
            "18",
            "21",
            "28",
            "33",
            "255",
            "256",
            "2147483647",
            "0x1F",
            "+7",
        ],
        Kind::Int8 => &[
            "-9223372036854775808",
            "-1",
            "0",
            "1",
            "9223372036854775807",
            "0x10",
        ],
        Kind::Double => &[
            "-1", "0", "0.05", "0.5", "1", "1e-5", "10", "99.5", "100", "101", "1e308", "inf",
            "nan", "0x1p3", "2,5",
        ],
        Kind::Text => &[
            "x",
            "yes",
            "no",
            "0",
            "6",
            "1-10",
            "7 qseqid sseqid",
            "20 64 1",
            "12 2.2 2.5",
            "",
        ],
        Kind::InputFile | Kind::OutputFile => &["x", "-", "in/q.fa"],
    };
    values.extend(more.iter().map(OsString::from));
    #[cfg(unix)]
    if matches!(kind, Kind::Text | Kind::InputFile | Kind::OutputFile) {
        use std::os::unix::ffi::OsStringExt;
        values.push(OsString::from_vec(b"a\xffb".to_vec()));
    }
    values
}

/// The words that a program's command line needs before the swept argument: makeblastdb's
/// mandatory `-dbtype` (`-in` takes its default).
fn base(program: &str, swept: &str) -> Vec<OsString> {
    match program {
        "makeblastdb" if swept != "dbtype" => words(&["-dbtype", "nucl"]),
        _ => Vec::new(),
    }
}

/// `words` (the swept argument) with the arguments that it requires in `tool`'s program
/// (`Dep::Requires`), each with the first value with which NCBI's reader takes the words; a
/// requirement without one keeps its first value (the reader then ends the words).
fn with_required(tool: &Tool, name: &str, mut words: Vec<OsString>) -> Vec<OsString> {
    for (first, dep, second) in crate::ncbi_app::args::dependencies(tool.program) {
        if first != name || dep != Dep::Requires {
            continue;
        }
        let Some(described) = described(tool).find(|arg| arg.name == second) else {
            continue;
        };
        words.push(OsString::from(format!("-{second}")));
        if described.kind != Kind::Flag {
            let values = values(described.kind, &described.choices);
            let taken = values.iter().find(|value| {
                let mut tried = words.clone();
                tried.push((*value).clone());
                read(tool, &tried).is_some()
            });
            words.push(taken.unwrap_or(&values[0]).clone());
        }
    }
    words
}

#[test]
fn the_options_are_clap_s_on_a_sweep_of_every_argument() {
    let mut counts = Counts::default();
    for tool in TOOLS.iter().copied() {
        let program = tool.program;
        // Not given: every option at its default.
        counts.add(compare(tool, &base(program, "")));
        let mut dry = base(program, "");
        dry.push("-dryrun".into());
        counts.add(compare(tool, &dry));
        let mut arguments: Vec<(&str, Kind, Vec<String>)> = described(tool)
            .map(|arg| (arg.name, arg.kind, arg.choices))
            .collect();
        if program == "makeblastdb" {
            arguments.push(("logfile", Kind::OutputFile, Vec::new()));
        }
        for (name, kind, choices) in arguments {
            let given: Vec<Option<OsString>> = match kind {
                Kind::Flag => vec![None],
                _ => values(kind, &choices).into_iter().map(Some).collect(),
            };
            for value in given {
                for dry_run in [false, true] {
                    let mut words = base(program, name);
                    words.push(format!("-{name}").into());
                    words.extend(value.clone());
                    let mut words = with_required(tool, name, words);
                    if dry_run {
                        words.push("-dryrun".into());
                    }
                    counts.add(compare(tool, &words));
                }
            }
        }
    }
    assert!(counts.same > 1000);
    counts.assert_same("sweep");
}

/// Every argument of each program that LOSAT does not refuse by name is read into its options:
/// the names read by `from_args` are its XML's (and makeblastdb's `-logfile`).
#[test]
fn every_argument_is_read_into_the_options() {
    for tool in TOOLS.iter().copied() {
        let (args, _) = read(tool, &base(tool.program, "")).expect("the base words read");
        options(tool, &args, false).expect("the options of the base words");
        let names = args.names_read();
        let mut expected: BTreeSet<&str> = described(tool).map(|arg| arg.name).collect();
        if tool.program == "makeblastdb" {
            expected.insert("logfile");
        }
        assert_eq!(names, expected, "{}", tool.program);
    }
}
