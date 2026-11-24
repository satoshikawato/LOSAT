//! Every argument of the ten programs on NCBI's reader that LOSAT does not refuse by name is read
//! into the program's options (`options`, each program's `from_args`; line RG item K4i-3c), and
//! the options read no other argument. The programs are one list. The options of each program
//! are frozen (`tests/fixtures/options/golden.txt`): at their defaults, and with each argument
//! given, so that a default that moves or two names that change places fail here.

use std::collections::BTreeSet;
use std::ffi::OsString;

use super::{not_converted, options, refusal_after_read, tool, TOOLS};
use crate::blastinput::ncbi_environment::{preparse, DiagTarget, Preparsed};
use crate::cli::Commands;
use crate::ncbi_app::args::{checked, described, parse, Args, Dep, Kind};
use crate::ncbi_app::Tool;

/// The programs whose options this reader reads (`TOOLS`), for the tests of other modules that go
/// over every one of them, so that a program added to the reader is in those tests too.
pub(crate) fn programs() -> impl Iterator<Item = &'static str> {
    TOOLS.iter().map(|tool| tool.program)
}

/// The names of the arguments that `tool`'s options read (`options`, `Args::names_read`), of a
/// command line with only what the program needs (makeblastdb's mandatory `-dbtype`; `-in` takes
/// its default): the tests that compare the arguments of the options with those of NCBI's XML,
/// LOSAT's descriptions and the search strategy's walk read them here.
pub(crate) fn names_read(tool: &Tool) -> BTreeSet<&'static str> {
    let args = parse(tool, &base(tool.program, "")).expect("the base words read");
    options(tool, &args, false).expect("the options of the base words");
    args.names_read()
}

#[test]
fn every_argument_is_read_into_the_options() {
    for tool in TOOLS.iter().copied() {
        let mut expected: BTreeSet<&str> = described(tool).map(|arg| arg.name).collect();
        if tool.program == "makeblastdb" {
            expected.insert("logfile");
        }
        assert_eq!(names_read(tool), expected, "{}", tool.program);
    }
}

/// The program of `command`'s variant: a new variant of `Commands` names its program here.
fn program_of(command: &Commands) -> &'static str {
    match command {
        Commands::Blastx(_) => "blastx",
        Commands::Blastn(_) => "blastn",
        Commands::Blastp(_) => "blastp",
        Commands::Tblastx(_) => "tblastx",
        Commands::Tblastn(_) => "tblastn",
        Commands::Psiblast(_) => "psiblast",
        Commands::Rpsblast(_) => "rpsblast",
        Commands::Rpstblastn(_) => "rpstblastn",
        Commands::Makedb(_) => "makeblastdb",
        Commands::Rpsbproc(_) => "rpsbproc",
    }
}

// The programs whose options are LOSAT's types are one list, the reader's table (`TOOLS`): the
// top level's help describes them in its order (`cli::PROGRAMS`), and each one reads its words into
// its own variant of `Commands` (`options`). A program added to one of them and not to the others
// fails here (line RG item K4i-3c, review L2).
#[test]
fn one_list_of_programs() {
    let described: Vec<&str> = crate::cli::PROGRAMS.iter().map(|(name, _)| *name).collect();
    let tools: Vec<&str> = TOOLS.iter().map(|tool| tool.program).collect();
    assert_eq!(described, tools);
    for tool in TOOLS.iter().copied() {
        assert!(crate::cli::is_program(tool.program), "{}", tool.program);
        let args = parse(tool, &base(tool.program, "")).expect("the base words read");
        let command = options(tool, &args, false).expect("the options of the base words");
        assert_eq!(program_of(&command), tool.program);
    }
}

// Line RG DECISIONS K4i-3-1 Q8: tblastn's and psiblast's `-query` is `-` (its default in the
// program's table) where a given `-in_pssm` or `-in_msa` excludes it, as before the options were
// read from NCBI's arguments.
#[test]
fn an_excluded_query_is_its_table_default() {
    for words in [
        ["losat", "tblastn", "-in_pssm", "x"],
        ["losat", "psiblast", "-in_pssm", "x"],
        ["losat", "psiblast", "-in_msa", "x"],
    ] {
        let query = match crate::cli::try_parse_from(words).unwrap().command {
            Commands::Tblastn(args) => args.query,
            Commands::Psiblast(args) => args.query,
            command => panic!("{words:?}: {command:?}"),
        };
        assert_eq!(query, std::path::Path::new("-"), "{words:?}");
    }
}

/// A word that is not UTF-8: `a`, a byte (Unix) or an unpaired surrogate (Windows), `b`.
fn not_utf8() -> OsString {
    #[cfg(unix)]
    let word = {
        use std::os::unix::ffi::OsStringExt;
        OsString::from_vec(b"a\xffb".to_vec())
    };
    #[cfg(windows)]
    let word = {
        use std::os::windows::ffi::OsStringExt;
        OsString::from_wide(&[0x61, 0xD800, 0x62])
    };
    assert!(word.to_str().is_none());
    word
}

// A dry run of the library (`-dryrun`, which ends before LOSAT's refusals) reads a value that
// LOSAT refuses for not being UTF-8 as its lossy text (`Args::with_lossy_text`); the run without
// `-dryrun` is refused.
#[cfg(any(unix, windows))]
#[test]
fn a_dry_run_reads_a_value_that_is_not_utf8_as_its_lossy_text() {
    let argv = |dry_run: bool| {
        let mut argv: Vec<OsString> = ["losat", "blastn", "-db"].map(OsString::from).to_vec();
        argv.push(not_utf8());
        argv.extend(dry_run.then(|| OsString::from("-dryrun")));
        crate::cli::try_parse_with_own_start(argv, &DiagTarget::Stderr)
    };
    let db = match argv(true).unwrap().command {
        Commands::Blastn(args) => args.db,
        command => panic!("{command:?}"),
    };
    assert_eq!(db.as_deref(), Some("a\u{FFFD}b"));
    assert!(matches!(
        argv(false),
        Err(crate::cli::ArgvError::Refused(_))
    ));
}

// The arguments that LOSAT Web's ABI v1 checked (`checked`): a value that its key's kind does not
// convert, or a name that is not the program's, gives none; a value that the option's type does
// not hold is LOSAT's error (`not_converted`), whose text is pinned here.
#[test]
fn checked_arguments_and_values_that_do_not_convert() {
    let blastx = tool("blastx").unwrap();
    let given = |name: &str, value: &str| vec![(name.to_string(), Some(OsString::from(value)))];
    assert!(checked(blastx, &given("evalue", "1e-5")).is_some());
    assert!(checked(blastx, &given("evalue", "x")).is_none());
    assert!(checked(blastx, &given("num_threads", "1.5")).is_none());
    assert!(checked(blastx, &given("nosuch", "1")).is_none());
    // `-num_descriptions -1` is an integer, which the option's `usize` does not hold.
    let args = checked(blastx, &given("num_descriptions", "-1")).unwrap();
    let text = "LOSAT's BLASTX does not convert the value of -num_descriptions that NCBI's \
                argument reader read";
    assert_eq!(not_converted(blastx, "num_descriptions"), text);
    assert_eq!(options(blastx, &args, false).unwrap_err(), text);
}

/// The arguments that NCBI's reader reads from `words` (after the toolkit's words that the
/// application start takes); `None` where the reader refuses them or prints a text.
fn read(tool: &Tool, words: &[OsString]) -> Option<Args> {
    match preparse(words) {
        Preparsed::Words { words, .. } => parse(tool, &words).ok(),
        Preparsed::Version(_) => None,
    }
}

/// The words that a program's command line needs besides the argument `given`: makeblastdb's
/// mandatory `-dbtype` (`-in` takes its default).
fn base(program: &str, given: &str) -> Vec<OsString> {
    match program {
        "makeblastdb" if given != "dbtype" => vec!["-dbtype".into(), "nucl".into()],
        _ => Vec::new(),
    }
}

/// The values of a key of `kind` that the golden tries, in their order: its choices, then values
/// of its kind.
fn values(kind: Kind, choices: &[String]) -> Vec<OsString> {
    let more: &[&str] = match kind {
        Kind::Flag => &[],
        Kind::Boolean => &["true", "false"],
        Kind::Integer => &[
            "2",
            "3",
            "1",
            "5",
            "11",
            "16",
            "28",
            "4",
            "6",
            "18",
            "21",
            "33",
            "0",
            "255",
            "256",
            "-1",
            "2147483647",
        ],
        Kind::Int8 => &["2", "1", "0", "-1", "9223372036854775807"],
        Kind::Double => &[
            "0.5", "2", "1e-5", "0.05", "1", "10", "99.5", "100", "0", "-1", "101", "1e308",
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
        Kind::InputFile | Kind::OutputFile => &["x", "in/q.fa", "-"],
    };
    choices
        .iter()
        .map(OsString::from)
        .chain(more.iter().map(OsString::from))
        .collect()
}

/// `words` (with the argument `name` given) and the arguments that `name` requires in `tool`'s
/// program (`Dep::Requires`), and those that they require, each with the first value with which
/// NCBI's reader takes the words.
fn with_required(tool: &Tool, name: &str, mut words: Vec<OsString>) -> Vec<OsString> {
    for (first, dep, second) in crate::ncbi_app::args::dependencies(tool.program) {
        let word = OsString::from(format!("-{second}"));
        if first != name || dep != Dep::Requires || words.contains(&word) {
            continue;
        }
        let Some(described) = described(tool).find(|arg| arg.name == second) else {
            continue;
        };
        let given: Vec<Option<OsString>> = match described.kind {
            Kind::Flag => vec![None],
            kind => values(kind, &described.choices)
                .into_iter()
                .map(Some)
                .collect(),
        };
        let with = |value: &Option<OsString>| {
            let mut tried = words.clone();
            tried.push(word.clone());
            tried.extend(value.clone());
            with_required(tool, second, tried)
        };
        words = given
            .iter()
            .map(with)
            .find(|tried| read(tool, tried).is_some())
            .unwrap_or_else(|| with(&given[0]));
    }
    words
}

/// The lines of `command`'s `Debug` (`{:#?}`): one per value, `path = value`, the path the names
/// of the fields (and the wrappers: `Some`) that hold it, below the program's type.
fn lines(command: &Commands) -> Vec<String> {
    let text = format!("{command:#?}");
    let mut path: Vec<String> = Vec::new();
    let mut lines = Vec::new();
    for line in text.lines() {
        let item = line.trim();
        if item.starts_with(['}', ')', ']']) {
            path.pop();
            continue;
        }
        let (name, value) = match item.split_once(": ") {
            Some((name, value)) if !name.starts_with('"') => (Some(name), value),
            _ => (None, item),
        };
        if let Some(wrapper) = value.strip_suffix(['{', '(', '[']) {
            path.push(name.unwrap_or(wrapper.trim_end()).to_string());
            continue;
        }
        let mut held: Vec<&str> = path.iter().skip(2).map(String::as_str).collect();
        held.extend(name);
        lines.push(format!(
            "{} = {}",
            held.join("."),
            value.strip_suffix(',').unwrap_or(value)
        ));
    }
    lines
}

/// The lines of `lines` that are not in `other` (each line as many times as it is there).
fn missing(lines: &[String], other: &[String]) -> Vec<String> {
    let mut other: Vec<&String> = other.iter().collect();
    lines
        .iter()
        .filter(|line| match other.iter().position(|seen| seen == line) {
            Some(at) => {
                other.remove(at);
                false
            }
            None => true,
        })
        .cloned()
        .collect()
}

/// A word as the golden writes it: in quotes when it is empty or has a space.
fn shown(word: &OsString) -> String {
    let word = word.to_string_lossy();
    match word.is_empty() || word.contains(' ') {
        true => format!("{word:?}"),
        false => word.into_owned(),
    }
}

/// The options that `convert` gives (`options`, or what gave them before), of `tool`'s base words
/// and of each argument given with them, as the golden writes them: `[program]` and the lines of
/// the options at their defaults, then for each argument `[program words]` and the lines that it
/// takes away (`-`) and adds (`+`). An argument's value is the first of `values` that NCBI's reader
/// and LOSAT take and that changes the options (or the first taken, else none: `(none taken)`).
pub(crate) fn golden(
    tool: &Tool,
    convert: impl Fn(&Tool, &Args) -> Result<Commands, String>,
) -> String {
    let program = tool.program;
    let of = |words: &[OsString]| -> Option<Vec<String>> {
        let args = read(tool, words)?;
        if refusal_after_read(tool, &args).is_some() {
            return None;
        }
        convert(tool, &args).ok().map(|command| lines(&command))
    };
    let defaults = of(&base(program, "")).expect("the base words read");
    let mut text = format!("[{program}]\n");
    for line in &defaults {
        text += &format!("{line}\n");
    }
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
        let mut taken: Option<(Vec<OsString>, Vec<String>)> = None;
        for value in given {
            let mut words = base(program, name);
            words.push(format!("-{name}").into());
            words.extend(value);
            let words = with_required(tool, name, words);
            let Some(options) = of(&words) else {
                continue;
            };
            let changes = options != defaults;
            if changes || taken.is_none() {
                taken = Some((words, options));
            }
            if changes {
                break;
            }
        }
        let Some((words, options)) = taken else {
            text += &format!("[{program} -{name}] (none taken)\n");
            continue;
        };
        let words: Vec<String> = words.iter().map(shown).collect();
        text += &format!("[{program} {}]\n", words.join(" "));
        for line in missing(&defaults, &options) {
            text += &format!("- {line}\n");
        }
        for line in missing(&options, &defaults) {
            text += &format!("+ {line}\n");
        }
    }
    text
}

/// The golden: the options of each program on NCBI's reader (`golden`), frozen from the options
/// that clap converted before line RG item K4i-3c (each program's `from_args` gave the same on
/// every case of its proof, and gives the same here).
const GOLDEN: &str = include_str!("../../tests/fixtures/options/golden.txt");

#[test]
fn the_options_of_each_program_are_frozen() {
    let mut expected = GOLDEN.lines().filter(|line| !line.starts_with('#'));
    for tool in TOOLS.iter().copied() {
        let text = golden(tool, |tool, args| options(tool, args, false));
        for (at, line) in text.lines().enumerate() {
            assert_eq!(
                Some(line),
                expected.next(),
                "{} line {at} of its options",
                tool.program
            );
        }
    }
    assert_eq!(expected.next(), None);
}

/// The golden's first lines.
const GOLDEN_HEAD: &str = "\
# The options of the ten programs on NCBI's reader (src/cli/options_tests.rs, `golden`): for each
# program its options at their defaults, `[program]`, one line per value (the path of fields that
# holds it = its Debug), then for each argument the words given and the lines that they take away
# (-) and add (+). Frozen from clap's conversion before line RG item K4i-3c; checked by
# `the_options_of_each_program_are_frozen`. After a deliberate change of an option type, rewrite it
# with `cargo test --lib write_the_options_golden -- --ignored` and review its diff.
";

#[test]
#[ignore = "rewrites tests/fixtures/options/golden.txt after a deliberate change of the options"]
fn write_the_options_golden() {
    let mut text = GOLDEN_HEAD.to_string();
    for tool in TOOLS.iter().copied() {
        text += &golden(tool, |tool, args| options(tool, args, false));
    }
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/options/golden.txt"
    );
    std::fs::write(path, text).unwrap();
}
