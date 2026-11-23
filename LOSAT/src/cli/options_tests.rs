//! Every argument of the ten programs on NCBI's reader that LOSAT does not refuse by name is read
//! into the program's options (`options`, each program's `from_args`; line RG item K4i-3c).

use std::collections::BTreeSet;
use std::ffi::OsString;

use super::{options, TOOLS};
use crate::blastinput::ncbi_environment::{preparse, Preparsed};
use crate::ncbi_app::args::{described, parse, Args};
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

fn words(words: &[&str]) -> Vec<OsString> {
    words.iter().map(OsString::from).collect()
}

/// The words that a program's command line needs before the swept argument: makeblastdb's
/// mandatory `-dbtype` (`-in` takes its default).
fn base(program: &str, swept: &str) -> Vec<OsString> {
    match program {
        "makeblastdb" if swept != "dbtype" => words(&["-dbtype", "nucl"]),
        _ => Vec::new(),
    }
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
