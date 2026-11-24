//! Every argument of the ten programs on NCBI's reader that LOSAT does not refuse by name is read
//! into the program's options (`options`, each program's `from_args`; line RG item K4i-3c), and
//! the options read no other argument.

use std::collections::BTreeSet;
use std::ffi::OsString;

use super::{options, TOOLS};
use crate::ncbi_app::args::{described, parse};
use crate::ncbi_app::Tool;

/// The names of the arguments that `tool`'s options read (`options`, `Args::names_read`), of a
/// command line with only what the program needs (makeblastdb's mandatory `-dbtype`; `-in` takes
/// its default): the tests that compare the arguments of the options with those of NCBI's XML,
/// LOSAT's descriptions and the search strategy's walk read them here.
pub(crate) fn names_read(tool: &Tool) -> BTreeSet<&'static str> {
    let words: Vec<OsString> = match tool.program {
        "makeblastdb" => vec!["-dbtype".into(), "nucl".into()],
        _ => Vec::new(),
    };
    let args = parse(tool, &words).expect("the base words read");
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
