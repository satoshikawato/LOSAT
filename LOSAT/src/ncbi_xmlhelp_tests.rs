//! `-xmlhelp` (line RG item E2x): NCBI's XML of each program without the arguments that its reader
//! refuses by name, against the readers themselves.

use super::*;
use crate::cli::{run_refusal, try_parse_from, xml_help_refusal, xml_help_text, Cli, RUN_REFUSALS};
use clap::error::ErrorKind;
use std::ffi::OsString;

/// The programs that `cli::try_parse_from` reads.
const CLAP_PROGRAMS: [&str; 10] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "makeblastdb",
    "rpsbproc",
];

/// NCBI's arguments that its application start reads before the arguments
/// (`ncbi_environment::preparse`), which the parser of the Web engine (`try_parse_from`: no
/// application start, as in these tests) refuses, and the help flags.
const READ_BEFORE_OR_HELP: [&str; 9] = [
    "version",
    "version-full",
    "version-full-xml",
    "version-full-json",
    "logfile",
    "conffile",
    "dryrun",
    "h",
    "help",
];

/// Whether `xml` has the element of the argument `name`.
fn has_argument(xml: &str, name: &str) -> bool {
    xml.split_inclusive('\n')
        .filter_map(argument_start)
        .any(|(_, found)| found == name)
}

/// The two names of each dependency of `xml` (`<arg1>`, `<arg2>`).
fn dependencies(xml: &str) -> Vec<(String, String)> {
    let lines: Vec<&str> = xml.split_inclusive('\n').collect();
    lines
        .windows(3)
        .filter(|window| {
            matches!(
                window[0],
                "<first_requires_second>\n" | "<first_excludes_second>\n"
            )
        })
        .map(|window| {
            (
                xml_line_value(window[1]).unwrap().1.to_string(),
                xml_line_value(window[2]).unwrap().1.to_string(),
            )
        })
        .collect()
}

// Each argument of NCBI's XML is read by the parser or refused by name (none is unknown), and
// the XML keeps exactly the ones read: the parser is asked with the argument (a value `1` for a
// key) before `-h`, whose help ends the reading when the argument passed.
#[test]
fn the_xml_keeps_the_arguments_that_the_parser_reads() {
    for program in CLAP_PROGRAMS {
        let xml = xml_help_text(program).unwrap();
        for (name, flag) in ncbi_argument_names(program) {
            let kept = has_argument(&xml, name);
            if READ_BEFORE_OR_HELP.contains(&name) {
                assert!(kept, "{program} -{name}");
                continue;
            }
            let word = format!("-{name}");
            let mut argv = vec!["losat", program, word.as_str()];
            if !flag {
                argv.push("1");
            }
            argv.push("-h");
            let refused = match try_parse_from::<Cli, _, _>(argv) {
                Ok(_) => false,
                Err(error) => {
                    assert_ne!(
                        error.kind(),
                        ErrorKind::UnknownArgument,
                        "{program} -{name}: NCBI's argument is neither read nor refused by name"
                    );
                    let text = error.to_string();
                    text.contains(&format!("option -{name} is not supported by LOSAT's"))
                        || text.contains(&format!("unsupported BLASTX option '-{name}'"))
                }
            };
            // The arguments whose run LOSAT refuses whatever their value pass the parser.
            let run_refused = RUN_REFUSALS
                .iter()
                .any(|refusal| refusal.0 == program && refusal.1 == name);
            assert!(!(refused && run_refused), "{program} -{name}");
            assert_eq!(kept, !refused && !run_refused, "{program} -{name}");
        }
    }
    // Each run refusal is an argument of NCBI's XML of its program.
    for (program, name, _) in RUN_REFUSALS {
        assert!(
            ncbi_argument_names(program)
                .iter()
                .any(|(found, _)| found == name),
            "{program} -{name}"
        );
    }
}

// The run refusals' texts (`run_refusal`), which the checks of the runs write.
#[test]
fn run_refusals_name_the_argument_and_the_program() {
    assert_eq!(
        run_refusal("BLASTP", "use_sw_tback"),
        "-use_sw_tback (the Smith-Waterman traceback) is not supported by LOSAT's BLASTP"
    );
    assert_eq!(
        run_refusal("TBLASTN", "ungapped"),
        "-ungapped (an ungapped search) is not supported by LOSAT's TBLASTN"
    );
}

// The XML is NCBI's lines in NCBI's order, of which only whole elements of the arguments refused
// and the dependencies that name one are missing.
#[test]
fn the_xml_is_ncbis_without_the_refused_arguments() {
    for program in CLAP_PROGRAMS {
        let ncbi = ncbi_xml(program).unwrap();
        let xml = xml_help_text(program).unwrap();
        let gone: Vec<&str> = ncbi_argument_names(program)
            .into_iter()
            .map(|(name, _)| name)
            .filter(|name| !has_argument(&xml, name))
            .collect();
        assert_eq!(
            xml,
            print(program, |name| gone.contains(&name)).unwrap(),
            "{program}"
        );
        // NCBI's lines, in order.
        let mut ncbi_lines = ncbi.split_inclusive('\n');
        for line in xml.split_inclusive('\n') {
            assert!(
                ncbi_lines.any(|ncbi_line| ncbi_line == line),
                "{program}: {line:?}"
            );
        }
        // The dependencies kept are NCBI's of two names kept.
        let expected: Vec<(String, String)> = dependencies(ncbi)
            .into_iter()
            .filter(|(a, b)| !gone.contains(&a.as_str()) && !gone.contains(&b.as_str()))
            .collect();
        assert_eq!(dependencies(&xml), expected, "{program}");
        assert!(
            !xml.contains("<dependencies>\n</dependencies>\n"),
            "{program}"
        );
        assert!(
            xml.ends_with("</arguments>\n</ncbi_application>\n"),
            "{program}"
        );
    }
}

// Without a refusal the XML is NCBI's bytes, for every program; with every argument refused, what
// is left is the program's part and an empty argument list (NCBI writes no `<dependencies>`
// without one).
#[test]
fn the_filter_leaves_whole_elements() {
    for program in programs() {
        let ncbi = ncbi_xml(program).unwrap();
        assert_eq!(print(program, |_| false).unwrap(), ncbi, "{program}");
        let none = print(program, |_| true).unwrap();
        let head = &ncbi[..ncbi.find("<arguments>\n").unwrap()];
        assert_eq!(
            none,
            format!("{head}<arguments>\n</arguments>\n</ncbi_application>\n"),
            "{program}"
        );
    }
    assert_eq!(print("windowmasker", |_| false), None);
}

// One refused name: its element and the dependencies naming it go, the rest stays (NCBI blastn:
// `-remote` has no element of its own besides its flag, and excludes or is required by others).
#[test]
fn a_refused_name_takes_its_dependencies() {
    let ncbi = ncbi_xml("blastn").unwrap();
    let xml = print("blastn", |name| name == "remote").unwrap();
    assert!(!has_argument(&xml, "remote"));
    assert!(!xml.contains("<arg1>remote</arg1>") && !xml.contains("<arg2>remote</arg2>"));
    let named = dependencies(ncbi)
        .iter()
        .filter(|(a, b)| a == "remote" || b == "remote")
        .count();
    assert!(named > 0);
    assert_eq!(dependencies(&xml).len(), dependencies(ncbi).len() - named);
    assert_eq!(
        ncbi_argument_names("blastn").len() - 1,
        xml.split_inclusive('\n').filter_map(argument_start).count()
    );
}

// The library API and the Web engine (no application start) refuse `-xmlhelp` by name where the
// reader takes it as a flag; as a value it is a value, and after `-h` it is not read.
#[test]
fn the_library_refuses_xmlhelp() {
    for program in CLAP_PROGRAMS {
        let error = try_parse_from::<Cli, _, _>(["losat", program, "-xmlhelp"]).unwrap_err();
        let shown = if program == "makeblastdb" {
            program.to_string()
        } else {
            program.to_ascii_uppercase()
        };
        assert_eq!(
            error.to_string().trim_end(),
            format!("error: {}", xml_help_refusal(&shown)),
            "{program}"
        );
        let short = try_parse_from::<Cli, _, _>(["losat", program, "-h", "-xmlhelp"]);
        assert!(
            !matches!(&short, Err(error) if error.to_string().contains("-xmlhelp")),
            "{program}"
        );
    }
    let value = try_parse_from::<Cli, _, _>(["losat", "blastn", "-out", "-xmlhelp"]);
    assert!(!matches!(&value, Err(error) if error.to_string().contains("-xmlhelp")));
    type RunWith = fn(
        &[OsString],
        &mut dyn std::io::Read,
        &mut dyn std::io::Write,
        &mut dyn std::io::Write,
    ) -> i32;
    let readers: [(&str, RunWith); 4] = [
        ("BLASTDBCMD", crate::blastdbcmd::run_with),
        ("DUSTMASKER", crate::dustmasker::run_with),
        ("SEGMASKER", crate::segmasker::run_with),
        ("CONVERT2BLASTMASK", crate::convert2blastmask::run_with),
    ];
    for (shown, run_with) in readers {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = run_with(
            &[OsString::from("-xmlhelp")],
            &mut std::io::empty(),
            &mut out,
            &mut err,
        );
        assert_eq!(status, 1, "{shown}");
        assert!(out.is_empty(), "{shown}");
        assert_eq!(
            String::from_utf8(err).unwrap(),
            format!("Error: {}\n", xml_help_refusal(shown)),
        );
    }
}

// The readers of blastdbcmd, dustmasker, segmasker and convert2blastmask (`ncbi_app`) read every
// argument of NCBI's XML (none is refused by name, so their XML is NCBI's whole): the argument
// with a value `1` for a key, then `-h`, writes the usage or an error of the value, never NCBI's
// `Unknown argument` or a refusal.
#[test]
fn the_ncbi_readers_read_every_argument_of_the_xml() {
    type RunWith = fn(
        &[OsString],
        &mut dyn std::io::Read,
        &mut dyn std::io::Write,
        &mut dyn std::io::Write,
    ) -> i32;
    let readers: [(&str, RunWith); 4] = [
        ("blastdbcmd", crate::blastdbcmd::run_with),
        ("dustmasker", crate::dustmasker::run_with),
        ("segmasker", crate::segmasker::run_with),
        ("convert2blastmask", crate::convert2blastmask::run_with),
    ];
    for (program, run_with) in readers {
        for (name, flag) in ncbi_argument_names(program) {
            if READ_BEFORE_OR_HELP.contains(&name) {
                continue;
            }
            let mut words = vec![OsString::from(format!("-{name}"))];
            if !flag {
                words.push("1".into());
            }
            words.push("-h".into());
            let (mut out, mut err) = (Vec::new(), Vec::new());
            run_with(&words, &mut std::io::empty(), &mut out, &mut err);
            let err = String::from_utf8_lossy(&err);
            assert!(
                !err.contains("Unknown argument") && !err.contains("not supported by LOSAT"),
                "{program} -{name}: {err}"
            );
        }
    }
}
