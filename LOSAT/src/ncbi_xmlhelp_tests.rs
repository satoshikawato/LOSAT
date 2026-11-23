//! `-xmlhelp` (line RG item E2x): NCBI's XML of each program without the arguments that its reader
//! refuses by name, against the readers themselves.

use super::*;
use crate::blastinput::ncbi_environment::DiagTarget;
use crate::cli::{run_refusal, try_parse_from, xml_help_refusal, xml_help_text, ArgvError};
use crate::ncbi_app::args::{refusal, Refusal};
use std::ffi::OsString;

/// The programs whose arguments NCBI's reader reads and clap converts (`cli::reader`).
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

// The arguments of NCBI's XML, which NCBI's reader reads (its tables are made from the XML), and
// those that clap converts are one table (line RG DECISIONS K4i-2b-2 (2)): each argument of the
// XML is converted by clap under its name, as a key or a flag as the XML says, or refused by name
// after the read (`Refusal::Read`, never converted); clap has no other argument. The XML that
// `-xmlhelp` prints keeps exactly the arguments not refused (by name or at the run, `Refusal::Run`).
#[test]
fn the_xml_names_are_the_arguments_that_clap_converts() {
    for program in CLAP_PROGRAMS {
        let xml = xml_help_text(program).unwrap();
        let command = crate::cli::command_of(program);
        let scope = command.find_subcommand(program).unwrap();
        let converted = |name: &str| {
            scope
                .get_arguments()
                .find(|arg| arg.get_long() == Some(name))
        };
        let names = ncbi_argument_names(program);
        for &(name, flag) in &names {
            let kept = has_argument(&xml, name);
            if READ_BEFORE_OR_HELP.contains(&name) || name == "xmlhelp" || name == "help-full" {
                assert!(kept, "{program} -{name}");
                continue;
            }
            let row = refusal(program, name);
            let by_name = matches!(row, Some(Refusal::Read(_)));
            let arg = converted(name);
            assert_eq!(arg.is_some(), !by_name, "{program} -{name}");
            if let Some(arg) = arg {
                assert_eq!(
                    arg.get_action().takes_values(),
                    !flag,
                    "{program} -{name}: a key of the XML is a key of clap"
                );
            }
            assert_eq!(kept, row.is_none(), "{program} -{name}");
        }
        // clap converts nothing that the reader does not read (the hooks of LOSAT Web are no
        // argument: `#[arg(skip)]`).
        for arg in scope.get_arguments() {
            let Some(long) = arg.get_long() else {
                continue;
            };
            assert!(
                names.iter().any(|&(name, _)| name == long),
                "{program} --{long} is no argument of NCBI's XML"
            );
        }
    }
}

// The run refusals' texts (`run_refusal`), which the checks of the runs write.
#[test]
fn run_refusals_name_the_argument_and_the_program() {
    assert_eq!(
        run_refusal("PSIBLAST", "use_sw_tback"),
        "-use_sw_tback (the Smith-Waterman traceback) is not supported by LOSAT's PSIBLAST"
    );
    assert_eq!(
        run_refusal("RPSTBLASTN", "ungapped"),
        "-ungapped (NCBI BLAST+ 2.17.0's rpstblastn crashes after an ungapped search) is not \
         supported by LOSAT's RPSTBLASTN"
    );
}

// The XML is NCBI's lines (in LOSAT's words) in NCBI's order, of which only whole elements of the arguments refused
// and the dependencies that name one are missing.
#[test]
fn the_xml_is_ncbis_without_the_refused_arguments() {
    for program in CLAP_PROGRAMS {
        // NCBI's XML in LOSAT's words (the program's description names LOSAT's version).
        let ncbi = crate::wording::losat_text(ncbi_xml(program).unwrap()).into_owned();
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
        let expected: Vec<(String, String)> = dependencies(&ncbi)
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

// Without a refusal the XML is NCBI's bytes in LOSAT's words, for every program; with every argument refused, what
// is left is the program's part and an empty argument list (NCBI writes no `<dependencies>`
// without one).
#[test]
fn the_filter_leaves_whole_elements() {
    for program in programs() {
        // NCBI's XML in LOSAT's words (the program's description names LOSAT's version).
        let ncbi = crate::wording::losat_text(ncbi_xml(program).unwrap()).into_owned();
        assert_eq!(print(program, |_| false).unwrap(), ncbi, "{program}");
        let none = print(program, |_| true).unwrap();
        let head = &ncbi[..ncbi.find("<arguments>\n").unwrap()];
        assert_eq!(
            none,
            format!("{head}<arguments>\n</arguments>\n</ncbi_application>\n"),
            "{program}"
        );
    }
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

// The library API and the Web engine's ABI v1 refuse `-xmlhelp` by name where the reader takes it as
// a flag; as a value it is a value, and after `-h` it is not read. The Web engine's ABI v2 prints the
// XML (`the_web_engine_prints_the_command_lines_help`).
#[test]
fn the_library_refuses_xmlhelp() {
    for program in CLAP_PROGRAMS {
        let error = try_parse_from(["losat", program, "-xmlhelp"]).unwrap_err();
        let shown = crate::cli::reader::tool(program).unwrap().refusal_name;
        assert_eq!(
            error,
            ArgvError::Refused(xml_help_refusal(shown)),
            "{program}"
        );
        // The printed `-h` may name `-xmlhelp` (rpsbproc's synopsis does); it is not refused.
        let short = try_parse_from(["losat", program, "-h", "-xmlhelp"]);
        assert!(
            !matches!(&short, Err(ArgvError::Refused(text)) if text.contains("-xmlhelp")),
            "{program}"
        );
    }
    let value = try_parse_from(["losat", "blastn", "-out", "-xmlhelp"]);
    assert!(!matches!(&value, Err(error) if error.to_string().contains("-xmlhelp")));
    type RunWith = fn(
        &[OsString],
        &mut dyn std::io::Read,
        &mut dyn std::io::Write,
        &mut dyn std::io::Write,
    ) -> i32;
    let readers: [(&str, RunWith); 8] = [
        ("BLASTDBCMD", crate::blastdbcmd::run_with),
        ("DUSTMASKER", crate::dustmasker::run_with),
        ("SEGMASKER", crate::segmasker::run_with),
        ("WINDOWMASKER", crate::windowmasker::run_with),
        ("CONVERT2BLASTMASK", crate::convert2blastmask::run_with),
        ("MAKEMBINDEX", crate::makembindex::run_with),
        ("MAKEPROFILEDB", crate::makeprofiledb::run_with),
        ("BLAST_FORMATTER", crate::blast_formatter::run_with),
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

// The readers of blastdbcmd, dustmasker, segmasker, windowmasker, convert2blastmask,
// makembindex, makeprofiledb and blast_formatter (`ncbi_app`) read every argument of NCBI's XML
// (none is refused by name; their XML is NCBI's whole but for makeprofiledb's `-binary`, which
// its run refuses, `Refusal::Run`, and blast_formatter's `-rid`, which its run refuses): the
// argument with a value `1` for a key, then `-h`, writes the usage or an error of the value,
// never NCBI's
// `Unknown argument` or a refusal.
#[test]
fn the_ncbi_readers_read_every_argument_of_the_xml() {
    type RunWith = fn(
        &[OsString],
        &mut dyn std::io::Read,
        &mut dyn std::io::Write,
        &mut dyn std::io::Write,
    ) -> i32;
    let readers: [(&str, RunWith); 8] = [
        ("blastdbcmd", crate::blastdbcmd::run_with),
        ("dustmasker", crate::dustmasker::run_with),
        ("segmasker", crate::segmasker::run_with),
        ("windowmasker", crate::windowmasker::run_with),
        ("convert2blastmask", crate::convert2blastmask::run_with),
        ("makembindex", crate::makembindex::run_with),
        ("makeprofiledb", crate::makeprofiledb::run_with),
        ("blast_formatter", crate::blast_formatter::run_with),
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

// Line WQA's D14 (2026-11-20): the reader of LOSAT Web's ABI v2 stops at the help words as the
// command line's does, `-xmlhelp` included, and the run prints the text that the command line's
// entry prints for them (`ncbi_app::read`, which `cli::main_entry` calls): NCBI's help, and the
// XML of `-xmlhelp`. The words before them are read first.
#[test]
fn the_web_engine_prints_the_command_lines_help() {
    use crate::cli::{try_parse_with_printing_start, WEB_V2_PROGRAMS};
    for program in WEB_V2_PROGRAMS {
        let web = |words: &[&str]| {
            let argv = ["losat", program].into_iter().chain(words.iter().copied());
            try_parse_with_printing_start(argv, &DiagTarget::Stderr).map(|_| ())
        };
        let xml = web(&["-evalue", "1", "-xmlhelp", "-bogus"]).unwrap_err();
        assert_eq!(
            xml,
            ArgvError::Printed(xml_help_text(program).unwrap()),
            "{program}"
        );
        let tool = crate::cli::reader::tool(program).unwrap();
        for (word, text) in [
            ("-h", crate::ncbi_app::short_usage(tool)),
            ("-help", crate::ncbi_app::help_text(tool.help).into_owned()),
            (
                "-help-full",
                crate::ncbi_app::help_text(tool.help).into_owned(),
            ),
        ] {
            let help = web(&["-dryrun", word]).unwrap_err();
            assert_eq!(help, ArgvError::Printed(text.clone()), "{program} {word}");
            let command_line = try_parse_from(["losat", program, word]).unwrap_err();
            assert_eq!(command_line, ArgvError::Printed(text), "{program} {word}");
        }
        let refused = web(&["-evalue", "x", "-xmlhelp"]).unwrap_err();
        assert!(matches!(refused, ArgvError::Report(_)), "{program}");
        // The library refuses `-xmlhelp`, naming the places that print it.
        assert!(xml_help_refusal(&program.to_ascii_uppercase())
            .contains("runs on the command line and in LOSAT Web of LOSAT's"));
    }
}
