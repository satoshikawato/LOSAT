//! A search: a program and its argv, checked by the command line's parser.
//!
//! The command line is the one definition of the options (NCBI's `CBlastAppArgs`); the library
//! API takes the same words and parses them with the same parser, so an option means, and fails,
//! exactly as on the command line (LOSAT Web ABI v2 §7 does the same).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,51,75
//! ```c++
//! const string kArgQuery("query");
//! const string kArgSubject("subject");
//! const string kArgNumThreads("num_threads");
//! ```
#![warn(missing_docs, unused)]

use super::error::Error;
use super::program::Program;
use crate::cli::{Cli, Commands};

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
// ```c
// arg_desc.AddDefaultKey(kArgOutput, "output_file",
//                "Output file name",
//                CArgDescriptions::eOutputFile, "-");
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
// ```c++
//     arg_desc.AddDefaultKey(kArgOutputFormat, "format",
//                            kOutputFormatDescription,
//                            CArgDescriptions::eString,
//                            NStr::IntToString(dft_outfmt));
// ```
/// The options that the library API sets itself and does not take in an argv: the inputs
/// (`Inputs` gives the records and their names), the output file and the output format
/// (arguments of `run_collect` and `run`).
pub const API_OWNED: &[&str] = &["query", "subject", "db", "out", "outfmt"];

/// Options that the library API rejects with a reason of its own.
const API_REJECTED: [(&str, &str); 1] = [(
    "parse_deflines",
    "records are read before the options are known (`read_fasta`), so their ids cannot be parsed",
)];

/// A program and an argv that the command line's parser accepted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Search {
    program: Program,
    argv: Vec<String>,
}

impl Search {
    /// The program.
    pub fn program(&self) -> Program {
        self.program
    }

    /// The argv as given: the options after the program name.
    pub fn argv(&self) -> &[String] {
        &self.argv
    }

    /// The parsed arguments of the search with the inputs named `query` and `subject` (the
    /// command line's `-query` and `-subject` values, shown in the reports). `validate`
    /// accepted the argv, so this parses.
    pub(crate) fn command(&self, query: &str, subject: &str) -> Commands {
        parse(self.program, &self.argv, query, subject).expect("validated argv")
    }
}

/// The key of an argv word: the option name, without a value joined by `=`.
fn key(word: &str) -> &str {
    word.split('=').next().unwrap_or("")
}

/// Parses the argv with the command line's parser, as `losat <program> -query <query>
/// -subject <subject> <argv...>`.
fn parse(program: Program, argv: &[String], query: &str, subject: &str) -> Result<Commands, Error> {
    let words = [
        "losat",
        program.name(),
        "-query",
        query,
        "-subject",
        subject,
    ];
    let cli: Cli =
        crate::cli::try_parse_from(words.iter().copied().chain(argv.iter().map(String::as_str)))
            .map_err(|error| Error::Argv(crate::cli::render_message(&error)))?;
    Ok(cli.command)
}

/// Checks an argv for `program`: the options after the program name, as on the command line
/// (`-evalue 1e-5 -num_threads 4`). Every word goes to the command line's parser, so an
/// option means, and fails, exactly as there; the error is the parser's message.
///
/// The library API sets the inputs and the outputs itself, so `-query`, `-subject`, `-db`,
/// `-out` and `-outfmt` are refused (`API_OWNED`), as is `-parse_deflines` (the records are
/// read before the options are known). The checks of the options that the command line makes
/// before it searches (ranges, scoring, the hit list size of each format) come with `run` and
/// `run_collect`, which know the formats and the inputs, and report them with the command
/// line's messages.
pub fn validate(program: Program, argv: &[&str]) -> Result<Search, Error> {
    for word in argv {
        let key = key(word);
        let Some(name) = key.strip_prefix('-') else {
            continue;
        };
        let name = name.strip_prefix('-').unwrap_or(name);
        if API_OWNED.contains(&name) {
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: the inputs are given as records and the outputs as formats"
            )));
        }
        if let Some((_, reason)) = API_REJECTED.iter().find(|(flag, _)| *flag == name) {
            return Err(Error::Argv(format!(
                "{key} is not supported by the library API yet: {reason}"
            )));
        }
    }
    let argv: Vec<String> = argv.iter().map(|word| word.to_string()).collect();
    parse(program, &argv, "query", "subject")?;
    Ok(Search { program, argv })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_argv_is_parsed_as_on_the_command_line() {
        let search = validate(Program::Blastp, &["-evalue", "0.1", "-num_threads", "2"]).unwrap();
        assert_eq!(search.program(), Program::Blastp);
        assert_eq!(search.argv(), ["-evalue", "0.1", "-num_threads", "2"]);
        match search.command("q.faa", "s.faa") {
            Commands::Blastp(args) => {
                assert_eq!(args.evalue, Some(0.1));
                assert_eq!(args.num_threads, 2);
                assert_eq!(args.query, std::path::PathBuf::from("q.faa"));
                assert_eq!(args.subject.as_deref(), Some(std::path::Path::new("s.faa")));
            }
            other => panic!("{other:?}"),
        }
        // An empty argv is the command line's defaults.
        validate(Program::Tblastx, &[]).unwrap();
        // A parser error is the command line's message.
        let error = validate(Program::Blastp, &["-evalue"]).unwrap_err();
        let cli = crate::cli::try_parse_from::<Cli, _, _>([
            "losat", "blastp", "-query", "query", "-subject", "subject", "-evalue",
        ])
        .unwrap_err();
        assert_eq!(error, Error::Argv(crate::cli::render_message(&cli)));
        let error = validate(Program::Blastn, &["-nosuchoption", "1"]).unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{error}");
    }

    #[test]
    fn the_options_the_api_owns_are_refused() {
        for word in [
            "-query",
            "--query",
            "-subject",
            "-db",
            "-out",
            "-outfmt",
            "-outfmt=6",
        ] {
            let error = validate(Program::Blastp, &[word, "x"]).unwrap_err();
            assert!(
                error.message().contains("not accepted by the library API"),
                "{word}: {error}"
            );
            assert!(error.message().starts_with(key(word)), "{word}: {error}");
        }
        let error = validate(Program::Blastp, &["-parse_deflines"]).unwrap_err();
        assert!(
            error
                .message()
                .contains("not supported by the library API yet"),
            "{error}"
        );
        // A value that looks like an option is a value.
        validate(Program::Blastp, &["-matrix", "BLOSUM62"]).unwrap();
    }
}
