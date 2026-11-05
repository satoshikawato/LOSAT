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
use super::input::Inputs;
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
/// (arguments of `run_collect`, `run` and `run_to_file`).
pub const API_OWNED: &[&str] = &["query", "subject", "db", "out", "outfmt"];

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

    /// The parsed arguments, with either kind of subject: `validate` accepted one of them.
    fn any_subject(&self) -> Commands {
        let names = ["query", "subject"];
        parse(self.program, &self.argv, names[0], SUBJECT, names[1])
            .or_else(|_| parse(self.program, &self.argv, names[0], DATABASE, names[1]))
            .expect("validated argv")
    }

    /// Whether the search's argv has `-parse_deflines` (`read_inputs`).
    pub(crate) fn parse_deflines(&self) -> bool {
        match self.any_subject() {
            Commands::Blastp(args) => args.parse_deflines,
            Commands::Tblastn(args) => args.parse_deflines,
            Commands::Blastn(args) => args.parse_deflines,
            Commands::Tblastx(args) => args.parse_deflines,
            Commands::Blastx(args) => args.parse_deflines,
            Commands::Psiblast(_) | Commands::Makedb(_) => false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:67-72
    // ```c++
    //     // Set an appropriate default for the strand
    //     if (m_Strand == eNa_strand_other) {
    //         m_Strand = (m_DLConfig.m_IsLoadingProteins)
    //             ? eNa_strand_unknown
    //             : eNa_strand_both;
    //     }
    // ```
    /// The `-strand` value of the search's argv, for the programs that have the option.
    pub(crate) fn strand(&self) -> Option<String> {
        match self.any_subject() {
            Commands::Blastn(args) => Some(args.strand),
            Commands::Blastx(args) => Some(args.strand),
            _ => None,
        }
    }

    /// The parsed arguments of the search of `inputs`: the command line's `-query` and
    /// `-subject` or `-db` values. An option can hold for one kind of subject only (`-subject_loc`
    /// excludes `-db`, `-seqidlist` excludes `-subject`), so this can fail with the parser's
    /// message where `validate`, which knows no inputs, could not.
    pub(crate) fn command(&self, inputs: Inputs<'_>) -> Result<Commands, Error> {
        let subject_flag = if inputs.database { DATABASE } else { SUBJECT };
        parse(
            self.program,
            &self.argv,
            inputs.query_name,
            subject_flag,
            inputs.subject_name,
        )
    }
}

/// The flags of the two kinds of subject.
const SUBJECT: &str = "-subject";
const DATABASE: &str = "-db";

/// The key of an argv word: the option name, without a value joined by `=`.
fn key(word: &str) -> &str {
    word.split('=').next().unwrap_or("")
}

/// Parses the argv with the command line's parser, as `losat <program> -query <query>
/// <subject_flag> <subject> <argv...>`.
fn parse(
    program: Program,
    argv: &[String],
    query: &str,
    subject_flag: &str,
    subject: &str,
) -> Result<Commands, Error> {
    let words = [
        "losat",
        program.name(),
        "-query",
        query,
        subject_flag,
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
/// `-out` and `-outfmt` are refused (`API_OWNED`). BLASTX's in-memory entry does not take
/// `-parse_deflines` yet. The checks of the options that the command line makes
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
    }
    let argv: Vec<String> = argv.iter().map(|word| word.to_string()).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2377-2381
    // ```c++
    //         ITERATE(vector<string>, dbarg, database_args) {
    //             arg_desc.SetDependency(kArgSubjectLocation,
    //                                    CArgDescriptions::eExcludes,
    //                                    *dbarg);
    //         }
    // ```
    // An argv is checked before the inputs are known: it holds if it holds for one kind of
    // subject, and the search parses it again for the kind it has (`Search::command`).
    if let Err(error) = parse(program, &argv, "query", SUBJECT, "subject") {
        parse(program, &argv, "query", DATABASE, "subject").map_err(|_| error)?;
    }
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
        let inputs = Inputs::new(&[], &[]).named("q.faa", "s.faa");
        match search.command(inputs).unwrap() {
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
        for program in [Program::Blastp, Program::Blastx] {
            assert!(validate(program, &["-parse_deflines"])
                .unwrap()
                .parse_deflines());
        }
        // A value that looks like an option is a value.
        validate(Program::Blastp, &["-matrix", "BLOSUM62"]).unwrap();
    }
}
