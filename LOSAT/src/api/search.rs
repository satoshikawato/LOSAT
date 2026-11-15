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
/// (`Inputs` gives the records and their names, a database's name, and the PSSM file of a PSSM
/// query, `-in_pssm` of psiblast and tblastn), the output file and the output format (arguments
/// of `run_collect`, `run` and `run_to_file`).
pub const API_OWNED: &[&str] = &["query", "subject", "db", "in_pssm", "out", "outfmt"];

// No NCBI counterpart: the search strategy file (`-export_search_strategy`,
// `report::archive::strategy`) is an output outside those the API owns too: a file the caller
// did not give. The command line writes it.
/// The options that write the search strategy file (line RG-F1).
const STRATEGY_FILES: &[&str] = &["export_search_strategy"];

/// The options of `program`'s command line that write a file the library API does not write (the
/// search strategy, `STRATEGY_FILES`, line RG item F1b-1): they are not in an argv that
/// `validate` accepts, and `describe` leaves them out. psiblast's PSSM files (`-out_pssm`,
/// `-out_ascii_pssm`, `-save_pssm_after_last_round`, `-save_each_pssm`) are written where the
/// command line writes them (series Q item Q6a). The searches of records of LOSAT Web take none
/// of them either, and tblastn's `-in_pssm` neither until line PS's PSW1
/// (`web_adapter::records_refused_options`).
pub(crate) fn file_options(program: Program) -> &'static [&'static str] {
    match program {
        Program::Blastn
        | Program::Blastp
        | Program::Blastx
        | Program::Tblastx
        | Program::Psiblast
        | Program::Tblastn => STRATEGY_FILES,
    }
}

/// A program and an argv that the command line's parser accepted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Search {
    program: Program,
    argv: Vec<String>,
    /// The `-out` of [`run_to_file`](super::run_to_file)'s searches (`Search::writing_to`).
    out: Option<std::path::PathBuf>,
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

    /// Whether the argv has `-dryrun`: every search of it runs the NCBI application start (it
    /// opens the file of `-logfile`, loads the file of `-conffile`) and checks the settings, then
    /// ends without a search, as the command line's dry run does. [`run`](super::run()) writes
    /// nothing, [`run_collect`](super::run_collect) gives an empty report per format,
    /// [`run_to_file`](super::run_to_file) makes no file and [`run_batches`](super::run_batches)
    /// one batch of no query. The records are not read then: the inputs may be empty.
    pub fn dry_run(&self) -> bool {
        use crate::blastinput::ncbi_environment::{preparse, Preparsed};
        matches!(preparse(&self.argv), Preparsed::Words { dryrun: true, .. })
    }

    /// Whether a search of this argv runs the NCBI application start at its beginning: the argv
    /// has `-logfile`, `-conffile` or `-dryrun`, which the start reads (`api::start`).
    pub(crate) fn starts(&self) -> bool {
        use crate::blastinput::ncbi_environment::{log_file, preparse, Preparsed};
        log_file(&self.argv).is_some()
            || matches!(
                preparse(&self.argv),
                Preparsed::Words {
                    conffile: Some(_),
                    ..
                } | Preparsed::Words { dryrun: true, .. }
            )
    }

    /// The parsed arguments, with either kind of subject: `validate` accepted one of them. The
    /// reader's options come from here before the search refuses what the API does not take
    /// (`parse`).
    fn any_subject(&self) -> Commands {
        let query = (QUERY, std::ffi::OsStr::new("query"));
        parse_words(self.program, &self.argv, query, [SUBJECT, "subject"], None)
            .or_else(|_| parse_words(self.program, &self.argv, query, [DATABASE, "subject"], None))
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
            Commands::Psiblast(args) => args.parse_deflines,
            Commands::Rpsblast(_)
            | Commands::Rpstblastn(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => false,
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

    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
    // (quoted at `blastinput::value_parsers::blastn_output_path`)
    /// The search of `run_to_file`, whose `-out` is `path`: the parser checks it as the command
    /// line's `-out` (a file name of 256 bytes or more is NCBI's `Illegal value`, dry run or not),
    /// and the arguments carry it to the search, which makes the file where the command line makes
    /// it (`ReportOutputs::open_out`).
    pub(crate) fn writing_to(&self, path: &std::path::Path) -> Search {
        Search {
            out: Some(path.to_path_buf()),
            ..self.clone()
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1240-1247
    // ```c++
    //     if (!m_IsDeltaBlast) {
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQuery);
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQueryLocation);
    //     }
    // ```
    /// The parsed arguments of the search of `inputs`: the command line's `-query` (or `-in_pssm`
    /// for a PSSM query, `Inputs::pssm`) and `-subject` or `-db` values, and `-out` for
    /// `run_to_file` (`Search::writing_to`). An option can hold for one kind of input only
    /// (`-subject_loc` excludes `-db`, `-seqidlist` excludes `-subject`, `-query_loc` excludes
    /// `-in_pssm`), so this can fail with the parser's message where `validate`, which knows no
    /// inputs, could not. A PSSM query is psiblast's and tblastn's only, and a file (not `-`).
    pub(crate) fn command(&self, inputs: Inputs<'_>) -> Result<Commands, Error> {
        let subject_flag = if inputs.database { DATABASE } else { SUBJECT };
        let query = match inputs.pssm {
            Some(path) => {
                if !matches!(self.program, Program::Psiblast | Program::Tblastn) {
                    return Err(Error::Argv(format!(
                        "a PSSM query (Inputs::pssm) is searched by psiblast and tblastn, not by {}",
                        self.program.name()
                    )));
                }
                // `-in_pssm -` reads the command line's standard input, which the library does
                // not read (`run_to_file` does not write it either).
                if path.as_os_str() == "-" {
                    return Err(Error::Argv(
                        "Inputs::pssm reads a file: `-` is the command line's standard input"
                            .to_string(),
                    ));
                }
                (PSSM, path.as_os_str())
            }
            None => (QUERY, std::ffi::OsStr::new(inputs.query_name)),
        };
        parse(
            self.program,
            &self.argv,
            query,
            [subject_flag, inputs.subject_name],
            self.out.as_deref(),
        )
    }
}

/// The flags of the two kinds of query and of subject.
const QUERY: &str = "-query";
const PSSM: &str = "-in_pssm";
const SUBJECT: &str = "-subject";
const DATABASE: &str = "-db";

/// The key of an argv word: the option name, without a value joined by `=`.
fn key(word: &str) -> &str {
    word.split('=').next().unwrap_or("")
}

/// The option that an argv key names, without its one or two dashes; `None` for a value.
fn option_name(key: &str) -> Option<&str> {
    let name = key.strip_prefix('-')?;
    Some(name.strip_prefix('-').unwrap_or(name))
}

/// `parse_words` after the refusal of the options that the library API does not take
/// (`file_options`). `validate` and every search (`Search::command`) come here.
fn parse(
    program: Program,
    argv: &[String],
    query: (&str, &std::ffi::OsStr),
    subject: [&str; 2],
    out: Option<&std::path::Path>,
) -> Result<Commands, Error> {
    for word in argv {
        let key = key(word);
        if let Some(name) = option_name(key).filter(|name| file_options(program).contains(name)) {
            debug_assert!(STRATEGY_FILES.contains(&name));
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: the library API writes no search strategy file"
            )));
        }
    }
    parse_words(program, argv, query, subject, out)
}

/// Parses the argv with the command line's parser, as `losat <program> <query_flag> <query>
/// <subject_flag> <subject> [-out=<out>] <argv...>` (`query`: `-query` and the query's name, or
/// `-in_pssm` and the PSSM file; `subject`: the subject flag and the subject). The `-out` value is
/// one word with its option, so that NCBI's reading of its toolkit's words
/// (`ncbi_environment::preparse`: `-dryrun`, the versions, `--`) takes no file name for one; the
/// parser keeps the bytes of a name that is not UTF-8, the PSSM file's too.
fn parse_words(
    program: Program,
    argv: &[String],
    (query_flag, query): (&str, &std::ffi::OsStr),
    [subject_flag, subject]: [&str; 2],
    out: Option<&std::path::Path>,
) -> Result<Commands, Error> {
    let words = [
        std::ffi::OsStr::new("losat"),
        std::ffi::OsStr::new(program.name()),
        std::ffi::OsStr::new(query_flag),
        query,
        std::ffi::OsStr::new(subject_flag),
        std::ffi::OsStr::new(subject),
    ];
    let out = out.map(|path| {
        let mut word = std::ffi::OsString::from("-out=");
        word.push(path);
        word
    });
    let cli: Cli = crate::cli::try_parse_with_own_start(
        words
            .iter()
            .map(|word| word.to_os_string())
            .chain(out)
            .chain(argv.iter().map(std::ffi::OsString::from)),
    )
    .map_err(|error| Error::Argv(crate::cli::render_message(&error)))?;
    Ok(cli.command)
}

/// Checks an argv for `program`: the options after the program name, as on the command line
/// (`-evalue 1e-5 -num_threads 4`). Every word goes to the command line's parser, so an
/// option means, and fails, exactly as there; the error is the parser's message.
///
/// The library API sets the inputs and the outputs itself, so `-query`, `-subject`, `-db`,
/// `-in_pssm`, `-out` and `-outfmt` are refused (`API_OWNED`; a PSSM query of psiblast and
/// tblastn is given as `Inputs::pssm`). It writes no search strategy file
/// (`-export_search_strategy` is refused, here and by every search). psiblast's PSSM files (`-out_pssm`,
/// `-out_ascii_pssm`, with `-save_pssm_after_last_round` and `-save_each_pssm`) are written by
/// every search where the command line writes them, with its names. The checks of the options that the command line makes before it
/// searches (ranges, scoring, the hit list size of each format) come with `run` and
/// `run_collect`, which know the formats and the inputs, and report them with the command line's
/// messages.
pub fn validate(program: Program, argv: &[&str]) -> Result<Search, Error> {
    for word in argv {
        let key = key(word);
        let Some(name) = option_name(key) else {
            continue;
        };
        if name == "in_pssm" {
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: a PSSM query is given as the inputs (Inputs::pssm)"
            )));
        }
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
    let query = (QUERY, std::ffi::OsStr::new("query"));
    if let Err(error) = parse(program, &argv, query, [SUBJECT, "subject"], None) {
        parse(program, &argv, query, [DATABASE, "subject"], None).map_err(|_| error)?;
    }
    Ok(Search {
        program,
        argv,
        out: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::Role;

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
            "-in_pssm",
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
        for program in [Program::Blastp, Program::Blastx, Program::Psiblast] {
            assert!(validate(program, &["-parse_deflines"])
                .unwrap()
                .parse_deflines());
        }
        // A value that looks like an option is a value.
        validate(Program::Blastp, &["-matrix", "BLOSUM62"]).unwrap();
    }

    #[test]
    fn the_search_strategy_file_is_refused_by_validate() {
        for program in [
            Program::Blastn,
            Program::Blastp,
            Program::Blastx,
            Program::Tblastn,
            Program::Tblastx,
            Program::Psiblast,
        ] {
            for argv in [
                &["-export_search_strategy", "s.asn"][..],
                &["-export_search_strategy=s.asn"],
                &["-evalue", "1", "--export_search_strategy", "s.asn"],
            ] {
                let key = argv.iter().find(|word| word.contains("export")).unwrap();
                assert_eq!(
                    validate(program, argv),
                    Err(Error::Argv(format!(
                        "{} is not accepted by the library API: the library API writes no search strategy file",
                        super::key(key)
                    ))),
                    "{program} {argv:?}"
                );
            }
        }
    }

    // Series Q item Q6a: psiblast's PSSM files go to the parser as on the command line, and the
    // search writes them (`tests/library_api.rs`).
    #[test]
    fn psiblast_pssm_files_are_the_command_lines_options() {
        for argv in [
            &["-out_ascii_pssm", "p.txt"][..],
            &["-out_ascii_pssm=p.txt"],
            &["-out_pssm", "p.asn"],
            &[
                "-evalue",
                "1",
                "-save_pssm_after_last_round",
                "-out_pssm",
                "p.asn",
            ],
            &["-save_each_pssm", "-out_ascii_pssm", "p.txt"],
        ] {
            let search = validate(Program::Psiblast, argv).unwrap();
            let query =
                crate::api::read_fasta(Program::Psiblast, Role::Query, b">q\nMKT\n").unwrap();
            match search.command(Inputs::new(&query, &query)).unwrap() {
                Commands::Psiblast(args) => assert!(
                    args.out_pssm.is_some() || args.out_ascii_pssm.is_some(),
                    "{argv:?}"
                ),
                other => panic!("{other:?}"),
            }
        }
        // NCBI's check of the options comes with the search, as on the command line
        // (`PsiblastPssmOptions::check`).
        validate(Program::Psiblast, &["-save_each_pssm"]).unwrap();
        // The parser takes one dash only.
        let error = validate(Program::Psiblast, &["--save_each_pssm"]).unwrap_err();
        assert!(!error.message().contains("library API"), "{error}");
        // blastp has no such option: the parser's message.
        let error = validate(Program::Blastp, &["-save_each_pssm"]).unwrap_err();
        assert!(!error.message().contains("PSSM"), "{error}");
    }

    // Series Q item Q6b: a PSSM query is an input (`Inputs::pssm`), which the search gives the
    // parser as `-in_pssm` in place of `-query`; in an argv the option is the API's.
    #[test]
    fn a_pssm_query_is_an_input_that_the_search_gives_the_parser() {
        let refused = |key: &str| {
            Error::Argv(format!(
                "{key} is not accepted by the library API: a PSSM query is given as the inputs (Inputs::pssm)"
            ))
        };
        for program in [Program::Psiblast, Program::Tblastn, Program::Blastp] {
            assert_eq!(
                validate(program, &["-in_pssm", "p.asn"]),
                Err(refused("-in_pssm"))
            );
            assert_eq!(
                validate(program, &["--in_pssm=p.asn"]),
                Err(refused("--in_pssm"))
            );
        }
        let path = std::path::Path::new("p.asn");
        for program in [Program::Psiblast, Program::Tblastn] {
            let search = validate(program, &["-evalue", "1"]).unwrap();
            for inputs in [Inputs::pssm(path, &[]), Inputs::pssm_database(path, "db")] {
                assert_eq!(inputs.query_count(), 1);
                let in_pssm = match search.command(inputs.named("q.faa", "s.faa")).unwrap() {
                    Commands::Psiblast(args) => args.in_pssm,
                    Commands::Tblastn(args) => args.in_pssm,
                    other => panic!("{other:?}"),
                };
                assert_eq!(in_pssm.as_deref(), Some(path), "{program}");
            }
            // `-in_pssm` excludes `-query_loc`: the parser's error of the search, as `-db` and
            // `-subject_loc`.
            let ranged = validate(program, &["-query_loc", "1-10"]).unwrap();
            let error = ranged.command(Inputs::pssm(path, &[])).unwrap_err();
            assert!(
                matches!(&error, Error::Argv(message) if message.contains("cannot be used with")),
                "{program}: {error}"
            );
        }
        let search = validate(Program::Psiblast, &[]).unwrap();
        let standard_input = std::path::Path::new("-");
        assert!(matches!(
            search.command(Inputs::pssm(standard_input, &[])),
            Err(Error::Argv(message)) if message.contains("standard input")
        ));
        let search = validate(Program::Blastp, &[]).unwrap();
        assert_eq!(
            search.command(Inputs::pssm(path, &[])).unwrap_err(),
            Error::Argv(
                "a PSSM query (Inputs::pssm) is searched by psiblast and tblastn, not by blastp"
                    .to_string()
            )
        );
    }

    // Line RG item Q5: psiblast's rounds (`-num_iterations` other than 1) run in the library too;
    // the records carry their round (`HspRecord::round`, `Outcome::rounds`).
    #[test]
    fn psiblast_rounds_are_accepted_by_validate_and_every_search() {
        let query = crate::api::read_fasta(
            Program::Psiblast,
            Role::Query,
            b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n",
        )
        .unwrap();
        let inputs = Inputs::new(&query, &query).named("q.faa", "s.faa");
        for (rounds, number) in [("0", 0), ("1", 1), ("3", 3)] {
            let search = validate(Program::Psiblast, &["-num_iterations", rounds]).unwrap();
            match search.command(inputs).unwrap() {
                Commands::Psiblast(args) => assert_eq!(args.num_iterations, number),
                other => panic!("{other:?}"),
            }
        }
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
    // (quoted at `blastinput::value_parsers::blastn_output_path`)
    /// `run_to_file`'s file is the parser's `-out` (`Search::writing_to`): any name, one that looks
    /// like a toolkit word (`-dryrun`, `--`) or is not UTF-8 among them, is the file and leaves the
    /// argv as it is; a name of 256 bytes or more is NCBI's error.
    #[test]
    fn the_file_of_run_to_file_is_the_parsers_out() {
        let search = validate(Program::Blastn, &["-ungapped"]).unwrap();
        let inputs = Inputs::new(&[], &[]).named("q.fna", "s.fna");
        let names = [
            std::ffi::OsString::from("-dryrun"),
            std::ffi::OsString::from("--"),
            std::ffi::OsString::from("-version"),
            std::ffi::OsString::from("dir/o.out"),
        ]
        .into_iter()
        // A name that is not UTF-8, which only Unix can make.
        .chain({
            #[cfg(unix)]
            let not_utf8 = Some(
                <std::ffi::OsStr as std::os::unix::ffi::OsStrExt>::from_bytes(b"o\xff.out")
                    .to_os_string(),
            );
            #[cfg(not(unix))]
            let not_utf8 = None;
            not_utf8
        });
        for name in names {
            let writing = search.writing_to(std::path::Path::new(&name));
            assert!(!writing.dry_run(), "{name:?}");
            match writing.command(inputs).unwrap() {
                Commands::Blastn(args) => {
                    assert_eq!(args.out.as_deref(), Some(std::path::Path::new(&name)));
                    assert!(args.ungapped, "{name:?}");
                }
                other => panic!("{other:?}"),
            }
        }
        let long = "a".repeat(256);
        let Err(Error::Argv(message)) = search
            .writing_to(std::path::Path::new(&long))
            .command(inputs)
        else {
            panic!("a name of 256 bytes");
        };
        assert!(
            message.contains("Illegal value, expected file name length < 256"),
            "{message}"
        );
    }
}
