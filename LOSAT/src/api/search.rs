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

// No NCBI counterpart: NCBI has no library API. LOSAT's library runs psiblast one query per call
// (`api::batches`, `psiblast_app.cpp:468`), so the numbers of `-save_each_pssm`'s files (`F.1`,
// `F.2`, ...; `CAutoOutputFileReset`) would start again at each call, and a PSSM file is an output
// outside those the API owns (`API_OWNED`): a file the caller did not give. The command line
// writes them (`algorithm::psiblast::pssm`).
/// psiblast's options that write a PSSM file, which the library API refuses, and `-in_pssm`, the
/// query as a PSSM file (its queries are records; line PS, DECISIONS 11, as 8 and 10); and its
/// search strategy file (`STRATEGY_FILES`, line RG item F1b-1).
const PSIBLAST_FILES: &[&str] = &[
    "out_pssm",
    "out_ascii_pssm",
    "save_pssm_after_last_round",
    "save_each_pssm",
    "in_pssm",
    "export_search_strategy",
];

// No NCBI counterpart: the search strategy file (`-export_search_strategy`,
// `report::archive::strategy`) is an output outside those the API owns too: a file the caller
// did not give. The command line writes it.
/// The options that write the search strategy file (line RG-F1).
const STRATEGY_FILES: &[&str] = &["export_search_strategy"];

/// tblastn's `-in_pssm` (PSI-TBLASTN, line PS item PS7), refused for the same reason as psiblast's,
/// and its search strategy file.
const TBLASTN_FILES: &[&str] = &["in_pssm", "export_search_strategy"];

/// The options of `program`'s command line that write a file the library API does not write (a
/// PSSM, `PSIBLAST_FILES`, or the search strategy, `STRATEGY_FILES`) or read the query from one: they
/// are not in an argv that `validate` accepts, and `describe` leaves them out.
/// The searches of records of LOSAT Web take none of them either (the browser refuses `-in_pssm`
/// until line PS's PSW1), and its `describe` leaves them out through the same list
/// (`web_adapter::records_refused_options`).
pub(crate) fn file_options(program: Program) -> &'static [&'static str] {
    match program {
        Program::Psiblast => PSIBLAST_FILES,
        Program::Tblastn => TBLASTN_FILES,
        Program::Blastn | Program::Blastp | Program::Blastx | Program::Tblastx => STRATEGY_FILES,
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
        parse_words(
            self.program,
            &self.argv,
            ["query", SUBJECT, "subject"],
            None,
        )
        .or_else(|_| {
            parse_words(
                self.program,
                &self.argv,
                ["query", DATABASE, "subject"],
                None,
            )
        })
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

    /// The parsed arguments of the search of `inputs`: the command line's `-query` and
    /// `-subject` or `-db` values, and `-out` for `run_to_file` (`Search::writing_to`). An option
    /// can hold for one kind of subject only (`-subject_loc` excludes `-db`, `-seqidlist` excludes
    /// `-subject`), so this can fail with the parser's message where `validate`, which knows no
    /// inputs, could not.
    pub(crate) fn command(&self, inputs: Inputs<'_>) -> Result<Commands, Error> {
        let subject_flag = if inputs.database { DATABASE } else { SUBJECT };
        parse(
            self.program,
            &self.argv,
            [inputs.query_name, subject_flag, inputs.subject_name],
            self.out.as_deref(),
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
    inputs: [&str; 3],
    out: Option<&std::path::Path>,
) -> Result<Commands, Error> {
    for word in argv {
        let key = key(word);
        if let Some(name) = option_name(key).filter(|name| file_options(program).contains(name)) {
            if name == "in_pssm" {
                return Err(Error::Argv(format!(
                    "{key} is not accepted by the library API: the queries are given as records"
                )));
            }
            if STRATEGY_FILES.contains(&name) {
                return Err(Error::Argv(format!(
                    "{key} is not accepted by the library API: the library API writes no search strategy file"
                )));
            }
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: the library API writes no PSSM file"
            )));
        }
    }
    let command = parse_words(program, argv, inputs, out)?;
    // No NCBI counterpart: NCBI has no library API. psiblast's rounds (`-num_iterations` other
    // than 1) run on the command line (`algorithm::psiblast::rounds`); the library's results
    // (`QueryHits`, the hit records and the observer's events) have no round, so the library
    // runs one round only (line PS, DECISIONS 15).
    if let Commands::Psiblast(args) = &command {
        if args.num_iterations != 1 {
            return Err(Error::Argv(format!(
                "-num_iterations {} is not accepted by the library API: its results have one round (the command line runs psiblast's rounds)",
                args.num_iterations
            )));
        }
    }
    Ok(command)
}

/// Parses the argv with the command line's parser, as `losat <program> -query <query>
/// <subject_flag> <subject> [-out <out>] <argv...>` (`inputs`: query, subject flag, subject).
fn parse_words(
    program: Program,
    argv: &[String],
    [query, subject_flag, subject]: [&str; 3],
    out: Option<&std::path::Path>,
) -> Result<Commands, Error> {
    use std::ffi::OsStr;
    let words = [
        "losat",
        program.name(),
        "-query",
        query,
        subject_flag,
        subject,
    ];
    let out = out.map(|path| [OsStr::new("-out"), path.as_os_str()]);
    let cli: Cli = crate::cli::try_parse_from_library(
        words
            .iter()
            .map(OsStr::new)
            .chain(out.iter().flatten().copied())
            .chain(argv.iter().map(OsStr::new)),
    )
    .map_err(|error| Error::Argv(crate::cli::render_message(&error)))?;
    Ok(cli.command)
}

/// Checks an argv for `program`: the options after the program name, as on the command line
/// (`-evalue 1e-5 -num_threads 4`). Every word goes to the command line's parser, so an
/// option means, and fails, exactly as there; the error is the parser's message.
///
/// The library API sets the inputs and the outputs itself, so `-query`, `-subject`, `-db`,
/// `-out` and `-outfmt` are refused (`API_OWNED`). It writes no PSSM file, so psiblast's
/// `-out_pssm`, `-out_ascii_pssm`, `-save_pssm_after_last_round` and `-save_each_pssm` are refused too, here
/// and by every search, and no search strategy file (`-export_search_strategy`). The checks of the options that the command line makes before it
/// searches (ranges, scoring, the hit list size of each format) come with `run` and
/// `run_collect`, which know the formats and the inputs, and report them with the command line's
/// messages.
pub fn validate(program: Program, argv: &[&str]) -> Result<Search, Error> {
    for word in argv {
        let key = key(word);
        let Some(name) = option_name(key) else {
            continue;
        };
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
    if let Err(error) = parse(program, &argv, ["query", SUBJECT, "subject"], None) {
        parse(program, &argv, ["query", DATABASE, "subject"], None).map_err(|_| error)?;
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

    /// The refusal of psiblast's PSSM files, with the word as the caller wrote it.
    fn writes_no_pssm_file(key: &str) -> Error {
        Error::Argv(format!(
            "{key} is not accepted by the library API: the library API writes no PSSM file"
        ))
    }

    #[test]
    fn the_search_strategy_file_is_refused_by_validate() {
        for program in [
            Program::Blastn,
            Program::Blastp,
            Program::Blastx,
            Program::Tblastn,
            Program::Tblastx,
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

    #[test]
    fn psiblast_pssm_files_are_refused_by_validate() {
        for argv in [
            &["-out_ascii_pssm", "p.txt"][..],
            &["-out_ascii_pssm=p.txt"],
            &["-out_pssm", "p.asn"],
            &["-evalue", "1", "-save_pssm_after_last_round"],
            &["-save_each_pssm"],
        ] {
            // The command line takes each of them.
            let words = ["losat", "psiblast", "-query", "q", "-subject", "s"];
            crate::cli::try_parse_from::<Cli, _, _>(words.iter().chain(argv)).unwrap();
            let refused = argv.iter().find(|word| word.contains("pssm")).unwrap();
            assert_eq!(
                validate(Program::Psiblast, argv),
                Err(writes_no_pssm_file(key(refused))),
                "{argv:?}"
            );
        }
        // The first of them is named.
        assert_eq!(
            validate(
                Program::Psiblast,
                &[
                    "-save_each_pssm",
                    "-out_ascii_pssm",
                    "-save_pssm_after_last_round"
                ]
            ),
            Err(writes_no_pssm_file("-save_each_pssm"))
        );
        // As the options the API owns, before the parser (which takes one dash only).
        assert_eq!(
            validate(Program::Psiblast, &["--save_each_pssm"]),
            Err(writes_no_pssm_file("--save_each_pssm"))
        );
        validate(Program::Psiblast, &["-pseudocount", "2"]).unwrap();
        // blastp has no such option: the parser's message.
        let error = validate(Program::Blastp, &["-save_each_pssm"]).unwrap_err();
        assert!(!error.message().contains("PSSM"), "{error}");
    }

    #[test]
    fn psiblast_pssm_files_are_refused_by_every_search_before_it_runs() {
        // `validate` makes every `Search`; this one is made without it, so the search's own
        // refusal (`Search::command`) is what stops it.
        let search = Search {
            program: Program::Psiblast,
            argv: vec![
                "-save_pssm_after_last_round".to_string(),
                "-out_ascii_pssm".to_string(),
                "p.txt".to_string(),
            ],
            out: None,
        };
        let refused = writes_no_pssm_file("-save_pssm_after_last_round");
        let query = crate::api::read_fasta(
            Program::Psiblast,
            Role::Query,
            b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n",
        )
        .unwrap();
        let inputs = Inputs::new(&query, &query).named("q.faa", "s.faa");
        assert_eq!(search.command(inputs).unwrap_err(), refused);
        assert_eq!(
            crate::api::run_collect(&search, inputs, &["6"]).unwrap_err(),
            refused
        );
        let mut report = Vec::new();
        let mut diagnostics = Vec::new();
        assert_eq!(
            crate::api::run(&search, inputs, &mut [("6", &mut report)], &mut diagnostics)
                .unwrap_err(),
            refused
        );
        assert!(report.is_empty() && diagnostics.is_empty());
        assert_eq!(
            crate::api::run_batches(&search, inputs).err(),
            Some(refused.clone())
        );
        let path =
            std::env::temp_dir().join(format!("losat-api-pssm-refusal-{}.out", std::process::id()));
        assert_eq!(
            crate::api::run_to_file(&search, inputs, "6", &path, &mut diagnostics).unwrap_err(),
            refused
        );
        assert!(!path.exists());
        assert!(!std::path::Path::new("p.txt").exists());
    }

    // Line PS, PS5a (DECISIONS 15): psiblast's rounds run on the command line only; the library's
    // results have one round, so `validate` and every search refuse `-num_iterations` other than 1.
    #[test]
    fn psiblast_rounds_are_refused_by_validate_and_every_search() {
        for rounds in ["0", "2"] {
            let error = validate(Program::Psiblast, &["-num_iterations", rounds]).unwrap_err();
            assert_eq!(
                error,
                Error::Argv(format!(
                    "-num_iterations {rounds} is not accepted by the library API: its results have one round (the command line runs psiblast's rounds)"
                ))
            );
        }
        validate(Program::Psiblast, &["-num_iterations", "1"]).unwrap();
        let search = Search {
            program: Program::Psiblast,
            argv: vec!["-num_iterations".to_string(), "3".to_string()],
            out: None,
        };
        let query = crate::api::read_fasta(
            Program::Psiblast,
            Role::Query,
            b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n",
        )
        .unwrap();
        let inputs = Inputs::new(&query, &query).named("q.faa", "s.faa");
        assert!(
            matches!(search.command(inputs), Err(Error::Argv(message)) if message.starts_with("-num_iterations 3 "))
        );
    }
}
