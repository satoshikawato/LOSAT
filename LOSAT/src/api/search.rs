//! A search: a program and its argv, checked by the command line's parser.
//!
//! The command line is the one definition of the options (NCBI's `CBlastAppArgs`); the library
//! API takes the same words and parses them with the same parser, so an option means, and fails,
//! exactly as on the command line (LOSAT Web ABI v2 §7 does the same).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,51,75
#![warn(missing_docs, unused)]

use super::error::Error;
use super::input::{Inputs, QueryFile};
use super::program::Program;
use crate::blastinput::ncbi_environment::DiagTarget;
use crate::cli::{Cli, Commands};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
/// The options that the library API sets itself and does not take in an argv: the inputs
/// (`Inputs` gives the records and their names, a database's name, and the PSSM file of a PSSM
/// query, `-in_pssm` of psiblast and tblastn, and the search strategy file of
/// `-import_search_strategy`, `Search::read_strategy` and `Inputs::strategy`), the output file and
/// the output format (arguments of `run_collect`, `run` and `run_to_file`). psiblast's `-in_msa` is
/// owned the same way for that program only (`Inputs::msa`; `alignment_owned`), and so is not in
/// this list.
pub const API_OWNED: &[&str] = &[
    "query",
    "subject",
    "db",
    "in_pssm",
    "import_search_strategy",
    "out",
    "outfmt",
];

// No NCBI counterpart: the search strategy file (`-export_search_strategy`,
// `report::archive::strategy`) is an output outside those the API owns too: a file the caller
// did not give. The command line writes it.
/// The options that write the search strategy file (line RG-F1).
const STRATEGY_FILES: &[&str] = &["export_search_strategy"];

/// The options of `program`'s command line that write a file the library API does not write (the
/// search strategy, `STRATEGY_FILES`, line RG item F1b-1; the same for every program): they are not in an argv that
/// `validate` accepts, and `describe` leaves them out. psiblast's PSSM files (`-out_pssm`,
/// `-out_ascii_pssm`, `-save_pssm_after_last_round`, `-save_each_pssm`) are written where the
/// command line writes them (series Q item Q6a). The searches of records of LOSAT Web take none
/// of them either, and tblastn's `-in_pssm` neither until line PS's PSW1
/// (`web_adapter::records_refused_options`).
pub(crate) fn file_options(_program: Program) -> &'static [&'static str] {
    STRATEGY_FILES
}

/// The option of `program`'s command line that gives the query as an alignment (psiblast's
/// `-in_msa`): the library API gives the alignment as the inputs (`Inputs::msa`), so `validate`
/// refuses the option and `describe` leaves it out. Its `-msa_master_idx` and
/// `-ignore_msa_master` are argv options (`ALIGNMENT_ARGUMENTS`), as on the command line.
pub(crate) fn alignment_owned(program: Program) -> &'static [&'static str] {
    match program {
        Program::Psiblast => &["in_msa"],
        _ => &[],
    }
}

/// The options that choose the master of the query alignment: the command line takes them with
/// `-in_msa` only (`-msa_master_idx` and `-ignore_msa_master`, psiblast's), which the parser
/// excludes from `-query`. An argv that holds one is checked with the alignment as the query.
const ALIGNMENT_ARGUMENTS: &[&str] = &["msa_master_idx", "ignore_msa_master"];

/// A program and an argv that the command line's parser accepted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Search {
    program: Program,
    argv: Vec<String>,
    /// The `-out` of [`run_to_file`](super::run_to_file)'s searches (`Search::writing_to`).
    out: Option<std::path::PathBuf>,
}

impl Search {
    /// A search of `argv` that `validate` has not checked: the tests of a call's reading of an argv
    /// that only its start can make fail (`start_tests`: an argument error with `-logfile`, which
    /// `validate` reports as without one).
    #[cfg(test)]
    pub(crate) fn unchecked(program: Program, argv: &[&str]) -> Search {
        Search {
            program,
            argv: argv.iter().map(|word| word.to_string()).collect(),
            out: None,
        }
    }

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

    /// The file of `run_to_file` (`writing_to`).
    pub(crate) fn out_path(&self) -> Option<&std::path::Path> {
        self.out.as_deref()
    }

    /// The parsed arguments of the search's argv with an `-import_search_strategy` (the file
    /// itself is the engine's import's, `Search::read_strategy`) and neither `-query`
    /// nor `-subject`: the command line's import without overrides. The posts of an argument
    /// error go to `target` (`Search::command`).
    pub(crate) fn strategy_command(
        &self,
        out: Option<&std::path::Path>,
        target: &DiagTarget,
    ) -> Result<Commands, Error> {
        let argv: Vec<String> = self
            .argv
            .iter()
            .cloned()
            .chain([
                "-import_search_strategy".to_string(),
                "strategy".to_string(),
            ])
            .collect();
        let query = (QUERY, std::ffi::OsStr::new("-"));
        parse_words(self.program, &argv, query, None, out, target)
    }

    /// The parsed arguments, with either kind of subject: `validate` accepted one of them. The
    /// reader's options come from here before the search refuses what the API does not take
    /// (`parse`).
    pub(super) fn any_subject(&self) -> Commands {
        parse_for_some_inputs(self.program, &self.argv).expect("validated argv")
    }

    /// Whether the search's argv has `-parse_deflines` (`read_inputs`).
    pub(crate) fn parse_deflines(&self) -> bool {
        super::input::parse_deflines_of(&self.any_subject())
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
    // (cited at `ncbi_app::constraint::Constraint::MaxFileNameLength`)
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1241-1248
    /// The parsed arguments of the search of `inputs`: the command line's `-query` (or `-in_pssm`
    /// for a PSSM query, `Inputs::pssm`, or `-in_msa` for an alignment, `Inputs::msa`) and
    /// `-subject` or `-db` values, and `-out` for `run_to_file` (`Search::writing_to`). An option
    /// can hold for one kind of input only (`-subject_loc` excludes `-db`, `-seqidlist` excludes
    /// `-subject`, `-query_loc` excludes `-in_pssm` and `-in_msa`, `-msa_master_idx` needs
    /// `-in_msa`), so this can fail with the parser's message where `validate`, which knows no
    /// inputs, could not. A PSSM query is psiblast's and tblastn's only, an alignment query
    /// psiblast's, and each a file (not `-`). The posts of an argument error go to `target`, the
    /// call's `-logfile` target (`CallStart::new`), as the command line's go to its own.
    pub(crate) fn command(
        &self,
        inputs: Inputs<'_>,
        target: &DiagTarget,
    ) -> Result<Commands, Error> {
        // A strategy's arguments are the import's (`Search::read_strategy`), set on the program's
        // defaults; the search's `-out` is the only word of its own that they take.
        if let Some(strategy) = inputs.strategy {
            if self.out.is_some() {
                // The parser checks the name of `-out` as for any search.
                self.strategy_command(self.out.as_deref(), target)?;
            }
            return strategy.command(self);
        }
        let subject_flag = if inputs.database { DATABASE } else { SUBJECT };
        let query = match inputs.query_file {
            Some(file) => {
                let (flag, programs, what) = match file {
                    QueryFile::Pssm(_) => {
                        (PSSM, "psiblast and tblastn", "a PSSM query (Inputs::pssm)")
                    }
                    QueryFile::Alignment(_) => (MSA, "psiblast", "a query alignment (Inputs::msa)"),
                };
                let takes = match file {
                    QueryFile::Pssm(_) => {
                        matches!(self.program, Program::Psiblast | Program::Tblastn)
                    }
                    QueryFile::Alignment(_) => self.program == Program::Psiblast,
                };
                if !takes {
                    return Err(Error::Argv(
                        format!(
                            "{what} is searched by {programs}, not by {}",
                            self.program.name()
                        )
                        .into(),
                    ));
                }
                // `-in_pssm -` and `-in_msa -` read the command line's standard input, which the
                // library does not read (`run_to_file` does not write it either).
                if file.path().as_os_str() == "-" {
                    return Err(Error::Argv(
                        format!(
                            "Inputs::{} reads a file: `-` is the command line's standard input",
                            match file {
                                QueryFile::Pssm(_) => "pssm",
                                QueryFile::Alignment(_) => "msa",
                            }
                        )
                        .into(),
                    ));
                }
                (flag, file.path().as_os_str())
            }
            None => (QUERY, std::ffi::OsStr::new(inputs.query_name)),
        };
        parse(
            self.program,
            &self.argv,
            query,
            [subject_flag, inputs.subject_name],
            self.out.as_deref(),
            target,
        )
    }
}

/// The flags of the two kinds of query and of subject.
const QUERY: &str = "-query";
const PSSM: &str = "-in_pssm";
const MSA: &str = "-in_msa";
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

/// The parsed arguments of an argv that `validate` checks or a search reads before it knows its
/// inputs: they hold if they hold for one kind of query (`-query`; psiblast's `-in_msa`, which
/// the options of the alignment need) and one kind of subject (`-subject` or `-db`). The error is
/// that of the kind the argv's options name (the alignment's, when it has `-msa_master_idx` or
/// `-ignore_msa_master`), else that of the records.
fn parse_for_some_inputs(program: Program, argv: &[String]) -> Result<Commands, Error> {
    let aligned = argv
        .iter()
        .any(|word| option_name(key(word)).is_some_and(|name| ALIGNMENT_ARGUMENTS.contains(&name)));
    let records = (QUERY, std::ffi::OsStr::new("query"));
    let alignment = (MSA, std::ffi::OsStr::new("alignment"));
    let queries = match (program, aligned) {
        (Program::Psiblast, true) => vec![alignment, records],
        (Program::Psiblast, false) => vec![records, alignment],
        _ => vec![records],
    };
    let mut first_error = None;
    for query in queries {
        for subject_flag in [SUBJECT, DATABASE] {
            match parse(
                program,
                argv,
                query,
                [subject_flag, "subject"],
                None,
                &DiagTarget::Stderr,
            ) {
                Ok(commands) => return Ok(commands),
                Err(error) => {
                    first_error.get_or_insert(error);
                }
            }
        }
    }
    Err(first_error.expect("at least one kind of inputs"))
}

/// `parse_words` after the refusal of the options that the library API does not take
/// (`file_options`). `validate` and every search (`Search::command`) come here.
fn parse(
    program: Program,
    argv: &[String],
    query: (&str, &std::ffi::OsStr),
    subject: [&str; 2],
    out: Option<&std::path::Path>,
    target: &DiagTarget,
) -> Result<Commands, Error> {
    for word in argv {
        let key = key(word);
        if let Some(name) = option_name(key).filter(|name| file_options(program).contains(name)) {
            debug_assert!(STRATEGY_FILES.contains(&name));
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: the library API writes no search strategy file"
            ).into()));
        }
    }
    parse_words(program, argv, query, Some(subject), out, target)
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
    subject: Option<[&str; 2]>,
    out: Option<&std::path::Path>,
    target: &DiagTarget,
) -> Result<Commands, Error> {
    let words = [
        std::ffi::OsStr::new("losat"),
        std::ffi::OsStr::new(program.name()),
        std::ffi::OsStr::new(query_flag),
        query,
    ];
    let subject = subject
        .into_iter()
        .flat_map(|[flag, name]| [std::ffi::OsStr::new(flag), std::ffi::OsStr::new(name)]);
    let out = out.map(|path| {
        let mut word = std::ffi::OsString::from("-out=");
        word.push(path);
        word
    });
    let cli: Cli = crate::cli::try_parse_with_own_start(
        words
            .iter()
            .copied()
            .chain(subject)
            .map(|word| word.to_os_string())
            .chain(out)
            .chain(argv.iter().map(std::ffi::OsString::from)),
        target,
    )
    .map_err(Error::of_argv)?;
    Ok(cli.command)
}

/// Checks an argv for `program`: the options after the program name, as on the command line
/// (`-evalue 1e-5 -num_threads 4`). Every word goes to the command line's parser, so an
/// option means, and fails, exactly as there; the error is the parser's message.
///
/// The library API sets the inputs and the outputs itself, so `-query`, `-subject`, `-db`,
/// `-in_pssm`, `-in_msa`, `-out` and `-outfmt` are refused (`API_OWNED`; a PSSM query of psiblast
/// and tblastn is given as `Inputs::pssm`, psiblast's alignment query as `Inputs::msa`, whose
/// `-msa_master_idx` and `-ignore_msa_master` are options of the argv). It writes no search strategy file
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
        // `-in_pssm` is psiblast's and tblastn's input (`Inputs::pssm`); the other programs' parser
        // knows no such option.
        if name == "in_pssm" {
            if !matches!(program, Program::Psiblast | Program::Tblastn) {
                continue;
            }
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: a PSSM query is given as the inputs (Inputs::pssm)"
            ).into()));
        }
        // `-in_msa` is the command line's alternative to a query (psiblast's alignment, series PS
        // item PS6), given as the inputs (`Inputs::msa`); the other programs' parser knows no
        // such option.
        if alignment_owned(program).contains(&name) {
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: a query alignment is given as the inputs (Inputs::msa)"
            ).into()));
        }
        if name == "import_search_strategy" {
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: a search strategy is given as bytes (Search::read_strategy, Inputs::strategy)"
            ).into()));
        }
        if API_OWNED.contains(&name) {
            return Err(Error::Argv(format!(
                "{key} is not accepted by the library API: the inputs are given as records and the outputs as formats"
            ).into()));
        }
    }
    let argv: Vec<String> = argv.iter().map(|word| word.to_string()).collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2377-2381
    // An argv is checked before the inputs are known: it holds if it holds for one kind of
    // query and subject, and the search parses it again for the kind it has (`Search::command`).
    parse_for_some_inputs(program, &argv)?;
    Ok(Search {
        program,
        argv,
        out: None,
    })
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
