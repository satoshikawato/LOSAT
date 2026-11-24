//! NCBI's argument reader (`ncbi_app::read`, line RG DECISIONS K4i-2) for the ten programs whose
//! arguments clap then converts into their types: blastn, blastp, blastx, tblastn, tblastx,
//! psiblast, rpsblast, rpstblastn, makeblastdb and rpsbproc. Each program is a `Tool` (its frozen
//! `-help` text; its argument XML is NCBI's), read in NCBI's order by every caller: the command
//! line's entry (`cli::main_entry`), the library API (`api::search`, `api::run_makeblastdb`) and
//! LOSAT Web's ABI v2 (`web_adapter`). NCBI's reader writes the help and the argument errors
//! (one producer of their bytes); then LOSAT refuses by name what it does not run
//! (`refusal_after_read`), and clap converts the arguments that the reader read
//! (`convert`): the argv that clap reads is made from them, so clap checks nothing that NCBI's
//! reader has not checked first.

use std::ffi::OsString;

use clap::FromArgMatches;

use super::{Cli, Commands};
use crate::blastinput::ncbi_environment::{preparse, DiagTarget, Preparsed};
use crate::ncbi_app::args::{is_file_argument, refusal, Args, Refusal};
use crate::ncbi_app::{Diag, ReadEnd, Tool};

/// A `Tool` of a program on clap: NCBI's frozen `-help` text, no diagnostic prefix, no
/// components.
macro_rules! tool {
    ($program:literal, $shown:literal, $help:expr) => {
        Tool {
            program: $program,
            refusal_name: $shown,
            help: $help,
            post_prefix: None,
            components: &[],
            help_full: None,
            short_description: None,
            write_failure: None,
            dry_run: None,
        }
    };
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:89-95 (each search program's
// `Init` describes its arguments; NCBI BLAST+ 2.17.0's `-help` of each is frozen)
static BLASTN: Tool = tool!(
    "blastn",
    "BLASTN",
    include_str!("../algorithm/blastn/help.txt")
);
static BLASTP: Tool = tool!(
    "blastp",
    "BLASTP",
    include_str!("../algorithm/blastp/help.txt")
);
static BLASTX: Tool = tool!(
    "blastx",
    "BLASTX",
    include_str!("../algorithm/blastx/help.txt")
);
static TBLASTN: Tool = tool!(
    "tblastn",
    "TBLASTN",
    include_str!("../algorithm/tblastn/help.txt")
);
static TBLASTX: Tool = tool!(
    "tblastx",
    "TBLASTX",
    include_str!("../algorithm/tblastx/help.txt")
);
static PSIBLAST: Tool = tool!(
    "psiblast",
    "PSIBLAST",
    include_str!("../algorithm/psiblast/help.txt")
);
static RPSBLAST: Tool = tool!(
    "rpsblast",
    "RPSBLAST",
    include_str!("../algorithm/rpsblast/help.txt")
);
static RPSTBLASTN: Tool = tool!(
    "rpstblastn",
    "RPSTBLASTN",
    crate::algorithm::rpstblastn::args::HELP
);

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:300-302
// `-randomize` is described `fHidden`: NCBI's `-help-full` lists it, `-help` does not (NCBI BLAST+
// 2.17.0, frozen).
static MAKEBLASTDB: Tool = Tool {
    help_full: Some(include_str!("../writedb/makeblastdb_help_full.txt")),
    ..tool!(
        "makeblastdb",
        "makeblastdb",
        include_str!("../writedb/makeblastdb_help.txt")
    )
};

// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:64-65
// rpsbproc's `Init` gives its usage context a short description, which `-h` writes, and a long one,
// which `-help` writes (RpsbProc 0.5.1, frozen). A failed write ends it with exit status 4
// (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 9).
static RPSBPROC: Tool = Tool {
    short_description: Some(
        "   Process output from offline rpsblast to annotate sequence with conserved\n   domain \
         information\n",
    ),
    write_failure: Some(crate::rpsbproc::E_OUTPUT_FILE_OPEN_ERROR),
    ..tool!("rpsbproc", "rpsbproc", include_str!("../rpsbproc/help.txt"))
};

/// The programs on this reader, in the order of `Commands`' variants.
static TOOLS: [&Tool; 10] = [
    &BLASTX,
    &BLASTN,
    &BLASTP,
    &TBLASTX,
    &TBLASTN,
    &PSIBLAST,
    &RPSBLAST,
    &RPSTBLASTN,
    &MAKEBLASTDB,
    &RPSBPROC,
];

/// The `Tool` of `program`, if it is one of the programs whose arguments clap converts.
pub(crate) fn tool(program: &str) -> Option<&'static Tool> {
    TOOLS.iter().copied().find(|tool| tool.program == program)
}

/// LOSAT's refusal of NCBI's argument `-name` of `program` that LOSAT does not run, by name
/// (`Refusal::Read`, the rows of `ncbi_app::args::refusal`): the program as LOSAT's refusals name
/// it (`Tool::refusal_name`) and why, if the row says; `None` for an argument that LOSAT runs.
pub(crate) fn named_refusal(tool: &Tool, name: &str) -> Option<String> {
    let Some(Refusal::Read(why)) = refusal(tool.program, name) else {
        return None;
    };
    let shown = match why {
        Some(why) => format!("{} ({why})", tool.refusal_name),
        None => tool.refusal_name.to_string(),
    };
    Some(format!(
        "the NCBI BLAST+ option -{name} is not supported by LOSAT's {shown}"
    ))
}

/// Whether LOSAT refuses a value of `tool`'s argument `name` that is not UTF-8: every value of
/// makeblastdb, whose arguments are text; every value of the search programs but the file names
/// that keep their bytes (`-query`, `-subject`, `-out`: `check_utf8_file_name`); none of rpsbproc,
/// whose file names keep their bytes and whose other keys take a value of their constraint. The
/// application start reads `-logfile` and `-conffile` by their bytes, before the reader.
fn refuses_bytes(tool: &Tool, name: &str) -> bool {
    match tool.program {
        "makeblastdb" => true,
        "rpsbproc" => false,
        _ => !is_file_argument(name) && !matches!(name, "query" | "subject" | "out"),
    }
}

// No NCBI counterpart: LOSAT's refusals of what it does not run (AGENTS.md rule 2), after NCBI's
// reader, so that NCBI's errors, help and `-dryrun` come first (line RG DECISIONS K4i-2 D4,
// K4i-2b-2 (1)).
/// LOSAT's refusal of the arguments `args` that `tool`'s reader read: the first given argument, in
/// the command line's order, that LOSAT refuses by name (`named_refusal`), or whose value is not
/// UTF-8 where LOSAT reads it as text (`cli::non_utf8_value_refusal`); `None` when it runs them.
/// The command line writes it as its refusals (`Error: ` and the text, exit status 1); the library
/// API's error is `Error::Argv` of the text.
pub(crate) fn refusal_after_read(tool: &Tool, args: &Args) -> Option<String> {
    args.given_values().find_map(|(name, _)| {
        named_refusal(tool, name).or_else(|| {
            (args.is_bytes(name) && refuses_bytes(tool, name))
                .then(|| super::non_utf8_value_refusal(tool.refusal_name, name))
        })
    })
}

// No NCBI counterpart: NCBI's `CArgs` is the program's arguments; LOSAT's programs read theirs from
// clap's types (line RG item K4i-3 takes them from `Args`).
/// The argv that clap converts for `tool`'s arguments `args` (`convert`): `losat <program>`, then
/// each argument given, in the command line's order, as `--name=` and its value's bytes, or
/// `--name` for a flag. The defaults are clap's (NCBI's reader adds them, and clap's types have
/// the same), the toolkit's file arguments are the application start's (only makeblastdb's clap
/// has its own `-logfile`), and an argument that LOSAT refuses is left out where the caller ends
/// before the refusals (a dry run, `refused`): its value is not read. A value that is not UTF-8,
/// which LOSAT refuses after a dry run's end, is clap's lossy text there.
fn clap_argv(tool: &Tool, args: &Args, command: &clap::Command, dry_run: bool) -> Vec<OsString> {
    let program = command
        .find_subcommand(tool.program)
        .expect("each program on the reader is a command");
    let has_long = |name: &str| {
        program
            .get_arguments()
            .any(|arg| arg.get_long() == Some(name))
    };
    let mut argv: Vec<OsString> = vec!["losat".into(), tool.program.into()];
    for (name, value) in args.given_values() {
        if is_file_argument(name) && !has_long(name) {
            continue;
        }
        if dry_run && named_refusal(tool, name).is_some() {
            continue;
        }
        match value {
            None => argv.push(format!("--{name}").into()),
            Some(value) => {
                let mut word = OsString::from(format!("--{name}="));
                match dry_run && args.is_bytes(name) && refuses_bytes(tool, name) {
                    true => word.push(value.to_string_lossy().as_ref()),
                    false => word.push(value),
                }
                argv.push(word);
            }
        }
    }
    argv
}

/// clap's conversion of `tool`'s arguments `args` into the program's types (`clap_argv`); `dry_run`
/// for a caller that ends before LOSAT's refusals (`clap_argv`). NCBI's reader has checked the
/// arguments, so clap's checks hold for every argv that it read (`ncbi_xmlhelp_tests`,
/// `cli_v2_tests`): a failure is LOSAT's error (`conversion_error`), never clap's text.
pub(crate) fn convert(tool: &Tool, args: &Args, dry_run: bool) -> Result<Commands, String> {
    let mut command = super::command_of(tool.program);
    let argv = clap_argv(tool, args, &command, dry_run);
    command
        .try_get_matches_from_mut(argv)
        .and_then(|mut matches| Cli::from_arg_matches_mut(&mut matches))
        .map(|cli| cli.command)
        .map_err(|error| conversion_error(tool, &error))
}

/// LOSAT's error where clap does not convert an argument that NCBI's reader read (`convert`): the
/// argument that clap names, without clap's words.
fn conversion_error(tool: &Tool, error: &clap::Error) -> String {
    use clap::error::{ContextKind, ContextValue};
    let name = match error.get(ContextKind::InvalidArg) {
        Some(ContextValue::String(arg)) => arg
            .split([' ', '='])
            .next()
            .map(|word| word.trim_start_matches('-').to_string()),
        _ => None,
    };
    match name {
        Some(name) => format!(
            "LOSAT's {} does not convert the value of -{name} that NCBI's argument reader read",
            tool.refusal_name
        ),
        None => format!(
            "LOSAT's {} does not convert the arguments that NCBI's argument reader read",
            tool.refusal_name
        ),
    }
}

/// Why a command line was not read into a program's arguments (`read_argv`, `try_parse_from`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArgvError {
    /// The reader stopped at a help word (`-h`, `-help`, `-help-full`, `-xmlhelp` where the
    /// caller prints the XML): the text that the command line prints to its standard output,
    /// exit status 0.
    Printed(String),
    /// NCBI's report of an argument error (the usage, the line of `=` and the posts), or the top
    /// level's error: the bytes that the command line writes to its standard error, as without a
    /// `-logfile` (exit status 1; the top level's, clap's, 2).
    Report(Vec<u8>),
    /// LOSAT's refusal: its text, which the command line writes after `Error: ` (exit status 1).
    Refused(String),
}

impl std::fmt::Display for ArgvError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ArgvError::Printed(text) | ArgvError::Refused(text) => f.write_str(text),
            ArgvError::Report(bytes) => f.write_str(&String::from_utf8_lossy(bytes)),
        }
    }
}

impl std::error::Error for ArgvError {}

/// Who reads a command line through `read_argv`, which decides what NCBI's toolkit words do there.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[expect(
    clippy::enum_variant_names,
    reason = "each variant names the NCBI application start that the caller runs"
)]
pub(crate) enum Caller {
    /// A caller that runs no NCBI application start (the engine's own argv: BLASTN's repeat
    /// search, the tests): `-logfile`, `-conffile` and `-dryrun` are refused
    /// (`toolkit_file_word_refusal`; makeblastdb keeps its own `-logfile`), and so are the
    /// versions (`version_word_refusal`) and `-xmlhelp` (`xml_help_refusal`).
    NoStart,
    /// The library API, whose calls run a start of their own (`api::start`): `-logfile`,
    /// `-conffile` and `-dryrun` are the start's; the versions and `-xmlhelp` are refused.
    OwnStart,
    /// LOSAT Web's ABI v2, whose runs make a start of their own (`web_adapter::start_run`) and
    /// print the command line's texts: the toolkit's words are the start's, the versions are
    /// printed before the reader (a version word that reaches it is refused), and `-xmlhelp`
    /// prints the XML, as on the command line.
    PrintingStart,
}

/// Reads `argv` (`losat <program> <words>`) as the command line reads it, for a caller other than
/// the command line's entry (`Caller`): the toolkit's words, NCBI's reader of the program (its
/// texts in memory: the help, or the report of an argument error, whose posts go to `target`, the
/// caller's start's `-logfile` target, as the command line's go to its own: `ArgvError::Report`
/// holds what the command line writes to its standard error, with NCBI's line that points to the
/// error logs when the target is not `cerr`), then, unless the argv is a dry run, LOSAT's refusals
/// (`refusal_after_read`), then clap's conversion (`convert`). A caller without a start
/// (`Caller::NoStart`) gives `cerr` (`DiagTarget::Stderr`): the posts are then the report's. An
/// argv that names no program on this reader is read by the top level (`top_level_error`).
pub(crate) fn read_argv(
    argv: &[OsString],
    caller: Caller,
    target: &DiagTarget,
) -> Result<Commands, ArgvError> {
    let Some(tool) = argv.get(1).and_then(|word| word.to_str()).and_then(tool) else {
        let error = super::top_level_error(argv);
        return Err(match error.use_stderr() {
            true => ArgvError::Report(super::top_level_message(&error).into_bytes()),
            false => ArgvError::Printed(super::top_level_message(&error)),
        });
    };
    let words = &argv[2..];
    let (words, dry_run) = match preparse(words) {
        Preparsed::Version(request) => {
            return Err(ArgvError::Refused(super::version_word_refusal(
                tool,
                request.word(),
            )))
        }
        Preparsed::Words {
            words: read,
            conffile,
            dryrun,
        } => {
            if caller == Caller::NoStart {
                if let Some(word) = super::toolkit_file_word(words, tool.program) {
                    return Err(ArgvError::Refused(super::toolkit_file_word_refusal(
                        tool, word,
                    )));
                }
            }
            debug_assert!(caller != Caller::NoStart || (conffile.is_none() && !dryrun));
            (read, dryrun)
        }
    };
    let mut out = Vec::new();
    let mut stderr = Vec::new();
    let mut diag = Diag::new(target, &mut stderr, tool);
    let xml_help = match caller {
        Caller::PrintingStart => Some(super::xml_help_text as fn(&str) -> Option<String>),
        Caller::NoStart | Caller::OwnStart => None,
    };
    let args = match crate::ncbi_app::read(
        tool,
        &words,
        &mut out,
        &mut diag,
        xml_help,
        &crate::report::run_files::OpenFiles::disk(),
    ) {
        Ok(args) => args,
        Err(ReadEnd::Exit(0)) => {
            return Err(ArgvError::Printed(
                String::from_utf8_lossy(&out).into_owned(),
            ))
        }
        Err(ReadEnd::Exit(_)) => {
            drop(diag);
            return Err(ArgvError::Report(stderr));
        }
        Err(ReadEnd::Refused(reason)) => return Err(ArgvError::Refused(reason)),
    };
    if !dry_run {
        if let Some(reason) = refusal_after_read(tool, &args) {
            return Err(ArgvError::Refused(reason));
        }
    }
    convert(tool, &args, dry_run).map_err(ArgvError::Refused)
}
