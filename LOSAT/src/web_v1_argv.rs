//! The words of a request of LOSAT Web's ABI v1 (frozen, plan TD-1), as its BLASTX and TBLASTN
//! read them: with the reader that LOSAT's command line had before it read its arguments with
//! NCBI's argument reader (`cli::reader`). v1 keeps its own copy of that reader and of its texts
//! (line RG DECISIONS K4i-2 D5), so that its frozen messages stay
//! (`tests/fixtures/web_v1_tblastn`); no other caller reads an argv here.
//!
//! No NCBI counterpart: NCBI has no Web engine.

use std::ffi::OsString;

use clap::error::ErrorKind;
use clap::{CommandFactory, FromArgMatches};

use crate::cli::{Cli, Commands};

/// The arguments of the request `argv` (`losat <program> <words>`, BLASTX or TBLASTN) read with
/// `command` (the program's command, `cli::command_of`, with v1's own value parsers), or the
/// message of their error as v1 gives it (`render`). ABI v1 runs no NCBI application start: the
/// versions and the toolkit's words that the start reads (`-logfile`, `-conffile`, `-dryrun`) are
/// refused; its callers refuse the help words before.
pub(crate) fn parse(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, String> {
    read(command, argv).map_err(|error| render(&error))
}

/// `parse` without the rendering of the error.
fn read(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, clap::Error> {
    use crate::blastinput::ncbi_environment::{preparse, Preparsed};
    let program = argv
        .get(1)
        .and_then(|word| word.to_str())
        .unwrap_or_default()
        .to_string();
    let words = argv.get(2..).unwrap_or(&[]);
    let words = match preparse(words) {
        Preparsed::Version(request) => {
            return Err(raw(crate::cli::version_word_refusal(
                &program,
                request.word(),
            )))
        }
        Preparsed::Words { words: read, .. } => {
            if let Some(word) = crate::cli::toolkit_file_word(words, &program) {
                return Err(raw(crate::cli::toolkit_file_word_refusal(&program, word)));
            }
            read
        }
    };
    let scope = command
        .find_subcommand(&program)
        .expect("v1 reads the words of one of its programs");
    let mut translated: Vec<OsString> = vec!["losat".into(), program.clone().into()];
    let mut input = words.into_iter();
    while let Some(token) = input.next() {
        let text = token.to_str().unwrap_or("").to_string();
        // A last `--` changes nothing in TBLASTN, as NCBI's reader; BLASTX's reader of v1 did not
        // take it.
        if text == "--" && input.len() == 0 && program == "tblastn" {
            continue;
        }
        let key = text.strip_prefix('-').filter(|key| !key.starts_with('-'));
        let (name, inline) = key
            .unwrap_or("")
            .split_once('=')
            .map_or((key.unwrap_or(""), None), |(k, v)| (k, Some(v)));
        let Some(arg) = scope
            .get_arguments()
            .find(|arg| arg.get_long() == Some(name))
        else {
            return Err(unknown_option(&program, name, &text));
        };
        if arg.get_action().takes_values() {
            let value = match inline {
                Some(value) => OsString::from(value),
                None => input.next().ok_or_else(|| {
                    clap::Error::raw(
                        ErrorKind::InvalidValue,
                        format!("-{name} requires one value"),
                    )
                })?,
            };
            let mut word = OsString::from(format!("--{name}="));
            word.push(value);
            translated.push(word);
        } else {
            if inline.is_some() {
                return Err(clap::Error::raw(
                    ErrorKind::TooManyValues,
                    format!("-{name} is a flag and does not take a value"),
                ));
            }
            translated.push(format!("--{name}").into());
        }
    }
    let mut matches = command.clone().try_get_matches_from(&translated)?;
    Cli::from_arg_matches_mut(&mut matches)
        .map(|cli| cli.command)
        .map_err(|error| error.format(&mut command.clone()))
}

/// A refusal of v1's reader, as its reader wrote them (a clap error of LOSAT's text).
fn raw(message: String) -> clap::Error {
    clap::Error::raw(ErrorKind::InvalidValue, message)
}

/// v1's rejection of an option that the program does not have: NCBI's options that LOSAT
/// refuses by name, and the rest.
fn unknown_option(program: &str, name: &str, text: &str) -> clap::Error {
    let tool = crate::cli::reader::tool(program).expect("v1's programs are on the reader");
    match crate::cli::reader::named_refusal(tool, name) {
        Some(refusal) => raw(refusal),
        None => clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!("unknown option or argument '{text}'; use -help to list the options"),
        ),
    }
}

/// The text of v1's `error`: clap's, with each program's long options written with one dash, as
/// v1's messages have them.
fn render(error: &clap::Error) -> String {
    let mut message = error
        .to_string()
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let command = Cli::command();
    let toolkit = ["logfile", "conffile"].into_iter();
    let longs = std::iter::once(&command)
        .chain(command.get_subcommands())
        .flat_map(|scope| scope.get_arguments().filter_map(clap::Arg::get_long))
        .map(str::to_string)
        .chain(toolkit.map(str::to_string));
    for name in longs.filter(|name| !matches!(name.as_str(), "help" | "version") && name.len() > 1)
    {
        message = message.replace(&format!("--{name}"), &format!("-{name}"));
    }
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
        .replace("see a summary with '-help'", "see a summary with '-h'")
        .replace("see more with '--help'", "see more with '-help'")
}
