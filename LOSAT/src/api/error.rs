//! The error of the library API.
//!
//! NCBI reports the three kinds of failure of a search separately: the arguments
//! (`CArgException`, `CInputException`), the input sequences (`CObjReaderParseException`) and
//! the search itself (`CBlastException`). The message of each kind is made from the text that the
//! command line writes to its standard error for the same failure, by the rules of [`Error`] (LOSAT
//! Web ABI v2 reports the engine's errors in the same shape; docs/web/abi_v2.md §4), as its bytes:
//! NCBI's diagnostic stream writes a message's `std::string` as it is, so a message that quotes an
//! input's bytes (a tool's diagnostics, `tools.rs`) holds a byte that is not UTF-8 where the
//! command line writes one.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
#![warn(missing_docs, unused)]

use std::borrow::Cow;
use std::fmt;

use super::id::DebugBytes;

/// Why a library call failed. The message is the text that the command line writes to its
/// standard error for the same failure, as bytes (also a byte that is not UTF-8, one that a tool's
/// diagnostics quote), by these rules:
///
/// - A message of NCBI's that a search reports (the command line's text for a failure of the
///   search that NCBI BLAST+ reports, such as the engine's error or `BLAST failed to write
///   output`, an [`Error::Search`]; and NCBI's text for an `-outfmt` value that its reading
///   refuses, an [`Error::Argv`]) is its bytes, with its trailing newline or newlines.
/// - LOSAT's own refusals, and the other errors that the command line writes after `Error: `, are
///   the text without that `Error: ` and without the newline after it. An error with causes, which
///   the command line writes under `Caused by:`, is the error and its causes joined by `: `
///   (anyhow's alternate form).
/// - [`Error::Argv`] of NCBI's argument reader (an argument error: NCBI's short usage, the line of
///   `=` and NCBI's posts) is the command line's standard error without its last newline. With
///   `-logfile`, a call's start opens the file and the reader posts there, as the command line's:
///   the error is then the usage, NCBI's line `Error in command-line arguments. See error logs for
///   more details.` and the line of `=`. [`validate`](super::validate), which makes no start and
///   opens no file, gives the report as the command line writes it without `-logfile`, its posts
///   included.
/// - The other messages, NCBI's too, are the command line's text without its last newlines:
///   [`Error::Input`] (the reader's, such as that of an input that does not start with a
///   definition line); the error of a tool ([`run_makeblastdb`](super::run_makeblastdb) and the
///   others), which is what the tool wrote to its standard error, held by the call's
///   `diagnostics` too; and that of NCBI's application start that ends the program
///   (`Critical: ...`: a `-conffile` that does not open or does not parse), which a search call's
///   `diagnostics` do not hold and a tool's do.
///
/// Text that a search writes to its standard error before the failure (warnings, notices) is not
/// in the message. It is in the diagnostics where the call returns them: [`run`](super::run()) and
/// [`run_to_file`](super::run_to_file) have written it to their writer;
/// [`run_collect`](super::run_collect) and [`run_batches`](super::run_batches) do not return it
/// with an error.
///
/// [`message`](Self::message) gives the bytes. As text, [`to_string_lossy`](Self::to_string_lossy)
/// and `Display` give one U+FFFD per maximal ill-formed subpart, as `String::from_utf8_lossy` and
/// Python's `decode("utf-8", "replace")` do: the text of a message that is UTF-8 is the message.
/// `Debug` writes the kind and the message as `str` does where it is UTF-8
/// (`Argv("-query is required")`), and otherwise as a quoted string with each byte that is not
/// UTF-8 as `\xNN` (`Search("Error: [blastdbcmd] Unrecognized format specification: '%\xC3'")`).
/// An error is made from text with `into`: `Error::Argv("-query is required".into())`.
#[derive(Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// The argv: NCBI's report of an argument error as the command line writes it to its standard
    /// error (the usage, a line of `=` and the posts), without its last newline; LOSAT's refusal
    /// of an argument (the text that the command line writes after `Error: `, without the
    /// newline); NCBI's message for an `-outfmt` value that NCBI's reading refuses (`Error:
    /// Formatting choice is out of range` and its newline), or an option or input that the
    /// library API does not take (`validate`, which has no command line counterpart). The
    /// argument error of a tool ([`run_makeblastdb`](super::run_makeblastdb)) is the text that the
    /// tool writes to its `diagnostics`, without its last newlines.
    Argv(Vec<u8>),
    /// A FASTA input: the reader's message (`read_fasta`, `read_inputs`, `read_input_file`)
    /// without its last newlines, and without the `Error: ` before LOSAT's own refusals (a file
    /// named `.gz`, `.bz2` or `.zst` that is not of that format). It may hold a NUL byte where
    /// NCBI writes one: a subject input that does not start with a defline ends with `LOSAT
    /// query/options error: Query input doesn't start with a defline or comment, line N`, a NUL,
    /// a newline and `Please refer to the LOSAT user manual.` (`read_inputs`).
    Input(Vec<u8>),
    /// The search, its reports, the start or a tool: a search's NCBI message with its trailing
    /// newlines, LOSAT's refusal without `Error: ` and the newline, or the text of a tool's
    /// failure or of a start that ends the program without its last newlines (the rules of
    /// [`Error`]).
    Search(Vec<u8>),
}

impl Error {
    /// The message, whatever the kind, as bytes (the rules of [`Error`]).
    pub fn message(&self) -> &[u8] {
        match self {
            Error::Argv(message) | Error::Input(message) | Error::Search(message) => message,
        }
    }

    /// The message as text, with one U+FFFD per maximal ill-formed subpart
    /// (`String::from_utf8_lossy`); borrowed where every byte is UTF-8. `Display` writes it.
    pub fn to_string_lossy(&self) -> Cow<'_, str> {
        String::from_utf8_lossy(self.message())
    }

    /// The engine's error as the API reports it (`chain_bytes`): a `NativeError`'s bytes as the
    /// command line writes them, other errors without the command line's `Error: ` and newline, the
    /// error and its causes joined by `: ` (anyhow's alternate form), the shape of LOSAT Web ABI v2
    /// (docs/web/abi_v2.md §4).
    pub(crate) fn search(error: &anyhow::Error) -> Self {
        Error::Search(chain_bytes(error))
    }

    /// The error of a run that failed after writing `text` to its standard error: that text,
    /// without its last line feeds (`Error::Search`).
    pub(crate) fn written(text: Vec<u8>) -> Self {
        Error::Search(without_line_feeds(text))
    }

    /// The parser's refusal of a tool's argv after the tool wrote `text`, the refusal's message
    /// last, to its standard error: that text without its last line feeds, as `written`'s
    /// (`Error::Argv`).
    pub(crate) fn written_argv(text: Vec<u8>) -> Self {
        Error::Argv(without_line_feeds(text))
    }

    /// The error of an argv that NCBI's reader and LOSAT's refusals end (`cli::ArgvError`): the
    /// report of an argument error as the command line writes it to its standard error with the
    /// reader's `-logfile` target (`written_argv`; `validate`'s target is `cerr`), LOSAT's refusal
    /// without `Error: ` and the newline, or the help text that the command line prints in place of
    /// a run.
    pub(crate) fn of_argv(error: crate::cli::ArgvError) -> Self {
        match error {
            crate::cli::ArgvError::Report(text) => Error::written_argv(text),
            crate::cli::ArgvError::Refused(text) | crate::cli::ArgvError::Printed(text) => {
                Error::Argv(text.into_bytes())
            }
        }
    }
}

/// The bytes of `error` in anyhow's alternate form (the error and its causes joined by `: `), with
/// the bytes of each `cli::NativeError` and `NcbiException` of the chain, the text the command
/// line writes: a name that is not UTF-8 keeps its bytes there (line LIB-B), where `{error:#}`
/// would write U+FFFD.
pub(crate) fn chain_bytes(error: &anyhow::Error) -> Vec<u8> {
    let mut bytes = Vec::new();
    for (at, cause) in error.chain().enumerate() {
        if at > 0 {
            bytes.extend_from_slice(b": ");
        }
        // A `NativeError`, or an exception of NCBI's that a report's writer raised mid-report
        // (the `-html` links' registry): `CATCH_ALL`'s report of it, as the command line writes
        // it, whichever program wrote the report (`cli::NativeError::caught`).
        if let Some(native) = crate::cli::NativeError::caught(cause) {
            bytes.extend_from_slice(&native.message);
        } else if let Some(exception) = cause.downcast_ref::<crate::ncbi_exception::NcbiException>()
        {
            bytes.extend_from_slice(&exception.what());
        } else {
            bytes.extend_from_slice(cause.to_string().as_bytes());
        }
    }
    bytes
}

/// `bytes` without its last line feeds.
fn without_line_feeds(mut bytes: Vec<u8>) -> Vec<u8> {
    while bytes.last() == Some(&b'\n') {
        bytes.pop();
    }
    bytes
}

/// [`chain_bytes`] without the last line feeds (an input's error, `Error::Input`).
pub(crate) fn input_bytes(error: &anyhow::Error) -> Vec<u8> {
    without_line_feeds(chain_bytes(error))
}

/// The text of the message (`to_string_lossy`).
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_string_lossy())
    }
}

/// The kind and the message, the message as `Id`'s `Debug` writes its bytes.
impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (kind, message) = match self {
            Error::Argv(message) => ("Argv", message),
            Error::Input(message) => ("Input", message),
            Error::Search(message) => ("Search", message),
        };
        f.debug_tuple(kind).field(&DebugBytes(message)).finish()
    }
}

impl std::error::Error for Error {}

#[cfg(test)]
#[path = "error_tests.rs"]
mod tests;
