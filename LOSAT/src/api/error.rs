//! The error of the library API.
//!
//! NCBI reports the three kinds of failure of a search separately: the arguments
//! (`CArgException`, `CInputException`), the input sequences (`CObjReaderParseException`) and
//! the search itself (`CBlastException`). The message of each kind is the command line's
//! (`losat::cli::render_message` for the parser, the error and its causes joined by `: ` for the
//! others, as LOSAT Web ABI v2 reports them; docs/web/abi_v2.md §4), as its bytes: NCBI's
//! diagnostic stream writes a message's `std::string` as it is, so a message that quotes an
//! input's bytes (a tool's diagnostics, `tools.rs`) holds a byte that is not UTF-8 where the
//! command line writes one.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
//! ```c++
//!     catch (const CArgException& e) {                                        \
//!         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//!         exit_code = BLAST_INPUT_ERROR;                                      \
//!     }                                                                       \
//!     catch (const CObjReaderParseException& e) {                             \
//!         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//!         exit_code = BLAST_INPUT_ERROR;                                      \
//!     }                                                                       \
//! ```
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
//! ```c++
//!     catch (const blast::CBlastException& e) {                               \
//!         const string& msg = e.GetMsg();                                     \
//!         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
//!             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
//!             exit_code = BLAST_INPUT_ERROR;                                  \
//! ```
#![warn(missing_docs, unused)]

use std::borrow::Cow;
use std::fmt;

use super::id::DebugBytes;

/// Why a library call failed. Every message is the command line's bytes for the same failure:
/// also a byte that is not UTF-8 (one that a tool's diagnostics quote) is kept.
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
    /// The argv: the parser's message (`cli::render_message`), or an option that the library
    /// API does not take (`validate`).
    Argv(Vec<u8>),
    /// A FASTA input: the reader's message (`read_fasta`), without its last newline. It may hold
    /// a NUL byte where NCBI writes one: a subject input that does not start with a defline ends
    /// with `BLAST query/options error: Query input doesn't start with a defline or comment, line
    /// N`, a NUL, a newline and `Please refer to the LOSAT user manual.` (`read_inputs`).
    Input(Vec<u8>),
    /// The search, its reports or a tool: the engine's error and its causes, joined by `: `, or
    /// what a tool wrote to its diagnostics, without the last newlines.
    Search(Vec<u8>),
}

impl Error {
    /// The message, whatever the kind: its bytes, as the command line writes them.
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

    /// The engine's error as the API reports it: the error and its causes joined by `: `
    /// (anyhow's alternate form), the shape of LOSAT Web ABI v2 (docs/web/abi_v2.md §4).
    pub(crate) fn search(error: &anyhow::Error) -> Self {
        Error::Search(format!("{error:#}").into_bytes())
    }

    /// The error of a run that failed after writing `text` to its standard error: that text,
    /// without its last line feeds (`Error::Search`).
    pub(crate) fn written(mut text: Vec<u8>) -> Self {
        while text.last() == Some(&b'\n') {
            text.pop();
        }
        Error::Search(text)
    }
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
mod tests {
    use super::Error;

    #[test]
    fn the_message_of_an_error_is_its_bytes_and_its_text_a_view_of_them() {
        let bytes = b"Error: [blastdbcmd] Unrecognized format specification: '%\xC3'".to_vec();
        let error = Error::Search(bytes.clone());
        assert_eq!(error.message(), &bytes[..]);
        let text = "Error: [blastdbcmd] Unrecognized format specification: '%\u{FFFD}'";
        assert_eq!(error.to_string_lossy(), text);
        assert_eq!(error.to_string(), text);
        // One U+FFFD per maximal ill-formed subpart: a truncated `€` (E2 82) is one subpart, a
        // lone continuation byte (80) and FF are one each.
        let error = Error::Input(b"a\xE2\x82b\x80\xFF".to_vec());
        assert_eq!(error.to_string(), "a\u{FFFD}b\u{FFFD}\u{FFFD}");
        // A message that is UTF-8 is its own text, borrowed.
        let error = Error::Argv("-query is required: é".into());
        assert_eq!(error.message(), "-query is required: é".as_bytes());
        assert!(matches!(
            error.to_string_lossy(),
            std::borrow::Cow::Borrowed("-query is required: é")
        ));
        assert_eq!(error.to_string(), "-query is required: é");
        assert_eq!(
            Error::Argv(String::from("x").into()),
            Error::Argv(b"x".to_vec())
        );
        // A tool's diagnostics lose their last line feeds only.
        assert_eq!(
            Error::written(b"a\n\xE9\n\n".to_vec()),
            Error::Search(b"a\n\xE9".to_vec())
        );
    }

    #[test]
    fn the_debug_form_of_an_error_is_its_kind_and_its_escaped_message() {
        assert_eq!(
            format!("{:?}", Error::Argv("-query is required".into())),
            r#"Argv("-query is required")"#
        );
        assert_eq!(
            format!("{:?}", Error::Search(b"q\xE9 \"t\"\n".to_vec())),
            r#"Search("q\xE9 \"t\"\n")"#
        );
        assert_eq!(
            format!("{:?}", Error::Input("é\0".into())),
            format!("Input({:?})", "é\0")
        );
        assert_eq!(
            format!("{:#?}", Error::Argv("x".into())),
            "Argv(\n    \"x\",\n)"
        );
    }
}
