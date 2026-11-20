//! The error of the library API.
//!
//! NCBI reports the three kinds of failure of a search separately: the arguments
//! (`CArgException`, `CInputException`), the input sequences (`CObjReaderParseException`) and
//! the search itself (`CBlastException`). The message of each kind is the command line's
//! (`losat::cli::render_message` for the parser, the error and its causes joined by `: ` for the
//! others, as LOSAT Web ABI v2 reports them; docs/web/abi_v2.md §4).
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

use std::fmt;

/// Why a library call failed. Every message is the command line's text for the same failure.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// The argv: the parser's message (`cli::render_message`), or an option that the library
    /// API does not take (`validate`).
    Argv(String),
    /// A FASTA input: the reader's message (`read_fasta`), without its last newline. It may hold
    /// a NUL byte where NCBI writes one: a subject input that does not start with a defline ends
    /// with `BLAST query/options error: Query input doesn't start with a defline or comment, line
    /// N`, a NUL, a newline and `Please refer to the BLAST+ user manual.` (`read_inputs`).
    Input(String),
    /// The search or its reports: the engine's error and its causes, joined by `: `.
    Search(String),
}

impl Error {
    /// The message, whatever the kind.
    pub fn message(&self) -> &str {
        match self {
            Error::Argv(message) | Error::Input(message) | Error::Search(message) => message,
        }
    }

    /// The engine's error as the API reports it: the error and its causes joined by `: `
    /// (anyhow's alternate form), the shape of LOSAT Web ABI v2 (docs/web/abi_v2.md §4).
    pub(crate) fn search(error: &anyhow::Error) -> Self {
        Error::Search(format!("{error:#}"))
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for Error {}
