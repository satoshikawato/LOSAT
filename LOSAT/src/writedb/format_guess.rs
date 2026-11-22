//! Port of the decision `makeblastdb` takes in `CMakeBlastDBApp::x_VerifyInputFilesType` (NCBI BLAST+ 2.17.0) when it prints
//! "Input file X does NOT appear to be FASTA (processing anyway)": whether `CFormatGuess::GuessFormat` returns `eFasta` for
//! the input stream. Only the decision is ported, not the warning.
//!
//! The rule, as read from the source:
//!
//! * `makeblastdb` prefers `eBinaryASN`, `eTextASN` and `eFasta` and disables every other format. `GuessFormat` then walks
//!   `sm_CheckOrder`, where `eFasta` comes before `eTextASN` and `eBinaryASN`, and returns the first preferred format whose
//!   test is true. `eFasta` is the first of the three, so the result is `eFasta` exactly when `TestFormatFasta` is true
//!   (a later format can only claim an input that `TestFormatFasta` rejected). The disabled formats are never tried.
//! * `GuessFormat` needs a test buffer first (`EnsureTestBuffer`). Without one it answers `eNewick` or `eUnknown`,
//!   never `eFasta`.
//! * `TestFormatFasta` works on that sample only: the first 8096 bytes of the stream, or more when the sample is all comment
//!   (see `ensure_test_buffer`).
//!
//! C `isalpha`, `isdigit`, `isspace` and `isprint` are the ones of the "C" locale (ASCII only). The C++ works on bytes, so
//! this port works on bytes.
//!
//! The guess is the one port of `CFormatGuess` (`utils::format_guess`), which the `-in_pssm` reader and makeblastdb's
//! mask files share. `TestFormatFasta`'s `NStr::Find` of `"|"` in `m_pTestBuffer` reads past the
//! sample when it holds 10 bytes or fewer without a NUL or `|`: this port searches the sample only
//! there (the answer is false). A file that cannot be opened is the caller's error, where LOSAT
//! stops and makeblastdb goes on with an empty stream and exits with status 0.

use std::io::Read;

use crate::utils::format_guess::FormatGuess;

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:841-850
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:452-491
/// `true` when makeblastdb takes the input for FASTA (no warning): `x_GuessFileType` returns `CFormatGuess::eFasta`.
///
/// The preferred formats are tried in the order of `sm_CheckOrder`, where `eFasta` precedes `eTextASN` and `eBinaryASN`,
/// and all other formats are disabled, so the answer is `TestFormatFasta` on the test buffer (false when there is none).
/// A read error counts as end of data.
pub fn appears_fasta<R: Read>(input: R) -> bool {
    let mut guess = FormatGuess::new(input);
    guess.ensure_test_buffer() && guess.test_format_fasta().unwrap_or(false)
}

#[cfg(test)]
#[path = "format_guess_tests.rs"]
mod tests;
