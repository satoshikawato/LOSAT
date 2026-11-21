//! The bytes of file names and environment values, as NCBI's `std::string` holds them.
//!
//! NCBI BLAST+ keeps a file name of the command line, of a directory listing, of a database
//! header or of an environment variable as the bytes the system gave it, and writes those bytes
//! back in its messages and reports. On Unix LOSAT does the same through these functions, so a
//! name with a byte that is not UTF-8 (`d\xE9` of Latin-1) opens, is listed and is written as in
//! NCBI (NCBI BLAST+ 2.17.0, probe R1-R8 of line LIB-B). Off Unix (Windows, WASI) a name is
//! Unicode text to the program and keeps the conversion LOSAT made before
//! (`to_string_lossy`): no NCBI run there tells such a name apart (line LIB-B).
//
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-621

use std::borrow::Cow;
use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

/// The bytes of `name` (on Unix the system's own; elsewhere its text, `to_string_lossy`).
#[cfg(unix)]
pub(crate) fn os_str_bytes(name: &OsStr) -> Cow<'_, [u8]> {
    use std::os::unix::ffi::OsStrExt;
    Cow::Borrowed(name.as_bytes())
}

/// The bytes of `name` (on Unix the system's own; elsewhere its text, `to_string_lossy`).
#[cfg(not(unix))]
pub(crate) fn os_str_bytes(name: &OsStr) -> Cow<'_, [u8]> {
    match name.to_string_lossy() {
        Cow::Borrowed(text) => Cow::Borrowed(text.as_bytes()),
        Cow::Owned(text) => Cow::Owned(text.into_bytes()),
    }
}

/// The bytes of `path` ([`os_str_bytes`]).
pub(crate) fn path_bytes(path: &Path) -> Cow<'_, [u8]> {
    os_str_bytes(path.as_os_str())
}

/// The name whose bytes are `bytes` (on Unix the bytes themselves; elsewhere their text,
/// `String::from_utf8_lossy`).
#[cfg(unix)]
pub(crate) fn os_string_from_bytes(bytes: Vec<u8>) -> OsString {
    use std::os::unix::ffi::OsStringExt;
    OsString::from_vec(bytes)
}

/// The name whose bytes are `bytes` (on Unix the bytes themselves; elsewhere their text,
/// `String::from_utf8_lossy`).
#[cfg(not(unix))]
pub(crate) fn os_string_from_bytes(bytes: Vec<u8>) -> OsString {
    OsString::from(String::from_utf8_lossy(&bytes).into_owned())
}

/// The path whose bytes are `bytes` ([`os_string_from_bytes`]).
pub(crate) fn path_from_bytes(bytes: impl Into<Vec<u8>>) -> PathBuf {
    PathBuf::from(os_string_from_bytes(bytes.into()))
}

#[cfg(all(test, unix))]
#[path = "os_bytes_tests.rs"]
mod tests;
