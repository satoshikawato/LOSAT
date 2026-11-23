//! The text of an operating system error as NCBI BLAST+ writes it in its messages.
//!
//! NCBI appends the C library's `strerror` text of `errno` to the messages of a file that it
//! cannot open or map (`... Cannot memory map file '<path>':Permission denied`). A native LOSAT
//! takes the same text from the system it runs on (std prints it before ` (os error N)`), as
//! NCBI does there. LOSAT Web (WASI) has no system of its own whose NCBI build it could match: it
//! gives the native CLI's bytes on the oracle platform, Linux with glibc, whose texts for the
//! codes that opening or reading a file can give are kept here. They are NCBI BLAST+ 2.17.0's
//! output for each code (line WQA D18 probe 2, 2026-11-21: `open()` of a database file made to
//! fail with the code), not text taken from a C library's source.

use std::borrow::Cow;

/// The text NCBI BLAST+ appends for the operating system error of `error`, or `None` where LOSAT
/// does not know it: an error without a system code (a read that ended early), a code outside
/// `ORACLE_TEXTS` in LOSAT Web, or a platform whose text NCBI takes from another call (Windows),
/// where the caller keeps LOSAT's own words.
pub(crate) fn ncbi_strerror(error: &std::io::Error) -> Option<Cow<'static, str>> {
    let code = error.raw_os_error()?;
    system_text(error, code)
}

/// The system's own text, before std's ` (os error N)`.
#[cfg(unix)]
fn system_text(error: &std::io::Error, code: i32) -> Option<Cow<'static, str>> {
    let text = error.to_string();
    text.strip_suffix(&format!(" (os error {code})"))
        .map(|strerror| Cow::Owned(strerror.to_owned()))
}

/// The oracle platform's text of the WASI code `code` (wasi-libc's own texts differ: `I/O error`).
#[cfg(target_os = "wasi")]
fn system_text(_error: &std::io::Error, code: i32) -> Option<Cow<'static, str>> {
    oracle_text(code).map(Cow::Borrowed)
}

/// Windows and the other targets: NCBI's text comes from another call there (`WIN_LAST_ERROR_STR`
/// on Windows), which LOSAT does not reproduce.
#[cfg(not(any(unix, target_os = "wasi")))]
fn system_text(_error: &std::io::Error, _code: i32) -> Option<Cow<'static, str>> {
    None
}

/// The texts of glibc on Linux, as NCBI BLAST+ 2.17.0 wrote them after `Cannot memory map file
/// '<path>':` when `open()` failed with the Linux code (EMFILE and ENFILE: NCBI writes its own
/// `Too many open files, please raise the open file limit` for a text that holds `Too many open
/// files`, which these do). The WASI values are those of the WASI preview 1 interface
/// (`wasi_snapshot_preview1.witx`), which wasi-libc's `errno` and so std's `raw_os_error` give.
/// Each entry: the WASI code, the Linux code (the tests check the text against the system's on
/// Linux), the text.
#[cfg(any(target_os = "wasi", test))]
pub(crate) const ORACLE_TEXTS: [(i32, i32, &str); 11] = [
    (63, 1, "Operation not permitted"),
    (44, 2, "No such file or directory"),
    (29, 5, "Input/output error"),
    (48, 12, "Cannot allocate memory"),
    (2, 13, "Permission denied"),
    (54, 20, "Not a directory"),
    (31, 21, "Is a directory"),
    (41, 23, "Too many open files in system"),
    (33, 24, "Too many open files"),
    (37, 36, "File name too long"),
    (32, 40, "Too many levels of symbolic links"),
];

/// The oracle platform's text of the WASI code `code`.
#[cfg(any(target_os = "wasi", test))]
pub(crate) fn oracle_text(code: i32) -> Option<&'static str> {
    ORACLE_TEXTS
        .iter()
        .find(|&&(wasi, _, _)| wasi == code)
        .map(|&(_, _, text)| text)
}

#[cfg(test)]
#[path = "os_error_text_tests.rs"]
mod tests;
