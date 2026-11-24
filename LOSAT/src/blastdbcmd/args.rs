//! The arguments of `losat blastdbcmd`: NCBI's argument reader (`ncbi_app::args`) reads the
//! descriptions of blastdbcmd 2.17.0 (`CBlastDBCmdApp::Init`, blastdbcmd.cpp:901-1167) from its
//! `-xmlhelp` XML (`ncbi_xmlhelp/blastdbcmd.xml`), so that every command line that NCBI's reader
//! refuses is refused with NCBI's text, also for the options that LOSAT's BLASTDBCMD does not serve
//! (the program refuses those after the reader accepts them); `args_tests.rs` tests that reading.
//! With it `string_to_int8_message`, `NStr::StringToInt8` (`ncbi_str`) for the values blastdbcmd
//! converts itself (`-mask_sequence_with`, the lines of its batch files).

/// `NStr::StringToInt8` of `value` (`ncbi_str::string_to_int8`): the number, or the message of
/// NCBI's exception as `GetMsg` gives it (without the position).
pub(super) fn string_to_int8_message(value: impl AsRef<[u8]>) -> Result<i64, String> {
    crate::ncbi_str::string_to_int8(value.as_ref(), 0)
        .map_err(|error| String::from_utf8_lossy(error.get_msg()).into_owned())
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
