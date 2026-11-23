//! The pattern file of `-phi_pattern` as NCBI's `CPhiBlastArgs` reads it: line by line through a
//! buffer of 4,096 bytes, each line cut at its first NUL, the text after the first four bytes of
//! the last line that starts with `PA` taken as the pattern (lines starting with `ID` are read the
//! same way and not used).

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1382-1411
/// The size of NCBI's line buffer: a line of more bytes than one less than it ends the reading.
const LINE_BUFFER: usize = 4096;

/// Why the file gives no pattern.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PatternFileError {
    /// No `PA` line, or the last one empty after its first four bytes: NCBI's input error
    /// `PHI pattern not read` (exit 1, before the prolog).
    NotRead,
    /// A line `ID` or `PA` of two or three bytes (`size`): `std::string::substr` past the end,
    /// which NCBI's application reports as a standard exception (exit 255).
    Substr { size: usize },
}

impl PatternFileError {
    /// The error as the run reports it.
    pub(crate) fn into_error(self) -> anyhow::Error {
        match self {
            Self::NotRead => crate::blastinput::app::options_error("PHI pattern not read"),
            Self::Substr { size } => crate::cli::NativeError::std_exception(&format!(
                "basic_string::substr: __pos (which is 4) > this->size() (which is {size})"
            ))
            .into(),
        }
    }
}

/// The lines `std::istream::getline` gives with NCBI's buffer: up to the next LF (not kept), a
/// last line without LF included; the reading ends at a line of `LINE_BUFFER - 1` bytes not
/// followed by LF (the stream fails and the line is not used).
fn buffer_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut rest = bytes;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let stored = rest.len().min(LINE_BUFFER - 1);
        let line;
        if let Some(end) = rest[..stored].iter().position(|&byte| byte == b'\n') {
            line = &rest[..end];
            rest = &rest[end + 1..];
        } else if rest.len() == stored {
            line = rest;
            rest = &[];
        } else if rest[stored] == b'\n' {
            line = &rest[..stored];
            rest = &rest[stored + 1..];
        } else {
            rest = &[];
            return None;
        }
        Some(line)
    })
}

/// The pattern of the file at `path` as text, for the search strategy's `PHIPattern`
/// (`OptionInputs::phi_pattern`); `None` when the file does not give one (the run then stops with
/// the error of `read_pattern` where NCBI reads the file, `BlastpArgs::check_options`).
pub(crate) fn pattern_text(path: &std::path::Path) -> Option<String> {
    let bytes = std::fs::read(path).ok()?;
    let pattern = read_pattern(&bytes).ok()?;
    Some(String::from_utf8_lossy(&pattern).into_owned())
}

/// The pattern of the file's `bytes`, as the bytes NCBI passes to the compiler.
pub(crate) fn read_pattern(bytes: &[u8]) -> Result<Vec<u8>, PatternFileError> {
    let mut pattern: &[u8] = &[];
    for buffer in buffer_lines(bytes) {
        let line = buffer
            .iter()
            .position(|&byte| byte == 0)
            .map_or(buffer, |end| &buffer[..end]);
        let kind = &line[..line.len().min(2)];
        if kind == b"ID" || kind == b"PA" {
            if line.len() < 4 {
                return Err(PatternFileError::Substr { size: line.len() });
            }
            if kind == b"PA" {
                pattern = &line[4..];
            }
        }
    }
    if pattern.is_empty() {
        return Err(PatternFileError::NotRead);
    }
    Ok(pattern.to_vec())
}

#[cfg(test)]
#[path = "file_tests.rs"]
mod tests;
