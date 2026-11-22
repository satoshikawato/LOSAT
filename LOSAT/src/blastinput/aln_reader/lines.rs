//! The peek-ahead line input of NCBI's alignment reader (`CPeekAheadStream`): the lines of the
//! stream as `getline` gives them, which the format guesser peeks at and the scanner then reads.

/// The lines of the input, `getline`'s: split at `\n` only, a last line without a newline is a
/// line, a `\r` stays in its line, and a file that ends with `\n` has no empty line after it.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_peek_ahead.hpp:64-94
/// The peeked lines and the lines read after them are one sequence of lines, so the stream is
/// kept as that sequence; `peeked` counts the lines taken by `peek_line` and `read` those of
/// `read_line`. The line numbers (from 1) only name lines in NCBI's messages, which the
/// alignment reader's caller never prints (every report throws and the throw is swallowed), but
/// the Clustal scanner tests whether `lineCount` is 1.
pub(crate) struct PeekAheadStream<'a> {
    lines: Vec<&'a [u8]>,
    /// Lines handed out by `peek_line` (they stay available to `read_line`).
    peeked: usize,
    /// Lines handed out by `read_line`.
    read: usize,
}

impl<'a> PeekAheadStream<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Self {
        let mut lines: Vec<&[u8]> = bytes.split(|&b| b == b'\n').collect();
        // `bytes.split` gives an empty piece after a final newline (and one piece for no bytes):
        // `getline` gives no line for either.
        if lines.last().is_some_and(|last| last.is_empty()) {
            lines.pop();
        }
        PeekAheadStream {
            lines,
            peeked: 0,
            read: 0,
        }
    }

    /// `PeekLine`: the next line that has not been peeked, `None` at the end of the stream.
    pub(crate) fn peek_line(&mut self) -> Option<&'a [u8]> {
        let line = self.lines.get(self.peeked).copied()?;
        self.peeked += 1;
        Some(line)
    }

    /// `ReadLine`: the next line with its number (from 1), `None` at the end of the stream.
    pub(crate) fn read_line(&mut self) -> Option<(&'a [u8], i32)> {
        let line = self.lines.get(self.read).copied()?;
        self.read += 1;
        Some((line, self.read as i32))
    }
}
