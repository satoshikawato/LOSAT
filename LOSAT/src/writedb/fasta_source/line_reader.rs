//! NCBI's `CBufferedLineReader` (`c++/src/util/line_reader.cpp`), the line reader that
//! `CFastaBioseqSource` gives to `CFastaReader`: lines end at `\n`, `\r` or `\r\n`, and
//! the input is read 32 KiB at a time.
//!
//! The buffer is kept as NCBI keeps it, because one result depends on its bytes: after a
//! `\r` that is the last byte of a buffer, when the next read brings only the `\n` that
//! follows it, the reader is not at the end of the input and `PeekChar` returns the byte
//! after that `\n` in the buffer, which is what the read before left there (the input's
//! byte 32,767 bytes before its last one). A `>` there ends makeblastdb with "Expected
//! defline around line N" (`tests::stale_peek_after_cr_lf_across_buffers`).

use std::io::Read;

/// `m_BufferSize`.
///
/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:433-446
pub(super) const BUFFER_SIZE: usize = 32 * 1024;

/// What `m_Line` points to.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Line {
    /// The null `CTempString`: before the first line, and after `operator++` at the end.
    Null,
    /// `m_Buffer[start..end]`.
    Buffer { start: usize, end: usize },
    /// `m_String`.
    String,
}

/// `CBufferedLineReader` over `CStreamReader`.
pub(super) struct BufferedLineReader<R: Read> {
    reader: R,
    /// The input returned no bytes once (`eRW_Eof`); it is not read again, like a C++
    /// stream whose `eofbit` is set.
    exhausted: bool,
    /// A failed read, kept for `FastaBioseqSource::next_bioseq`.
    read_error: Option<std::io::Error>,
    /// `m_Eof`.
    eof: bool,
    /// `m_UngetLine`.
    unget_line: bool,
    /// `m_Buffer`.
    buffer: Box<[u8]>,
    /// `m_Pos` and `m_End`, as indexes of `buffer`.
    pos: usize,
    end: usize,
    /// `m_Line`.
    line: Line,
    /// `m_String`.
    string: Vec<u8>,
    /// `m_LineNumber`.
    line_number: u64,
}

/// The first `\n` or `\r` of `bytes`.
fn find_eol(bytes: &[u8]) -> Option<usize> {
    bytes
        .iter()
        .position(|&byte| byte == b'\n' || byte == b'\r')
}

impl<R: Read> BufferedLineReader<R> {
    pub(super) fn new(reader: R) -> Self {
        let mut lines = Self {
            reader,
            exhausted: false,
            read_error: None,
            eof: false,
            unget_line: false,
            buffer: vec![0; BUFFER_SIZE].into_boxed_slice(),
            pos: 0,
            end: 0,
            line: Line::Null,
            string: Vec::new(),
            line_number: 0,
        };
        lines.read_buffer();
        lines
    }

    /// The error of a failed read, once.
    pub(super) fn take_read_error(&mut self) -> Option<std::io::Error> {
        self.read_error.take()
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:469-472
    pub(super) fn at_eof(&self) -> bool {
        self.eof && !self.unget_line
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:475-495
    /// `*m_Pos` with `m_Pos` at `m_End` (the case of the module comment) is the byte that an
    /// earlier read left in the buffer; the buffer is 32 KiB, so the byte exists. The reader
    /// never calls this at the end of the input.
    pub(super) fn peek_char(&self) -> u8 {
        if self.unget_line {
            return self.line().first().copied().unwrap_or(0);
        }
        match self.buffer.get(self.pos) {
            Some(b'\n' | b'\r') | None => 0,
            Some(&byte) => byte,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:498-506
    pub(super) fn unget_line(&mut self) {
        if self.unget_line || self.line == Line::Null {
            return;
        }
        self.line_number -= 1;
        self.unget_line = true;
    }

    /// `operator++`: reads the next line.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:509-568
    /// `m_LastReadSize` serves `GetPosition` only, which the FASTA reader does not call.
    pub(super) fn next_line(&mut self) {
        if self.at_eof() {
            self.line = Line::Null;
            return;
        }
        self.line_number += 1;
        if self.unget_line {
            self.unget_line = false;
            return;
        }
        let (start, end) = (self.pos, self.end);
        let Some(offset) = find_eol(&self.buffer[start..end]) else {
            self.load_long();
            return;
        };
        let mut p = start + offset;
        let terminator = self.buffer[p];
        self.line = Line::Buffer { start, end: p };
        p += 1;
        self.pos = p;
        if terminator == b'\n' {
            if p == end {
                self.keep_line();
                self.read_buffer();
            }
            return;
        }
        if p == end {
            self.keep_line();
            if self.read_buffer() && self.buffer[self.pos] == b'\n' {
                self.pos += 1;
            }
            return;
        }
        if self.buffer[p] != b'\n' {
            return;
        }
        p += 1;
        self.pos = p;
        if p == end {
            self.keep_line();
            self.read_buffer();
        }
    }

    /// Copies the line into `m_String` and makes `m_Line` refer to it, before the buffer is read
    /// again.
    fn keep_line(&mut self) {
        if let Line::Buffer { start, end } = self.line {
            self.string.clear();
            self.string.extend_from_slice(&self.buffer[start..end]);
            self.line = Line::String;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:571-616
    fn load_long(&mut self) {
        self.string.clear();
        self.string
            .extend_from_slice(&self.buffer[self.pos..self.end]);
        self.line = Line::String;
        while self.read_buffer() {
            let (start, end) = (self.pos, self.end);
            let Some(offset) = find_eol(&self.buffer[start..end]) else {
                self.string.extend_from_slice(&self.buffer[start..end]);
                continue;
            };
            let mut p = start + offset;
            let terminator = self.buffer[p];
            self.string.extend_from_slice(&self.buffer[start..p]);
            p += 1;
            if p == end {
                if self.read_buffer()
                    && self.pos < self.end
                    && terminator == b'\r'
                    && self.buffer[self.pos] == b'\n'
                {
                    self.pos += 1;
                }
            } else {
                if terminator == b'\r' && self.buffer[p] == b'\n' {
                    p += 1;
                    if p == end {
                        self.read_buffer();
                        p = self.pos;
                    }
                }
                self.pos = p;
            }
            return;
        }
    }

    /// Fills the buffer from its start; false at the end of the input.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:619-653
    /// `CStreamReader::Read` takes the bytes with `sgetn`, which returns fewer than asked
    /// only at the end of the input, and answers `eRW_Eof` when it got none:
    ///
    /// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:604-622
    /// A read that fails ends the input here and is reported by the source (LOSAT's: a
    /// Rust reader returns the error, where `sgetn` returns a short count or throws).
    fn read_buffer(&mut self) -> bool {
        if self.eof {
            return false;
        }
        self.pos = 0;
        self.end = 0;
        let mut size = 0;
        while !self.exhausted && size < BUFFER_SIZE {
            match self.reader.read(&mut self.buffer[size..]) {
                Ok(0) => self.exhausted = true,
                Ok(count) => size += count,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(error) => {
                    self.read_error = Some(error);
                    self.exhausted = true;
                }
            }
        }
        if size == 0 {
            self.eof = true;
            return false;
        }
        self.end = size;
        true
    }

    /// `operator*`: the current line, without its end.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:656-669
    /// The reader takes the line only after `next_line` (never after `unget_line`).
    pub(super) fn line(&self) -> &[u8] {
        match self.line {
            Line::Null => &[],
            Line::Buffer { start, end } => &self.buffer[start..end],
            Line::String => &self.string,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:682-690
    pub(super) fn line_number(&self) -> u64 {
        self.line_number
    }
}
