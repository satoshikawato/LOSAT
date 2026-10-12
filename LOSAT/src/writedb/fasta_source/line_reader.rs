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
/// ```c++
/// CBufferedLineReader::CBufferedLineReader(CNcbiIstream& is,
///                                          EOwnership ownership)
///     : m_Reader(new CStreamReader(is, ownership)),
///       m_Eof(false),
///       m_UngetLine(false),
///       m_BufferSize(32*1024),
///       m_Buffer(new char[m_BufferSize]),
///       m_Pos(m_Buffer.get()),
///       m_End(m_Pos),
///       m_InputPos(0),
///       m_LineNumber(0)
/// {
///     x_ReadBuffer();
/// }
/// ```
pub(super) const BUFFER_SIZE: usize = 32 * 1024;

/// What `m_Line` points to.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Line {
    /// `CTempString(NULL)`: before the first line, and after `operator++` at the end.
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
    /// ```c++
    /// bool CBufferedLineReader::AtEOF(void) const
    /// {
    ///     return m_Eof && !m_UngetLine;
    /// }
    /// ```
    pub(super) fn at_eof(&self) -> bool {
        self.eof && !self.unget_line
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:475-495
    /// ```c++
    /// char CBufferedLineReader::PeekChar(void) const
    /// {
    ///     _ASSERT(!AtEOF());
    ///     /* If at EOF - undefined behavior */
    ///     if (AtEOF()) {
    ///         return *m_Pos;
    ///     }
    ///     /* If line was ungot - return its first symbol */
    ///     if (m_UngetLine) {
    ///         /* if line is empty - return 0 */
    ///         if (m_Line.empty()) {
    ///             return 0;
    ///         }
    ///         return *m_Line.begin();
    ///     }
    ///     /* If line is empty - return 0 */
    ///     if (*m_Pos == '\n' || *m_Pos == '\r') {
    ///         return 0;
    ///     }
    ///     return *m_Pos;
    /// }
    /// ```
    /// `*m_Pos` with `m_Pos == m_End` (the case of the module comment) is the byte that an
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
    /// ```c++
    /// void CBufferedLineReader::UngetLine(void)
    /// {
    ///     _ASSERT(!m_UngetLine && m_Line.begin());
    ///     /* If after UngetLine() or after constructor - noop */
    ///     if (m_UngetLine || m_Line.begin() == NULL)
    ///         return;
    ///     --m_LineNumber;
    ///     m_UngetLine = true;
    /// }
    /// ```
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
    /// ```c++
    /// CBufferedLineReader& CBufferedLineReader::operator++(void)
    /// {
    ///     /* If at EOF - noop */
    ///     if (AtEOF()) {
    ///         m_Line = CTempString(NULL);
    ///         return *this;
    ///     }
    ///     ++m_LineNumber;
    ///     if ( m_UngetLine ) {
    ///         _ASSERT(m_Line.begin());
    ///         m_UngetLine = false;
    ///         return *this;
    ///     }
    ///     // check if we are at the buffer end
    ///     const char* start = m_Pos;
    ///     const char* end = m_End;
    ///     for ( const char* p = start; p < end; ++p ) {
    ///         if ( *p == '\n' ) {
    ///             m_Line = CTempString(start, p - start);
    ///             m_LastReadSize = p + 1 - start;
    ///             m_Pos = ++p;
    ///             if ( p == end ) {
    ///                 m_String = m_Line;
    ///                 m_Line = m_String;
    ///                 x_ReadBuffer();
    ///             }
    ///             return *this;
    ///         }
    ///         else if ( *p == '\r' ) {
    ///             m_Line = CTempString(start, p - start);
    ///             m_LastReadSize = p + 1 - start;
    ///             m_Pos = ++p;
    ///             if ( p == end ) {
    ///                 m_String = m_Line;
    ///                 m_Line = m_String;
    ///                 if ( x_ReadBuffer() ) {
    ///                     p = m_Pos;
    ///                     if ( *p == '\n' ) {
    ///                         m_Pos = p+1;
    ///                         ++m_LastReadSize;
    ///                     }
    ///                 }
    ///                 return *this;
    ///             }
    ///             if ( *p != '\n' ) {
    ///                 return *this;
    ///             }
    ///             ++m_LastReadSize;
    ///             m_Pos = ++p;
    ///             if ( p == end ) {
    ///                 m_String = m_Line;
    ///                 m_Line = m_String;
    ///                 x_ReadBuffer();
    ///             }
    ///             return *this;
    ///         }
    ///     }
    ///     x_LoadLong();
    ///     return *this;
    /// }
    /// ```
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

    /// `m_String = m_Line; m_Line = m_String;` before the buffer is read again.
    fn keep_line(&mut self) {
        if let Line::Buffer { start, end } = self.line {
            self.string.clear();
            self.string.extend_from_slice(&self.buffer[start..end]);
            self.line = Line::String;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:571-616
    /// ```c++
    /// void CBufferedLineReader::x_LoadLong(void)
    /// {
    ///     const char* start = m_Pos;
    ///     const char* end = m_End;
    ///     m_String.assign(start, end);
    ///     while ( x_ReadBuffer() ) {
    ///         start = m_Pos;
    ///         end = m_End;
    ///         for ( const char* p = start; p < end; ++p ) {
    ///             char c = *p;
    ///             if ( c == '\r' || c == '\n' ) {
    ///                 m_String.append(start, p - start);
    ///                 m_Line = m_String;
    ///                 m_LastReadSize = m_Line.size() + 1;
    ///                 if ( ++p == end ) {
    ///                     m_String = m_Line;
    ///                     m_Line = m_String;
    ///                     if ( x_ReadBuffer() ) {
    ///                         p = m_Pos;
    ///                         end = m_End;
    ///                         if ( p < end && c == '\r' && *p == '\n' ) {
    ///                             ++p;
    ///                             m_Pos = p;
    ///                             ++m_LastReadSize;
    ///                         }
    ///                     }
    ///                 }
    ///                 else {
    ///                     if ( c == '\r' && *p == '\n' ) {
    ///                         if ( ++p == end ) {
    ///                             x_ReadBuffer();
    ///                             p = m_Pos;
    ///                             ++m_LastReadSize;
    ///                         }
    ///                     }
    ///                     m_Pos = p;
    ///                 }
    ///                 return;
    ///             }
    ///         }
    ///         m_String.append(start, end - start);
    ///     }
    ///     m_Line = m_String;
    ///     m_LastReadSize = m_Line.size();
    ///     return;
    /// }
    /// ```
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
    /// ```c++
    /// bool CBufferedLineReader::x_ReadBuffer()
    /// {
    ///     _ASSERT(m_Reader);
    ///
    ///     if ( m_Eof ) {
    ///         return false;
    ///     }
    ///
    ///     m_InputPos += CT_OFF_TYPE(m_End - m_Buffer.get());
    ///     m_Pos = m_End = m_Buffer.get();
    ///     for (bool flag = true; flag; ) {
    ///         size_t size;
    ///         ERW_Result result =
    ///             m_Reader->Read(m_Buffer.get(), m_BufferSize, &size);
    ///         switch (result) {
    ///         case eRW_NotImplemented:
    ///         case eRW_Error:
    ///             NCBI_THROW(CIOException, eRead, "Read error");
    ///             /*NOTREACHED*/
    ///             break;
    ///         case eRW_Timeout:
    ///             // keep spinning around
    ///             break;
    ///         case eRW_Eof:
    ///             m_Eof = true;
    ///             // fall through
    ///         case eRW_Success:
    ///             m_End = m_Pos + size;
    ///             return (result == eRW_Success  ||  size > 0);
    ///         default:
    ///             _ASSERT(0);
    ///         }
    ///     } // for
    ///     return false;
    /// }
    /// ```
    /// `CStreamReader::Read` takes the bytes with `sgetn`, which returns fewer than asked
    /// only at the end of the input, and answers `eRW_Eof` when it got none:
    ///
    /// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:604-622
    /// ```c++
    ///     streambuf* sb = m_Stream->rdbuf();
    ///     bool       ok = sb  &&  m_Stream->good();
    ///     streamsize r  = ok ? sb->sgetn(static_cast<char*>(buf), count) : 0;
    ///     ...
    ///     ERW_Result result = eRW_Success;
    ///     if ( bytes_read ) {
    ///         *bytes_read = (size_t) r;
    ///     } else if ((size_t) r < count) {
    ///         result = eRW_Error;
    ///     }
    ///     if (!r) {
    ///         m_Stream->setstate(!sb  ? NcbiBadbit
    ///                            : ok ? NcbiEofbit : NcbiFailbit);
    ///         result = ok ? eRW_Eof : eRW_Error;
    ///     }
    /// ```
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
    /// ```c++
    /// CTempString CBufferedLineReader::operator*(void) const
    /// {
    ///     _ASSERT(!m_UngetLine);
    ///     /* After UngetLine() - undefined behavior */
    ///     if (m_UngetLine) {
    ///         return CTempString(NULL);
    ///     }
    ///     /* Right after constructor (m_LineNumber is 0 and UngetLine() was not run,
    ///     the latter was already checked) - returns NULL */
    ///     if (m_Line.begin() == NULL) {
    ///         return CTempString(NULL);
    ///     }
    ///     return m_Line;
    /// }
    /// ```
    /// The reader takes the line only after `next_line` (never after `unget_line`).
    pub(super) fn line(&self) -> &[u8] {
        match self.line {
            Line::Null => &[],
            Line::Buffer { start, end } => &self.buffer[start..end],
            Line::String => &self.string,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:682-690
    /// ```c++
    /// Uint8 CBufferedLineReader::GetLineNumber(void) const
    /// {
    ///     /* Right after constructor (m_LineNumber is 0 and UngetLine() was not run)
    ///        - 0 */
    ///     /* If at EOF - returns the number of the last string */
    ///     /* After UngetLine() - number of the previous string */
    ///     /* Not at EOF, not after UngetLine() - number of the current string */
    ///     return m_LineNumber;
    /// }
    /// ```
    pub(super) fn line_number(&self) -> u64 {
        self.line_number
    }
}
