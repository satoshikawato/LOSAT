//! NCBI's `CStreamLineReader` over the opened input stream (`IsIStreamEmpty`, the
//! end-of-line detection, `UngetLine`, `PeekChar` and the line number).
//!
//! The stream and its pushback buffers are copied from BLASTX's port of the same NCBI
//! code (`algorithm/blastx/input.rs`, `FastaStream`), which was checked against an
//! independent C++ oracle of the pinned stream behaviour
//! (`LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`); BLASTX keeps its own copy
//! until session SX.
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
pub(crate) fn stream_is_empty<R: Read + Seek>(stream: &mut FastaStream<R>) -> bool {
    if stream.no_position {
        return false;
    }
    let input = &mut stream.input;
    let position = match input.stream_position() {
        Ok(position) => position,
        Err(_) => return false,
    };
    let mut byte = [0];
    loop {
        match input.read(&mut byte) {
            Ok(0) | Err(_) => return true,
            Ok(_) if crate::ncbi_str::is_space(byte[0]) => {}
            Ok(_) => {
                // IsIStreamEmpty restores position only after successful extraction.
                // NCBI unconditionally clears seekg failure and restores orig_state.
                let _ = input.seek(SeekFrom::Start(position));
                return false;
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:166-174
#[derive(Clone, Copy, PartialEq)]
enum EolStyle {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
struct PushbackBuffer {
    bytes: Vec<u8>,
    position: usize,
    capacity: usize,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
pub(crate) struct FastaStream<R: Read> {
    input: BufReader<R>,
    failed: bool,
    eol: EolStyle,
    pushback: Vec<PushbackBuffer>,
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    backend_available: Option<fn(&mut R) -> std::io::Result<usize>>,
    /// Whether the input's buffered bytes count as readable (`in_avail` of a `filebuf`'s get
    /// area); false for `cin`, whose `stdio_sync_filebuf` has no get area.
    get_area: bool,
    /// Whether lines are copied from the buffered input at once while no pushback buffer
    /// is active (LOSAT's speed). Tests turn it off to compare with one byte at a time.
    pub(crate) bulk: bool,
    /// Whether the input has no stream position, as a pipe has none (a negative `tellg()` in
    /// `IsIStreamEmpty`; `input_files::has_stream_position`): then it is never empty.
    no_position: bool,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
    pub(crate) fn new(input: R) -> Self {
        Self {
            // NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
            // The pinned CNcbiIfstream exposes 8191 readable bytes per default
            // window (independent C++ file_window oracle; tested at EOF/pushback
            // boundaries). Preserve that window because reuse retains EOF state.
            input: BufReader::with_capacity(8191, input),
            failed: false,
            eol: EolStyle::Unknown,
            pushback: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
            backend_available: None,
            get_area: true,
            bulk: true,
            no_position: false,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-358
    fn raw_peek(&mut self) -> Option<u8> {
        while let Some(buffer) = self.pushback.last() {
            if let Some(&byte) = buffer.bytes.get(buffer.position) {
                return Some(byte);
            }
            if self.pushback.len() > 1 {
                let below = self.pushback.remove(self.pushback.len() - 2);
                if below.position < below.bytes.len() {
                    *self.pushback.last_mut().unwrap() = below;
                }
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-338
            let buffered = if self.get_area {
                self.input.buffer().len()
            } else {
                0
            };
            let available = if buffered != 0 {
                buffered
            } else if let Some(in_avail) = self.backend_available {
                in_avail(self.input.get_mut()).unwrap_or(0)
            } else {
                0
            };
            let capacity = self.pushback.last().unwrap().capacity.max(4096);
            let requested = capacity.min(available.max(1));
            let mut bytes = vec![0; requested];
            let mut length = 0;
            while length < requested {
                match self.input.read(&mut bytes[length..]) {
                    Ok(0) | Err(_) => break,
                    Ok(n) => length += n,
                }
            }
            if length == 0 {
                return None;
            }
            bytes.truncate(length);
            let buffer = self.pushback.last_mut().unwrap();
            buffer.capacity = capacity;
            buffer.bytes = bytes;
            buffer.position = 0;
        }
        match self.input.fill_buf() {
            Ok(bytes) => bytes.first().copied(),
            Err(_) => None,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-98
    fn raw_take(&mut self) -> Option<u8> {
        let byte = self.raw_peek()?;
        if let Some(buffer) = self.pushback.last_mut() {
            buffer.position += 1;
        } else {
            self.input.consume(1);
        }
        Some(byte)
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-443
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
    fn push_back(&mut self, bytes: &[u8]) {
        let mut remaining = bytes.len();
        if let Some(buffer) = self.pushback.last_mut() {
            if bytes.len() <= (4096 >> 4) {
                let take = buffer.position.min(remaining);
                buffer.position -= take;
                remaining -= take;
                buffer.bytes[buffer.position..buffer.position + take]
                    .copy_from_slice(&bytes[remaining..]);
            }
        }
        if remaining != 0 {
            self.pushback.push(PushbackBuffer {
                bytes: bytes[..remaining].to_vec(),
                position: 0,
                capacity: remaining,
            });
            // Installing a new streambuf clears state; recycling the old one does not.
            self.failed = false;
        }
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    pub(crate) fn at_eof(&mut self) -> bool {
        if self.failed {
            return true;
        }
        if self.raw_peek().is_none() {
            self.failed = true;
        }
        self.failed
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-166
    fn next_byte(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.raw_take();
        if byte.is_none() {
            self.failed = true;
        }
        byte
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-105
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-174
    fn getline(&mut self, delimiters: &[u8], line: &mut Vec<u8>) -> Option<u8> {
        line.clear();
        loop {
            // Without pushback buffers the bytes come from the buffered input in order, so
            // the bytes before a delimiter are copied at once and the delimiter is taken
            // from the buffer (LOSAT's speed; the bytes, the delimiter taken and the
            // stream state are those of one byte at a time).
            if self.bulk && self.pushback.is_empty() && (delimiters.len() > 1 || !self.failed) {
                let buffered = match self.input.fill_buf() {
                    Ok(bytes) => bytes,
                    Err(_) => &[],
                };
                if !buffered.is_empty() {
                    match find_delimiter(buffered, delimiters) {
                        Some(end) => {
                            let delimiter = buffered[end];
                            line.extend_from_slice(&buffered[..end]);
                            self.input.consume(end + 1);
                            return self.end_of_getline(delimiters, delimiter);
                        }
                        None => {
                            line.extend_from_slice(buffered);
                            let length = buffered.len();
                            self.input.consume(length);
                            continue;
                        }
                    }
                }
            }
            let byte = if delimiters.len() == 1 {
                self.next_byte()
            } else {
                self.raw_take()
            };
            let Some(byte) = byte else {
                self.failed = true;
                return line.last().copied();
            };
            if delimiters.contains(&byte) {
                return self.end_of_getline(delimiters, byte);
            }
            line.push(byte);
        }
    }

    /// `getline` with the one delimiter `eol`, which also gives the first position of
    /// `watch` in the line (the `find` of `alt_eol` in `m_Line` in `x_AdvanceEOLSimple`, found while
    /// the line is read: LOSAT's speed only).
    fn getline_watching(&mut self, eol: u8, watch: u8, line: &mut Vec<u8>) -> Option<usize> {
        line.clear();
        let mut watched = None;
        loop {
            if self.bulk && self.pushback.is_empty() && !self.failed {
                let buffered = match self.input.fill_buf() {
                    Ok(bytes) => bytes,
                    Err(_) => &[],
                };
                if !buffered.is_empty() {
                    let found = if watched.is_none() {
                        find_delimiter(buffered, &[eol, watch])
                    } else {
                        find_byte(buffered, eol)
                    };
                    match found {
                        Some(end) if buffered[end] == eol => {
                            line.extend_from_slice(&buffered[..end]);
                            self.input.consume(end + 1);
                            return watched;
                        }
                        Some(end) => {
                            watched = Some(line.len() + end);
                            line.extend_from_slice(&buffered[..=end]);
                            self.input.consume(end + 1);
                        }
                        None => {
                            line.extend_from_slice(buffered);
                            let length = buffered.len();
                            self.input.consume(length);
                        }
                    }
                    continue;
                }
            }
            let Some(byte) = self.next_byte() else {
                self.failed = true;
                return watched;
            };
            if byte == eol {
                return watched;
            }
            if byte == watch && watched.is_none() {
                watched = Some(line.len());
            }
            line.push(byte);
        }
    }

    /// The end of `getline` at `delimiter`: with two delimiters, a second one that comes
    /// after it in `delimiters` is taken too (CR LF is one end of line, LF CR two). The
    /// last byte taken.
    fn end_of_getline(&mut self, delimiters: &[u8], delimiter: u8) -> Option<u8> {
        let mut last = Some(delimiter);
        if delimiters.len() > 1 {
            let position = delimiters
                .iter()
                .position(|&other| other == delimiter)
                .unwrap_or(delimiters.len());
            if let Some(next) = self.raw_peek() {
                if delimiters[position + 1..].contains(&next) {
                    last = self.raw_take();
                }
            }
        }
        last
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-267
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
    fn advance_simple(&mut self, eol: u8, alternate: u8, line: &mut Vec<u8>) -> EolStyle {
        if let Some(position) = self.getline_watching(eol, alternate, line) {
            let position = position + 1;
            if eol != b'\n' || position != line.len() {
                let rest = line[position..].to_vec();
                self.push_back(&rest);
                self.eol = EolStyle::Mixed;
            }
            line.truncate(position - 1);
            return if self.eol == EolStyle::Mixed {
                EolStyle::Mixed
            } else {
                EolStyle::CrLf
            };
        }
        if eol == b'\r' && !self.failed && self.raw_peek() == Some(alternate) {
            self.raw_take();
            return EolStyle::CrLf;
        }
        if eol == b'\r' {
            EolStyle::Cr
        } else {
            EolStyle::Lf
        }
    }

    /// `PeekChar`'s `m_Stream->peek()`: the next byte, or `None` at the end (which sets
    /// the stream's end-of-file state, as `peek` does).
    pub(crate) fn peek(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.raw_peek();
        if byte.is_none() {
            self.failed = true;
        }
        byte
    }
}
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
// The registered std::ifstream backend reports buffered availability first,
// then remaining regular-file bytes without forcing the next 8191-byte window.
// Independent instrumented pinned CPushback_Streambuf calibrates this hint;
// sgetn may span windows when the preserved allocation is larger than 8191.
fn file_available(file: &mut crate::blastinput::input_files::InputFile) -> std::io::Result<usize> {
    let metadata = file.metadata()?;
    // Windows reports a pipe or a console as a file whose position can be read
    // (`input_files::has_stream_position`); elsewhere `stream_position` below fails on one.
    if metadata.is_file()
        && (cfg!(not(windows)) || crate::blastinput::input_files::has_stream_position(file))
    {
        let position = file.stream_position()?;
        Ok(metadata.len().saturating_sub(position) as usize)
    } else {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // The registered native CNcbiIfstream exposes queued pipe bytes via
        // FIONREAD before forcing a file-buffer refill (independent C++ oracle).
        #[cfg(unix)]
        {
            rustix::io::ioctl_fionread(&**file)
                .map(|n| n as usize)
                .map_err(Into::into)
        }
        #[cfg(not(unix))]
        {
            Ok(0)
        }
    }
}
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
// Bytes in memory are read as a regular file with the same bytes (`file_available`): the
// readable-byte hint is the number of bytes not yet taken from them.
fn bytes_available(rest: &mut &[u8]) -> std::io::Result<usize> {
    Ok(rest.len())
}

impl<'a> FastaStream<&'a [u8]> {
    /// A stream over bytes in memory that reads them as `from_file` reads a regular file
    /// with the same bytes (the same refills of the pushback buffers, so the same lines).
    pub(crate) fn from_bytes(input: &'a [u8]) -> Self {
        let mut stream = Self::new(input);
        stream.backend_available = Some(bytes_available);
        stream
    }
}

impl FastaStream<crate::blastinput::input_files::InputFile> {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    /// A stream over an opened input, whose refills read as NCBI's stream of it answers them
    /// (`InputBackend::of`: `cin` for `-`, a `CNcbiIfstream` for a path).
    pub(crate) fn from_file(input: impl Into<crate::blastinput::input_files::InputFile>) -> Self {
        let input = input.into();
        let no_position = !crate::blastinput::input_files::has_stream_position(&input);
        let backend = InputBackend::of(&input);
        let mut stream = Self::new(input);
        stream.no_position = no_position;
        stream.backend_available = Some(backend.available);
        stream.get_area = backend.get_area;
        stream
    }
}

/// How NCBI's stream over an opened input answers the refill of a pushback buffer, which asks
/// for `in_avail()` bytes, at least one (`CPushback_Streambuf::underflow`, `x_FillBuffer`, quoted
/// at `raw_peek`). The one place that tells standard input from a file for every FASTA reader
/// of an opened input (`FastaStream::from_file` here and BLASTX's and RPSTBLASTN's
/// `algorithm::blastx::input::FastaStream::from_file`).
#[derive(Clone, Copy)]
pub(crate) struct InputBackend {
    /// `in_avail()` once the get area is used up (`backend_available`).
    pub(crate) available:
        fn(&mut crate::blastinput::input_files::InputFile) -> std::io::Result<usize>,
    /// Whether the stream's buffered bytes count as readable (`get_area`).
    pub(crate) get_area: bool,
}

impl InputBackend {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1342-1344
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    /// The backend of `input`. Standard input (`-`, `InputFile::is_standard_input`) is `cin`,
    /// which the BLAST applications keep synchronised with stdio: its streambuf (libstdc++'s
    /// `stdio_sync_filebuf`) has no get area and the default `showmanyc` of 0, so `in_avail` is
    /// always 0 and each refill reads one byte (`x_FillBuffer` raises 0 to 1), whether standard
    /// input is a pipe or a regular file. A file named by its path is a `CNcbiIfstream`, whose
    /// buffered bytes and then the bytes left in the file (a pipe's queued bytes) are readable
    /// (`file_available`). The two read the same bytes into different lines where a refill ends
    /// inside a line end (fixture rows `stdin_q_lost_g55_pipe` and `stdin_q_lost_g55_file` keep
    /// the tail that the same file named by its path loses).
    pub(crate) fn of(input: &crate::blastinput::input_files::InputFile) -> Self {
        if input.is_standard_input() {
            Self {
                available: |_| Ok(0),
                get_area: false,
            }
        } else {
            Self {
                available: file_available,
                get_area: true,
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-241
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    /// The switch of `operator++` on the end-of-line style (the caller has checked
    /// `AtEOF` and the ungot line): the next line into `line`, without its end of line.
    pub(crate) fn advance(&mut self, line: &mut Vec<u8>) {
        match self.eol {
            EolStyle::Unknown => {
                let last = self.getline(b"\r\n", line);
                // Successful unget clears eofbit before get rereads the last byte.
                if self.failed && !line.is_empty() {
                    self.failed = false;
                }
                match last {
                    Some(b'\r') => self.eol = EolStyle::Cr,
                    Some(b'\n') => self.eol = EolStyle::CrLf,
                    _ => {}
                }
            }
            EolStyle::Cr => {
                self.advance_simple(b'\r', b'\n', line);
            }
            EolStyle::Lf => {
                self.advance_simple(b'\n', b'\r', line);
            }
            EolStyle::CrLf => {
                let style = self.advance_simple(b'\n', b'\r', line);
                if style == EolStyle::Mixed {
                    self.eol = EolStyle::Cr;
                } else if style != EolStyle::CrLf {
                    self.eol = EolStyle::Lf;
                }
            }
            EolStyle::Mixed => {
                self.getline(b"\r\n", line);
            }
        }
    }
}
/// The bytes of a word that are `byte`, as the high bit of each byte of the word. The
/// lowest set bit marks the first such byte exactly; bits above it can be set for bytes
/// that are not `byte`, so only the lowest one is used.
fn equal_bytes(word: u64, byte: u8) -> u64 {
    const LOW: u64 = 0x0101_0101_0101_0101;
    const HIGH: u64 = 0x8080_8080_8080_8080;
    let difference = word ^ (LOW * u64::from(byte));
    difference.wrapping_sub(LOW) & !difference & HIGH
}

/// The first position of `needle` in `haystack` (`memchr`, a word at a time: LOSAT's
/// speed only).
pub(crate) fn find_byte(haystack: &[u8], needle: u8) -> Option<usize> {
    let mut words = haystack.chunks_exact(8);
    let mut offset = 0;
    for word in &mut words {
        let found = equal_bytes(u64::from_le_bytes(word.try_into().unwrap()), needle);
        if found != 0 {
            return Some(offset + (found.trailing_zeros() / 8) as usize);
        }
        offset += 8;
    }
    words
        .remainder()
        .iter()
        .position(|&byte| byte == needle)
        .map(|position| offset + position)
}

/// The first position in `haystack` of one of the one or two `delimiters`.
fn find_delimiter(haystack: &[u8], delimiters: &[u8]) -> Option<usize> {
    match *delimiters {
        [only] => find_byte(haystack, only),
        [first, second] => {
            let mut words = haystack.chunks_exact(8);
            let mut offset = 0;
            for word in &mut words {
                let word = u64::from_le_bytes(word.try_into().unwrap());
                let found = equal_bytes(word, first) | equal_bytes(word, second);
                if found != 0 {
                    return Some(offset + (found.trailing_zeros() / 8) as usize);
                }
                offset += 8;
            }
            words
                .remainder()
                .iter()
                .position(|&byte| byte == first || byte == second)
                .map(|position| offset + position)
        }
        _ => haystack.iter().position(|byte| delimiters.contains(byte)),
    }
}

/// NCBI's `CStreamLineReader` over a `FastaStream`: the current line, the line number
/// and the ungot line.
///
/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:86-92
pub(crate) struct LineReader<R: Read> {
    stream: FastaStream<R>,
    line: Vec<u8>,
    line_number: u64,
    unget_line: bool,
}

impl<R: Read + Seek> LineReader<R> {
    /// `IsIStreamEmpty` on the stream under the line reader, before its first line is
    /// read (`stream_is_empty`).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:845-847
    pub(crate) fn stream_is_empty(&mut self) -> bool {
        debug_assert!(self.line_number == 0 && !self.unget_line);
        stream_is_empty(&mut self.stream)
    }
}

impl<R: Read> LineReader<R> {
    pub(crate) fn new(stream: FastaStream<R>) -> Self {
        Self {
            stream,
            line: Vec::new(),
            line_number: 0,
            unget_line: false,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    pub(crate) fn at_eof(&mut self) -> bool {
        !self.unget_line && self.stream.at_eof()
    }

    /// The first byte of the next line, `0` for an empty line (`None` at the end).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:107-140
    pub(crate) fn peek_char(&mut self) -> Option<u8> {
        if self.at_eof() {
            return self.stream.peek();
        }
        if self.unget_line {
            return Some(self.line.first().copied().unwrap_or(0));
        }
        match self.stream.peek() {
            Some(b'\n' | b'\r') => Some(0),
            other => other,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:142-151
    pub(crate) fn unget_line(&mut self) {
        if self.unget_line || self.line_number == 0 {
            return;
        }
        self.line_number -= 1;
        self.unget_line = true;
    }

    /// An increment of the reader followed by its dereference: the next line (the ungot line again after
    /// `UngetLine`, an empty line at the end).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-175
    pub(crate) fn next_line(&mut self) -> &[u8] {
        if self.at_eof() {
            self.line.clear();
            return &self.line;
        }
        self.line_number += 1;
        if self.unget_line {
            self.unget_line = false;
            return &self.line;
        }
        self.stream.advance(&mut self.line);
        &self.line
    }

    /// Takes the current line out (`*reader`), to be read while the reader goes on with
    /// other work; `put_line` gives it back before the next line, `UngetLine` or
    /// `PeekChar` (LOSAT's speed: no copy of the line).
    pub(crate) fn take_line(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.line)
    }

    /// Gives back the line of `take_line`.
    pub(crate) fn put_line(&mut self, line: Vec<u8>) {
        self.line = line;
    }

    /// `FastaStream::bulk`.
    pub(crate) fn bulk(&self) -> bool {
        self.stream.bulk
    }

    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:208-216
    pub(crate) fn line_number(&self) -> u64 {
        self.line_number
    }

    /// `FastaReader::x_start_at`: a line reader that has read nothing counts its lines from `lines`.
    /// No NCBI counterpart: `CStreamLineReader` starts with `m_LineNumber` at 0 (line_reader.cpp:86-92,
    /// `LineReader`); a reader of a part of the file starts where the file's one reader is.
    #[cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]
    pub(crate) fn x_start_at(&mut self, lines: u64) {
        debug_assert!(self.line_number == 0 && !self.unget_line);
        self.line_number = lines;
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
