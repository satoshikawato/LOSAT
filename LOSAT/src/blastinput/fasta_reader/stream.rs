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
// ```c++
// 	char c;
// 	CNcbiStreampos orig_p = in.tellg();
// 	// Piped input
// 	if(orig_p < 0)
// 		return false;
//
// 	IOS_BASE::iostate orig_state = in.rdstate();
// 	IOS_BASE::fmtflags orig_flags = in.setf(ios::skipws);
//
// 	if(! (in >> c))
// 		return true;
//
// 	in.seekg(orig_p);
// 	in.flags(orig_flags);
// 	in.clear();
// 	in.setstate(orig_state);
//
// 	return false;
// ```
pub(crate) fn stream_is_empty<R: Read + Seek>(stream: &mut FastaStream<R>) -> bool {
    let input = &mut stream.input;
    let position = match input.stream_position() {
        Ok(position) => position,
        Err(_) => return false,
    };
    let mut byte = [0];
    loop {
        match input.read(&mut byte) {
            Ok(0) | Err(_) => return true,
            Ok(_) if input_space(byte[0]) => {}
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
// ```c++
//
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
//     return *this;
// ```
#[derive(Clone, Copy, PartialEq)]
enum EolStyle {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
// ```c++
//     if (!del_ptr  &&  how != ePushback_NoCopy) {
//         del_ptr = new CT_CHAR_TYPE[buf_size];
//         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
//     }
//
//     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
// ```c++
// CPushback_Streambuf::CPushback_Streambuf(istream&      is,
//                                          CT_CHAR_TYPE* buf,
//                                          size_t        buf_size,
//                                          void*         del_ptr)
//     : m_Is(is), m_Next(0), m_Buf(buf), m_BufSize(buf_size), m_DelPtr(del_ptr)
// {
//     _ASSERT(m_Buf  &&  m_BufSize);
//     setp(0, 0);  // unbuffered output at this level of streambuf's hierarchy
//     setg(m_Buf, m_Buf, m_Buf + m_BufSize);
//     m_Sb = m_Is.rdbuf(this);
// ```
struct PushbackBuffer {
    bytes: Vec<u8>,
    position: usize,
    capacity: usize,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
pub(crate) struct FastaStream<R: Read> {
    input: BufReader<R>,
    failed: bool,
    eol: EolStyle,
    pushback: Vec<PushbackBuffer>,
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    // ```c++
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    backend_available: Option<fn(&mut R) -> std::io::Result<usize>>,
    /// Whether the input's buffered bytes count as readable (`in_avail` of a `filebuf`'s get
    /// area); false for `cin`, whose `stdio_sync_filebuf` has no get area.
    get_area: bool,
    /// Whether lines are copied from the buffered input at once while no pushback buffer
    /// is active (LOSAT's speed). Tests turn it off to compare with one byte at a time.
    pub(crate) bulk: bool,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// ```c++
//       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
//       m_EOLStyle(eol_style)
// {
// }
//
//
// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
//                                      EOwnership ownership)
//     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
// ```
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
    // ```c++
    //       m_UngetLine(false), m_AutoEOL(eol_style == eEOL_unknown),
    //       m_EOLStyle(eol_style)
    // {
    // }
    //
    //
    // CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
    //                                      EOwnership ownership)
    //     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
    // ```
    pub(crate) fn new(input: R) -> Self {
        Self {
            // NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
            // ```c++
            // #else
            // /// Portable alias for ifstream.
            // typedef IO_PREFIX::ifstream      CNcbiIfstream;
            // #endif
            // ```
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
            // ```c++
            // CT_INT_TYPE CPushback_Streambuf::underflow(void)
            // {
            //     // we are here because there is no more data in the pushback buffer
            //     _ASSERT(gptr()  &&  gptr() >= egptr());
            //
            // #ifdef NCBI_COMPILER_MIPSPRO
            //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
            //         return CT_EOF;
            //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
            // #endif //NCBI_COMPILER_MIPSPRO
            //
            //     x_FillBuffer((size_t) m_Sb->in_avail());
            //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
            // ```
            // The pinned CNcbiIfstream exposes 8191 readable bytes per default
            // window (independent C++ file_window oracle; tested at EOF/pushback
            // boundaries). Preserve that window because reuse retains EOF state.
            input: BufReader::with_capacity(8191, input),
            failed: false,
            eol: EolStyle::Unknown,
            pushback: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
            // ```c++
            //     x_FillBuffer((size_t) m_Sb->in_avail());
            //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
            // ```
            backend_available: None,
            get_area: true,
            bulk: true,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // ```c++
    // CT_INT_TYPE CPushback_Streambuf::underflow(void)
    // {
    //     // we are here because there is no more data in the pushback buffer
    //     _ASSERT(gptr()  &&  gptr() >= egptr());
    //
    // #ifdef NCBI_COMPILER_MIPSPRO
    //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
    //         return CT_EOF;
    //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
    // #endif //NCBI_COMPILER_MIPSPRO
    //
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-358
    // ```c++
    // void CPushback_Streambuf::x_FillBuffer(size_t max_size)
    // {
    //     _ASSERT(m_Sb);
    //     if ( !max_size ) {
    //         ++max_size;
    //     }
    //
    //     CPushback_Streambuf* sb = dynamic_cast<CPushback_Streambuf*> (m_Sb);
    //     if ( !sb ) {
    //         CT_CHAR_TYPE* bp = 0;
    //         size_t buf_size = m_DelPtr
    //             ? (size_t)(m_Buf - (CT_CHAR_TYPE*) m_DelPtr) + m_BufSize : 0;
    //         if (buf_size < kMinBufSize) {
    //             buf_size = kMinBufSize;
    //             bp = new CT_CHAR_TYPE[buf_size];
    //         }
    //         streamsize r = (streamsize)(buf_size < max_size ? buf_size : max_size);
    //         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
    //         if (n <= 0) {
    //             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
    //             delete[] bp;
    //             return;
    //         }
    //         if (bp) {
    //             delete[] (CT_CHAR_TYPE*) m_DelPtr;
    //             m_DelPtr = bp;
    //         }
    //         m_Buf = (CT_CHAR_TYPE*) m_DelPtr;
    //         m_BufSize = buf_size;
    //         setg(m_Buf, m_Buf, m_Buf + n);
    //         return;
    //     }
    //
    //     _ASSERT(&m_Is  == &sb->m_Is);
    //     _ASSERT(m_Next == sb);
    //     m_Sb       = sb->m_Sb;
    //     m_Next     = sb->m_Next;
    //     sb->m_Sb   = 0;
    //     sb->m_Next = 0;
    //     if (sb->gptr() >= sb->egptr()) {
    //         delete sb;
    //         x_FillBuffer(max_size);
    //         return;
    //     }
    //     delete[] (CT_CHAR_TYPE*) m_DelPtr;
    //     m_Buf        = sb->m_Buf;
    //     m_BufSize    = sb->m_BufSize;
    //     m_DelPtr     = sb->m_DelPtr;
    //     sb->m_DelPtr = 0;
    //     setg(sb->gptr(), sb->gptr(), sb->egptr());
    //     delete sb;
    // }
    //
    // ```
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
            // ```c++
            // void CPushback_Streambuf::x_FillBuffer(size_t max_size)
            // {
            //     _ASSERT(m_Sb);
            //     if ( !max_size ) {
            //         ++max_size;
            //     }
            //
            //     CPushback_Streambuf* sb = dynamic_cast<CPushback_Streambuf*> (m_Sb);
            //     if ( !sb ) {
            //         CT_CHAR_TYPE* bp = 0;
            //         size_t buf_size = m_DelPtr
            //             ? (size_t)(m_Buf - (CT_CHAR_TYPE*) m_DelPtr) + m_BufSize : 0;
            //         if (buf_size < kMinBufSize) {
            //             buf_size = kMinBufSize;
            //             bp = new CT_CHAR_TYPE[buf_size];
            //         }
            //         streamsize r = (streamsize)(buf_size < max_size ? buf_size : max_size);
            //         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
            //         if (n <= 0) {
            //             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
            //             delete[] bp;
            //             return;
            //         }
            //         if (bp) {
            //             delete[] (CT_CHAR_TYPE*) m_DelPtr;
            //             m_DelPtr = bp;
            //         }
            //         m_Buf = (CT_CHAR_TYPE*) m_DelPtr;
            //         m_BufSize = buf_size;
            //         setg(m_Buf, m_Buf, m_Buf + n);
            //         return;
            //     }
            //
            // ```
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
    // ```c++
    //             iostate = NcbiEofbit;
    //             break;
    //         }
    //         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
    //         if (delim_pos != NPOS) {
    //             // Special case -- if two different delimiters are back to
    //             // back and in the same order as in delims, treat them as
    //             // a single delimiter (necessary for correct handling of
    //             // DOS/MAC-style CR/LF endings).
    //             ch = is.rdbuf()->sgetc();
    //             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
    //                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
    // ```
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
    // ```c++
    //         if (how == ePushback_Stepback
    //             ||  (how == ePushback_Copy
    //                  &&  buf_size <= (del_ptr
    //                                   ? CPushback_Streambuf::kMinBufSize
    //                                   : CPushback_Streambuf::kMinBufSize >> 4))) {
    //             CT_CHAR_TYPE* bp = sb->gptr();
    //             size_t avail = bp - sb->m_Buf;
    //             size_t take  = avail < buf_size ? avail : buf_size;
    //             if (take) {
    //                 bp -= take;
    //                 buf_size -= take;
    //                 if (how != ePushback_Stepback  &&  bp != buf + buf_size) {
    //                     memmove(bp, buf + buf_size, take);
    //                 }
    //                 sb->setg(bp, bp, sb->egptr());
    //             }
    //         }
    //     }
    //
    //     if ( !buf_size ) {
    //         delete[] (CT_CHAR_TYPE*) del_ptr;
    //         return;
    //     }
    //
    //     if (!del_ptr  &&  how != ePushback_NoCopy) {
    //         del_ptr = new CT_CHAR_TYPE[buf_size];
    //         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
    //     }
    //
    //     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
    // ```c++
    // CPushback_Streambuf::CPushback_Streambuf(istream&      is,
    //                                          CT_CHAR_TYPE* buf,
    //                                          size_t        buf_size,
    //                                          void*         del_ptr)
    //     : m_Is(is), m_Next(0), m_Buf(buf), m_BufSize(buf_size), m_DelPtr(del_ptr)
    // {
    //     _ASSERT(m_Buf  &&  m_BufSize);
    //     setp(0, 0);  // unbuffered output at this level of streambuf's hierarchy
    //     setg(m_Buf, m_Buf, m_Buf + m_BufSize);
    //     m_Sb = m_Is.rdbuf(this);
    // ```
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
    // ```c++
    // bool CStreamLineReader::AtEOF(void) const
    // {
    //     return !m_UngetLine &&
    //         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    // }
    // ```
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
    // ```c++
    //     SIZE_TYPE size = 0;
    //     SIZE_TYPE max_size = str.max_size();
    //     do {
    //         CT_INT_TYPE nextc = is.get();
    //         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
    //             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
    //             ++size;
    //             break;
    //         }
    //         if ( !is.unget() )
    //             break;
    //         if (size == max_size) {
    //             is.clear(NcbiFailbit);
    //             break;
    //         }
    //         SIZE_TYPE n = max_size - size;
    // ```
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
    // ```c++
    //             iostate = NcbiEofbit;
    //             break;
    //         }
    //         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
    //         if (delim_pos != NPOS) {
    //             // Special case -- if two different delimiters are back to
    //             // back and in the same order as in delims, treat them as
    //             // a single delimiter (necessary for correct handling of
    //             // DOS/MAC-style CR/LF endings).
    //             ch = is.rdbuf()->sgetc();
    //             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
    //                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
    //                 is.rdbuf()->sbumpc();
    //                 delim_count = 2;
    //             } else {
    //                 delim_count = 1;
    //             }
    //             break;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-174
    // ```c++
    //     SIZE_TYPE size = 0;
    //     SIZE_TYPE max_size = str.max_size();
    //     do {
    //         CT_INT_TYPE nextc = is.get();
    //         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
    //             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
    //             ++size;
    //             break;
    //         }
    //         if ( !is.unget() )
    //             break;
    //         if (size == max_size) {
    //             is.clear(NcbiFailbit);
    //             break;
    //         }
    //         SIZE_TYPE n = max_size - size;
    //         is.get(buf, n < sizeof(buf) ? n : sizeof(buf), delim);
    //         n = (size_t) is.gcount();
    //         str.append(buf, n);
    //         size += n;
    //         _ASSERT(size == str.length());
    //     } while ( is.good() );
    // #endif
    //
    // ```
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
    /// `watch` in the line (the `m_Line.find(alt_eol)` of `x_AdvanceEOLSimple`, found while
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
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLSimple(char eol,
    //                                                                    char alt_eol)
    // {
    //     SIZE_TYPE pos;
    //     NcbiGetline(*m_Stream, m_Line, eol, &m_LastReadSize);
    //     if (m_AutoEOL  &&  (pos = m_Line.find(alt_eol)) != NPOS) {
    //         ++pos;
    //         if (eol != '\n'  ||  pos != m_Line.size()) {
    //             // an *immediately* preceding CR is quite all right
    //             CStreamUtils::Pushback(*m_Stream, m_Line.data() + pos,
    //                                    m_Line.size() - pos);
    //             m_EOLStyle = eEOL_mixed;
    //         }
    //         m_Line.resize(pos - 1);
    //         m_LastReadSize = pos;
    //         return (m_EOLStyle == eEOL_mixed) ? m_EOLStyle : eEOL_crlf;
    //     } else if (m_AutoEOL  &&  eol == '\r'  &&
    //                CT_EQ_INT_TYPE(m_Stream->peek(), CT_TO_INT_TYPE(alt_eol))) {
    //         m_Stream->get();
    //         ++m_LastReadSize;
    //         return eEOL_crlf;
    //     }
    //     return (eol == '\r') ? eEOL_cr : eEOL_lf;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
    // ```c++
    //     if (!del_ptr  &&  how != ePushback_NoCopy) {
    //         del_ptr = new CT_CHAR_TYPE[buf_size];
    //         buf = (CT_CHAR_TYPE*) memcpy(del_ptr, buf, buf_size);
    //     }
    //
    //     (void) new CPushback_Streambuf(is, buf, buf_size, del_ptr);
    // }
    // ```
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
// ```c++
// #else
// /// Portable alias for ifstream.
// typedef IO_PREFIX::ifstream      CNcbiIfstream;
// #endif
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
// ```c++
// CT_INT_TYPE CPushback_Streambuf::underflow(void)
// {
//     // we are here because there is no more data in the pushback buffer
//     _ASSERT(gptr()  &&  gptr() >= egptr());
//
// #ifdef NCBI_COMPILER_MIPSPRO
//     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
//         return CT_EOF;
//     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
// #endif //NCBI_COMPILER_MIPSPRO
//
//     x_FillBuffer((size_t) m_Sb->in_avail());
//     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
// ```
// The registered std::ifstream backend reports buffered availability first,
// then remaining regular-file bytes without forcing the next 8191-byte window.
// Independent instrumented pinned CPushback_Streambuf calibrates this hint;
// sgetn may span windows when the preserved allocation is larger than 8191.
fn file_available(file: &mut std::fs::File) -> std::io::Result<usize> {
    let metadata = file.metadata()?;
    if metadata.is_file() {
        let position = file.stream_position()?;
        Ok(metadata.len().saturating_sub(position) as usize)
    } else {
        // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
        // ```c++
        //     x_FillBuffer((size_t) m_Sb->in_avail());
        //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
        // ```
        // The registered native CNcbiIfstream exposes queued pipe bytes via
        // FIONREAD before forcing a file-buffer refill (independent C++ oracle).
        #[cfg(unix)]
        {
            rustix::io::ioctl_fionread(file)
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
// ```c++
// #else
// /// Portable alias for ifstream.
// typedef IO_PREFIX::ifstream      CNcbiIfstream;
// #endif
// ```
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
// ```c++
//     x_FillBuffer((size_t) m_Sb->in_avail());
//     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
// ```
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

impl FastaStream<std::fs::File> {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // ```c++
    // CT_INT_TYPE CPushback_Streambuf::underflow(void)
    // {
    //     // we are here because there is no more data in the pushback buffer
    //     _ASSERT(gptr()  &&  gptr() >= egptr());
    //
    // #ifdef NCBI_COMPILER_MIPSPRO
    //     if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
    //         return CT_EOF;
    //     m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
    // #endif //NCBI_COMPILER_MIPSPRO
    //
    //     x_FillBuffer((size_t) m_Sb->in_avail());
    //     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    // ```
    pub(crate) fn from_file(input: std::fs::File) -> Self {
        let mut stream = Self::new(input);
        stream.backend_available = Some(file_available);
        stream
    }

    /// A stream over standard input (`-`), which NCBI reads through `cin`.
    ///
    /// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
    /// ```c++
    ///     if (AsString() == "-") {
    /// #if defined(NCBI_OS_MSWIN)
    ///         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
    /// #endif
    ///         m_Ios  = &cin;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1342-1344
    /// ```c++
    ///     if ((m_StdioFlags & fNoSyncWithStdio) != 0) {
    ///         IOS_BASE::sync_with_stdio(false);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    /// ```c++
    ///     x_FillBuffer((size_t) m_Sb->in_avail());
    ///     return gptr() < egptr() ? CT_TO_INT_TYPE(*gptr()) : CT_EOF;
    /// ```
    /// The BLAST applications keep `cin` synchronised with stdio, so its streambuf (libstdc++'s
    /// `stdio_sync_filebuf`) has no get area and the default `showmanyc` of 0: `in_avail` is
    /// always 0, and each refill of a pushback buffer reads one byte (`x_FillBuffer` raises 0
    /// to 1), whether standard input is a pipe or a regular file (fixture rows
    /// `stdin_q_lost_g55_pipe` and `stdin_q_lost_g55_file` keep the tail that a file named by
    /// its path loses).
    pub(crate) fn from_standard_input(input: std::fs::File) -> Self {
        let mut stream = Self::new(input);
        stream.backend_available = Some(|_| Ok(0));
        stream.get_area = false;
        stream
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
// ```c++
// CStreamLineReader& CStreamLineReader::operator++(void)
// {
//     /* If at EOF - noop */
//     if (AtEOF()) {
//         m_Line = string();
//         return *this;
//     }
//     ++m_LineNumber;
//     if ( m_UngetLine ) {
//         m_UngetLine = false;
//         return *this;
//     }
//
//     switch (m_EOLStyle) {
//     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
//     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
//     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
//     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
//     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
//     }
//     return *this;
// ```
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-241
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLUnknown(void)
    // {
    //     _ASSERT(m_AutoEOL);
    //     NcbiGetline(*m_Stream, m_Line, "\r\n", &m_LastReadSize);
    //     m_Stream->unget();
    //     CT_INT_TYPE eol = m_Stream->get();
    //     if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\r'))) {
    //         m_EOLStyle = eEOL_cr;
    //     } else if (CT_EQ_INT_TYPE(eol, CT_TO_INT_TYPE('\n'))) {
    //         // NcbiGetline doesn't yield enough information to determine
    //         // whether eEOL_lf or eEOL_crlf is more appropriate, and not
    //         // all streams allow tellg() (which could otherwise resolve
    //         // matters), so defer further analysis to x_AdvanceEOLCRLF,
    //         // which will be responsible for reading the next line and
    //         // supports switching to eEOL_lf as appropriate.
    //         //
    //         // An alternative approach would have been to pass \n\r rather
    //         // than \r\n, and then check for an immediately following \n
    //         // if eol turned out to be \r, but that would miscount an
    //         // actual(!) \n\r sequence as a single line break.
    //         m_EOLStyle = eEOL_crlf;
    //     }
    //     return m_EOLStyle;
    // ```
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    // ```c++
    // CStreamLineReader::EEOLStyle CStreamLineReader::x_AdvanceEOLCRLF(void)
    // {
    //     if (m_AutoEOL) {
    //         EEOLStyle style = x_AdvanceEOLSimple('\n', '\r');
    //         if (style == eEOL_mixed) {
    //             // found an embedded CR
    //             m_EOLStyle = eEOL_cr;
    //         } else if (style != eEOL_crlf) {
    //             m_EOLStyle = eEOL_lf;
    //         }
    //     } else {
    // ```
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

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:930-931
// ```c++
//         case '\t': case '\n': case '\v': case '\f': case '\r': case ' ':
//             continue;
// ```
pub(crate) const fn input_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3153-3161
// ```c++
//     SIZE_TYPE beg = 0;
//     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
//         _ASSERT(beg < length);
//         while ( isspace((unsigned char) str[beg]) ) {
//             if (++beg == length) {
//                 return empty_str;
//             }
//         }
//     }
// ```
pub(crate) fn trim_input_start(bytes: &[u8]) -> &[u8] {
    let start = bytes
        .iter()
        .position(|&byte| !input_space(byte))
        .unwrap_or(bytes.len());
    &bytes[start..]
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3162-3172
// ```c++
//     SIZE_TYPE end = length;
//     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
//         _ASSERT(beg < end);
//         while (isspace((unsigned char) str[--end])) {
//             if (beg == end) {
//                 return empty_str;
//             }
//         }
//         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
//         ++end;
//     }
// ```
pub(crate) fn trim_input_end(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .rposition(|&byte| !input_space(byte))
        .map_or(0, |p| p + 1);
    &bytes[..end]
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3187-3190
// ```c++
// CTempString NStr::TruncateSpaces_Unsafe(const CTempString str, ETrunc where)
// {
//     return s_TruncateSpaces(str, where, CTempString());
// }
// ```
pub(crate) fn trim_input_space(bytes: &[u8]) -> &[u8] {
    trim_input_end(trim_input_start(bytes))
}

/// NCBI's `CStreamLineReader` over a `FastaStream`: the current line, the line number
/// and the ungot line.
///
/// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:86-92
/// ```c++
/// CStreamLineReader::CStreamLineReader(CNcbiIstream& is,
///                                      EOwnership ownership)
///     : m_Stream(&is, ownership), m_LineNumber(0), m_LastReadSize(0),
///       m_UngetLine(false), m_AutoEOL(true), m_EOLStyle(eEOL_unknown)
/// {
/// }
/// ```
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
    /// ```c++
    /// bool
    /// IsIStreamEmpty(CNcbiIstream & in)
    /// {
    /// ```
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
    /// ```c++
    /// bool CStreamLineReader::AtEOF(void) const
    /// {
    ///     return !m_UngetLine &&
    ///         (m_Stream->eof()  ||  CT_EQ_INT_TYPE(m_Stream->peek(), CT_EOF));
    /// }
    /// ```
    pub(crate) fn at_eof(&mut self) -> bool {
        !self.unget_line && self.stream.at_eof()
    }

    /// The first byte of the next line, `0` for an empty line (`None` at the end).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:107-140
    /// ```c++
    /// char CStreamLineReader::PeekChar(void) const
    /// {
    ///     _ASSERT(!AtEOF());
    ///     /* If at EOF - undefined behavior, return m_Stream->peek() */
    ///     if (AtEOF()) {
    ///         return (char)m_Stream->peek();
    ///     }
    ///     /* If right after constructor (line number is 0 and line was not ungot) -
    ///        return the first character of the first line */
    ///     if (m_LineNumber == 0 && !m_UngetLine) {
    ///         char c = (char)m_Stream->peek();
    ///         /* If there are delimiters right from the start (which means the
    ///         first line is empty), return 0 */
    ///         if (c == '\n' || c == '\r') {
    ///             return 0;
    ///         }
    ///         return c;
    ///     }
    ///     /* If line was ungot - return its first symbol */
    ///     if (m_UngetLine) {
    ///         /* if line is empty - return 0 */
    ///         if (m_Line.empty()) {
    ///             return 0;
    ///         }
    ///         return *m_Line.begin();
    ///     }
    ///     char c = (char)m_Stream->peek();
    ///     /* If line is empty - return 0 */
    ///     if (c == '\n' || c == '\r') {
    ///         return 0;
    ///     }
    ///     return c;
    /// }
    /// ```
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
    /// ```c++
    /// void CStreamLineReader::UngetLine(void)
    /// {
    ///     _ASSERT(!m_UngetLine && m_LineNumber != 0);
    ///     /* If after UngetLine() or after constructor - noop */
    ///     if (m_UngetLine || m_LineNumber == 0) {
    ///         return;
    ///     }
    ///     --m_LineNumber;
    ///     m_UngetLine = true;
    /// }
    /// ```
    pub(crate) fn unget_line(&mut self) {
        if self.unget_line || self.line_number == 0 {
            return;
        }
        self.line_number -= 1;
        self.unget_line = true;
    }

    /// `++reader` followed by `*reader`: the next line (the ungot line again after
    /// `UngetLine`, an empty line at the end).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-175
    /// ```c++
    /// CStreamLineReader& CStreamLineReader::operator++(void)
    /// {
    ///     /* If at EOF - noop */
    ///     if (AtEOF()) {
    ///         m_Line = string();
    ///         return *this;
    ///     }
    ///     ++m_LineNumber;
    ///     if ( m_UngetLine ) {
    ///         m_UngetLine = false;
    ///         return *this;
    ///     }
    ///
    ///     switch (m_EOLStyle) {
    ///     case eEOL_unknown: x_AdvanceEOLUnknown();                   break;
    ///     case eEOL_cr:      x_AdvanceEOLSimple('\r', '\n');          break;
    ///     case eEOL_lf:      x_AdvanceEOLSimple('\n', '\r');          break;
    ///     case eEOL_crlf:    x_AdvanceEOLCRLF();                      break;
    ///     case eEOL_mixed:   NcbiGetline(*m_Stream, m_Line, "\r\n");  break;
    ///     }
    ///     return *this;
    /// }
    /// ```
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
    /// ```c++
    /// Uint8 CStreamLineReader::GetLineNumber(void) const
    /// {
    ///     /* Right after constructor (m_LineNumber is 0 and UngetLine() was not run)
    ///        - 0 */
    ///     /* If at EOF - returns the number of the last string */
    ///     /* After UngetLine() - number of the previous string */
    ///     /* Not at EOF, not after UngetLine() - number of the current string */
    ///     return m_LineNumber;
    /// }
    /// ```
    pub(crate) fn line_number(&self) -> u64 {
        self.line_number
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
