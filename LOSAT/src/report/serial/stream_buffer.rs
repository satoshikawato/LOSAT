//! NCBI's output buffer of an object stream (`COStreamBuffer`, c++/src/util/strbuffer.cpp), as
//! far as the order of a stream's bytes and NCBI's messages goes: the bytes of a document that
//! the buffer has handed to its `CNcbiOstream` when NCBI posts a message (a replaced VisibleString
//! byte, `visible::BadChars`). NCBI's diagnostics go to `cerr`, which is tied to `cout`, so a
//! post first writes what `cout` holds: on a stream that both go to (`2>&1`) the message follows
//! the bytes this buffer has handed out, and precedes the rest of the document.
//!
//! This is not `cli::NcbiBuffer`: that one models the C library's and libstdc++'s buffer under
//! `cout` or an `ofstream` (glibc's `_IO_new_file_xsputn`, libstdc++'s `basic_filebuf::xsputn`),
//! the layer below this one, which the tie to `cerr` empties at every post. This buffer is
//! NCBI's own, one per object stream (`CObjectOStream::m_Output`), between the stream's writes
//! and the `CNcbiOstream`: 4096 bytes at first, flushed when a reservation does not fit, keeping
//! the stream's back limit (ASN.1 text 80, XML 1, ASN.1 binary 0).
//!
//! The encoders keep it only in their NCBI form (`asn_text::AsnTextEncoder::ncbi`,
//! `asn_ber::AsnBinaryEncoder::ncbi`), which `encoder::Traced` writes through when a document
//! is written again for its messages; their plain form keeps nothing.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:46
// ```c++
// static const size_t KInitialBufferSize = 4096;
// ```
// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:64-72 (ASN.1 text: 80)
// ```c++
// CObjectOStreamAsn::CObjectOStreamAsn(CNcbiOstream& out,
//                                      EFixNonPrint how)
//     : CObjectOStream(eSerial_AsnText, out)
// {
//     FixNonPrint(how);
//     m_Output.SetBackLimit(80);
// ```
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:783-791 (ASN.1 binary keeps the
// constructor's 0: objostrasnb.cpp sets no back limit)
// ```c++
// COStreamBuffer::COStreamBuffer(CNcbiOstream& out, bool deleteOut)
//     THROWS1((bad_alloc))
//     : m_Output(out), m_DeleteOutput(deleteOut), m_Closed(false), m_Error(0),
//       m_IndentLevel(0), m_BufferPos(0),
//       m_Buffer(new char[KInitialBufferSize]),
//       m_CurrentPos(m_Buffer),
//       m_BufferEnd(m_Buffer + KInitialBufferSize),
//       m_Line(1), m_LineLength(0),
//       m_BackLimit(0), m_UseIndentation(true), m_UseEol(true),
// ```

/// `KInitialBufferSize`.
pub(crate) const INITIAL_SIZE: usize = 4096;

/// The back limit of the ASN.1 text stream (`SetBackLimit(80)`).
pub(crate) const ASN_TEXT_BACK_LIMIT: usize = 80;

/// The back limit of the ASN.1 binary stream (`COStreamBuffer`'s own 0).
pub(crate) const ASN_BINARY_BACK_LIMIT: usize = 0;

/// The bytes that one of NCBI's object streams has written and handed to its output.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct StreamBuffer {
    /// `GetUsedSpace()`: the bytes in the buffer.
    used: usize,
    /// `m_BufferEnd - m_Buffer`.
    size: usize,
    /// `m_BufferPos`: the bytes handed to the output.
    flushed: usize,
    back_limit: usize,
}

impl StreamBuffer {
    /// The empty buffer of a new stream whose back limit is `back_limit`.
    pub(crate) const fn new(back_limit: usize) -> Self {
        StreamBuffer {
            used: 0,
            size: INITIAL_SIZE,
            flushed: 0,
            back_limit,
        }
    }

    /// The bytes of the document handed to the output (`m_BufferPos`).
    pub(crate) fn flushed(&self) -> usize {
        self.flushed
    }

    /// `Skip(count)` (and `Reserve(count)` followed by `count` bytes, as `PutChar`, `PutEol` and
    /// `PutEolAtWordEnd`'s inserted end of line): when the bytes do not fit, `DoReserve` hands every
    /// byte but the back limit to the output and grows the buffer if it is still too small.
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:283-296
    // ```c++
    // char* COStreamBuffer::Skip(size_t count)
    // {
    //     char* pos = m_CurrentPos;
    //     char* end = pos + count;
    //     if ( end <= m_BufferEnd ) {
    //         // enough space in buffer
    //         m_CurrentPos = end;
    //         m_LineLength += count;
    //         return pos;
    //     }
    //     else {
    //         return DoSkip(count);
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-942
    // ```c++
    // void COStreamBuffer::FlushBuffer(bool fullBuffer)
    // ...
    //     size_t used = GetUsedSpace();
    //     size_t count;
    //     size_t leave;
    //     if ( fullBuffer ) {
    //         count = used;
    //         leave = 0;
    //     }
    //     else {
    //         leave = m_BackLimit;
    //         if ( used < leave )
    //             return; // none to flush
    //         count = used - leave;
    //     }
    //     if ( count != 0 ) {
    //         if ( !m_Output.write(m_Buffer, count) ) {
    // ...
    // char* COStreamBuffer::DoReserve(size_t count)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     FlushBuffer(false);
    //     size_t usedSize = m_CurrentPos - m_Buffer;
    //     size_t needSize = usedSize + count;
    //     size_t bufferSize = m_BufferEnd - m_Buffer;
    //     if ( bufferSize < needSize ) {
    //         // realloc too small buffer
    //         do {
    //             bufferSize = BiggerBufferSize(bufferSize);
    //         } while ( bufferSize < needSize );
    // ```
    #[inline]
    pub(crate) fn skip(&mut self, count: usize) {
        if self.used + count > self.size {
            self.do_reserve(count);
        }
        self.used += count;
    }

    /// `DoReserve`: the flush and the growth of `skip`.
    #[cold]
    #[inline(never)]
    fn do_reserve(&mut self, count: usize) {
        self.flush_but_back_limit();
        while self.size < self.used + count {
            self.size = bigger_buffer_size(self.size);
        }
    }

    /// `FlushBuffer(false)`.
    fn flush_but_back_limit(&mut self) {
        if self.used >= self.back_limit {
            self.flushed += self.used - self.back_limit;
            self.used = self.back_limit;
        }
    }

    /// `PutString(str, length)`: `Skip` below 1024 bytes, else `Write`.
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:329-339
    // ```c++
    // void COStreamBuffer::PutString(const char* str, size_t length)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     if ( length < 1024 ) {
    //         memcpy(Skip(length), str, length);
    //     }
    //     else {
    //         Write(str, length);
    //     }
    // }
    // ```
    #[inline]
    pub(crate) fn put_string(&mut self, length: usize) {
        if length < 1024 {
            self.skip(length);
        } else {
            self.write(length);
        }
    }

    /// `Write(data, dataLength)`: the buffer filled to its end and flushed (but the back limit)
    /// as often as the bytes need, the rest kept.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1143-1161
    // ```c++
    // void COStreamBuffer::Write(const char* data, size_t dataLength)
    //     THROWS1((CIOException, bad_alloc))
    // {
    //     while ( dataLength > 0 ) {
    //         size_t available = GetAvailableSpace();
    //         if ( available == 0 ) {
    //             FlushBuffer(false);
    //             available = GetAvailableSpace();
    //         }
    //         if ( available >= dataLength )
    //             break; // current chunk will fit in buffer
    //         memcpy(m_CurrentPos, data, available);
    //         m_CurrentPos += available;
    //         data += available;
    //         dataLength -= available;
    //     }
    //     memcpy(m_CurrentPos, data, dataLength);
    //     m_CurrentPos += dataLength;
    // }
    // ```
    #[cold]
    #[inline(never)]
    fn write(&mut self, mut length: usize) {
        while length > 0 {
            let mut available = self.size - self.used;
            if available == 0 {
                self.flush_but_back_limit();
                available = self.size - self.used;
            }
            if available >= length {
                break;
            }
            self.used += available;
            length -= available;
        }
        self.used += length;
    }

    /// `FlushBuffer()` (`fullBuffer` true): every byte to the output, as the ASN.1 text stream's
    /// separator after the object does (`WriteSeparator`, objostrasn.cpp:941-945).
    pub(crate) fn flush_all(&mut self) {
        self.flushed += self.used;
        self.used = 0;
    }
}

/// `BiggerBufferSize`: twice the size.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:48-52
// ```c++
// static inline
// size_t BiggerBufferSize(size_t size) THROWS1_NONE
// {
//     return size * 2;
// }
// ```
fn bigger_buffer_size(size: usize) -> usize {
    size * 2
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ASN.1 text stream: the reservation that does not fit hands all but the last 80 bytes
    /// out (4,015 when 4,095 are held and two more come, as the archive of fasta_input
    /// `blastp.arch_title_s.clean.o11.merged`).
    #[test]
    fn a_reservation_that_does_not_fit_hands_all_but_the_back_limit_out() {
        let mut buffer = StreamBuffer::new(ASN_TEXT_BACK_LIMIT);
        buffer.skip(4095);
        assert_eq!(buffer.flushed(), 0);
        buffer.skip(2);
        assert_eq!(buffer.flushed(), 4015);
        buffer.flush_all();
        assert_eq!(buffer.flushed(), 4097);
    }

    /// The binary stream (back limit 0): bytes one by one fill the buffer to its end; a string of
    /// 1024 bytes or more fills it and flushes it whole as often as it needs.
    #[test]
    fn the_binary_stream_flushes_everything_and_writes_long_strings_in_pieces() {
        let mut buffer = StreamBuffer::new(ASN_BINARY_BACK_LIMIT);
        for _ in 0..4096 {
            buffer.skip(1);
        }
        assert_eq!(buffer.flushed(), 0);
        buffer.skip(1);
        assert_eq!(buffer.flushed(), 4096);
        // 1 held + 9000: fills 4095, flushes 4096, fills 4096, flushes 4096, keeps 809.
        buffer.put_string(9000);
        assert_eq!(buffer.flushed(), 4096 * 3);
        assert_eq!(buffer.used, 809);
        // 809 + 1000 fit; a 4000-byte reservation after them flushes 1809.
        buffer.put_string(1000);
        assert_eq!(buffer.flushed(), 4096 * 3);
        buffer.skip(3000);
        assert_eq!(buffer.flushed(), 4096 * 3 + 1809);
    }
}
