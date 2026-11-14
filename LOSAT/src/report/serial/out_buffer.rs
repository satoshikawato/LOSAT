//! The output buffer of NCBI's object streams (`COStreamBuffer`) as `encoder::Traced` follows it:
//! which bytes of a document NCBI has written to its output stream when it posts a byte that it
//! replaced (`encoder::ReplacedChar`). NCBI's `cerr` is tied to `cout`, so with both in one file
//! (`2>&1`) each post follows exactly those bytes: `Traced` hands its sink the same bytes before
//! each post (`OutBuffer::written`).
//!
//! The buffer holds 4096 bytes at first. A call that puts `count` bytes (`Skip`, `Reserve`) and
//! does not fit first writes all the buffer but its last `back limit` bytes to the output
//! (`FlushBuffer(false)`), then doubles the buffer until the call fits; a string of 1024 bytes or
//! more goes through `Write`, which fills the buffer up and flushes it as often as it needs
//! without growing it; the end of a document writes everything (`EndOfWrite`). The back limit
//! is 80 for ASN.1 text (the string wrap inserts its end of line up to a line back), 1 for XML
//! (`BackChar` takes back the `>` of a tag) and 0 for ASN.1 binary.
//!
//! NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:46-52
//! ```c++
//! static const size_t KInitialBufferSize = 4096;
//!
//! static inline
//! size_t BiggerBufferSize(size_t size) THROWS1_NONE
//! {
//!     return size * 2;
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-947
//! ```c++
//! void COStreamBuffer::FlushBuffer(bool fullBuffer)
//!     THROWS1((CIOException))
//! {
//!     ...
//!     size_t used = GetUsedSpace();
//!     size_t count;
//!     size_t leave;
//!     if ( fullBuffer ) {
//!         count = used;
//!         leave = 0;
//!     }
//!     else {
//!         leave = m_BackLimit;
//!         if ( used < leave )
//!             return; // none to flush
//!         count = used - leave;
//!     }
//!     if ( count != 0 ) {
//!         if ( !m_Output.write(m_Buffer, count) ) {
//!     ...
//! char* COStreamBuffer::DoReserve(size_t count)
//!     THROWS1((CIOException, bad_alloc))
//! {
//!     FlushBuffer(false);
//!     size_t usedSize = m_CurrentPos - m_Buffer;
//!     size_t needSize = usedSize + count;
//!     size_t bufferSize = m_BufferEnd - m_Buffer;
//!     if ( bufferSize < needSize ) {
//!         // realloc too small buffer
//!         do {
//!             bufferSize = BiggerBufferSize(bufferSize);
//!         } while ( bufferSize < needSize );
//!     ...
//! ```
//! NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:281-339
//! ```c++
//! inline
//! char* COStreamBuffer::Skip(size_t count)
//!     THROWS1((CIOException, bad_alloc))
//! {
//!     char* pos = m_CurrentPos;
//!     char* end = pos + count;
//!     if ( end <= m_BufferEnd ) {
//!         // enough space in buffer
//!         m_CurrentPos = end;
//!         m_LineLength += count;
//!         return pos;
//!     }
//!     else {
//!         return DoSkip(count);
//!     }
//! }
//! ...
//! inline
//! void COStreamBuffer::PutString(const char* str, size_t length)
//!     THROWS1((CIOException, bad_alloc))
//! {
//!     if ( length < 1024 ) {
//!         memcpy(Skip(length), str, length);
//!     }
//!     else {
//!         Write(str, length);
//!     }
//! }
//! ```
//! (`PutChar` is `Skip(1)`; `PutEol` is `Reserve(1)` and the end of line, then `PutIndent`, a
//! `Skip` of the indentation; `PutInt4` and `PutInt8` are one `Skip` of the digits.)
//!
//! Each `Encoder` call of `Traced` gives the bytes that it wrote and what it wrote (`Piece`);
//! `OutBuffer::put` splits them into the buffer calls that NCBI's stream makes for them.

/// NCBI's object streams whose buffer `Traced` follows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stream {
    /// `CObjectOStreamAsn`.
    AsnText,
    /// `CObjectOStreamAsnBinary`.
    AsnBinary,
    /// `CObjectOStreamXml`.
    Xml,
}

impl Stream {
    /// `COStreamBuffer::m_BackLimit` of the stream.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:64-69
    /// ```c++
    /// CObjectOStreamAsn::CObjectOStreamAsn(CNcbiOstream& out,
    ///                                      EFixNonPrint how)
    ///     : CObjectOStream(eSerial_AsnText, out)
    /// {
    ///     FixNonPrint(how);
    ///     m_Output.SetBackLimit(80);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:84-87
    /// ```c++
    ///       m_SkipIndent( false ), m_SkipNextTag( false )
    /// {
    ///     m_Output.SetBackLimit(1);
    /// }
    /// ```
    /// (`CObjectOStreamAsnBinary` sets none: `m_BackLimit(0)` of the `COStreamBuffer`
    /// constructor, util/strbuffer.cpp:783-796.)
    fn back_limit(self) -> usize {
        match self {
            Stream::AsnText => 80,
            Stream::Xml => 1,
            Stream::AsnBinary => 0,
        }
    }

    /// Whether the stream writes `#` for the byte `byte` of a VisibleString and posts it
    /// (`ReplaceVisibleChar` with `eFNP_ReplaceAndWarn`, `encoder::ReplacedChar`): the one rule
    /// that the encoders (`write_string` of `asn_text` and `asn_ber`; `xml`'s match arms, which
    /// a test holds to it) and the archive's check of its strings
    /// (`archive::ArchiveEncoding::writes_text`) read. Text and binary: a byte outside
    /// `' '..='~'` (`GoodVisibleChar`); XML: a control byte other than NUL, TAB, LF and CR
    /// (`BAD_CHAR`, with NCBI's signed `char`: no byte from 0x80 on).
    ///
    /// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:455-459
    /// ```c++
    /// inline
    /// bool GoodVisibleChar(char c)
    /// {
    ///     return c >= ' ' && c <= '~';
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:548-550
    /// ```c++
    /// inline bool BAD_CHAR(char x) {
    ///     return (x < 0x20 && x > 0x0 && x != 0x9 && x != 0xA && x != 0xD);
    /// }
    /// ```
    #[inline(always)]
    pub const fn replaces(self, byte: u8) -> bool {
        match self {
            Stream::AsnText | Stream::AsnBinary => !matches!(byte, 0x20..=0x7e),
            Stream::Xml => matches!(byte, 0x01..=0x08 | 0x0b | 0x0c | 0x0e..=0x1f),
        }
    }
}

/// What one `Encoder` call wrote, which decides NCBI's buffer calls for its bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Piece {
    /// `WriteFileHeader`.
    Header,
    /// Brackets, separators, member and variant names, tags, ends of line and indentation.
    Markup,
    /// `WriteInt4` / `WriteInt8`, after the markup that opens it.
    Number,
    /// `WriteDouble`.
    Real,
    /// `WriteBool`.
    Bool,
    /// `WriteEnum` (an ENUMERATED value or an INTEGER with named values).
    Enum,
    /// `WriteString` / `WriteStringStore`.
    Str,
    /// An OCTET STRING (`BeginBytes`, `WriteBytes`, `EndBytes`).
    Octets,
    /// `EndOfWrite` and the separator.
    End,
}

/// NCBI's buffer of one document (each document has a stream of its own,
/// `serialobject.cpp` `WriteObject`).
#[derive(Clone, Debug)]
pub(super) struct OutBuffer {
    stream: Stream,
    /// `m_BufferEnd - m_Buffer`.
    capacity: usize,
    /// `GetUsedSpace()`.
    used: usize,
    /// The bytes of the document written to the output (`m_BufferPos`).
    written: u64,
}

impl OutBuffer {
    pub(super) fn new(stream: Stream) -> Self {
        Self {
            stream,
            capacity: 4096,
            used: 0,
            written: 0,
        }
    }

    /// The bytes of the document that NCBI has written to its output.
    pub(super) fn written(&self) -> u64 {
        self.written
    }

    /// `FlushBuffer(false)`.
    fn flush_buffer(&mut self) {
        let leave = self.stream.back_limit();
        if self.used < leave {
            return;
        }
        self.written += (self.used - leave) as u64;
        self.used = leave;
    }

    /// `Skip(count)`, or `Reserve(count)` and `count` bytes.
    fn skip(&mut self, count: usize) {
        if self.used + count > self.capacity {
            // `DoReserve`.
            self.flush_buffer();
            while self.capacity < self.used + count {
                self.capacity *= 2;
            }
        }
        self.used += count;
    }

    /// `count` calls of `Skip(1)` (`PutChar`).
    fn skip_ones(&mut self, mut count: usize) {
        while count > 0 {
            if self.used == self.capacity {
                self.skip(1);
                count -= 1;
                continue;
            }
            let room = (self.capacity - self.used).min(count);
            self.used += room;
            count -= room;
        }
    }

    /// `PutString(str, length)`: `Skip` below 1024 bytes, `Write` from 1024.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1143-1161
    /// ```c++
    /// void COStreamBuffer::Write(const char* data, size_t dataLength)
    ///     THROWS1((CIOException, bad_alloc))
    /// {
    ///     while ( dataLength > 0 ) {
    ///         size_t available = GetAvailableSpace();
    ///         if ( available == 0 ) {
    ///             FlushBuffer(false);
    ///             available = GetAvailableSpace();
    ///         }
    ///         if ( available >= dataLength )
    ///             break; // current chunk will fit in buffer
    ///         memcpy(m_CurrentPos, data, available);
    ///         m_CurrentPos += available;
    ///         data += available;
    ///         dataLength -= available;
    ///     }
    ///     memcpy(m_CurrentPos, data, dataLength);
    ///     m_CurrentPos += dataLength;
    /// }
    /// ```
    fn put_string(&mut self, mut length: usize) {
        if length < 1024 {
            self.skip(length);
            return;
        }
        while length > 0 {
            if self.capacity == self.used {
                self.flush_buffer();
            }
            let available = self.capacity - self.used;
            if available >= length {
                break;
            }
            self.used += available;
            length -= available;
        }
        self.used += length;
    }

    /// `BackChar`: the stream takes its last byte back (XML's `OpenTagEndBack`).
    fn back_char(&mut self) {
        self.used -= 1;
    }

    /// The buffer calls of one `Encoder` call that wrote `bytes` (`piece`). `posts` are the
    /// offsets in `bytes` at which NCBI posts a replaced byte (before the byte's `#`); `written`
    /// gets `written()` at each of them. The end of a document writes everything (the ASN.1
    /// text separator after it, objostrasn.cpp:941-945).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:559-564
    /// ```c++
    /// void CObjectOStream::EndOfWrite(void)
    /// {
    ///     FlushBuffer();
    ///     if ( m_Objects )
    ///         m_Objects->Clear();
    /// }
    /// ```
    pub(super) fn put(
        &mut self,
        piece: Piece,
        bytes: &[u8],
        posts: &[usize],
        written: &mut Vec<u64>,
    ) {
        if piece == Piece::End {
            self.written += (self.used + bytes.len()) as u64;
            self.used = 0;
            return;
        }
        match self.stream {
            Stream::AsnText => self.asn_text(piece, bytes, posts, written),
            Stream::AsnBinary => self.asn_binary(piece, bytes, posts, written),
            Stream::Xml => self.xml(piece, bytes, posts, written),
        }
    }

    /// The calls of `CObjectOStreamAsn` (objostrasn.cpp).
    fn asn_text(&mut self, piece: Piece, bytes: &[u8], posts: &[usize], written: &mut Vec<u64>) {
        match piece {
            // `WriteId(type name)` and `PutString(" ::= ")` (objostrasn.cpp:105-111).
            Piece::Header => {
                let name = bytes.len().saturating_sub(5);
                self.skip(name);
                self.skip(bytes.len() - name);
            }
            Piece::Markup => self.asn_text_markup(bytes),
            // `PutInt4`, `PutInt8`, `PutString("TRUE")`.
            Piece::Number | Piece::Bool => self.skip(bytes.len()),
            // `PutChar(tolower(name[0]))` and `PutString(name + 1)` (objostrasn.cpp:113-127).
            Piece::Enum => {
                self.skip(1);
                self.skip(bytes.len() - 1);
            }
            Piece::Real => self.asn_text_real(bytes),
            // `WriteString` puts each byte with `PutChar`, after `WrapAt(78, true)` (an end of
            // line inserted with `Reserve(1)`); a doubled `"` is a `PutChar` of its own
            // (objostrasn.cpp:413-444). Every call puts one byte, so the encoder's bytes before a
            // post count NCBI's calls before it.
            Piece::Str => self.ones_with_posts(bytes.len(), posts, written),
            // `PutChar('\'')`, `WrapAt(78, false)` and two `PutChar`s per byte, then `WrapAt` and
            // `PutString("'H")` (objostrasn.cpp:868-895).
            Piece::Octets => {
                let end = bytes.len().saturating_sub(2);
                self.skip_ones(end);
                self.skip(bytes.len() - end);
            }
            Piece::End => {}
        }
    }

    /// The markup of ASN.1 text: `PutChar` of `,` `{` `}` and of the space after a name,
    /// `PutEol` (an end of line, then a `Skip` of the indentation), `PutString` of a member or
    /// variant name (`NextElement`, `StartBlock`, `EndBlock`, `WriteMemberId`,
    /// objostrasn.cpp:539-561, 646-667).
    fn asn_text_markup(&mut self, bytes: &[u8]) {
        let mut i = 0;
        while i < bytes.len() {
            match bytes[i] {
                b'\n' => {
                    self.skip(1);
                    i += 1;
                    let indent = bytes[i..].iter().take_while(|&&b| b == b' ').count();
                    self.skip(indent);
                    i += indent;
                }
                b',' | b'{' | b'}' | b' ' => {
                    self.skip(1);
                    i += 1;
                }
                _ => {
                    let name = bytes[i..]
                        .iter()
                        .take_while(|&&b| !matches!(b, b' ' | b'\n' | b',' | b'{' | b'}'))
                        .count();
                    self.skip(name);
                    i += name;
                }
            }
        }
    }

    /// `CObjectOStreamAsn::WriteDouble2` with `m_FastWriteDouble` (the default): `PutString`
    /// of `{ `, of the sign, of the digits and of `, 10, `, `PutInt4` of the exponent,
    /// `PutString(" }")`; zero, a NaN and an infinity are one `PutString`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:198-222
    /// ```c++
    ///     if ( data == 0.0 ) {
    ///         double zero = 0.;
    ///         if (memcmp(&data, &zero, sizeof(double)) == 0) {
    ///             m_Output.PutString("{ 0, 10, 0 }");
    ///         } else {
    ///             m_Output.PutString("{ -0, 10, 0 }");
    ///         }
    ///         return;
    ///     }
    ///
    ///     char buffer[128];
    ///     if (m_FastWriteDouble) {
    ///         int dec, sign;
    ///         size_t len = NStr::DoubleToString_Ecvt(
    ///             data, digits, buffer, sizeof(buffer), &dec, &sign);
    ///         _ASSERT(len > 0);
    ///         m_Output.PutString("{ ");
    ///         if (sign < 0) {
    ///             m_Output.PutString("-");
    ///         }
    ///         m_Output.PutString(buffer,len);
    ///         m_Output.PutString(", 10, ");
    ///         m_Output.PutInt4(dec - (int)(len-1));
    ///         m_Output.PutString(" }");
    /// ```
    fn asn_text_real(&mut self, bytes: &[u8]) {
        if !bytes.starts_with(b"{ ") || bytes == b"{ 0, 10, 0 }" || bytes == b"{ -0, 10, 0 }" {
            self.skip(bytes.len());
            return;
        }
        self.skip(2);
        let mut i = 2;
        if bytes.get(i) == Some(&b'-') {
            self.skip(1);
            i += 1;
        }
        let digits = bytes[i..].iter().take_while(|&&b| b != b',').count();
        self.skip(digits);
        i += digits;
        self.skip(6);
        i += 6;
        let exponent = bytes[i..].iter().take_while(|&&b| b != b' ').count();
        self.skip(exponent);
        self.skip(2);
    }

    /// `count` calls of `PutChar`, with the posts after `posts[k]` of them.
    fn ones_with_posts(&mut self, count: usize, posts: &[usize], written: &mut Vec<u64>) {
        let mut done = 0;
        for &at in posts {
            self.skip_ones(at - done);
            done = at;
            written.push(self.written);
        }
        self.skip_ones(count - done);
    }

    /// The calls of `CObjectOStreamAsnBinary` (objostrasnb.cpp): each tag, length and number
    /// byte is a `WriteByte` (`PutChar`), the bytes of a string, an OCTET STRING and a REAL's
    /// digits one `WriteBytes` (`PutString`).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:831-872
    /// ```c++
    /// void CObjectOStreamAsnBinary::WriteString(const string& str, EStringType type)
    /// {
    ///     size_t length = str.size(), fixed = 0;
    ///     WriteStringTag(type);
    ///     ...
    ///     WriteLength(length - fixed);
    ///     if ( type == eStringTypeVisible && x_FixCharsMethod() != eFNP_Allow ) {
    ///         size_t done = 0;
    ///         for ( size_t i = 0; i < length; ++i ) {
    ///             char c = str[i];
    ///             if ( !GoodVisibleChar(c) ) {
    ///                 ...
    ///                 if ( i > done ) {
    ///                     WriteBytes(str.data() + done, i - done);
    ///                 }
    ///                 c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, str, x_FixCharsSubst());
    ///                 if (c != 0) {
    ///                     WriteByte(c);
    ///                 }
    ///                 done = i + 1;
    ///             }
    ///         }
    ///         if ( done < length ) {
    ///             WriteBytes(str.data() + done, length - done);
    ///         }
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:798-813
    /// ```c++
    ///     WriteLength(width + 1);
    ///     if (width) {
    ///         ...
    ///         WriteByte(type);
    ///         WriteBytes(buffer, width);
    /// ```
    /// (`WriteBytes` puts nothing for 0 bytes, objostrasnb.cpp:230-247.)
    fn asn_binary(&mut self, piece: Piece, bytes: &[u8], posts: &[usize], written: &mut Vec<u64>) {
        match piece {
            Piece::Str | Piece::Octets => {
                // The tag, the length (`0x8n` and n bytes in the long form), the bytes.
                let length_bytes = match bytes.get(1) {
                    Some(&first) if first >= 0x80 => 1 + usize::from(first & 0x7f),
                    _ => 1,
                };
                let head = (1 + length_bytes).min(bytes.len());
                self.skip_ones(head);
                let mut done = head;
                for &at in posts {
                    if at > done {
                        self.put_string(at - done);
                    }
                    written.push(self.written);
                    self.skip(1);
                    done = at + 1;
                }
                if bytes.len() > done {
                    self.put_string(bytes.len() - done);
                }
            }
            Piece::Real if bytes.len() > 3 => {
                self.skip_ones(3);
                self.put_string(bytes.len() - 3);
            }
            _ => self.skip_ones(bytes.len()),
        }
    }

    /// The calls of `CObjectOStreamXml` (objostrxml.cpp): the markup of tags first, then the
    /// value.
    fn xml(&mut self, piece: Piece, bytes: &[u8], posts: &[usize], written: &mut Vec<u64>) {
        if piece == Piece::Header {
            // At the start of the buffer: its calls do not reach the end of the buffer.
            self.skip(bytes.len());
            return;
        }
        let (mut i, after_name) = self.xml_markup(bytes);
        let rest = bytes.len() - i;
        match piece {
            Piece::Number | Piece::Real => {
                // `WriteDouble2` of an infinity: `PutChar('-')` and `PutString("INF", 3)`.
                if bytes[i..].starts_with(b"-INF") {
                    self.skip(1);
                    self.skip(3);
                } else if rest > 0 {
                    self.skip(rest);
                }
            }
            Piece::Str => {
                let mut posts = posts.iter().peekable();
                while i < bytes.len() {
                    while posts.next_if(|&&at| at <= i).is_some() {
                        written.push(self.written);
                    }
                    i += self.xml_char(&bytes[i..]);
                }
            }
            Piece::Octets => self.skip_ones(rest),
            Piece::Bool | Piece::Enum => {
                // `OpenTagEnd` put the `>` of the tag (here unless an earlier call put it), and
                // `OpenTagEndBack` takes it back (objostrxml.cpp:1029-1050).
                if after_name {
                    self.skip(1);
                }
                self.back_char();
                if piece == Piece::Bool {
                    // `PutString(" value=\"true\"")`, `SelfCloseTagEnd` (`PutString("/>")`).
                    self.skip(rest - 2);
                    self.skip(2);
                } else {
                    // `PutString(" value=\"")`, the name, `PutChar('"')`, then
                    // `SelfCloseTagEnd` or, for an INTEGER with named values, `OpenTagEnd` and
                    // `PutInt4` (objostrxml.cpp:427-487).
                    self.skip(8);
                    let name = bytes[i + 8..].iter().take_while(|&&b| b != b'"').count();
                    self.skip(name);
                    self.skip(1);
                    i += 8 + name + 1;
                    if bytes[i..].starts_with(b"/>") {
                        self.skip(2);
                    } else {
                        self.skip(1);
                        if bytes.len() > i + 1 {
                            self.skip(bytes.len() - i - 1);
                        }
                    }
                }
            }
            Piece::Header | Piece::Markup | Piece::End => {}
        }
    }

    /// The tags at the start of `bytes`: `OpenTagStart` (`PutEol(false)`, `PutIndent`,
    /// `PutChar('<')`) or `CloseTagStart` (`PutEol(false)` and `PutIndent` after a closing tag,
    /// `PutString("</")`), the name, `OpenTagEnd` / `CloseTagEnd` (`PutChar('>')`). A name is
    /// put as `PrintTagName` puts it: the type's name, then for each member or variant
    /// `PutChar('_')` and `PutString(name)`, for a list element `PutString("_E")`.
    /// Returns where the tags end and whether the last tag's name ends there without its `>`
    /// (a `value` attribute follows).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1089-1132
    /// ```c++
    ///     case TFrame::eFrameClassMember:
    ///     case TFrame::eFrameChoiceVariant:
    ///         {
    ///             bool attlist = m_Attlist;
    ///             if (!x_IsStdXml()) {
    ///                 PrintTagName(level + 1);
    ///                 m_Output.PutChar('_');
    ///                 m_Attlist = true;
    ///             }
    ///             WriteTag(frame.GetMemberId().GetName());
    ///             ...
    ///     case TFrame::eFrameArrayElement:
    ///         {
    ///             PrintTagName(level + 1);
    ///             if (!x_IsStdXml()) {
    ///                 m_Output.PutString("_E");
    ///             }
    ///             return;
    ///         }
    /// ```
    /// (ASN.1 member names start with a lower-case letter, so a piece `E` is an element.)
    fn xml_markup(&mut self, bytes: &[u8]) -> (usize, bool) {
        let mut i = 0;
        while i < bytes.len() {
            match bytes[i] {
                b'\n' => {
                    self.skip(1);
                    i += 1;
                    let indent = bytes[i..].iter().take_while(|&&b| b == b' ').count();
                    self.skip(indent);
                    i += indent;
                }
                b'<' => {
                    if bytes.get(i + 1) == Some(&b'/') {
                        self.skip(2);
                        i += 2;
                    } else {
                        self.skip(1);
                        i += 1;
                    }
                    let name = bytes[i..]
                        .iter()
                        .take_while(|&&b| !matches!(b, b'>' | b' '))
                        .count();
                    for (k, part) in bytes[i..i + name].split(|&b| b == b'_').enumerate() {
                        match (k, part) {
                            (0, _) => self.skip(part.len()),
                            (_, b"E") => self.skip(2),
                            _ => {
                                self.skip(1);
                                self.skip(part.len());
                            }
                        }
                    }
                    i += name;
                    if bytes.get(i) == Some(&b'>') {
                        self.skip(1);
                        i += 1;
                    } else {
                        return (i, true);
                    }
                }
                _ => break,
            }
        }
        (i, false)
    }

    /// `WriteEscapedChar` of the byte at the start of `bytes` (`PutString` of an entity, for a
    /// control byte `PutString("&#x")`, a `PutChar` per hex digit and `PutChar(';')`, else
    /// `PutChar`); the bytes that it wrote.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:502-538
    /// ```c++
    ///     case '&':
    ///         m_Output.PutString("&amp;");
    ///         break;
    ///     ...
    ///     default:
    ///         if ((unsigned int)c < 0x20) {
    ///             m_Output.PutString("&#x");
    ///             Uint1 ch = c;
    ///             unsigned hi = ch >> 4;
    ///             unsigned lo = ch & 0xF;
    ///             if ( hi ) {
    ///                 m_Output.PutChar("0123456789abcdef"[hi]);
    ///             }
    ///             m_Output.PutChar("0123456789abcdef"[lo]);
    ///             m_Output.PutChar(';');
    ///         } else {
    ///             m_Output.PutChar(c);
    ///         }
    /// ```
    fn xml_char(&mut self, bytes: &[u8]) -> usize {
        if bytes[0] != b'&' {
            self.skip(1);
            return 1;
        }
        let end = bytes
            .iter()
            .position(|&b| b == b';')
            .map_or(bytes.len(), |at| at + 1);
        if bytes.starts_with(b"&#x") {
            self.skip(3);
            self.skip_ones(end - 3);
        } else {
            self.skip(end);
        }
        end
    }
}

#[cfg(test)]
mod tests {
    use super::{OutBuffer, Piece, Stream};

    /// What NCBI had written at each post of `bytes` (one `piece`) after `before` bytes of the
    /// same kind of call.
    fn written_at(
        stream: Stream,
        before: &[(Piece, &[u8])],
        piece: Piece,
        bytes: &[u8],
        posts: &[usize],
    ) -> Vec<u64> {
        let mut buffer = OutBuffer::new(stream);
        for (piece, bytes) in before {
            buffer.put(*piece, bytes, &[], &mut Vec::new());
        }
        let mut written = Vec::new();
        buffer.put(piece, bytes, posts, &mut written);
        written
    }

    #[test]
    fn a_full_buffer_keeps_its_back_limit() {
        // 4096 bytes of `PutChar` fill the buffer; the next one flushes it first: ASN.1 text
        // keeps 80 bytes, XML 1, binary none. A post before that byte comes before the flush.
        let chars = vec![b'a'; 4098];
        assert_eq!(
            written_at(Stream::AsnText, &[], Piece::Str, &chars, &[4096, 4097]),
            [0, 4016]
        );
        for (stream, written) in [(Stream::Xml, 4095), (Stream::AsnBinary, 4096)] {
            let mut buffer = OutBuffer::new(stream);
            buffer.put(Piece::Octets, &chars[..4097], &[], &mut Vec::new());
            assert_eq!(buffer.written(), written);
        }
        // The next flush comes when the buffer is full again.
        let chars = vec![b'a'; 8200];
        assert_eq!(
            written_at(Stream::AsnText, &[], Piece::Str, &chars, &[8112, 8113]),
            [4016, 8032]
        );
    }

    #[test]
    fn a_call_that_does_not_fit_flushes_before_it() {
        // 4090 bytes, then a name of 10 bytes: the buffer is flushed before the name (4010 of
        // 4090 bytes written, 80 kept), not at 4096.
        let mut buffer = OutBuffer::new(Stream::AsnText);
        buffer.put(Piece::Str, &[b'a'; 4090], &[], &mut Vec::new());
        buffer.put(Piece::Number, b"1234567890", &[], &mut Vec::new());
        assert_eq!(buffer.written(), 4010);
        assert_eq!(buffer.used, 90);
    }

    #[test]
    fn a_long_string_fills_the_buffer_up() {
        // ASN.1 binary: a run of 1100 bytes of a string after 3004 bytes goes through `Write`,
        // which fills the buffer (4096 written) without growing it; a run of 1000 bytes that
        // does not fit is one `Skip`, which flushes first.
        let mut string = vec![0x1a, 0x82, 0x05, 0xdc];
        string.extend(std::iter::repeat_n(b'L', 1500));
        string[4 + 1100] = b'#';
        let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[0u8; 3000])];
        assert_eq!(
            written_at(Stream::AsnBinary, before, Piece::Str, &string, &[4 + 1100]),
            [4096]
        );
        string[4 + 1100] = b'L';
        string[4 + 1000] = b'#';
        let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[0u8; 3500])];
        assert_eq!(
            written_at(Stream::AsnBinary, before, Piece::Str, &string, &[4 + 1000]),
            [3504]
        );
    }

    /// The buffer after a number of `before` digits (one `Skip` in ASN.1 text and XML) and the
    /// calls of `piece`.
    fn after(stream: Stream, before: usize, piece: Piece, bytes: &[u8]) -> OutBuffer {
        let mut buffer = OutBuffer::new(stream);
        buffer.put(Piece::Number, &vec![b'1'; before], &[], &mut Vec::new());
        buffer.put(piece, bytes, &[], &mut Vec::new());
        buffer
    }

    #[test]
    fn asn_text_header_enum_and_real_are_split_as_ncbi_calls_them() {
        // `WriteFileHeader`: `WriteId(name)` (14 bytes fit at 4080) and `PutString(" ::= ")`,
        // which flushes 4094 - 80 bytes; one call of 19 bytes would flush at 4080.
        let buffer = after(Stream::AsnText, 4080, Piece::Header, b"Blast4-archive ::= ");
        assert_eq!((buffer.written(), buffer.used), (4014, 85));
        // `WriteEnum`: `PutChar` of the first letter (the 4096th byte), then `PutString` of the
        // rest, which flushes.
        let buffer = after(Stream::AsnText, 4095, Piece::Enum, b"plus");
        assert_eq!((buffer.written(), buffer.used), (4016, 83));
        // `WriteDouble2`: `{ ` fits at 4090, the 9 digits flush 4092 - 80 bytes; then `, 10, `,
        // the exponent and ` }`.
        let buffer = after(Stream::AsnText, 4090, Piece::Real, b"{ 123456789, 10, -9 }");
        assert_eq!((buffer.written(), buffer.used), (4012, 99));
        // A negative mantissa puts its sign alone.
        let buffer = after(Stream::AsnText, 4093, Piece::Real, b"{ -12, 10, 3 }");
        assert_eq!((buffer.written(), buffer.used), (4016, 80 + 2 + 6 + 1 + 2));
        // Zero is one `PutString`, which flushes before it.
        let buffer = after(Stream::AsnText, 4090, Piece::Real, b"{ 0, 10, 0 }");
        assert_eq!((buffer.written(), buffer.used), (4010, 92));
    }

    #[test]
    fn asn_text_octets_wrap_with_one_byte_calls() {
        // `BeginBytes` (`'`), two `PutChar`s per byte with `WrapAt(78, false)`'s end of line
        // (`PutEol(false)`: `Reserve(1)`, no indentation) between them, then `PutString("'H")`:
        // the 97th of the one-byte calls flushes 4096 - 80 bytes.
        let mut octets = b"'".to_vec();
        octets.extend(b"AB".repeat(39));
        octets.push(b'\n');
        octets.extend(b"CD".repeat(10));
        octets.extend(b"'H");
        let buffer = after(Stream::AsnText, 4000, Piece::Octets, &octets);
        assert_eq!((buffer.written(), buffer.used), (4016, 80 + 4 + 2));
    }

    #[test]
    fn a_text_post_after_a_wrap_of_its_string_counts_the_wrap() {
        // The string's wrap (`PutEolAtWordEnd`) inserted one end of line, `Reserve(1)`, before
        // the replaced byte: the post after 112 one-byte calls follows the flush at the 107th.
        let mut string = b"\"".to_vec();
        string.extend([b'a'; 90]);
        string.push(b'\n');
        string.extend([b'b'; 20]);
        string.extend(b"#\"");
        let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[b'x'; 3990])];
        assert_eq!(
            written_at(Stream::AsnText, before, Piece::Str, &string, &[112]),
            [4016]
        );
        // Without the flush: the same string at the start of the buffer.
        assert_eq!(
            written_at(Stream::AsnText, &[], Piece::Str, &string, &[112]),
            [0]
        );
    }

    #[test]
    fn an_xml_character_reference_is_a_put_string_and_one_byte_calls() {
        // `&#x9;`: `PutString("&#x")` fits at 4092, `PutChar('9')` ends the buffer and
        // `PutChar(';')` flushes 4095 bytes (XML keeps 1); one call of 5 bytes would flush at
        // 4092.
        let buffer = after(Stream::Xml, 4092, Piece::Str, b"&#x9;");
        assert_eq!((buffer.written(), buffer.used), (4095, 2));
        // An entity of a markup character is one `PutString`.
        let buffer = after(Stream::Xml, 4093, Piece::Str, b"&amp;");
        assert_eq!((buffer.written(), buffer.used), (4092, 6));
    }

    #[test]
    fn a_call_larger_than_the_buffer_doubles_it() {
        // `DoReserve`: a `Skip` of 5000 bytes flushes and doubles the buffer to 8192; the next
        // call that does not fit flushes again and keeps that size.
        let mut buffer = after(Stream::AsnText, 100, Piece::Markup, &[b'n'; 5000]);
        assert_eq!(
            (buffer.written(), buffer.used, buffer.capacity),
            (20, 5080, 8192)
        );
        buffer.put(Piece::Markup, &[b'n'; 4000], &[], &mut Vec::new());
        assert_eq!(
            (buffer.written(), buffer.used, buffer.capacity),
            (5020, 4080, 8192)
        );
        // Two doublings for a call of more than twice the buffer.
        let buffer = after(Stream::AsnText, 0, Piece::Markup, &[b'n'; 9000]);
        assert_eq!(
            (buffer.written(), buffer.used, buffer.capacity),
            (0, 9000, 16384)
        );
    }

    #[test]
    fn an_xml_value_attribute_takes_the_tag_end_back() {
        // `<T_b>` ends the buffer; `value` replaces its `>`: one byte of the buffer is given
        // back before the attribute.
        let mut buffer = OutBuffer::new(Stream::Xml);
        buffer.put(Piece::Markup, b"\n<T_b>", &[], &mut Vec::new());
        buffer.put(Piece::Bool, b" value=\"true\"/>", &[], &mut Vec::new());
        assert_eq!(buffer.used, 6 - 1 + 15);
        // An element whose tag opens in the same call.
        let mut buffer = OutBuffer::new(Stream::Xml);
        buffer.put(
            Piece::Enum,
            b"\n  <T_E value=\"plus\"/>",
            &[],
            &mut Vec::new(),
        );
        assert_eq!(buffer.used, 22);
    }
}
