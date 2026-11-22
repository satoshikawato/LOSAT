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
//! NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-947
//! NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:281-339
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:84-87
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:548-550
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
    /// `WriteBool` and `WriteNull`: `PutString` of `TRUE`, `FALSE` or `NULL`.
    Bool,
    /// `WriteEnum` (an ENUMERATED value or an INTEGER with named values).
    Enum,
    /// `WriteString` / `WriteStringStore`.
    Str,
    /// An OCTET STRING (`BeginBytes`, `WriteBytes`, `EndBytes`).
    Octets,
    /// The text of a delay buffer (`Encoder::delayed`), copied with
    /// `COStreamBuffer::Write(CByteSourceReader&)`.
    Raw,
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
    /// The bytes that the buffer held when it first wrote to the output (`first_write`).
    first_write: Option<usize>,
}

impl OutBuffer {
    pub(super) fn new(stream: Stream) -> Self {
        Self {
            stream,
            capacity: 4096,
            used: 0,
            written: 0,
            first_write: None,
        }
    }

    /// The bytes of the document that NCBI has written to its output.
    pub(super) fn written(&self) -> u64 {
        self.written
    }

    /// The bytes at the start of the document that the buffer held when it first wrote to its
    /// output (a `FlushBuffer` with bytes to write, or the end of the document), `None` before:
    /// all that an output stream that has failed gets of the document (the object manager's
    /// messages after its diagnostic stream failed, `pssm_search::annot_index`). The write to
    /// such a stream throws before the buffer changes, and the object stream's destructor writes
    /// what the buffer holds with the stream's state cleared (`Close`, `DefaultFlush`, `Flush`).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:893-897
    /// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:317-336
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:910-918
    pub(super) fn first_write(&self) -> Option<usize> {
        self.first_write
    }

    /// `FlushBuffer(false)`.
    fn flush_buffer(&mut self) {
        let leave = self.stream.back_limit();
        if self.used < leave {
            return;
        }
        if self.used > leave && self.first_write.is_none() {
            self.first_write = Some(self.used);
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

    /// `Write(CByteSourceReader&)` of `length` bytes: the buffer filled up and flushed whenever it
    /// is full, without growing it, and flushed once more when the bytes end exactly at its end
    /// (the reader's last `Read` finds no room first).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1164-1182
    /// (The reader of a delay buffer reads at most one of its chunks per call: the buffer fills
    /// and flushes at the same bytes.)
    fn write_reader(&mut self, mut length: usize) {
        loop {
            if self.used == self.capacity {
                self.flush_buffer();
            }
            let count = (self.capacity - self.used).min(length);
            if count == 0 {
                return;
            }
            self.used += count;
            length -= count;
        }
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
    pub(super) fn put(
        &mut self,
        piece: Piece,
        bytes: &[u8],
        posts: &[usize],
        written: &mut Vec<u64>,
    ) {
        if piece == Piece::End {
            // `EndOfWrite`'s `FlushBuffer`, then the separator.
            if self.used > 0 && self.first_write.is_none() {
                self.first_write = Some(self.used);
            }
            self.written += (self.used + bytes.len()) as u64;
            self.used = 0;
            return;
        }
        // A call that wrote nothing (a value that the stream refused, `SerialError`) made no
        // buffer call.
        if bytes.is_empty() {
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
            // `PutInt4`, `PutInt8`, `PutString("TRUE")`, `PutString("NULL")`.
            Piece::Number | Piece::Bool => self.skip(bytes.len()),
            Piece::Raw => self.write_reader(bytes.len()),
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:798-815
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
            Piece::Header | Piece::Markup | Piece::Raw | Piece::End => {}
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

/// The bytes of `document`, a document of NCBI's XML stream in its standard-XML mode without a
/// file header (a `BlastOutput2` of outfmt 16), that NCBI's buffer had written to its output at
/// each offset of `posts` (ascending): the stream's calls taken back from the document's bytes
/// as `encoder::Traced` would give them, a run of tags as `Piece::Markup` and the text of an
/// element as `Piece::Str` where `is_string` says so of the element's tag, else as a number. The
/// XML writers of outfmt 16 make their reports without `Traced` (`report::blast_xml2::writer`),
/// and only a report with a replaced control character asks this.
pub(crate) fn xml_written_at(
    document: &[u8],
    posts: &[usize],
    is_string: impl Fn(&[u8]) -> bool,
) -> Vec<u64> {
    let mut buffer = OutBuffer::new(Stream::Xml);
    let mut written = Vec::with_capacity(posts.len());
    let mut next_post = 0;
    let mut i = 0;
    let mut tag: &[u8] = &[];
    while i < document.len() && written.len() < posts.len() {
        let start = i;
        // The tags (ends of line with their indentation, opening and closing tags).
        loop {
            match document.get(i) {
                Some(b'\n') => {
                    i += 1;
                    i += document[i..].iter().take_while(|&&b| b == b' ').count();
                }
                Some(b'<') => {
                    let end = document[i..]
                        .iter()
                        .position(|&b| b == b'>')
                        .map_or(document.len(), |at| i + at + 1);
                    if document.get(i + 1) != Some(&b'/') {
                        tag = &document[i + 1..end.saturating_sub(1).max(i + 1)];
                    }
                    i = end;
                }
                _ => break,
            }
        }
        let text_end = document[i..]
            .iter()
            .position(|&b| b == b'<')
            .map_or(document.len(), |at| i + at);
        let piece = if text_end == i {
            Piece::Markup
        } else if is_string(tag) {
            Piece::Str
        } else {
            Piece::Number
        };
        let mut here = Vec::new();
        while next_post < posts.len() && posts[next_post] < text_end.max(i + 1) {
            here.push(posts[next_post] - start);
            next_post += 1;
        }
        buffer.put(piece, &document[start..text_end], &here, &mut written);
        i = text_end.max(i);
        if i == start {
            break;
        }
    }
    written.resize(posts.len(), buffer.written());
    written
}

#[cfg(test)]
#[path = "out_buffer_tests.rs"]
mod tests;
