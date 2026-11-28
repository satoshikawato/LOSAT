//! NCBI's ASN.1 text object stream (`CObjectOStreamAsn`, outfmt 8; the Seq-align part of
//! outfmt 11) over its output buffer (`COStreamBuffer`), for the calls of `Encoder`.
//!
//! Each document is written by a new stream, so the indentation and the line length start
//! at 0 for every document.
//!
//! NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:816-827
//! NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:64-72
//! The back limit keeps at least the last 80 bytes in the buffer when it is flushed; the
//! string wrap below looks back at most one line (79 columns), so it always finds the bytes
//! in the buffer. The encoder's `out` holds the document from its start, or under
//! `encoder::Traced` from where NCBI's buffer was last flushed (`Traced` hands on only what
//! NCBI wrote, `out_buffer::OutBuffer`, which keeps those 80 bytes): either way it holds every
//! byte that the wrap looks at.

use super::encoder::{Encoder, Replaced, ReplacedChar, SerialError, Stream, StreamPosition};
use super::real::write_asn_text_double;

/// The column after which a VisibleString is wrapped (`WrapAt` at 78 with `true`).
const STRING_WRAP: usize = 78;

/// The line length from which an OCTET STRING starts a new line (`WrapAt` at 78 with `false`).
const OCTET_WRAP: usize = 78;

/// The ASN.1 text encoder, appending one document to `out`.
pub struct AsnTextEncoder<'b> {
    out: &'b mut Vec<u8>,
    /// Where the document starts in `out` (the start of NCBI's buffer, `m_Buffer`).
    start: usize,
    /// `COStreamBuffer::m_IndentLevel`.
    indent: usize,
    /// `COStreamBuffer::m_LineLength`: the bytes written since the last end of line.
    line_length: usize,
    /// `CObjectOStreamAsn::m_BlockStart`.
    block_start: bool,
    error: Option<SerialError>,
    /// The bytes that `visible_str_replacing` replaced (`replace`).
    replaced: Vec<Replaced>,
    /// `line`'s count of the ends of line: those before `lines_to` (the start of the line it
    /// was on at its last call) in `out`.
    lines_to: usize,
    lines: u64,
    /// The bytes of the document after which NCBI's write to a failed stream threw
    /// (`failed_after`); `usize::MAX` for a stream that writes.
    failed_after: usize,
}

impl<'b> AsnTextEncoder<'b> {
    /// An encoder that appends a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        let start = out.len();
        Self {
            out,
            start,
            indent: 0,
            line_length: 0,
            block_start: false,
            error: None,
            replaced: Vec::new(),
            lines_to: start,
            lines: 0,
            failed_after: usize::MAX,
        }
    }

    /// An encoder of a document whose stream failed once NCBI's buffer held `held` bytes of it
    /// (`asn_types::write::failed_stream_text`): the bytes from there on are not NCBI's (the
    /// caller cuts them), and the string wrap that would come after inserts no end of line before
    /// them (`put_eol_at_word_end`'s `Reserve(1)` threw first).
    pub(crate) fn failed_after(out: &'b mut Vec<u8>, held: usize) -> Self {
        let mut enc = Self::new(out);
        enc.failed_after = held;
        enc
    }

    /// `ReplaceVisibleChar` of the byte `byte`, at the line before the byte's wrap (`Replaced`,
    /// `encoder::Traced`); `write_string` then writes `#` as it writes a byte.
    #[cold]
    #[inline(never)]
    fn replace(&mut self, byte: u8) {
        let line = self.line();
        self.replaced.push(Replaced {
            byte,
            position: StreamPosition::Line(line),
            with_string: true,
            at: self.out.len() - self.start,
        });
    }

    /// `m_Output.GetLine()`: 1 + the ends of line that the encoder holds (`encoder::Traced` adds
    /// those that `flush_to` handed on, which starts the count again). The count goes on from the
    /// start of the line of the last call: the string wrap inserts its end of line in the current
    /// line only, and the lines before it are final.
    fn line(&mut self) -> u64 {
        if self.lines_to < self.start || self.lines_to > self.out.len() {
            self.lines_to = self.start;
            self.lines = 0;
        }
        let tail = &self.out[self.lines_to..];
        if let Some(last) = tail.iter().rposition(|&byte| byte == b'\n') {
            self.lines += tail[..=last].iter().filter(|&&byte| byte == b'\n').count() as u64;
            self.lines_to += last + 1;
        }
        1 + self.lines
    }

    /// `COStreamBuffer::PutString` / `PutChar`: bytes without an end of line.
    fn put(&mut self, bytes: &[u8]) {
        self.out.extend_from_slice(bytes);
        self.line_length += bytes.len();
    }

    /// `COStreamBuffer::PutEol` with its flag true: an end of line and the indentation.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:356-379
    fn eol(&mut self) {
        self.out.push(b'\n');
        self.line_length = 0;
        let indent = self.indent;
        self.out.resize(self.out.len() + indent, b' ');
        self.line_length += indent;
    }

    /// `CObjectOStreamAsn::StartBlock`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:539-544
    /// (`IncIndentLevel`, whose `step` defaults to 2, include/util/strbuffer.hpp.)
    fn start_block(&mut self) {
        self.put(b"{");
        self.indent += 2;
        self.block_start = true;
    }

    /// `CObjectOStreamAsn::EndBlock`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:546-552
    fn end_block(&mut self) {
        self.indent -= 2;
        self.eol();
        self.put(b"}");
        self.block_start = false;
    }

    /// `CObjectOStreamAsn::NextElement`.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:554-561
    fn next_element(&mut self) {
        if self.block_start {
            self.block_start = false;
        } else {
            self.put(b",");
        }
        self.eol();
    }

    /// `COStreamBuffer::WrapAt` of `lineLength`, flag false: an end of line without indentation
    /// when the line holds `line_length` bytes or more.
    fn wrap_at_line_length(&mut self, line_length: usize) {
        if self.line_length >= line_length {
            self.out.push(b'\n');
            self.line_length = 0;
        }
    }

    /// `COStreamBuffer::PutEolAtWordEnd` of `lineLength`: inserts an end of line (without
    /// indentation) before the last white space or `'` of the current line at a column up to
    /// `line_length`, looking back no further than a `"` or an end of line; without one, at
    /// column `line_length`, moved back before any `"` that precedes it (a doubled `""` is
    /// not split). The bytes after the inserted end of line start the next line.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:1084-1139
    /// `isspace` is that of the C locale; the bytes of the stream are visible characters,
    /// spaces and ends of line.
    // Inlined into both forms of `write_string`, as into its one caller before the replacing form.
    #[inline(always)]
    fn put_eol_at_word_end(&mut self, line_length: usize) {
        if self.out.len() - self.start >= self.failed_after {
            return;
        }
        let buffer = self.start;
        let mut line_pos = self.line_length;
        let mut pos = self.out.len();
        let mut good_place = false;
        while pos > buffer && line_pos > 0 {
            pos -= 1;
            line_pos -= 1;
            let byte = self.out[pos];
            if line_pos <= line_length
                && (matches!(byte, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
                    || byte == b'\'')
            {
                good_place = true;
                break;
            } else if byte == b'\n' || byte == b'"' {
                break;
            }
        }
        if pos > buffer && self.out[pos - 1] == b'\n' {
            good_place = false;
        }
        if !good_place {
            if line_pos < line_length {
                pos += line_length - line_pos;
            }
            while pos > buffer && self.out[pos - 1] == b'"' {
                pos -= 1;
            }
            if pos == buffer {
                while pos < self.out.len() && self.out[pos] == b'"' {
                    pos += 1;
                }
            }
        }
        let count = self.out.len() - pos;
        self.out.insert(pos, b'\n');
        self.line_length = count;
    }
}

impl Encoder for AsnTextEncoder<'_> {
    const STREAM: Option<super::encoder::Stream> = Some(super::encoder::Stream::AsnText);

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:105-111
    /// (`WriteId` brackets a name with a space, `<` or `:`, objostrasn.cpp:489-504; a type
    /// name has none.)
    fn begin_document(&mut self, type_name: &str) {
        self.put(type_name.as_bytes());
        self.put(b" ::= ");
    }

    /// The separator `"\n"` that the auto separator writes after the object
    /// (`CObjectOStream::Write`, objostr.cpp:623-624).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:941-945
    fn end_document(&mut self) {
        self.put(b"\n");
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:669-672
    fn begin_sequence(&mut self, _type_name: &str) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:674-677
    fn end_sequence(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:679-684
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:646-667
    /// (The member names of the ASN.1 modules start with a lower-case letter.)
    fn member(&mut self, name: &str, _tag: u8) {
        self.next_element();
        self.put(name.as_bytes());
        self.put(b" ");
    }

    fn end_member(&mut self) {}

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:563-566
    fn begin_list(&mut self, _set: bool) {
        self.start_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:573-576
    fn element(&mut self) {
        self.next_element();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:568-571
    fn end_list(&mut self) {
        self.end_block();
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:834-868
    /// The variant's name and a space, on the line of what precedes it.
    fn variant(&mut self, _type_name: &str, name: &str, _tag: u8) {
        self.block_start = false;
        self.put(name.as_bytes());
        self.put(b" ");
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:842-848
    fn end_variant(&mut self) {
        self.block_start = false;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:155-158
    fn int(&mut self, value: i32) {
        let mut digits = itoa_i32(value);
        let text = digits.text();
        self.put(text);
    }

    /// `WriteInt8`: the decimal digits, as `WriteInt4`.
    fn big_int(&mut self, value: i64) {
        self.put(value.to_string().as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:113-127
    fn enumerated(&mut self, name: &str, _value: i32) {
        self.put(name.as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:127-131
    fn named_int(&mut self, name: &str, _value: i32) {
        self.put(name.as_bytes());
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:140-146
    fn boolean(&mut self, value: bool) {
        self.put(if value { b"TRUE" } else { b"FALSE" });
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:381-393
    /// A byte outside 0x20..0x7E (`GoodVisibleChar`, include/serial/impl/objistr.inl:456-459)
    /// is NCBI's replacement with an error message: here an error (`write_string`).
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string::<CHECK>(value);
    }

    /// `visible_str` that writes NCBI's `#` for a byte outside 0x20..0x7E and records it
    /// (`replace`).
    fn visible_str_replacing(&mut self, value: &[u8]) {
        self.write_string::<REPLACE>(value);
    }

    /// `string_store` (`visible_str`) that writes NCBI's `#` (`visible_str_replacing`).
    fn string_store_replacing(&mut self, value: &[u8]) {
        self.write_string::<REPLACE>(value);
    }

    /// `WriteDouble2` of `data` with `DBL_DIG` (`real::write_asn_text_double`).
    fn real(&mut self, value: f64) {
        let before = self.out.len();
        write_asn_text_double(value, self.out);
        self.line_length += self.out.len() - before;
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:868-895
    /// NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:381-393
    fn octets(&mut self, value: &[u8]) {
        const HEX: &[u8; 16] = b"0123456789ABCDEF";
        self.put(b"'");
        for &byte in value {
            self.wrap_at_line_length(OCTET_WRAP);
            self.put(&[HEX[usize::from(byte >> 4)], HEX[usize::from(byte & 0xf)]]);
        }
        self.wrap_at_line_length(OCTET_WRAP);
        self.put(b"'H");
    }

    /// `WriteStringStore` writes the bytes as `WriteString` does.
    fn string_store(&mut self, value: &[u8]) {
        self.visible_str(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:315-318
    fn null(&mut self) {
        self.put(b"NULL");
    }

    /// `WriteString` with `eFNP_Allow` (`write_string`).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:459-467
    fn utf8_str(&mut self, value: &[u8]) {
        self.write_string::<ALLOW>(value);
    }

    /// The member's name was written (`member`); the bytes follow as they were read, and the line
    /// length stays what it was (`COStreamBuffer::Write` counts no byte of the line).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:716-733
    fn delayed(&mut self, text: &[u8]) {
        self.out.extend_from_slice(text);
    }

    fn unread_delayed(&mut self) {
        if self.error.is_none() {
            self.error = Some(SerialError::Unwritten(
                "a member read into a delay buffer of another format whose read LOSAT does not report \
                 as NCBI's or that posts a message other than a replaced character",
            ));
        }
    }

    /// The document's bytes stop where NCBI's exception leaves the stream, before the member
    /// (`encoder::Traced`): the encoder's error.
    fn delayed_read_failed(&mut self, _error: &crate::report::serial::encoder::DelayedReadError) {
        if self.error.is_none() {
            self.error = Some(SerialError::DelayedRead);
        }
    }

    fn delayed_posts(&mut self, posts: &[ReplacedChar]) {
        if let (None, Some(first)) = (&self.error, posts.first()) {
            self.error = Some(SerialError::NonVisibleChar(first.byte));
        }
    }

    fn buffered(&self) -> usize {
        self.out.len() - self.start
    }

    /// Keeps the last end of line and the line after it: `put_eol_at_word_end` looks back to
    /// the start of the current line and to the byte before it (an end of line stops it, and
    /// an end of line right before a place stops it from adding a second), and the buffer
    /// start of NCBI's stream only matters when the line is the first one.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        if all {
            sink.write_all(&self.out[self.start..])?;
            self.out.truncate(self.start);
            return Ok(());
        }
        let Some(last_eol) = self.out[self.start..]
            .iter()
            .rposition(|&byte| byte == b'\n')
        else {
            return Ok(());
        };
        let keep_from = self.start + last_eol;
        if keep_from == self.start {
            return Ok(());
        }
        sink.write_all(&self.out[self.start..keep_from])?;
        self.out.drain(self.start..keep_from);
        // `line` counts again from the bytes held (`encoder::Traced` adds those handed on).
        (self.lines_to, self.lines) = (self.start, 0);
        Ok(())
    }

    fn held(&self) -> &[u8] {
        &self.out[self.start..]
    }

    fn hand_on(&mut self, sink: &mut dyn std::io::Write, count: usize) -> std::io::Result<()> {
        let end = self.start + count;
        sink.write_all(&self.out[self.start..end])?;
        self.out.drain(self.start..end);
        // `line` counts again from the bytes held (`encoder::Traced` adds those handed on).
        (self.lines_to, self.lines) = (self.start, 0);
        Ok(())
    }

    /// The document's bytes stop where NCBI's stream throws (`encoder::NullElement`): the
    /// encoder's error; what it writes after is not NCBI's (`encoder::Traced` writes nothing).
    fn null_element(&mut self, _element_type: &str) {
        if self.error.is_none() {
            self.error = Some(SerialError::NullElement);
        }
    }

    /// `line N` of the bytes held (`line`).
    fn position(&mut self) -> Option<StreamPosition> {
        Some(StreamPosition::Line(self.line()))
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }

    fn take_replaced(&mut self) -> Vec<Replaced> {
        std::mem::take(&mut self.replaced)
    }
}

/// The line `line` (its bytes from the start of the line to where a VisibleString starts) and the
/// string `value` after it, as the encoder writes them (`write_string`, with the wraps of
/// `put_eol_at_word_end`, which may also break the line before the opening quote): for
/// `crate::wording::losat_output_bytes`, which writes an archive's message again in LOSAT's words.
/// `value` holds visible bytes only (it was read from such a string).
pub(crate) fn visible_string_line(line: &[u8], value: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(line.len() + value.len() + 8);
    let mut enc = AsnTextEncoder::new(&mut out);
    // The end of the line before: the wrap looks back to it, as in the document.
    enc.out.push(b'\n');
    enc.out.extend_from_slice(line);
    enc.line_length = line.len();
    enc.write_string::<CHECK>(value);
    debug_assert!(enc.error.is_none());
    out.split_off(1)
}

/// `write_string`'s handling of a byte that NCBI's stream replaces in a VisibleString: the error
/// `NonVisibleChar` (`visible_str`), `#` and its record (`visible_str_replacing`), or the byte
/// as it is (`utf8_str`, `eFNP_Allow`).
const CHECK: u8 = 0;
const REPLACE: u8 = 1;
const ALLOW: u8 = 2;

impl AsnTextEncoder<'_> {
    /// The loop of `WriteString` (cited at `visible_str`): with `CHECK` a byte outside
    /// 0x20..0x7E is the error `NonVisibleChar` and is not written, with `REPLACE` the byte is
    /// recorded at its line (before its wrap, `replace`) and `#` is written in its place, with
    /// `ALLOW` it is written as it is (a zero byte left out).
    #[inline(always)]
    fn write_string<const MODE: u8>(&mut self, value: &[u8]) {
        self.put(b"\"");
        for &byte in value {
            let byte = if MODE == ALLOW {
                // A NUL byte is skipped (objostrasn.cpp:436-438).
                if byte == 0 {
                    continue;
                }
                byte
            } else if !Stream::AsnText.replaces(byte) {
                byte
            } else if MODE == REPLACE {
                self.replace(byte);
                b'#'
            } else {
                if self.error.is_none() {
                    self.error = Some(SerialError::NonVisibleChar(byte));
                }
                continue;
            };
            if self.line_length > STRING_WRAP {
                self.put_eol_at_word_end(STRING_WRAP);
            }
            self.put(&[byte]);
            if byte == b'"' {
                self.put(b"\"");
            }
        }
        self.put(b"\"");
    }
}

/// The decimal digits of an `Int4` (`COStreamBuffer::PutInt4`, NStr's `IntToString`).
pub(crate) struct DecimalI32 {
    buffer: [u8; 11],
    start: usize,
}

impl DecimalI32 {
    pub(crate) fn text(&mut self) -> &[u8] {
        &self.buffer[self.start..]
    }
}

/// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:951-975
/// `i32::MIN` is written as its digits. NCBI's `-n` overflows there and its build writes `-8`
/// (the archive's `WordThreshold` of a `-threshold` past `Int4`, NCBI BLAST+ 2.17.0); LOSAT does
/// not reproduce that overflow. Positions, lengths, counts and scores of BLAST do not reach it.
pub(crate) fn itoa_i32(value: i32) -> DecimalI32 {
    let mut buffer = [0u8; 11];
    let mut pos = buffer.len();
    let mut n = value.unsigned_abs();
    loop {
        pos -= 1;
        buffer[pos] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if value < 0 {
        pos -= 1;
        buffer[pos] = b'-';
    }
    DecimalI32 { buffer, start: pos }
}
