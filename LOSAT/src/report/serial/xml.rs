//! NCBI's XML object stream (`CObjectOStreamXml`, outfmt 5; later the XML2 formats 14 and 16
//! in its standard-XML mode) for the calls of `Encoder`.
//!
//! BLAST's XML writer (`BlastXML_FormatReport`, algo/blast/format/blastxml_format.cpp)
//! serializes a document in pieces and joins them, so this encoder has three kinds of start:
//! `document` (the XML declaration and the DOCTYPE, then the root), `object` (one object
//! written with `WriteObject`: the first tag after an end of line, no declaration, no final
//! end of line) and `fragment` (one object written with `Write` whose declaration line was
//! erased: the first tag without an end of line, and the end of line of `EndOfWrite` after the
//! last). The caller cuts and joins the pieces (`report::blast_xml::writer`) and may add text
//! between the calls with `raw`.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1009-1042
//! ```c++
//! void CObjectOStreamXml::OpenTagStart(void)
//! {
//!     ...
//!         if (m_SkipIndent) {
//!             m_SkipIndent = false;
//!         } else {
//!             m_Output.PutEol(false);
//!             m_Output.PutIndent();
//!         }
//!         m_Output.PutChar('<');
//!     ...
//!     m_EndTag = false;
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1069-1079
//! ```c++
//! void CObjectOStreamXml::CloseTagStart(void)
//! {
//!     _ASSERT(m_LastTagAction != eTagSelfClosed);
//!     m_Output.DecIndentLevel();
//!     if (m_EndTag && !m_SkipIndent) {
//!         m_Output.PutEol(false);
//!         m_Output.PutIndent();
//!     }
//!     m_Output.PutString("</");
//! ```
//! (The indentation step is 2: `COStreamBuffer::IncIndentLevel`, util/strbuffer.inl:252-264.)
//! So an opening tag is an end of line, the indentation and `<tag>`; a closing tag has its
//! own line only when the last thing written was a closing tag (`m_EndTag`): a primitive
//! member is one line, `<T_m>value</T_m>`.

use super::encoder::{Encoder, SerialError};
use super::real::write_xml_double;

/// What the stream writes at the start and end of a document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Piece {
    /// `Write` of a root object: the declaration, the DOCTYPE and, at the end, `EndOfWrite`'s
    /// end of line.
    Document,
    /// `WriteObject`: neither.
    Object,
    /// `Write` with its first line erased.
    Fragment,
}

/// The XML encoder, appending to `out`.
pub struct XmlEncoder<'b> {
    out: &'b mut Vec<u8>,
    piece: Piece,
    /// The module of the root type, for the DOCTYPE of a document (`NCBI-BlastOutput`).
    module: &'static str,
    /// `COStreamBuffer::m_IndentLevel`.
    indent: usize,
    /// `CObjectOStreamXml::m_EndTag`: the last thing written was a closing tag.
    end_tag: bool,
    /// `CObjectOStreamXml::m_SkipIndent`: the next opening tag is not preceded by an end of line.
    skip_indent: bool,
    /// The names of the open tags, innermost last: a class value's type name (`begin_sequence`)
    /// or a member's tag (`member`), each one a slice of `names` between two `marks`.
    names: Vec<u8>,
    marks: Vec<usize>,
    error: Option<SerialError>,
}

impl<'b> XmlEncoder<'b> {
    fn new(out: &'b mut Vec<u8>, piece: Piece, module: &'static str) -> Self {
        Self {
            out,
            piece,
            module,
            indent: 0,
            end_tag: false,
            skip_indent: piece == Piece::Fragment,
            names: Vec::new(),
            marks: Vec::new(),
            error: None,
        }
    }

    /// A document of a root type of the ASN.1 module `module`: `begin_document` writes the
    /// declaration and the DOCTYPE.
    pub fn document(out: &'b mut Vec<u8>, module: &'static str) -> Self {
        Self::new(out, Piece::Document, module)
    }

    /// One object written by `WriteObject`: the root tag follows an end of line.
    pub fn object(out: &'b mut Vec<u8>) -> Self {
        Self::new(out, Piece::Object, "")
    }

    /// One object written by `Write` with the declaration line erased: the root tag without an
    /// end of line before it, an end of line after it.
    pub fn fragment(out: &'b mut Vec<u8>) -> Self {
        Self::new(out, Piece::Fragment, "")
    }

    /// Text of the caller between two calls (the joints of the pieces).
    pub fn raw(&mut self, bytes: &[u8]) {
        self.out.extend_from_slice(bytes);
    }

    fn push_name(&mut self, parts: &[&[u8]]) {
        self.marks.push(self.names.len());
        for part in parts {
            self.names.extend_from_slice(part);
        }
    }

    fn pop_name(&mut self) -> usize {
        self.marks.pop().expect("a close call matches an open call")
    }

    /// `OpenTagStart`, the name, `OpenTagEnd`.
    fn open_tag(&mut self, name_start: usize) {
        if self.skip_indent {
            self.skip_indent = false;
        } else {
            self.out.push(b'\n');
            self.out.resize(self.out.len() + self.indent, b' ');
        }
        self.out.push(b'<');
        self.out.extend_from_slice(&self.names[name_start..]);
        self.out.push(b'>');
        self.indent += 2;
        self.end_tag = false;
    }

    /// `CloseTagStart`, the name, `CloseTagEnd`; pops the name.
    fn close_tag(&mut self) {
        let name_start = self.pop_name();
        self.indent -= 2;
        if self.end_tag && !self.skip_indent {
            self.out.push(b'\n');
            self.out.resize(self.out.len() + self.indent, b' ');
        }
        self.out.extend_from_slice(b"</");
        self.out.extend_from_slice(&self.names[name_start..]);
        self.out.push(b'>');
        self.names.truncate(name_start);
        self.end_tag = true;
        self.skip_indent = false;
    }

    fn fail(&mut self, error: SerialError) {
        if self.error.is_none() {
            self.error = Some(error);
        }
    }
}

impl Encoder for XmlEncoder<'_> {
    /// `WriteFileHeader`: the declaration and the DOCTYPE of a document.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:232-276
    /// ```c++
    /// void CObjectOStreamXml::WriteFileHeader(TTypeInfo type)
    /// {
    ///     m_Output.PutString("<?xml version=\"1.0\"");
    ///     if (m_Encoding != eEncoding_Unknown) { ... m_Output.PutString(" encoding=\"") ... }
    ///     m_Output.PutString("?>");
    ///     ...
    ///             m_Output.PutString("\n<!DOCTYPE ");
    ///             ... " PUBLIC \"" ... "-//NCBI//" ... module name with '-' as ' ' ... "/EN\"" ...
    /// ```
    /// (`SetEncoding(eEncoding_Ascii)`, blastxml_format.cpp:786: `US-ASCII`; the DTD prefix
    /// `SetDefaultDTDFilePrefix`, :789.)
    fn begin_document(&mut self, type_name: &str) {
        if self.piece != Piece::Document {
            return;
        }
        self.out
            .extend_from_slice(b"<?xml version=\"1.0\" encoding=\"US-ASCII\"?>");
        self.out.extend_from_slice(b"\n<!DOCTYPE ");
        self.out.extend_from_slice(type_name.as_bytes());
        self.out.extend_from_slice(b" PUBLIC \"-//NCBI//");
        self.out
            .extend_from_slice(self.module.replace('-', " ").as_bytes());
        self.out
            .extend_from_slice(b"/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/");
        self.out
            .extend_from_slice(self.module.replace('-', "_").as_bytes());
        self.out.extend_from_slice(b".dtd\">");
    }

    /// `EndOfWrite`: an end of line after a document or a fragment.
    fn end_document(&mut self) {
        if self.piece != Piece::Object {
            self.out.push(b'\n');
        }
    }

    /// `BeginClass`: `<Type>`.
    fn begin_sequence(&mut self, type_name: &str) {
        self.push_name(&[type_name.as_bytes()]);
        let start = *self.marks.last().expect("just pushed");
        self.open_tag(start);
    }

    /// `EndClass`: `</Type>`.
    fn end_sequence(&mut self) {
        self.close_tag();
    }

    /// `BeginClassMember`: `<Type_member>` (`PrintTagName`, objostrxml.cpp:1089-1138: the
    /// class's name, an underscore and the member's name, outside the standard-XML mode).
    fn member(&mut self, name: &str, _tag: u8) {
        let class_start = self.marks.last().copied().unwrap_or(0);
        let class_end = self.names.len();
        let start = class_end;
        self.marks.push(start);
        self.names.extend_from_within(class_start..class_end);
        self.names.push(b'_');
        self.names.extend_from_slice(name.as_bytes());
        self.open_tag(start);
    }

    /// `EndClassMember`.
    fn end_member(&mut self) {
        self.close_tag();
    }

    /// `BeginContainer`: the member's tag holds the elements; a container writes nothing of
    /// its own.
    fn begin_list(&mut self, _set: bool) {}

    /// `BeginContainerElement`: a class element writes its own `<Type>` tag (`begin_sequence`).
    fn element(&mut self) {}

    /// `EndContainer`.
    fn end_list(&mut self) {}

    fn variant(&mut self, _name: &str, _tag: u8) {
        self.fail(SerialError::XmlUnsupported("CHOICE"));
    }

    fn end_variant(&mut self) {}

    /// `WriteInt4`: `PutInt4`, the decimal digits.
    fn int(&mut self, value: i32) {
        let mut buffer = itoa_i64(i64::from(value));
        self.out.extend_from_slice(buffer.text());
    }

    /// `WriteInt8`: `PutInt8`, the decimal digits.
    fn big_int(&mut self, value: i64) {
        let mut buffer = itoa_i64(value);
        self.out.extend_from_slice(buffer.text());
    }

    fn enumerated(&mut self, _name: &str, _value: i32) {
        self.fail(SerialError::XmlUnsupported("ENUMERATED value"));
    }

    fn boolean(&mut self, _value: bool) {
        self.fail(SerialError::XmlUnsupported("BOOLEAN"));
    }

    /// `WriteString`: the escaped text.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:502-549
    /// ```c++
    /// void CObjectOStreamXml::WriteEscapedChar(char c)
    /// {
    ///     switch ( c ) {
    ///     case '&':  m_Output.PutString("&amp;");  break;
    ///     case '<':  m_Output.PutString("&lt;");   break;
    ///     case '>':  m_Output.PutString("&gt;");   break;
    ///     case '\'': m_Output.PutString("&apos;"); break;
    ///     case '"':  m_Output.PutString("&quot;"); break;
    ///     default:
    ///         if ((unsigned int)c < 0x20) {
    ///             m_Output.PutString("&#x");
    ///             ... hi (if not 0) and lo in lower-case hex ... ';'
    /// ```
    /// A byte of 0x01-0x1F other than TAB, LF and CR is replaced by `#` with an error
    /// message (`x_VerifyChar`, `BAD_CHAR`), and a byte >= 0x80 makes the US-ASCII stream
    /// throw (`CUtf8::SymbolToChar`, corelib/ncbistr.cpp:6900-6903); both are errors here (the
    /// callers refuse such strings first). A NUL byte is dropped (`if (s != '\0')`).
    fn visible_str(&mut self, value: &[u8]) {
        for &byte in value {
            match byte {
                b'&' => self.out.extend_from_slice(b"&amp;"),
                b'<' => self.out.extend_from_slice(b"&lt;"),
                b'>' => self.out.extend_from_slice(b"&gt;"),
                b'\'' => self.out.extend_from_slice(b"&apos;"),
                b'"' => self.out.extend_from_slice(b"&quot;"),
                0 => {}
                0x09 | 0x0a | 0x0d => {
                    self.out.extend_from_slice(b"&#x");
                    self.out.push(b"0123456789abcdef"[usize::from(byte & 0xf)]);
                    self.out.push(b';');
                }
                0x01..=0x1f => self.fail(SerialError::NonVisibleChar(byte)),
                0x80..=0xff => self.fail(SerialError::NonAsciiChar(byte)),
                _ => self.out.push(byte),
            }
        }
    }

    /// `WriteDouble`: `WriteDouble2(data, DBL_DIG)`, `NStr::DoubleToStringPosix(data, 15)`
    /// (`real::write_xml_double`).
    fn real(&mut self, value: f64) {
        write_xml_double(value, self.out);
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }
}

/// The decimal digits of an `Int8`.
struct DecimalI64 {
    buffer: [u8; 20],
    start: usize,
}

impl DecimalI64 {
    fn text(&mut self) -> &[u8] {
        &self.buffer[self.start..]
    }
}

fn itoa_i64(value: i64) -> DecimalI64 {
    let mut buffer = [0u8; 20];
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
    DecimalI64 { buffer, start: pos }
}
