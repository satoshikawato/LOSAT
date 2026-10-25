//! NCBI's XML object stream (`CObjectOStreamXml`, outfmt 5, and the XML2 formats 14 and 16 in
//! its standard-XML mode) for the calls of `Encoder`.
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
//! The XML2 formats (`CBlastOutput2`, objects/blastxml2/blastxml2.asn) use the stream in its
//! standard-XML mode (`SetEnforcedStdXml`, blastxml2_format.cpp:723, :738, :789), whose tags
//! differ from outfmt 5's: a member's tag is the member name alone, a CHOICE writes its type tag and
//! then the variant as a member, a SEQUENCE OF VisibleString repeats the member's tag for each
//! string. A 14 side file is a `xml2_document` (the declaration, the root with its namespace
//! attributes) and a 16 report an `xml2_report` (the root after an end of line, a final end of
//! line, no declaration).
//!
//! The archive of outfmt 11 (`ARCHIVE_FORMAT=xml`, `Blast4-archive`) is the same stream used on
//! NCBI's own objects (`XmlEncoder::archive`, the object mode): what outfmt 5 never needs is
//! written as well. The type info tells NCBI's stream which values are of a named type, which
//! the `Encoder` calls say only for a SEQUENCE (`begin_sequence`) and a CHOICE (`variant`), so
//! a value of any other named type is wrapped in `begin_named` / `end_named`. In that mode:
//! a named CHOICE writes its own tag and the variant's (`<Seq-id><Seq-id_local>`), an anonymous
//! one (its name is the open member's tag, `Score_value`) only the variant's
//! (`<Score_value_int>`); an ENUMERATED or BOOLEAN value is an attribute that takes the place
//! of the `>` of the innermost open tag (`<Seq-align_type value="partial"/>`, the closing tag
//! is not written); a primitive element of a list is wrapped in `<Parent_E>`, an element that
//! writes its own tag (a SEQUENCE, a CHOICE or a named type) is not; an OCTET STRING is
//! unwrapped upper-case hex.
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
    /// `Write` of a `Blast4-archive` (the object mode): the declaration (without an encoding)
    /// and the DOCTYPE of the module `NCBI-Blast4` without a DTD URL, and `EndOfWrite`'s end of
    /// line.
    NcbiObject,
    /// `Write` of a `BlastOutput2` root (the 14 side file): the declaration without a DOCTYPE,
    /// the root tag with the namespace attributes, and `EndOfWrite`'s end of line.
    Xml2Document,
    /// A `BlastOutput2` of the 16 document (`CBlastOStreamXml`, blastxml2_format.cpp:704-714,
    /// whose file header is empty): the root after an end of line, and `EndOfWrite`'s end of
    /// line after it.
    Xml2Report,
}

/// The attributes of the root tag of a 14 side file (`x_WriteClassNamespace`, objostrxml.cpp:
/// 284-342): the namespace, the schema instance namespace and the schema location, each on a
/// line of its own, and a line break before the `>`.
pub const XML2_ROOT_ATTRIBUTES: &str = "\n    xmlns=\"http://www.ncbi.nlm.nih.gov\"\n    xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n    xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n";

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
    /// The standard-XML mode (`SetEnforcedStdXml`): the XML2 formats.
    std_xml: bool,
    /// The root tag still to write has the namespace attributes (a 14 side file).
    root_attributes: bool,
    /// The object mode (`Piece::NcbiObject`).
    object: bool,
    /// Where this document starts in `out` (what `flush_to` hands over starts here).
    start: usize,
    /// The object mode: for each open `begin_sequence` / `variant`, how many tags it opened
    /// (an anonymous type opens none or only the variant's).
    opened: Vec<u8>,
    /// The object mode: `element` was called and the element's first call has not come yet.
    pending_elem: bool,
    /// The object mode: a `<Parent_E>` tag of a primitive element is open.
    open_e: bool,
    /// The object mode: the innermost open tag was written as `<tag value="x"/>`, so that its
    /// closing call writes nothing.
    self_closed: bool,
    /// The names of the open tags, innermost last: a class value's type name (`begin_sequence`)
    /// or a member's tag (`member`), each one a slice of `names` between two `marks`.
    names: Vec<u8>,
    marks: Vec<usize>,
    error: Option<SerialError>,
}

impl<'b> XmlEncoder<'b> {
    fn new(out: &'b mut Vec<u8>, piece: Piece, module: &'static str) -> Self {
        let start = out.len();
        Self {
            start,
            out,
            piece,
            module,
            indent: 0,
            end_tag: false,
            skip_indent: piece == Piece::Fragment,
            std_xml: matches!(piece, Piece::Xml2Document | Piece::Xml2Report),
            root_attributes: piece == Piece::Xml2Document,
            object: piece == Piece::NcbiObject,
            opened: Vec::new(),
            pending_elem: false,
            open_e: false,
            self_closed: false,
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

    /// One `Blast4-archive` document of outfmt 11 with `ARCHIVE_FORMAT=xml`, in NCBI's object
    /// mode; `Encoder::flush_to` hands its bytes to a sink element by element.
    pub fn archive(out: &'b mut Vec<u8>) -> Self {
        Self::new(out, Piece::NcbiObject, "NCBI-Blast4")
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

    /// A 14 side file: the declaration, the `BlastOutput2` root with the namespace attributes,
    /// the standard-XML mode.
    pub fn xml2_document(out: &'b mut Vec<u8>) -> Self {
        Self::new(out, Piece::Xml2Document, "")
    }

    /// One `BlastOutput2` of a 16 document, in the standard-XML mode.
    pub fn xml2_report(out: &'b mut Vec<u8>) -> Self {
        Self::new(out, Piece::Xml2Report, "")
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
        if std::mem::take(&mut self.root_attributes) {
            self.out.extend_from_slice(XML2_ROOT_ATTRIBUTES.as_bytes());
        }
        self.out.push(b'>');
        self.indent += 2;
        self.end_tag = false;
    }

    /// `CloseTagStart`, the name, `CloseTagEnd`; pops the name.
    fn close_tag(&mut self) {
        let name_start = self.pop_name();
        self.indent -= 2;
        if self.self_closed {
            // The tag was written as `<tag value="x"/>`: nothing closes it, and the next
            // closing tag is on a line of its own (`m_EndTag`).
            self.self_closed = false;
            self.names.truncate(name_start);
            self.end_tag = true;
            return;
        }
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

    /// The tag `Class_name` of a member or a variant (`PrintTagName`): the innermost open
    /// name, an underscore and `name`.
    fn child_tag(&mut self, name: &str) {
        let class_start = self.marks.last().copied().unwrap_or(0);
        let class_end = self.names.len();
        self.marks.push(class_end);
        self.names.extend_from_within(class_start..class_end);
        self.names.push(b'_');
        self.names.extend_from_slice(name.as_bytes());
        self.open_tag(class_end);
    }

    /// The first call of a primitive element of a list opens its `<Parent_E>` tag (the object
    /// mode); `element` closes it.
    #[cold]
    fn open_element(&mut self) {
        self.pending_elem = false;
        let parent_start = self.marks.last().copied().unwrap_or(0);
        let parent_end = self.names.len();
        self.marks.push(parent_end);
        self.names.extend_from_within(parent_start..parent_end);
        self.names.extend_from_slice(b"_E");
        self.open_tag(parent_end);
        self.open_e = true;
    }

    /// Closes the `<Parent_E>` tag of the previous primitive element.
    fn close_element(&mut self) {
        if self.open_e {
            self.open_e = false;
            self.close_tag();
        }
    }

    /// `WriteEnum` and `WriteBool` of the object mode: ` value="x"/>` in place of the `>`
    /// that closes the innermost open tag (NCBI backs over it, `OpenTagEndBack`).
    fn value_attribute(&mut self, value: &str) {
        if self.pending_elem {
            self.open_element();
        }
        debug_assert_eq!(self.out.last(), Some(&b'>'));
        self.out.pop();
        self.out.extend_from_slice(b" value=\"");
        self.out.extend_from_slice(value.as_bytes());
        self.out.extend_from_slice(b"\"/>");
        self.self_closed = true;
    }

    /// `begin_sequence` of the object mode: `<Type>`, unless the type is anonymous (its name,
    /// the open member's tag, has an underscore), which writes nothing.
    fn begin_object_sequence(&mut self, type_name: &str) {
        if type_name.contains('_') {
            debug_assert!(!self.pending_elem);
            debug_assert_eq!(
                &self.names[self.marks.last().copied().unwrap_or(0)..],
                type_name.as_bytes()
            );
            self.opened.push(0);
            return;
        }
        self.pending_elem = false;
        self.opened.push(1);
        self.push_name(&[type_name.as_bytes()]);
        let start = *self.marks.last().expect("just pushed");
        self.open_tag(start);
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
        if self.piece == Piece::Xml2Document {
            // `WriteFileHeader` of the standard-XML mode: no DOCTYPE (objostrxml.cpp:232-276,
            // `m_UseSchemaRef` / `m_EnforcedStdXml`).
            self.out
                .extend_from_slice(b"<?xml version=\"1.0\" encoding=\"US-ASCII\"?>");
            return;
        }
        if self.piece == Piece::NcbiObject {
            // `WriteFileHeader` of the stream without an encoding: no `encoding`, and a DOCTYPE
            // whose DTD is the bare file name (the stream has no default DTD prefix).
            self.out
                .extend_from_slice(b"<?xml version=\"1.0\"?>\n<!DOCTYPE ");
            self.out.extend_from_slice(type_name.as_bytes());
            self.out.extend_from_slice(b" PUBLIC \"-//NCBI//");
            self.out
                .extend_from_slice(self.module.replace('-', " ").as_bytes());
            self.out.extend_from_slice(b"/EN\" \"");
            self.out
                .extend_from_slice(self.module.replace('-', "_").as_bytes());
            self.out.extend_from_slice(b".dtd\">");
            return;
        }
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
        if self.object {
            self.begin_object_sequence(type_name);
            return;
        }
        self.push_name(&[type_name.as_bytes()]);
        let start = *self.marks.last().expect("just pushed");
        self.open_tag(start);
    }

    /// `EndClass`: `</Type>`.
    fn end_sequence(&mut self) {
        if self.object && self.opened.pop() == Some(0) {
            return;
        }
        self.close_tag();
    }

    /// `BeginClassMember`: `<Type_member>` (`PrintTagName`, objostrxml.cpp:1089-1138: the
    /// class's name, an underscore and the member's name, outside the standard-XML mode; in it
    /// the member's name alone).
    fn member(&mut self, name: &str, _tag: u8) {
        if self.std_xml {
            self.push_name(&[name.as_bytes()]);
            let start = *self.marks.last().expect("just pushed");
            self.open_tag(start);
            return;
        }
        self.child_tag(name);
    }

    /// `EndClassMember`.
    fn end_member(&mut self) {
        self.close_tag();
    }

    /// `BeginNamedType`: `<Type>` (the object mode).
    fn begin_named(&mut self, type_name: &str) {
        if self.object {
            self.pending_elem = false;
            self.push_name(&[type_name.as_bytes()]);
            let start = *self.marks.last().expect("just pushed");
            self.open_tag(start);
        }
    }

    /// `EndNamedType`: `</Type>`.
    fn end_named(&mut self) {
        if self.object {
            self.close_tag();
        }
    }

    /// `BeginContainer`: the member's tag holds the elements; a container writes nothing of
    /// its own.
    fn begin_list(&mut self, _set: bool) {}

    /// `BeginContainerElement`: a class element writes its own `<Type>` tag (`begin_sequence`).
    /// In the object mode the element's first call decides: a primitive is wrapped in
    /// `<Parent_E>`, which this call (or `end_list`) closes.
    fn element(&mut self) {
        if self.object {
            self.close_element();
            self.pending_elem = true;
        }
    }

    /// `EndContainer`.
    fn end_list(&mut self) {
        if self.object {
            self.close_element();
            self.pending_elem = false;
        }
    }

    /// `BeginChoice` and `BeginChoiceVariant`: in the standard-XML mode the CHOICE's type tag
    /// and the variant's name as a member tag (`<Target><db>`); outside it NCBI's stream
    /// writes an element LOSAT's outfmt 5 never needs.
    fn variant(&mut self, type_name: &str, name: &str, tag: u8) {
        if self.object {
            if type_name.contains('_') {
                // An anonymous CHOICE: the member's tag is already open.
                debug_assert!(!self.pending_elem);
                debug_assert_eq!(
                    &self.names[self.marks.last().copied().unwrap_or(0)..],
                    type_name.as_bytes()
                );
                self.opened.push(1);
            } else {
                self.pending_elem = false;
                self.opened.push(2);
                self.push_name(&[type_name.as_bytes()]);
                let start = *self.marks.last().expect("just pushed");
                self.open_tag(start);
            }
            self.child_tag(name);
            return;
        }
        if !self.std_xml {
            self.fail(SerialError::XmlUnsupported("CHOICE"));
            return;
        }
        self.begin_sequence(type_name);
        self.member(name, tag);
    }

    /// `EndChoiceVariant` and `EndChoice`.
    fn end_variant(&mut self) {
        if self.object {
            self.close_tag();
            if self.opened.pop() == Some(2) {
                self.close_tag();
            }
        } else if self.std_xml {
            self.end_member();
            self.end_sequence();
        }
    }

    /// A SEQUENCE OF VisibleString variant in the standard-XML mode: the member's tag is
    /// repeated for each string, with no wrapper around them
    /// (`<Target><subjects>Subject_1</subjects><subjects>Subject_2</subjects></Target>`).
    fn variant_string_list<I>(&mut self, type_name: &str, name: &str, tag: u8, items: I)
    where
        I: IntoIterator,
        I::Item: AsRef<[u8]>,
    {
        if !self.std_xml {
            self.fail(SerialError::XmlUnsupported("CHOICE"));
            return;
        }
        self.begin_sequence(type_name);
        for item in items {
            self.member(name, tag);
            self.visible_str(item.as_ref());
            self.end_member();
        }
        self.end_sequence();
    }

    /// `WriteInt4`: `PutInt4`, the decimal digits.
    fn int(&mut self, value: i32) {
        if self.pending_elem {
            self.open_element();
        }
        let mut buffer = itoa_i64(i64::from(value));
        self.out.extend_from_slice(buffer.text());
    }

    /// `WriteInt8`: `PutInt8`, the decimal digits.
    fn big_int(&mut self, value: i64) {
        if self.pending_elem {
            self.open_element();
        }
        let mut buffer = itoa_i64(value);
        self.out.extend_from_slice(buffer.text());
    }

    /// `WriteEnum` of the object mode: ` value="name"/>` in place of the `>` of the open tag.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:473-484
    /// ```c++
    ///                 OpenTagEndBack();
    ///                 m_Output.PutString(" value=\"");
    ///                 m_Output.PutString(valueName);
    ///                 m_Output.PutChar('"');
    ///                 if ( values.IsInteger() ) {
    ///                     OpenTagEnd();
    ///                     m_Output.PutInt4(value);
    ///                 }
    ///                 else {
    ///                     SelfCloseTagEnd();
    ///                 }
    /// ```
    fn enumerated(&mut self, name: &str, _value: i32) {
        if self.object {
            self.value_attribute(name);
        } else {
            self.fail(SerialError::XmlUnsupported("ENUMERATED value"));
        }
    }

    /// `WriteBool` of the object mode: ` value="true"/>` or ` value="false"/>`.
    fn boolean(&mut self, value: bool) {
        if self.object {
            self.value_attribute(if value { "true" } else { "false" });
        } else {
            self.fail(SerialError::XmlUnsupported("BOOLEAN"));
        }
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
        if self.pending_elem {
            self.open_element();
        }
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

    /// `WriteBytes` of the object mode: the bytes as upper-case hex digits, unwrapped
    /// (`WriteBytes`, objostrxml.cpp:1702-1709); the other pieces write no OCTET STRING.
    fn octets(&mut self, value: &[u8]) {
        if !self.object {
            self.fail(SerialError::OctetString);
            return;
        }
        if self.pending_elem {
            self.open_element();
        }
        self.out.reserve(2 * value.len());
        for &byte in value {
            self.out.push(b"0123456789ABCDEF"[usize::from(byte >> 4)]);
            self.out.push(b"0123456789ABCDEF"[usize::from(byte & 0xf)]);
        }
    }

    /// `WriteDouble`: `WriteDouble2(data, DBL_DIG)`, `NStr::DoubleToStringPosix(data, 15)`
    /// (`real::write_xml_double`).
    fn real(&mut self, value: f64) {
        if self.pending_elem {
            self.open_element();
        }
        write_xml_double(value, self.out);
    }

    fn buffered(&self) -> usize {
        if self.object {
            self.out.len() - self.start
        } else {
            0
        }
    }

    /// The object mode keeps the last byte until the document is complete: it may be the `>` of
    /// an open tag that `enumerated` or `boolean` replaces.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        if !self.object {
            return Ok(());
        }
        let end = if all {
            self.out.len()
        } else {
            self.out.len().saturating_sub(1).max(self.start)
        };
        if end > self.start {
            sink.write_all(&self.out[self.start..end])?;
            self.out.drain(self.start..end);
        }
        Ok(())
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
