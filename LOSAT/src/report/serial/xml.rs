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
//! NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1069-1079
//! (The indentation step is 2: `COStreamBuffer::IncIndentLevel`, util/strbuffer.inl:252-264.)
//! So an opening tag is an end of line, the indentation and `<tag>`; a closing tag has its
//! own line only when the last thing written was a closing tag (`m_EndTag`): a primitive
//! member is one line, `<T_m>value</T_m>`.

use super::encoder::{Encoder, Replaced, ReplacedChar, SerialError, StreamPosition};
use super::real::write_xml_double;
use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::utf8::{decode_to_ascii, Utf8Error};

/// Whether `type_name`, the name that a writer gives `begin_sequence` or `variant`, stands for an
/// anonymous ASN.1 type (NCBI's `GetName()` is empty for it): the typed writers pass the open
/// member's tag (`Score_value`), which has an underscore, and the writer of the type table's
/// values (`asn_types::write`) the type's empty name, where a named type's name (`Seq-align`)
/// has neither. The object mode writes no tag of its own for it (`begin_object_sequence`,
/// `variant`), and `encoder::Traced` keeps its `.E` frame (`WillHaveName`).
#[inline(always)]
pub(crate) fn anonymous_type(type_name: &str) -> bool {
    type_name.is_empty() || type_name.contains('_')
}

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
///
/// `OBJECT` is the object mode of the archive (`XmlEncoder::archive`, `Piece::NcbiObject`): a
/// const parameter, so that the instance of outfmt 5 / 14 / 16 (`OBJECT = false`, the default
/// of the type) compiles without any test of the object mode and runs the instructions it ran
/// before the mode existed (the archive's pending element, `value=` attribute and `_E` state
/// are only read under `OBJECT`).
pub struct XmlEncoder<'b, const OBJECT: bool = false> {
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
    /// The control characters that `visible_str_replacing` replaced (`replace`), in the object
    /// mode (`encoder::Traced` adds their stack trace).
    replaced: Vec<Replaced>,
    /// `replace`'s count of the ends of line: those before `lines_to` in `out`.
    lines_to: usize,
    lines: u64,
    /// The standard-XML mode: a bit for each depth of `marks` whose tag is a CHOICE's type tag
    /// (`variant`), for the frames of a post or an exception (`frames`).
    choices: u64,
    /// The control characters that NCBI's stream of outfmt 5 / 14 / 16 writes as `#` and posts,
    /// with their stack trace (`replace`, `XmlPost`).
    posts: Vec<XmlPost>,
    /// The exception that the stream of outfmt 5 / 14 / 16 threw (`throw`, `XmlThrow`).
    thrown: Option<Box<XmlThrow>>,
}

/// A control character that NCBI's XML stream of outfmt 5 / 14 / 16 writes as `#` and posts
/// (`ReplaceVisibleChar` with `eFNP_ReplaceAndWarn`, error 803.7): the post's text
/// (`encoder::ReplacedChar`, with the stack trace of the open tags and the stream's line) and
/// where the byte's `#` goes in the encoder's document (`at`, from its start).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XmlPost {
    pub at: usize,
    pub replaced: ReplacedChar,
}

/// A frame of NCBI's object stack that the stream's exception passes on its way out
/// (`END_OBJECT_FRAME`, `NCBI_RETHROW_SAME` with `CObjectStackFrame::GetFrameInfo`).
#[derive(Clone, Debug, PartialEq, Eq)]
enum ThrowFrame {
    /// `WriteClassMember` (objostrxml.cpp:1530-1540): `eFrameClassMember`, the member's name.
    Member(String),
    /// `WriteClass` (objostrxml.cpp:1503-1520): `eFrameClass`, the class's type.
    Class(String),
    /// `WriteChoiceDefault`'s variant frame (choice.cpp:452-459): `eFrameChoiceVariant`.
    Variant(String),
    /// `WriteChoiceDefault`'s choice frame (choice.cpp:463): `eFrameChoice`, the CHOICE's type.
    Choice(String),
    /// `CObjectOStream::Write` (objostr.cpp:610-625): `eFrameNamed`, the root's type; an object
    /// written with `WriteObject` (outfmt 5's `Iteration`) has none.
    Named(String),
}

impl ThrowFrame {
    /// The place of the frame's rethrow and its message (`GetFrameInfo`).
    ///
    /// NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:209-213
    /// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:399-410
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1528-1540
    /// NCBI reference (598d8ae6): c++/src/serial/choice.cpp:452-463
    /// (`WriteClass`'s frame ends at objostrxml.cpp:1519, `Write`'s at objostr.cpp:621.)
    fn rethrow(&self, exception: NcbiException) -> NcbiException {
        let (file, line, function, info) = match self {
            ThrowFrame::Member(name) => (
                "src/serial/objostrxml.cpp",
                1538,
                "ncbi::CObjectOStreamXml::WriteClassMember",
                format!(" Frame type= eFrameClassMember, Member name= {name}"),
            ),
            ThrowFrame::Class(name) => (
                "src/serial/objostrxml.cpp",
                1519,
                "ncbi::CObjectOStreamXml::WriteClass",
                format!(" Frame type= eFrameClass, Object type= {name}"),
            ),
            ThrowFrame::Variant(name) => (
                "src/serial/choice.cpp",
                459,
                "ncbi::CChoiceTypeInfoFunctions::WriteChoiceDefault",
                format!(" Frame type= eFrameChoiceVariant, Member name= {name}"),
            ),
            ThrowFrame::Choice(name) => (
                "src/serial/choice.cpp",
                463,
                "ncbi::CChoiceTypeInfoFunctions::WriteChoiceDefault",
                format!(" Frame type= eFrameChoice, Object type= {name}"),
            ),
            ThrowFrame::Named(name) => (
                "src/serial/objostr.cpp",
                621,
                "ncbi::CObjectOStream::Write",
                format!(" Frame type= eFrameNamed, Object type= {name}"),
            ),
        };
        exception.rethrow_same(file, line, function, info)
    }
}

/// The exception that NCBI's XML stream of outfmt 5 / 14 / 16 throws at a byte >= 0x80 that its
/// US-ASCII encoding cannot write (`CUtf8`, `ncbi_str::utf8`): which one, where the document is
/// cut (`at`, the bytes from the document's start that NCBI's stream wrote before it), and the
/// frames of NCBI's object stack that it passes, innermost first.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct XmlThrow {
    pub error: Utf8Error,
    pub at: usize,
    frames: Vec<ThrowFrame>,
}

impl XmlThrow {
    /// The exception as it leaves the stream: thrown by `CUtf8`, rethrown by each frame.
    pub fn exception(&self) -> NcbiException {
        self.frames
            .iter()
            .fold(self.error.exception(), |exception, frame| {
                frame.rethrow(exception)
            })
    }
}

/// Whether NCBI's XML stream of outfmt 5 / 14 / 16 throws where it writes the string `value`
/// (`XmlThrow`: a byte >= 0x80 that does not decode to a US-ASCII character before the string's
/// first NUL): the stream's own `write_string` on a scratch stream, for a question asked before
/// the document is written (`report::blast_xml::names::iteration_cut_before_hits`).
pub(crate) fn us_ascii_stream_throws(value: &[u8]) -> bool {
    if value.is_ascii() {
        return false;
    }
    let mut out = Vec::new();
    let mut stream = XmlEncoder::fragment(&mut out);
    stream.write_string::<true>(value);
    stream.thrown.is_some()
}

impl<'b> XmlEncoder<'b, false> {
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

    /// The length of the output written so far, where the next `raw` bytes start.
    pub fn written(&self) -> usize {
        self.out.len()
    }
}

impl<'b> XmlEncoder<'b, true> {
    /// One `Blast4-archive` document of outfmt 11 with `ARCHIVE_FORMAT=xml`, in NCBI's object
    /// mode; `Encoder::flush_to` hands its bytes to a sink element by element.
    pub fn archive(out: &'b mut Vec<u8>) -> Self {
        Self::ncbi_object(out, "NCBI-Blast4")
    }

    /// One document of a root type of the ASN.1 module `module` written as the archive is
    /// (`os` takes `MSerial_Xml` and then `object`): dustmasker's `Seq-loc` (`NCBI-Seqloc`) and
    /// `Blast-db-mask-info` (`NCBI-BlastDL`).
    pub fn ncbi_object(out: &'b mut Vec<u8>, module: &'static str) -> Self {
        Self::new(out, Piece::NcbiObject, module)
    }
}

impl<'b, const OBJECT: bool> XmlEncoder<'b, OBJECT> {
    fn new(out: &'b mut Vec<u8>, piece: Piece, module: &'static str) -> Self {
        debug_assert_eq!(OBJECT, piece == Piece::NcbiObject);
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
            opened: Vec::new(),
            pending_elem: false,
            open_e: false,
            self_closed: false,
            names: Vec::new(),
            marks: Vec::new(),
            error: None,
            replaced: Vec::new(),
            lines_to: start,
            lines: 0,
            choices: 0,
            posts: Vec::new(),
            thrown: None,
        }
    }

    /// Whether the stream of outfmt 5 / 14 / 16 replaced a control character (`take_posts`).
    #[inline]
    pub fn has_posts(&self) -> bool {
        !self.posts.is_empty()
    }

    /// The control characters that the stream of outfmt 5 / 14 / 16 replaced and posts
    /// (`XmlPost`), which the encoder forgets.
    pub fn take_posts(&mut self) -> Vec<XmlPost> {
        std::mem::take(&mut self.posts)
    }

    /// The exception that the stream of outfmt 5 / 14 / 16 threw, if any (`XmlThrow`; the error
    /// is then `SerialError::Utf8`).
    pub fn take_thrown(&mut self) -> Option<Box<XmlThrow>> {
        self.thrown.take()
    }

    /// The name of the open tag `index` (0 the outermost) in `names`.
    fn open_tag_name(&self, index: usize) -> &[u8] {
        let end = self
            .marks
            .get(index + 1)
            .copied()
            .unwrap_or(self.names.len());
        &self.names[self.marks[index]..end]
    }

    /// The name of the member or variant of the open tag `index` (odd: the tags of a class or
    /// CHOICE value and of its member or variant alternate): outside the standard-XML mode the
    /// tag is the class's name, an underscore and the member's (`PrintTagName`).
    fn open_member_name(&self, index: usize) -> String {
        let tag = self.open_tag_name(index);
        let name = if self.std_xml {
            tag
        } else {
            &tag[self.open_tag_name(index - 1).len() + 1..]
        };
        String::from_utf8_lossy(name).into_owned()
    }

    /// `CObjectStack::GetStackTrace` of the stream (`encoder::Traced::stack_trace`): the root's
    /// type, then `.name` of each open member or variant (NCBI's XML stream pushes no element
    /// frame for an element of a named type, so no `.E`: every list of these documents but
    /// `Target.subjects`, whose strings are never replaced, holds a named type).
    #[cold]
    fn asn_path(&self) -> String {
        let mut path = String::from_utf8_lossy(self.open_tag_name(0)).into_owned();
        for index in (1..self.marks.len()).step_by(2) {
            path.push('.');
            path.push_str(&self.open_member_name(index));
        }
        path
    }

    /// The frames of NCBI's object stack at the current byte, innermost first (`ThrowFrame`): the
    /// open tags alternate between a class or CHOICE value (its type tag) and its member or
    /// variant, and a document written with `Write` (all but `WriteObject`'s `object`) has the
    /// root's named frame outside them.
    #[cold]
    fn frames(&self) -> Vec<ThrowFrame> {
        let mut frames = Vec::with_capacity(self.marks.len() + 1);
        for index in (0..self.marks.len()).rev() {
            let choice = |depth: usize| self.choices & (1u64 << depth) != 0;
            frames.push(if index % 2 == 1 {
                let name = self.open_member_name(index);
                if choice(index - 1) {
                    ThrowFrame::Variant(name)
                } else {
                    ThrowFrame::Member(name)
                }
            } else {
                let name = String::from_utf8_lossy(self.open_tag_name(index)).into_owned();
                if choice(index) {
                    ThrowFrame::Choice(name)
                } else {
                    ThrowFrame::Class(name)
                }
            });
        }
        if self.piece != Piece::Object && !self.marks.is_empty() {
            frames.push(ThrowFrame::Named(
                String::from_utf8_lossy(self.open_tag_name(0)).into_owned(),
            ));
        }
        frames
    }

    /// The stream throws `error` (`XmlThrow`): the document is cut at the bytes written so far,
    /// and nothing that follows is posted. Only the first exception counts.
    #[cold]
    #[inline(never)]
    fn throw(&mut self, error: Utf8Error) {
        if self.thrown.is_none() {
            self.thrown = Some(Box::new(XmlThrow {
                error,
                at: self.out.len() - self.start,
                frames: self.frames(),
            }));
        }
        self.fail(SerialError::Utf8(error));
    }

    /// The branch of `WriteEncodedChar` for a byte >= 0x80 of the US-ASCII stream of outfmt 5 /
    /// 14 / 16: the character that the UTF-8 sequence at the start of `src` decodes to
    /// (`CUtf8::Decode`, `CUtf8::SymbolToChar`), verified and written as an ASCII byte is, or the
    /// stream's exception (`throw`).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:556-572
    /// (The strings are UTF-8, `m_StringEncoding`'s default, objostrxml.cpp:81, and the stream
    /// US-ASCII, `SetEncoding` with `eEncoding_Ascii`, blastxml_format.cpp:502, :786,
    /// blastxml2_format.cpp:721, :733.) A decoded NUL (`C0 80`) writes nothing and the string
    /// goes on.
    ///
    /// `src` is the rest of the string from such a byte on: each UTF-8 sequence is written so and
    /// the ASCII bytes between them as `write_string` writes them, until the string's end, a NUL
    /// or the throw.
    #[cold]
    #[inline(never)]
    fn encoded_tail(&mut self, mut src: &[u8]) {
        loop {
            match decode_to_ascii(src) {
                Ok((character, used)) => {
                    self.write_string::<true>(&[character]);
                    src = &src[used..];
                }
                Err(error) => {
                    self.throw(error);
                    return;
                }
            }
            let ascii = src
                .iter()
                .position(|&byte| byte >= 0x80 || byte == 0)
                .unwrap_or(src.len());
            self.write_string::<true>(&src[..ascii]);
            match src.get(ascii) {
                Some(&byte) if byte >= 0x80 => src = &src[ascii..],
                _ => return,
            }
        }
    }

    /// `x_VerifyChar` of the control character `byte`: `ReplaceVisibleChar` without the
    /// string, reported at the line (`Replaced`, `encoder::Traced`); `write_string` then writes
    /// `#`. The line is 1 + the ends of line that the encoder holds (`Traced` adds those that
    /// `flush_to` handed on, which starts the count again), counted on from the last call: the XML
    /// stream only appends (the `>` that `value_attribute` takes back is no end of line).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:548-554
    ///
    /// Outside the object mode the encoder posts by itself (`XmlPost`, with the stack trace of its
    /// open tags, `asn_path`); a `fragment`'s stream wrote the declaration's line that BLAST then
    /// erased, so its lines count one more. Nothing is posted after the stream threw (`throw`).
    #[cold]
    #[inline(never)]
    fn replace(&mut self, byte: u8) {
        if !OBJECT && self.thrown.is_some() {
            return;
        }
        if self.lines_to < self.start {
            self.lines_to = self.start;
            self.lines = 0;
        }
        self.lines_to = self.lines_to.min(self.out.len());
        self.lines += self.out[self.lines_to..]
            .iter()
            .filter(|&&byte| byte == b'\n')
            .count() as u64;
        self.lines_to = self.out.len();
        let at = self.out.len() - self.start;
        if !OBJECT {
            let erased = u64::from(self.piece == Piece::Fragment);
            self.posts.push(XmlPost {
                at,
                replaced: ReplacedChar {
                    trace: self.asn_path(),
                    byte,
                    position: StreamPosition::Line(1 + erased + self.lines),
                    string: std::sync::Arc::from(&[][..]),
                },
            });
            return;
        }
        self.replaced.push(Replaced {
            byte,
            position: StreamPosition::Line(1 + self.lines),
            with_string: false,
            at,
        });
    }

    /// The loop of `WriteString` (cited at `visible_str`): the string ends at its first NUL. In
    /// the object mode a control character is the error `NonVisibleChar` without `REPLACE`, and
    /// with it the character is recorded (`replace`) and `#` written in its place; a byte >= 0x80
    /// is written as it is (the stream has no encoding). Outside it (outfmt 5 / 14 / 16) a
    /// control character is always replaced and posted, and a byte >= 0x80 starts a UTF-8
    /// sequence (`encoded_tail`), whose failure ends the string.
    #[inline(always)]
    fn write_string<const REPLACE: bool>(&mut self, value: &[u8]) {
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        let mut bytes = value.iter();
        while let Some(&byte) = bytes.next() {
            match byte {
                b'&' => self.out.extend_from_slice(b"&amp;"),
                b'<' => self.out.extend_from_slice(b"&lt;"),
                b'>' => self.out.extend_from_slice(b"&gt;"),
                b'\'' => self.out.extend_from_slice(b"&apos;"),
                b'"' => self.out.extend_from_slice(b"&quot;"),
                0 => return,
                0x09 | 0x0a | 0x0d => {
                    self.out.extend_from_slice(b"&#x");
                    self.out.push(b"0123456789abcdef"[usize::from(byte & 0xf)]);
                    self.out.push(b';');
                }
                0x01..=0x1f if REPLACE || !OBJECT => {
                    self.replace(byte);
                    self.out.push(b'#');
                }
                0x01..=0x1f => self.fail(SerialError::NonVisibleChar(byte)),
                0x80..=0xff if OBJECT => self.out.push(byte),
                0x80..=0xff => {
                    // The rest of the string goes through the cold path, which leaves the loop of
                    // an ASCII string as it was.
                    let rest = bytes.as_slice();
                    self.encoded_tail(&value[value.len() - rest.len() - 1..]);
                    return;
                }
                _ => self.out.push(byte),
            }
        }
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
        if OBJECT && self.self_closed {
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
    #[inline(always)]
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
    /// that closes the innermost open tag (NCBI backs over it, `OpenTagEndBack`); `WriteNull`:
    /// `/>` there (`None`).
    fn value_attribute(&mut self, value: Option<&str>) {
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        debug_assert_eq!(self.out.last(), Some(&b'>'));
        self.out.pop();
        if let Some(value) = value {
            self.out.extend_from_slice(b" value=\"");
            self.out.extend_from_slice(value.as_bytes());
            self.out.push(b'"');
        }
        self.out.extend_from_slice(b"/>");
        self.self_closed = true;
    }

    /// `begin_sequence` of the object mode: `<Type>`, unless the type is anonymous (its name is
    /// the open member's tag, which has an underscore, or empty), which writes nothing.
    fn begin_object_sequence(&mut self, type_name: &str) {
        if anonymous_type(type_name) {
            debug_assert!(!self.pending_elem);
            debug_assert!(
                type_name.is_empty()
                    || &self.names[self.marks.last().copied().unwrap_or(0)..]
                        == type_name.as_bytes()
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

impl<const OBJECT: bool> Encoder for XmlEncoder<'_, OBJECT> {
    const STREAM: Option<super::encoder::Stream> = Some(super::encoder::Stream::Xml);

    /// `WriteFileHeader`: the declaration and the DOCTYPE of a document.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:232-276
    /// (`SetEncoding` with `eEncoding_Ascii`, blastxml_format.cpp:786: `US-ASCII`; the DTD prefix
    /// `SetDefaultDTDFilePrefix`, :789.)
    fn begin_document(&mut self, type_name: &str) {
        if self.piece == Piece::Xml2Document {
            // `WriteFileHeader` of the standard-XML mode: no DOCTYPE (objostrxml.cpp:232-276,
            // `m_UseSchemaRef` / `m_EnforcedStdXml`).
            self.out
                .extend_from_slice(b"<?xml version=\"1.0\" encoding=\"US-ASCII\"?>");
            return;
        }
        if OBJECT {
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
        if OBJECT {
            self.begin_object_sequence(type_name);
            return;
        }
        self.push_name(&[type_name.as_bytes()]);
        let start = *self.marks.last().expect("just pushed");
        self.open_tag(start);
    }

    /// `EndClass`: `</Type>`.
    fn end_sequence(&mut self) {
        if OBJECT && self.opened.pop() == Some(0) {
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
        if OBJECT {
            self.pending_elem = false;
            self.push_name(&[type_name.as_bytes()]);
            let start = *self.marks.last().expect("just pushed");
            self.open_tag(start);
        }
    }

    /// `EndNamedType`: `</Type>`.
    fn end_named(&mut self) {
        if OBJECT {
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
        if OBJECT {
            self.close_element();
            self.pending_elem = true;
        }
    }

    /// `EndContainer`.
    fn end_list(&mut self) {
        if OBJECT {
            self.close_element();
            self.pending_elem = false;
        }
    }

    /// `BeginChoice` and `BeginChoiceVariant`: in the standard-XML mode the CHOICE's type tag
    /// and the variant's name as a member tag (`<Target><db>`); outside it NCBI's stream
    /// writes an element LOSAT's outfmt 5 never needs.
    fn variant(&mut self, type_name: &str, name: &str, tag: u8) {
        if OBJECT {
            if anonymous_type(type_name) {
                // An anonymous CHOICE: the member's tag is already open.
                debug_assert!(!self.pending_elem);
                debug_assert!(
                    type_name.is_empty()
                        || &self.names[self.marks.last().copied().unwrap_or(0)..]
                            == type_name.as_bytes()
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
        self.choices |= 1u64 << self.marks.len();
        self.begin_sequence(type_name);
        self.member(name, tag);
    }

    /// `EndChoiceVariant` and `EndChoice`.
    fn end_variant(&mut self) {
        if OBJECT {
            self.close_tag();
            if self.opened.pop() == Some(2) {
                self.close_tag();
            }
        } else if self.std_xml {
            self.end_member();
            self.end_sequence();
            self.choices &= !(1u64 << self.marks.len());
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
        if OBJECT {
            // The object mode writes a CHOICE and a list as the trait's default calls them.
            self.variant(type_name, name, tag);
            self.begin_list(false);
            for item in items {
                self.element();
                self.visible_str(item.as_ref());
            }
            self.end_list();
            self.end_variant();
            return;
        }
        if !self.std_xml {
            self.fail(SerialError::XmlUnsupported("CHOICE"));
            return;
        }
        let depth = self.marks.len();
        self.choices |= 1u64 << depth;
        self.begin_sequence(type_name);
        for item in items {
            self.member(name, tag);
            self.visible_str(item.as_ref());
            self.end_member();
        }
        self.end_sequence();
        self.choices &= !(1u64 << depth);
    }

    /// `WriteInt4`: `PutInt4`, the decimal digits.
    fn int(&mut self, value: i32) {
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        let mut buffer = itoa_i64(i64::from(value));
        self.out.extend_from_slice(buffer.text());
    }

    /// `WriteInt8`: `PutInt8`, the decimal digits.
    fn big_int(&mut self, value: i64) {
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        let mut buffer = itoa_i64(value);
        self.out.extend_from_slice(buffer.text());
    }

    /// `WriteEnum` of the object mode: ` value="name"/>` in place of the `>` of the open tag.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:473-484
    fn enumerated(&mut self, name: &str, _value: i32) {
        if OBJECT {
            self.value_attribute(Some(name));
        } else {
            self.fail(SerialError::XmlUnsupported("ENUMERATED value"));
        }
    }

    /// `WriteEnum` of an INTEGER with named values in the object mode: ` value="name">`, the
    /// number, and the closing tag that `end_named` writes.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:442-456
    fn named_int(&mut self, name: &str, value: i32) {
        if OBJECT {
            if self.pending_elem {
                self.open_element();
            }
            debug_assert_eq!(self.out.last(), Some(&b'>'));
            self.out.pop();
            self.out.extend_from_slice(b" value=\"");
            self.out.extend_from_slice(name.as_bytes());
            self.out.extend_from_slice(b"\">");
            let mut buffer = itoa_i64(i64::from(value));
            self.out.extend_from_slice(buffer.text());
        } else {
            self.fail(SerialError::XmlUnsupported("INTEGER with named values"));
        }
    }

    /// `WriteBool` of the object mode: ` value="true"/>` or ` value="false"/>`.
    fn boolean(&mut self, value: bool) {
        if OBJECT {
            self.value_attribute(Some(if value { "true" } else { "false" }));
        } else {
            self.fail(SerialError::XmlUnsupported("BOOLEAN"));
        }
    }

    /// `WriteString`: the escaped text.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:502-549
    /// A byte of 0x01-0x1F other than TAB, LF and CR is replaced by `#` with an error
    /// message (`x_VerifyChar`, `BAD_CHAR`; `replace`): in the object mode an error here, `#`
    /// and a record in `visible_str_replacing`; in outfmt 5 / 14 / 16 always `#` and a post
    /// (`XmlPost`). A byte >= 0x80 is decoded as UTF-8 by the US-ASCII stream of outfmt 5 / 14 /
    /// 16, which throws where it cannot write the character (`encoded_tail`, `XmlThrow`); the
    /// archive's stream has no encoding, so `WriteEncodedChar` takes its first branch
    /// (objostrxml.cpp:559-565) and the byte passes `x_VerifyChar` and is written as it is (the
    /// object mode). A NUL byte ends the string.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:927-935
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string::<false>(value);
    }

    /// `visible_str` that writes NCBI's `#` for a control character and records it
    /// (`replace`).
    fn visible_str_replacing(&mut self, value: &[u8]) {
        self.write_string::<true>(value);
    }

    /// `WriteNull` of the object mode (a NULL of a PSSM query's parts, `asn_types::write`):
    /// `/>` in place of the `>` that closes the innermost open tag (`<Seq-loc_null/>`). The
    /// reports' XML (outfmt 5, 14, 16) has none.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:720-724
    fn null(&mut self) {
        if OBJECT {
            self.value_attribute(None);
        } else {
            self.fail(SerialError::Unwritten("a NULL value"));
        }
    }

    fn null_element(&mut self) {
        self.fail(SerialError::NullElement);
    }

    /// `WriteString` of a UTF8String (a PSSM query's `User-field` `str`, `asn_types::write`) in
    /// the object mode: the archive's stream has no encoding, so the string's encoding is the
    /// stream's (UTF-8) and NCBI checks and writes each byte as a VisibleString's
    /// (`visible_str`; NCBI 2.17.0 probe X of line RG-F3 PR 2, `é` and the byte 0x01: `é#` and
    /// the post 803.7). The reports' XML (outfmt 5, 14, 16) has none.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:556-565
    fn utf8_str(&mut self, value: &[u8]) {
        if OBJECT {
            self.write_string::<false>(value);
        } else {
            self.fail(SerialError::Unwritten("a UTF8String"));
        }
    }

    /// `utf8_str` that writes NCBI's `#` for a control character and records it, as
    /// `visible_str_replacing`.
    fn utf8_str_replacing(&mut self, value: &[u8]) {
        if OBJECT {
            self.write_string::<true>(value);
        } else {
            self.fail(SerialError::Unwritten("a UTF8String"));
        }
    }

    /// The XML documents have no delayed member (`Auth-list.names`, a PSSM file's member read
    /// into a delay buffer, `asn_types::write`), which NCBI parses before it writes it in XML:
    /// LOSAT's psiblast refuses `-outfmt 11` and its tblastn `-in_pssm` refuses the archive of a
    /// search, so no archive of a PSSM query is written as XML.
    fn delayed(&mut self, _text: &[u8]) {
        self.fail(SerialError::Unwritten("a member read into a delay buffer"));
    }

    /// `WriteBytes` of the object mode: the bytes as upper-case hex digits, unwrapped
    /// (`WriteBytes`, objostrxml.cpp:1702-1709); the other pieces write no OCTET STRING.
    fn octets(&mut self, value: &[u8]) {
        if !OBJECT {
            self.fail(SerialError::OctetString);
            return;
        }
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        self.out.reserve(2 * value.len());
        for &byte in value {
            self.out.push(b"0123456789ABCDEF"[usize::from(byte >> 4)]);
            self.out.push(b"0123456789ABCDEF"[usize::from(byte & 0xf)]);
        }
    }

    /// `WriteDouble`: `WriteDouble2` of `data` with `DBL_DIG`, `NStr::DoubleToStringPosix` of
    /// `data` with 15 (`real::write_xml_double`).
    fn real(&mut self, value: f64) {
        if OBJECT && self.pending_elem {
            self.open_element();
        }
        write_xml_double(value, self.out);
    }

    fn buffered(&self) -> usize {
        if OBJECT {
            self.out.len() - self.start
        } else {
            0
        }
    }

    /// The object mode keeps the last byte until the document is complete: it may be the `>` of
    /// an open tag that `enumerated` or `boolean` replaces.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        if !OBJECT {
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
            // `replace` counts again from the bytes held (`encoder::Traced` adds those handed on).
            (self.lines_to, self.lines) = (self.start, 0);
        }
        Ok(())
    }

    fn held(&self) -> &[u8] {
        &self.out[self.start..]
    }

    fn hand_on(&mut self, sink: &mut dyn std::io::Write, count: usize) -> std::io::Result<()> {
        let end = self.start + count;
        sink.write_all(&self.out[self.start..end])?;
        self.out.drain(self.start..end);
        // `replace` counts again from the bytes held (`encoder::Traced` adds those handed on).
        (self.lines_to, self.lines) = (self.start, 0);
        Ok(())
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }

    fn take_replaced(&mut self) -> Vec<Replaced> {
        std::mem::take(&mut self.replaced)
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
