//! The object-stream calls that NCBI's serial writers make, shared by the encoders of the
//! serialized formats (ASN.1 text `asn_text`, JSON `json`, ASN.1 binary `asn_ber`), so that
//! every type writer (`seqalign`, `report::seqid`) is written once.
//!
//! NCBI's generated type info writes an object as a sequence of `CObjectOStream` calls:
//! `WriteFileHeader`, `BeginClass` / `BeginClassMember` / `EndClassMember` / `EndClass` for a
//! SEQUENCE, `BeginContainer` / `BeginContainerElement` / `EndContainer` for a SEQUENCE OF or
//! SET OF, `BeginChoice` / `BeginChoiceVariant` / `EndChoiceVariant` / `EndChoice` for a
//! CHOICE, and `WriteInt4`, `WriteEnum`, `WriteBool`, `WriteString`, `WriteDouble` for the
//! primitive values, then `EndOfWrite`.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:610-625
//!
//! The calls below are those pairs, one method each; a type writer calls them in the member
//! order of the ASN.1 module and skips an absent OPTIONAL member. `tag` is the member's or
//! variant's 0-based position in its ASN.1 SEQUENCE or CHOICE (absent OPTIONAL members still
//! count), which only the binary encoder writes (its context tag `[tag]`).

pub use super::out_buffer::{BufferPiece, Stream};
use super::out_buffer::{OutBuffer, Piece};

/// An error that NCBI's object stream raises while it writes a value, which LOSAT reports
/// instead of writing the document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerialError {
    /// `CObjectOStreamJson::WriteDouble2` on a NaN or an infinity
    /// (`real::JsonDoubleError`).
    JsonDouble(super::real::JsonDoubleError),
    /// A VisibleString byte outside 0x20..0x7E in the ASN.1 text or binary stream (in the XML
    /// object stream of the archive a control character other than TAB, LF and CR), which NCBI
    /// replaces by `#` with an error message (`eFNP_ReplaceAndWarn`, `ReplacedChar`):
    /// `visible_str` reports it here; a writer that posts NCBI's messages writes the document
    /// again through `Traced` (`visible_str_replacing`: the archive, the search strategy, the
    /// Seq-annots of outfmt 8 and 9, the mask writers, psiblast's PSSM file); for the others it
    /// is the writer's error. (The XML streams of outfmt 5 / 14 / 16 replace and post it
    /// themselves, `xml::XmlPost`.)
    NonVisibleChar(u8),
    /// The exception of NCBI's US-ASCII XML stream of outfmt 5 / 14 / 16 at a byte >= 0x80 that it
    /// cannot write (`CUtf8::DecodeFirst`, `DecodeNext` or `SymbolToChar`): the encoder keeps the
    /// rest of what NCBI reports (`xml::XmlThrow`: where the document is cut, the frames of the
    /// trace), which the document writers report as NCBI does at their level
    /// (`report::blast_xml::writer`, `report::blast_xml2::writer`). The archive's stream has no
    /// encoding and writes the byte.
    Utf8(crate::ncbi_str::utf8::Utf8Error),
    /// A call that the XML stream does not write outside the archive's object mode (a BOOLEAN,
    /// an ENUMERATED value or a CHOICE of outfmt 5).
    XmlUnsupported(&'static str),
    /// An OCTET STRING in a stream that does not write one (the archive of outfmt 11 is the
    /// only document with one; JSON and the XML of outfmt 5 / 14 / 16 do not write it).
    OctetString,
    /// A value of a PSSM file's members that only its ASN.1 text writes (`Encoder::null`,
    /// `utf8_str`, `delayed`; `report::serial::asn_types::write`) in another stream: no BLAST
    /// document that LOSAT writes in ASN.1 binary, XML or JSON holds one (psiblast and tblastn
    /// `-in_pssm` refuse `-outfmt 11`).
    Unwritten(&'static str),
    /// An element of a list of `CRef<>`s that holds no object (`Encoder::null_element`), where
    /// NCBI's stream throws: a writer that reports NCBI's error writes the document again through
    /// `Traced` (`NullElement`, the ASN.1 text stream of psiblast's PSSM file and of the search
    /// strategy); for the others it is the writer's error.
    NullElement,
    /// A member that NCBI's reader kept in a delay buffer of another format, whose read before
    /// the member is written fails (`Encoder::delayed_read_failed`), where NCBI's exception
    /// leaves its writer: a writer that reports NCBI's error writes the document again through
    /// `Traced` (`Thrown::DelayedRead`); for the others it is the writer's error.
    DelayedRead,
}

impl std::fmt::Display for SerialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SerialError::JsonDouble(error) => f.write_str(error.message()),
            SerialError::NonVisibleChar(byte) => write!(
                f,
                "a character 0x{byte:02x} outside the visible ASCII range in an ASN.1 VisibleString is not supported by LOSAT"
            ),
            SerialError::Utf8(error) => f.write_str(&String::from_utf8_lossy(error.exception().get_msg())),
            SerialError::XmlUnsupported(what) => {
                write!(f, "an XML {what} is not supported by LOSAT")
            }
            SerialError::OctetString => {
                f.write_str("an OCTET STRING in this serialized format is not supported by LOSAT")
            }
            SerialError::Unwritten(what) => {
                write!(f, "{what} in this serialized format is not supported by LOSAT")
            }
            SerialError::NullElement => f.write_str(
                "an element of a list stated as NULL (where NCBI BLAST+'s writer stops with a C++ exception) in this serialized format is not supported by LOSAT",
            ),
            SerialError::DelayedRead => f.write_str(
                "a member read into a delay buffer of another format whose value NCBI reads again with an error (where NCBI BLAST+'s writer stops with a C++ exception) in this serialized format is not supported by LOSAT",
            ),
        }
    }
}

impl std::error::Error for SerialError {}

/// Where NCBI's object stream stands when it replaces a byte (`GetPosition`): the line of the
/// ASN.1 text or XML stream (`m_Output.GetLine()`, 1 at the stream's start and one more for each
/// end of line written), the offset of the binary stream (`GetStreamPos()`, 0 at its start).
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:100-103
/// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:510-513
/// (`CObjectOStreamXml::GetPosition`, objostrxml.cpp:192-195, is the text stream's.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamPosition {
    Line(u64),
    Byte(u64),
}

/// `GetPosition()`'s text: `line N`, `byte N` (`CObjectIStream::ThrowError1`'s prefix, the
/// position of `ReplacedChar::message`).
impl std::fmt::Display for StreamPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StreamPosition::Line(line) => write!(f, "line {line}"),
            StreamPosition::Byte(byte) => write!(f, "byte {byte}"),
        }
    }
}

/// A byte of a VisibleString that NCBI's object stream writes as `#` and reports with
/// `ERR_POST_X` of subcode 7 with a message, from `Serial_IStream` (error 803.7: the function is in
/// objistr.cpp).
///
/// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897
/// The text and binary streams pass the whole string as `str`, the XML stream `kEmptyStr`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReplacedChar {
    /// `GetStackTrace()` (`Traced::stack_trace`).
    pub trace: String,
    pub byte: u8,
    pub position: StreamPosition,
    /// The string that holds the byte (empty in XML), one copy for the bytes of a string.
    pub string: std::sync::Arc<[u8]>,
}

/// A byte that an encoder replaced by `#`, as the encoder knows it: `Traced` adds the stack
/// trace and the string (`ReplacedChar`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Replaced {
    pub byte: u8,
    /// The position in the bytes that the encoder holds.
    pub position: StreamPosition,
    /// Whether NCBI's message quotes the string (the text and binary streams, not XML).
    pub with_string: bool,
    /// The length of `Encoder::held` when NCBI posts the byte (before its wrap and its `#`).
    pub at: usize,
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:578-605
// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1225-1250
// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1255-1300
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:522-555
/// NCBI's throw at an element of a list of `CRef<>`s that holds no object (`Traced::null_element`)
/// in the ASN.1 text, ASN.1 binary or XML stream: `WriteContainer` (XML:
/// `WriteContainerContents`) throws `fUnassigned` there, before the element's separator; the throw
/// flushes what the stream's buffer holds (so the output ends with the bytes before the element),
/// posts the stream's first failure (`SetFailFlags`: error 802.5, `Serial_OStream`) and throws
/// `CUnassignedMember::eWrite`, to which each frame it leaves adds its name
/// (`CSerialException::AddFrameInfo`). The message names the container's type, which is empty:
/// the SET OF or SEQUENCE OF that NCBI's generated class of a named list type holds (`Seq-descr`)
/// has no name (psiblast and tblastn 2.17.0: `... while writing container `, probe D3c-3b
/// `probe4`; psiblast `-outfmt 11` in the three streams, probe D3c-3c `probeA`, `probeA2`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NullElement {
    /// The stream that threw (`site`).
    pub stream: Stream,
    /// Whether the stream has a frame for the element (always but in XML, whose elements of a
    /// named type have none: `Traced::name_element`), which decides XML's branch (`site`).
    pub element_frame: bool,
    /// `GetStackTrace()` at the element (`Traced::stack_trace`, ending with `.E`).
    pub trace: String,
    /// `GetPosition()` after the flush.
    pub position: StreamPosition,
    /// `GetFrameName()` of each frame that the exception leaves, from the document's root to the
    /// element (`object_frame::frame_name`).
    pub frames: Vec<String>,
}

/// The `ThrowError` at an element that holds no object of `file` at `line`, in `function`.
const fn write_container(
    file: &'static str,
    line: u32,
    function: &'static str,
) -> crate::ncbi_exception::ThrowSite {
    crate::ncbi_exception::ThrowSite {
        file,
        line,
        kind: "CUnassignedMember::eWrite",
        function: Some(function),
    }
}

impl NullElement {
    /// The error code of the post (`Serial_OStream`, include/serial/error_codes.hpp:46).
    pub const ERR_CODE: (u32, u32) = (802, 5);

    /// The message of the throw (`ThrowError1` of `fUnassigned`).
    const MESSAGE: &'static str =
        "cannot write unassigned member NULL element while writing container ";

    /// The post of the stream's failure on the diagnostics stream: `Error: (802.5) `, the
    /// application's prefix in brackets when it has one, `CObjectOStream: error at <position>:
    /// <stack trace>: <message>` and an end of line (`ncbi_exception::error_post_line`).
    pub fn post(&self, prefix: Option<&str>) -> Vec<u8> {
        let (code, subcode) = Self::ERR_CODE;
        crate::ncbi_exception::error_post_line(
            Some(&format!("{code}.{subcode}")),
            prefix,
            format!(
                "CObjectOStream: error at {}: {}: {}",
                self.position,
                self.trace,
                Self::MESSAGE
            )
            .as_bytes(),
        )
    }

    /// Where the stream throws: its `WriteContainer` (the ASN.1 text and binary streams) or, in
    /// XML, the branch of `WriteContainerContents` for the elements of a named type, or the other
    /// one, which has a frame for the element (`element_frame`).
    fn site(&self) -> crate::ncbi_exception::ThrowSite {
        match self.stream {
            Stream::AsnText => write_container(
                "src/serial/objostrasn.cpp",
                597,
                "ncbi::CObjectOStreamAsn::WriteContainer",
            ),
            Stream::AsnBinary => write_container(
                "src/serial/objostrasnb.cpp",
                1244,
                "ncbi::CObjectOStreamAsnBinary::WriteContainer",
            ),
            Stream::Xml => write_container(
                "src/serial/objostrxml.cpp",
                if self.element_frame { 1296 } else { 1272 },
                "ncbi::CObjectOStreamXml::WriteContainerContents",
            ),
            // `Traced::null_element` makes no throw of the JSON stream, whose documents have no
            // NULL element (`json::JsonEncoder::null_element`).
            Stream::Json => unreachable!("the JSON stream writes no NULL element"),
        }
    }

    /// The post (`post`) at the start of `message`, the text of `cli::NativeError::thrown`,
    /// and the rest (`CATCH_ALL`'s report); no post in another message.
    pub fn split_post(message: &[u8]) -> (Option<&[u8]>, &[u8]) {
        let (code, subcode) = Self::ERR_CODE;
        let head = format!("Error: ({code}.{subcode}) ");
        match message.iter().position(|&byte| byte == b'\n') {
            Some(end) if message.starts_with(head.as_bytes()) => {
                (Some(&message[..end]), &message[end + 1..])
            }
            _ => (None, message),
        }
    }

    /// The exception: the position and the message, with the frames' names (the innermost
    /// added first, as the frames are left).
    pub fn exception(&self) -> crate::ncbi_exception::NcbiException {
        let mut exception = crate::ncbi_exception::NcbiException::new(
            self.site(),
            format!("{}: {}", self.position, Self::MESSAGE),
        );
        let names: Vec<&str> = self.frames.iter().rev().map(String::as_str).collect();
        exception.add_frame_infos(&names);
        exception
    }
}

/// The exception of NCBI's read of a delay buffer of another format that failed
/// (`CDelayBuffer::Update`, `asn_types::write`), with the frames of the read, which leaves NCBI's
/// writer through each frame the writer has open (`Traced::delayed_read_failed`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DelayedReadError {
    /// A `CSerialException`: each frame adds its name to the exception's frame stack
    /// (`CSerialException::AddFrameInfo`).
    Serial(crate::ncbi_exception::NcbiException),
    /// Another `CException` (the `CUtilException` of a number past its storage): each frame
    /// rethrows it with the frame's text (`NCBI_RETHROW_SAME`).
    Other(crate::ncbi_exception::NcbiException),
}

/// What NCBI's object stream threw while it wrote a document (`Traced::take_thrown`): the bytes
/// held are those NCBI flushed, and its exception ends the run (`cli::NativeError::thrown`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Thrown {
    /// An element that holds no object (`NullElement`).
    NullElement(NullElement),
    /// A delay buffer whose read failed (`DelayedReadError`): the exception with the writer's
    /// frames.
    DelayedRead(crate::ncbi_exception::NcbiException),
}

impl Thrown {
    /// The post of the stream's failure before the exception's report: `NullElement::post`; none
    /// for a failed read, whose exception passes through the writer without setting its fail
    /// flags.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
    pub fn post(&self, prefix: Option<&str>) -> Vec<u8> {
        match self {
            Thrown::NullElement(null) => null.post(prefix),
            Thrown::DelayedRead(_) => Vec::new(),
        }
    }

    /// The exception that NCBI's `CATCH_ALL` reports.
    pub fn exception(&self) -> crate::ncbi_exception::NcbiException {
        match self {
            Thrown::NullElement(null) => null.exception(),
            Thrown::DelayedRead(exception) => exception.clone(),
        }
    }
}

/// The place of the `END_OBJECT_FRAME` of a frame of NCBI's ASN.1 text or binary output stream,
/// where `NCBI_RETHROW_SAME` adds the frame's line to an exception other than a
/// `CSerialException` that leaves it (`Traced::delayed_read_failed`): the stream's `WriteClass`,
/// `WriteClassMember` and `WriteContainer`, `WriteChoiceDefault`'s two frames, `Write`'s root
/// frame (NCBI BLAST+ 2.17.0, line RG-F3 item O-2 probe 1: psiblast's and tblastn's strategy,
/// psiblast's `-out_pssm` and archives).
// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:579-614
// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:687-713
// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1206-1264
// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1381-1455
// NCBI reference (598d8ae6): c++/src/serial/choice.cpp:428-463
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:610-621
fn write_frame_site(
    stream: Stream,
    kind: super::object_frame::FrameKind,
) -> Option<(&'static str, u32, &'static str)> {
    use super::object_frame::FrameKind;
    const CHOICE: &str = "src/serial/choice.cpp";
    const WRITE_CHOICE: &str = "ncbi::CChoiceTypeInfoFunctions::WriteChoiceDefault";
    // The file, then the line and function of `WriteClass`, `WriteClassMember`, and
    // `WriteContainer`'s element and container frames.
    let (file, [class, member, element, array], [write_class, write_member, write_container]) =
        match stream {
            Stream::AsnText => (
                "src/serial/objostrasn.cpp",
                [698, 713, 610, 614],
                [
                    "ncbi::CObjectOStreamAsn::WriteClass",
                    "ncbi::CObjectOStreamAsn::WriteClassMember",
                    "ncbi::CObjectOStreamAsn::WriteContainer",
                ],
            ),
            Stream::AsnBinary => (
                "src/serial/objostrasnb.cpp",
                [1417, 1455, 1255, 1264],
                [
                    "ncbi::CObjectOStreamAsnBinary::WriteClass",
                    "ncbi::CObjectOStreamAsnBinary::WriteClassMember",
                    "ncbi::CObjectOStreamAsnBinary::WriteContainer",
                ],
            ),
            // The XML stream's frames differ (a class or a container frame only for a named
            // type, its own choice functions); a buffer of another format whose read fails with
            // an exception other than a `CSerialException` does not reach it (an XML buffer is
            // copied, the binary reader throws `CSerialException`s, the text reader's skip reads
            // the numbers).
            Stream::Xml | Stream::Json => return None,
        };
    Some(match kind {
        FrameKind::Named | FrameKind::NamedType => {
            ("src/serial/objostr.cpp", 621, "ncbi::CObjectOStream::Write")
        }
        FrameKind::Class => (file, class, write_class),
        FrameKind::ClassMember => (file, member, write_member),
        FrameKind::Choice => (CHOICE, 463, WRITE_CHOICE),
        FrameKind::ChoiceVariant => (CHOICE, 459, WRITE_CHOICE),
        FrameKind::ArrayElement => (file, element, write_container),
        FrameKind::Array => (file, array, write_container),
    })
}

/// A replaced byte (`ReplacedChar`) and the bytes of its document that NCBI's stream had
/// written to its output when NCBI posted it (`out_buffer::OutBuffer::written`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Posted {
    pub written: u64,
    pub replaced: ReplacedChar,
}

impl ReplacedChar {
    /// The error code of the post (`Serial_IStream`, include/serial/error_codes.hpp:47).
    pub const ERR_CODE: (u32, u32) = (803, 7);

    /// The message of the post: the stack trace, `Bad char [0x<hex>] in string at <position>`
    /// (the hex digits in upper case without padding, `NumericToString` of `c` in base 16) and the
    /// string when there is one.
    pub fn message(&self) -> Vec<u8> {
        use std::io::Write;
        // One buffer of the message's length (the text between the trace and the string, with a
        // position of 20 digits, is under 64 bytes).
        let mut message = Vec::with_capacity(self.trace.len() + 64 + self.string.len());
        write!(
            message,
            "{}\nBad char [0x{:X}] in string at {}",
            self.trace, self.byte, self.position
        )
        .expect("a write to a Vec");
        if !self.string.is_empty() {
            message.push(b'\n');
            message.extend_from_slice(&self.string);
        }
        message
    }

    /// The text of the post on the diagnostics stream: `Error: (803.7) `, the application's
    /// prefix (`SetDiagPostPrefix`) in brackets when it has one, the message and an end of line
    /// (`ncbi_exception::error_post_line`, `SDiagMessage::x_OldWrite`).
    pub fn post(&self, prefix: Option<&str>) -> Vec<u8> {
        let (code, subcode) = Self::ERR_CODE;
        crate::ncbi_exception::error_post_line(
            Some(&format!("{code}.{subcode}")),
            prefix,
            &self.message(),
        )
    }
}

/// The calls of NCBI's `CObjectOStream` that the serialized BLAST reports make.
pub trait Encoder {
    /// NCBI's object stream that the encoder writes, whose output buffer `Traced` follows
    /// (`out_buffer`); `None` for an encoder that `Traced` does not wrap.
    const STREAM: Option<Stream> = None;
    /// `IsWritingDefaultValuesEnforced`: whether a member with a DEFAULT that was not set is
    /// written anyway, with its default value (`report::seqid` quotes the rule). Only NCBI's
    /// JSON stream enforces it (`m_EnforceWritingDefaults`, `json`).
    fn writes_defaults(&self) -> bool {
        false
    }
    /// `GetDataFormat()`: the encoding of the delay buffers that the stream copies as they were
    /// read (`CObjectOStream::WriteClassMember` of a `CDelayBuffer` that `HaveFormat` it;
    /// `asn_types::write`), the one of its `STREAM`; `None` for a stream of no ASN.1 encoding.
    fn data_format(&self) -> Option<crate::report::serial::asn_input::Encoding> {
        use crate::report::serial::asn_input::Encoding;
        Self::STREAM.and_then(|stream| match stream {
            Stream::AsnText => Some(Encoding::Text),
            Stream::AsnBinary => Some(Encoding::Binary),
            Stream::Xml => Some(Encoding::Xml),
            // The JSON stream copies no delay buffer (`json::JsonEncoder::delayed`).
            Stream::Json => None,
        })
    }
    /// `WriteFileHeader` of a type: the start of a document whose top-level object is of the
    /// ASN.1 type `type_name` (`Seq-annot`).
    fn begin_document(&mut self, type_name: &str);
    /// `EndOfWrite` and the separator after the top-level object.
    fn end_document(&mut self);
    /// `BeginClass`: a SEQUENCE value of the ASN.1 type `type_name` starts (only the XML
    /// stream reads the name: it writes `<Type>` and `<Type_member>` tags).
    fn begin_sequence(&mut self, type_name: &str);
    /// `EndClass`.
    fn end_sequence(&mut self);
    /// `BeginNamedType` (`CObjectOStream::WriteNamedType`): a value of the named ASN.1 type
    /// `type_name` that is neither a SEQUENCE nor a CHOICE starts (a `Seq-descr`, a `Na-strand`
    /// or a `NCBI2na`: a SEQUENCE OF / SET OF alias, an ENUMERATED type, an OCTET STRING or
    /// StringStore alias, any other alias), around the value's own calls. Only the XML object
    /// stream of the archive (`XmlEncoder::archive`) writes anything for it (`<Type>`); the
    /// ASN.1 text, ASN.1 binary and JSON streams write nothing for a named type, so the
    /// default does nothing. The name is one of the writers' constants, which `Traced` keeps for
    /// its frames without a copy.
    fn begin_named(&mut self, type_name: &'static str) {
        let _ = type_name;
    }
    /// `EndNamedType`.
    fn end_named(&mut self) {}
    /// `BeginClassMember` of an `id`: member `name` (position `tag`) of the current SEQUENCE.
    fn member(&mut self, name: &str, tag: u8);
    /// `EndClassMember`.
    fn end_member(&mut self);
    /// `BeginContainer`: a SEQUENCE OF (`set` false) or SET OF (`set` true) value starts.
    fn begin_list(&mut self, set: bool);
    /// `BeginContainerElement`: the next element of the current container.
    fn element(&mut self);
    /// `EndContainer`.
    fn end_list(&mut self);
    /// `BeginChoice` and `BeginChoiceVariant` of an `id`: a CHOICE value of the ASN.1 type
    /// `type_name` with variant `name` (position `tag`). Only the XML stream reads the type
    /// name (in its standard-XML mode it writes `<Type>` and the variant as a member).
    fn variant(&mut self, type_name: &str, name: &str, tag: u8);
    /// `EndChoiceVariant` and `EndChoice`.
    fn end_variant(&mut self);
    /// `WriteInt4`.
    fn int(&mut self, value: i32);
    /// `WriteInt8`: an INTEGER written as a BigInt (XML and the text formats print it
    /// decimal).
    fn big_int(&mut self, value: i64);
    /// `WriteEnum` of an ENUMERATED value: its name and number.
    fn enumerated(&mut self, name: &str, value: i32);
    /// `WriteEnum` of an INTEGER with named values (`CEnumeratedTypeValues::IsInteger`, as
    /// `Blast-filter-program`): the name `name` of the value `value` in ASN.1 text and JSON,
    /// an INTEGER in BER, `<Type value="name">value</Type>` in XML (inside `begin_named`).
    fn named_int(&mut self, name: &str, value: i32);
    /// `WriteBool`.
    fn boolean(&mut self, value: bool);
    /// `WriteString` of a VisibleString.
    fn visible_str(&mut self, value: &[u8]);
    /// `visible_str` that writes NCBI's `#` for each byte that NCBI's stream replaces and
    /// records it (`Replaced`, `take_replaced`) instead of the error `NonVisibleChar`: what
    /// `Traced` calls. The default is `visible_str` (an encoder that replaces nothing).
    fn visible_str_replacing(&mut self, value: &[u8]) {
        self.visible_str(value);
    }
    /// `WriteDouble` of a REAL.
    fn real(&mut self, value: f64);
    /// `WriteBytes` of an OCTET STRING (`Seq-data` `ncbi2na` and `ncbi4na`): in ASN.1 text
    /// `'` then two upper-case hex digits per byte, an end of line without indentation before
    /// a digit pair whenever the line holds 78 bytes or more, then `'H`; in ASN.1 binary the
    /// tag 0x04, the length and the bytes; in the XML object stream of the archive two upper-case
    /// hex digits per byte, unwrapped.
    fn octets(&mut self, value: &[u8]);
    /// `WriteStringStore` of a StringStore (`Seq-data` `ncbieaa`): a VisibleString in ASN.1
    /// text; the application tag 0x41, the length and the bytes in ASN.1 binary.
    fn string_store(&mut self, value: &[u8]) {
        self.visible_str(value);
    }
    /// `string_store` with `visible_str_replacing`'s replacement where the stream checks a
    /// StringStore (ASN.1 text); the default is `string_store`.
    fn string_store_replacing(&mut self, value: &[u8]) {
        self.string_store(value);
    }
    /// `WriteNull` of a NULL value (`Seq-loc` `null`, `Seq-loc` `empty`'s neighbours in a PSSM's
    /// annotations; `report::serial::asn_types::write`): `NULL` in ASN.1 text, the tag 0x05 and the
    /// length 0 in ASN.1 binary.
    fn null(&mut self);
    /// `WriteString` of a string with `eStringTypeUTF8`, for a UTF8String: NCBI's stream writes it
    /// with `eFNP_Allow` (no byte is replaced or reported; in ASN.1 text a zero byte is left out, a
    /// `"` doubled and the line wrapped as a VisibleString's).
    fn utf8_str(&mut self, value: &[u8]);
    /// `utf8_str` with `visible_str_replacing`'s replacement where the stream checks a
    /// UTF8String's bytes (NCBI's XML stream); the default is `utf8_str`.
    fn utf8_str_replacing(&mut self, value: &[u8]) {
        self.utf8_str(value);
    }
    /// `WriteClassMember` of an `id` and a `CDelayBuffer` of the stream's own format
    /// (`data_format`), for a member that NCBI's reader kept in a delay buffer (`Auth-list.names`,
    /// `asn_types::DELAYED`): the bytes read, copied as they are between the member's start and
    /// end (`Write` of the buffer's source; in ASN.1 text `COStreamBuffer::Write` from a
    /// `CByteSourceReader`, which leaves the line length as it was). A buffer of another format is
    /// read first and written by its values (`CDelayBuffer::Update`, `asn_types::write`).
    fn delayed(&mut self, bytes: &[u8]);
    /// A member that NCBI's reader kept in a delay buffer of another format and that NCBI's
    /// `CDelayBuffer::Update` reads again with a message that the encoder does not post
    /// (`delayed_posts`), or with an error that LOSAT's reader does not report as NCBI's
    /// (`asn_types::read::ReadError`), before it writes it (`asn_types::write`): LOSAT does not
    /// write it (`SerialError::Unwritten`).
    fn unread_delayed(&mut self);
    /// A member that NCBI's reader kept in a delay buffer of another format, whose read by
    /// `CDelayBuffer::Update` fails with `error` before the member starts (`WriteLongMember`
    /// reads the buffer before `WriteClassMember`; `asn_types::write`): NCBI's exception leaves
    /// the writer there, after the bytes written before the member. `Traced` makes NCBI's
    /// exception (`Thrown::DelayedRead`) and writes nothing after it; the encoders' error is
    /// `SerialError::DelayedRead`. The default writes no such member (`unread_delayed`).
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:1084-1103
    fn delayed_read_failed(&mut self, error: &DelayedReadError) {
        let _ = error;
        self.unread_delayed();
    }
    /// The posts of the bytes that `CDelayBuffer::Update`'s read of a delayed member of another
    /// format replaced (803.7, `ReplacedChar`: the trace of that read's own stack, the position in
    /// the buffer), which NCBI makes before it writes the member (`WriteLongMember` calls `Update`
    /// before `WriteClassMember`; `asn_types::write` then writes the value read). `Traced` posts
    /// them as its own (`Posted`, at what NCBI's buffer had written before the member); an encoder
    /// that keeps no stack reports them as a byte it replaced (`SerialError::NonVisibleChar`,
    /// after which its callers write the document again through `Traced`); the default writes
    /// no such member (`unread_delayed`).
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:1084-1103
    fn delayed_posts(&mut self, posts: &[ReplacedChar]) {
        let _ = posts;
        self.unread_delayed();
    }
    /// The bytes that `flush_to` would hand over (0 for an encoder that cannot flush).
    fn buffered(&self) -> usize {
        0
    }
    /// Hands the bytes that are final to `sink` and drops them from the encoder's buffer, so
    /// that a document of any size is written element by element: the ASN.1 text encoder keeps
    /// the line that it is writing and the end of line before it (the string wrap looks back
    /// no further), the binary encoder keeps nothing. With `all` (the document is complete)
    /// every byte goes. The default keeps everything.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        let _ = (sink, all);
        Ok(())
    }
    /// The bytes of the document that the encoder holds (that `flush_to` has not handed on), for
    /// `Traced`; the binary encoder's buffer also holds what preceded the document until its
    /// first hand-on. The default holds nothing.
    fn held(&self) -> &[u8] {
        &[]
    }
    /// Hands the first `count` bytes of `held` to `sink` and drops them (`Traced` hands on the
    /// bytes that NCBI's buffer has written to its output). The default hands nothing.
    fn hand_on(&mut self, sink: &mut dyn std::io::Write, count: usize) -> std::io::Result<()> {
        let _ = (sink, count);
        Ok(())
    }
    /// A CHOICE value of the type `type_name` whose variant `name` (position `tag`) is a
    /// SEQUENCE OF VisibleString: the strings of `items` in order. The default is `variant`,
    /// a list of string elements and `end_variant`; the XML stream in its standard-XML mode
    /// writes the type's tag around one member tag for each string, with no list around
    /// them (`<Target><subjects>a</subjects><subjects>b</subjects></Target>`).
    fn variant_string_list<I>(&mut self, type_name: &str, name: &str, tag: u8, items: I)
    where
        I: IntoIterator,
        I::Item: AsRef<[u8]>,
    {
        self.variant(type_name, name, tag);
        self.begin_list(false);
        for item in items {
            self.element();
            self.visible_str(item.as_ref());
        }
        self.end_list();
        self.end_variant();
    }
    /// An element of a list of `CRef<>`s that holds no object (an empty reference that the
    /// reader kept, `pssm::Seqdesc::Empty`), in place of `element` and its value: NCBI's stream
    /// throws there (`NullElement`), so the encoder's error is `SerialError::NullElement`, and
    /// `Traced` writes nothing after it. `element_type` is the name of the list's element type
    /// (empty for an anonymous type), which decides whether NCBI's XML stream has a frame for the
    /// element (`Traced::name_element`).
    fn null_element(&mut self, element_type: &str);
    /// `GetPosition()` where the encoder stands, in the bytes it holds (`Traced` adds what it
    /// handed on); `None` for an encoder whose position `Traced` does not report.
    fn position(&mut self) -> Option<StreamPosition> {
        None
    }
    /// The first error raised while writing (the document is then not written), if any.
    fn error(&self) -> Option<SerialError>;
    /// The bytes that `visible_str_replacing` and `string_store_replacing` replaced so far
    /// (`Replaced`), which the encoder forgets. `Traced` takes them after each string and
    /// returns `ReplacedChar`s from its `take_replaced_chars`.
    fn take_replaced(&mut self) -> Vec<Replaced> {
        Vec::new()
    }
}

/// A frame of NCBI's object stack that `GetStackTraceASN` reads.
#[derive(Clone, Debug)]
enum Frame {
    /// A member: `.name`.
    Named(String),
    /// A SEQUENCE value of the type (`begin_sequence`, `eFrameClass`): nothing in the trace.
    Class(String),
    /// A CHOICE value of the type `choice` and its variant `name` (`variant`, `eFrameChoice` and
    /// `eFrameChoiceVariant`): `.name`.
    Variant { choice: String, name: String },
    /// A container: nothing.
    List,
    /// A container element: `.E`, unless `shown` is false (an element of a named type in
    /// NCBI's XML stream, `Traced::name_element`).
    Element { shown: bool },
    /// A named type that is neither a SEQUENCE nor a CHOICE (`Encoder::begin_named`): nothing in
    /// the trace; NCBI's XML stream alone pushes a frame for it (`eFrameNamed`, whose name the
    /// exception's frames show: `Seq-descr`), the ASN.1 text and binary streams none.
    Type(&'static str),
}

/// The encoder `E` (ASN.1 text, ASN.1 binary, or XML in the archive's object mode, whose
/// `variant_string_list` is the trait's sequence of calls; not the standard-XML mode of outfmt
/// 14 / 16, which writes it otherwise) with NCBI's object stack, for the posts of the bytes that
/// it replaces (`ReplacedChar`). A caller that posts NCBI's messages writes a document through
/// it when one of its strings has such a byte (the archive and the search strategy decide that
/// before they write, `report::archive`; the mask writers, the Seq-annots of outfmt 8 and 9
/// (`seqalign::write_query_to`) and psiblast's PSSM file (`pssm::asn_text_posting`) write a
/// document again when the first writing reported `SerialError::NonVisibleChar`); it writes the
/// strings with `visible_str_replacing`. The encoders by themselves keep no stack and replace nothing, so that
/// the documents without such a byte are written as before.
///
/// The positions count from the document's start: `flush_to` counts the bytes that it hands on
/// (the bytes for the binary stream, the ends of line for the text and XML streams), and `base`
/// is the length of the encoder's buffer before the document.
///
/// NCBI posts each byte when it writes it, after the bytes that its output buffer had written to
/// the output by then (`out_buffer`): `Traced` follows that buffer from the bytes of each call
/// (`Encoder::held`) and hands its sink no more than NCBI had written (`flush_to`, which hands
/// what NCBI has written so far; at the end, everything). A `Traced` with a post function
/// (`posting`) posts each byte there, after handing on exactly the bytes that NCBI had written
/// before it, so that the output and the posts in one file (`2>&1`) follow NCBI's order; without
/// one the posts wait for `take_replaced_chars` / `take_posted` (`Posted` holds those bytes).
///
/// NCBI's XML stream pushes no container element frame for an element of a named type
/// (`WillHaveName`, objostrxml.cpp:1186-1209), so its trace has no `.E` there: `Traced::xml`
/// drops the frame of an element whose value starts with a named type
/// (`begin_sequence`, `begin_named`, `variant`).
pub struct Traced<E, P = fn(&ReplacedChar) -> std::io::Result<()>> {
    inner: E,
    base: u64,
    /// Whether the stream is NCBI's XML stream (`name_element`).
    xml: bool,
    /// The bytes and the ends of line of the encoder's buffer that `flush_to` handed on.
    flushed_bytes: u64,
    flushed_lines: u64,
    root: String,
    frames: Vec<Frame>,
    /// NCBI's output buffer of the document (`E::STREAM`).
    buffer: Option<OutBuffer>,
    /// The bytes of the document handed on.
    handed: u64,
    /// The bytes at the start of `inner.held()` that precede the document (`Encoder::held`).
    foreign: usize,
    posted: Vec<Posted>,
    post: Option<P>,
    /// NCBI's throw (`Thrown`: at an element that holds no object, `null_element`, or at a delay
    /// buffer whose read failed, `delayed_read_failed`), after which nothing is written.
    thrown: Option<Box<Thrown>>,
}

/// The empty string of the XML stream's messages.
fn no_string() -> std::sync::Arc<[u8]> {
    std::sync::Arc::from(&[][..])
}

/// A sink that counts what it hands on (`Traced::flush_to`).
struct Counted<'s> {
    sink: &'s mut dyn std::io::Write,
    bytes: u64,
    lines: u64,
}

impl std::io::Write for Counted<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.sink.write(buf)?;
        self.bytes += written as u64;
        self.lines += buf[..written].iter().filter(|&&byte| byte == b'\n').count() as u64;
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.sink.flush()
    }
}

impl<E: Encoder> Traced<E> {
    /// `inner` (ASN.1 text or binary), which appends a document to a buffer of `base` bytes.
    pub fn new(inner: E, base: usize) -> Self {
        Self {
            foreign: inner.held().len(),
            inner,
            base: base as u64,
            xml: false,
            flushed_bytes: 0,
            flushed_lines: 0,
            root: String::new(),
            frames: Vec::new(),
            buffer: E::STREAM.map(OutBuffer::new),
            handed: 0,
            posted: Vec::new(),
            post: None,
            thrown: None,
        }
    }

    /// `new` of NCBI's XML object stream (`name_element`).
    pub fn xml(inner: E, base: usize) -> Self {
        Self {
            xml: true,
            ..Self::new(inner, base)
        }
    }

    /// The same encoder, which posts each replaced byte with `post` as NCBI does
    /// (`flush_to`).
    pub fn posting<Q>(self, post: Q) -> Traced<E, Q>
    where
        Q: FnMut(&ReplacedChar) -> std::io::Result<()>,
    {
        Traced {
            inner: self.inner,
            base: self.base,
            xml: self.xml,
            flushed_bytes: self.flushed_bytes,
            flushed_lines: self.flushed_lines,
            root: self.root,
            frames: self.frames,
            buffer: self.buffer,
            handed: self.handed,
            foreign: self.foreign,
            posted: self.posted,
            post: Some(post),
            thrown: self.thrown,
        }
    }
}

/// A frame of NCBI's object stack that shows in its path (`object_stack_trace`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StackStep<'a> {
    /// `eFrameClassMember` or `eFrameChoiceVariant` with its member's name.
    Member(&'a str),
    /// `eFrameArrayElement`.
    Element,
}

/// `CObjectStack::GetStackTraceASN`, the one writer of it for NCBI's output and input streams
/// (`Traced` and the typed reader): the document's type, then `.name` for each member or
/// variant frame and `.E` for each container element.
///
/// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-120
/// (A member without a name, `.[tag]`, is not in LOSAT's types.)
pub(crate) fn object_stack_trace<'a>(
    root: &str,
    steps: impl IntoIterator<Item = StackStep<'a>>,
) -> String {
    let mut stack = root.to_string();
    for step in steps {
        match step {
            StackStep::Member(name) => {
                stack.push('.');
                stack.push_str(name);
            }
            StackStep::Element => stack.push_str(".E"),
        }
    }
    stack
}

impl<E: Encoder, P> Traced<E, P> {
    /// `CObjectStack::GetStackTraceASN` of the frames (`object_stack_trace`).
    pub fn stack_trace(&self) -> String {
        object_stack_trace(
            &self.root,
            self.frames.iter().filter_map(|frame| match frame {
                Frame::Named(name) | Frame::Variant { name, .. } => Some(StackStep::Member(name)),
                Frame::Element { shown: true } => Some(StackStep::Element),
                Frame::Element { shown: false }
                | Frame::List
                | Frame::Type(_)
                | Frame::Class(_) => None,
            }),
        )
    }

    /// The writes of NCBI's buffer of the document to its output so far (`OutBuffer::pieces`):
    /// NCBI's insertions of it into its output stream with the bytes left behind in the buffer;
    /// none for a stream without a model of its buffer.
    pub fn pieces(&self) -> &[BufferPiece] {
        self.buffer.as_ref().map_or(&[], OutBuffer::pieces)
    }

    /// The byte sizes of [`Self::pieces`], in order (`OutBuffer::piece_sizes`).
    pub fn piece_sizes(&self) -> Vec<usize> {
        self.buffer
            .as_ref()
            .map_or_else(Vec::new, OutBuffer::piece_sizes)
    }

    /// The stream's position at the first of those writes (`OutBuffer::first_write`).
    pub fn first_write(&self) -> Option<u64> {
        self.buffer.as_ref().and_then(OutBuffer::first_write)
    }

    /// The bytes of the document that `inner` holds (`Encoder::held`, which `Traced` itself does
    /// not forward: what it hands on is `flush_to`'s).
    pub(crate) fn held_len(&self) -> usize {
        self.inner.held().len()
    }

    /// NCBI's throw (`Thrown`), if the document met one.
    pub fn take_thrown(&mut self) -> Option<Box<Thrown>> {
        self.thrown.take()
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:366-381
    /// `error`, an exception other than a `CSerialException`, as it leaves NCBI's writer of the
    /// stream `stream`: `END_OBJECT_FRAME`'s `NCBI_RETHROW_SAME` with `GetFrameInfo()` at each
    /// frame open, innermost first, then the root's (`write_frame_site`); `None` for a stream
    /// whose frames LOSAT does not place. A container's frames name no type: the container of a
    /// member is anonymous, and its elements on the way to a delayed member are objects held by
    /// `CRef<>`s, whose pointer type has no name (NCBI BLAST+ 2.17.0, line RG-F3 item O-2 probe 1).
    fn rethrown(
        &self,
        stream: Stream,
        mut error: crate::ncbi_exception::NcbiException,
    ) -> Option<crate::ncbi_exception::NcbiException> {
        use super::object_frame::{frame_info, FrameKind};
        let through = |error: crate::ncbi_exception::NcbiException,
                       kind: FrameKind,
                       type_name: Option<&str>,
                       member: Option<&str>| {
            write_frame_site(stream, kind).map(|(file, line, function)| {
                error.rethrow_same(file, line, function, frame_info(kind, type_name, member))
            })
        };
        for frame in self.frames.iter().rev() {
            error = match frame {
                Frame::Named(name) => through(error, FrameKind::ClassMember, None, Some(name))?,
                Frame::Class(name) => through(error, FrameKind::Class, Some(name), None)?,
                Frame::Variant { choice, name } => {
                    let error = through(error, FrameKind::ChoiceVariant, None, Some(name))?;
                    through(error, FrameKind::Choice, Some(choice), None)?
                }
                Frame::List => through(error, FrameKind::Array, Some(""), None)?,
                Frame::Element { .. } => through(error, FrameKind::ArrayElement, Some(""), None)?,
                Frame::Type(_) => error,
            };
        }
        through(error, FrameKind::Named, Some(&self.root), None)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1211-1227
    // NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1366-1376
    /// `GetFrameName()` of the root's frame and of each frame open, outermost first
    /// (`NullElement::frames`). NCBI's XML stream pushes no frame for a container of LOSAT's types
    /// (an anonymous SET OF or SEQUENCE OF: its named lists are the types of `begin_named`) nor
    /// for an element of a named type, and a frame for each named type of `begin_named`; the
    /// ASN.1 text and binary streams the reverse.
    fn frame_names(&self) -> Vec<String> {
        use super::object_frame::{frame_name, FrameKind};
        let xml = self.xml;
        std::iter::once(frame_name(FrameKind::Named, Some(&self.root), None))
            .chain(self.frames.iter().filter_map(|frame| match frame {
                Frame::Named(name) | Frame::Variant { name, .. } => {
                    Some(frame_name(FrameKind::ClassMember, None, Some(name)))
                }
                Frame::Class(_) => None,
                Frame::List => (!xml).then(|| frame_name(FrameKind::Array, None, None)),
                Frame::Element { shown } => {
                    (!xml || *shown).then(|| frame_name(FrameKind::ArrayElement, None, None))
                }
                Frame::Type(name) => xml.then(|| frame_name(FrameKind::Named, Some(name), None)),
            }))
            .map(str::to_string)
            .collect()
    }

    /// The bytes replaced so far (`ReplacedChar`) and not posted, which it forgets.
    pub fn take_replaced_chars(&mut self) -> Vec<ReplacedChar> {
        self.take_posted()
            .into_iter()
            .map(|posted| posted.replaced)
            .collect()
    }

    /// `take_replaced_chars` with the bytes that NCBI had written before each post (`Posted`).
    pub fn take_posted(&mut self) -> Vec<Posted> {
        std::mem::take(&mut self.posted)
    }

    /// `EndContainerElement` of the open element, if any.
    fn end_element(&mut self) {
        if matches!(self.frames.last(), Some(Frame::Element { .. })) {
            self.frames.pop();
        }
    }

    /// The value of the element just opened starts with the named type `type_name`: NCBI's
    /// XML stream has no element frame for it.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:1186-1202
    /// (`WriteContainerContents`, 1255-1316, pushes `eFrameArrayElement` under the same test.)
    /// An anonymous type (an empty name, or the open member's tag that the writers pass for one,
    /// `xml::anonymous_type`) has no name, so its element keeps its frame.
    fn name_element(&mut self, type_name: &str) {
        if self.xml && !type_name.is_empty() && !super::xml::anonymous_type(type_name) {
            if let Some(Frame::Element { shown }) = self.frames.last_mut() {
                *shown = false;
            }
        }
    }
}

impl<E: Encoder, P> Traced<E, P> {
    /// `write` (one call of `inner`, which writes `piece`), followed in NCBI's buffer; nothing
    /// after NCBI's stream threw (`null_element`).
    fn call(&mut self, piece: Piece, write: impl FnOnce(&mut E)) {
        if self.thrown.is_some() {
            return;
        }
        let before = self.inner.held().len();
        let last = self.inner.held().last().copied();
        write(&mut self.inner);
        self.follow(piece, before, last, &[]);
    }

    /// The buffer calls of the bytes that `inner` wrote after `before` held bytes, whose last
    /// was `last` (`out_buffer::OutBuffer::put`); what NCBI had written at each of `posts` (the
    /// offsets of `held` where NCBI posts).
    fn follow(
        &mut self,
        piece: Piece,
        before: usize,
        last: Option<u8>,
        posts: &[usize],
    ) -> Vec<u64> {
        let mut written = Vec::new();
        let Some(buffer) = self.buffer.as_mut() else {
            written.resize(posts.len(), 0);
            return written;
        };
        let held = self.inner.held();
        // XML's `value` attribute replaces the `>` that ended the bytes before the call
        // (`OpenTagEndBack`); a string's wrap inserts its end of line after them.
        let start = if before > 0 && held.get(before - 1).copied() != last {
            before - 1
        } else {
            before
        };
        let offsets: Vec<usize> = posts.iter().map(|&at| at - start).collect();
        buffer.put(piece, &held[start..], &offsets, &mut written);
        let byte_block_from = buffer.byte_block_from();
        let has_byte_block = byte_block_from < buffer.pieces().len();
        if has_byte_block {
            let trace = self.stack_trace();
            if let Some(buffer) = self.buffer.as_mut() {
                for bp in &mut buffer.pieces_mut()[byte_block_from..] {
                    bp.byte_block_trace = Some(trace.clone());
                }
            }
        }
        written
    }

    /// A string written by `write`: the bytes that `inner` replaced in it, with the stack trace,
    /// one shared copy of the string and what NCBI had written when it posted each.
    fn string(&mut self, string: &[u8], write: impl FnOnce(&mut E)) {
        if self.thrown.is_some() {
            return;
        }
        let before = self.inner.held().len();
        let last = self.inner.held().last().copied();
        write(&mut self.inner);
        let replaced = self.inner.take_replaced();
        let posts: Vec<usize> = replaced.iter().map(|byte| byte.at).collect();
        let written = self.follow(Piece::Str, before, last, &posts);
        if replaced.is_empty() {
            return;
        }
        let trace = self.stack_trace();
        let shared: std::sync::Arc<[u8]> = std::sync::Arc::from(string);
        for (byte, written) in replaced.into_iter().zip(written) {
            let position = match byte.position {
                StreamPosition::Byte(offset) => {
                    StreamPosition::Byte(self.flushed_bytes + offset - self.base)
                }
                StreamPosition::Line(line) => StreamPosition::Line(self.flushed_lines + line),
            };
            self.posted.push(Posted {
                written,
                replaced: ReplacedChar {
                    trace: trace.clone(),
                    byte: byte.byte,
                    position,
                    string: if byte.with_string {
                        shared.clone()
                    } else {
                        no_string()
                    },
                },
            });
        }
    }

    /// Hands on the bytes of the document up to `target` (bytes from its start).
    /// Each of NCBI's writes of its buffer (`pieces`) is handed on by itself: one insertion of
    /// NCBI's into the report's stream (`ostream::ReportStream`).
    fn hand_to(&mut self, sink: &mut dyn std::io::Write, target: u64) -> std::io::Result<()> {
        let mut end = 0u64;
        let ends: Vec<u64> = self
            .pieces()
            .iter()
            .map(|piece| {
                end += piece.size as u64;
                end
            })
            .filter(|&end| end > self.handed && end < target)
            .collect();
        for end in ends {
            self.hand_piece(sink, end)?;
        }
        self.hand_piece(sink, target)
    }

    /// Hands on the bytes of the document up to `target`, in one hand-on.
    fn hand_piece(&mut self, sink: &mut dyn std::io::Write, target: u64) -> std::io::Result<()> {
        if target > self.handed {
            let count = (target - self.handed) as usize + std::mem::take(&mut self.foreign);
            // NCBI's buffer never writes more than the encoder holds: a larger count is a model
            // (`out_buffer`) that drifted from the encoder.
            debug_assert!(
                count <= self.inner.held().len(),
                "Traced hands on {count} bytes of the {} held",
                self.inner.held().len()
            );
            let count = count.min(self.inner.held().len());
            self.inner.hand_on(sink, count)?;
            self.handed = target;
        }
        Ok(())
    }
}

impl<E: Encoder, P: FnMut(&ReplacedChar) -> std::io::Result<()>> Encoder for Traced<E, P> {
    fn writes_defaults(&self) -> bool {
        self.inner.writes_defaults()
    }
    fn data_format(&self) -> Option<crate::report::serial::asn_input::Encoding> {
        self.inner.data_format()
    }
    fn begin_document(&mut self, type_name: &str) {
        self.root = type_name.to_string();
        self.call(Piece::Header, |inner| inner.begin_document(type_name));
    }
    fn end_document(&mut self) {
        self.call(Piece::End, |inner| inner.end_document());
    }
    fn begin_sequence(&mut self, type_name: &str) {
        self.name_element(type_name);
        self.frames.push(Frame::Class(type_name.to_string()));
        self.call(Piece::Markup, |inner| inner.begin_sequence(type_name));
    }
    fn end_sequence(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_sequence());
        self.frames.pop();
    }
    fn begin_named(&mut self, type_name: &'static str) {
        self.name_element(type_name);
        self.frames.push(Frame::Type(type_name));
        self.call(Piece::Markup, |inner| inner.begin_named(type_name));
    }
    fn end_named(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_named());
        self.frames.pop();
    }
    fn member(&mut self, name: &str, tag: u8) {
        self.frames.push(Frame::Named(name.to_string()));
        self.call(Piece::Markup, |inner| inner.member(name, tag));
    }
    fn end_member(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_member());
        self.frames.pop();
    }
    fn begin_list(&mut self, set: bool) {
        self.frames.push(Frame::List);
        self.call(Piece::Markup, |inner| inner.begin_list(set));
    }
    fn element(&mut self) {
        self.end_element();
        self.frames.push(Frame::Element { shown: true });
        self.call(Piece::Markup, |inner| inner.element());
    }
    fn end_list(&mut self) {
        self.end_element();
        self.frames.pop();
        self.call(Piece::Markup, |inner| inner.end_list());
    }
    fn variant(&mut self, type_name: &str, name: &str, tag: u8) {
        self.name_element(type_name);
        self.frames.push(Frame::Variant {
            choice: type_name.to_string(),
            name: name.to_string(),
        });
        self.call(Piece::Markup, |inner| inner.variant(type_name, name, tag));
    }
    fn end_variant(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_variant());
        self.frames.pop();
    }
    fn int(&mut self, value: i32) {
        self.call(Piece::Number, |inner| inner.int(value));
    }
    fn big_int(&mut self, value: i64) {
        self.call(Piece::Number, |inner| inner.big_int(value));
    }
    fn enumerated(&mut self, name: &str, value: i32) {
        self.call(Piece::Enum, |inner| inner.enumerated(name, value));
    }
    fn named_int(&mut self, name: &str, value: i32) {
        self.call(Piece::Enum, |inner| inner.named_int(name, value));
    }
    fn boolean(&mut self, value: bool) {
        self.call(Piece::Bool, |inner| inner.boolean(value));
    }
    fn visible_str(&mut self, value: &[u8]) {
        self.string(value, |inner| inner.visible_str_replacing(value));
    }
    fn real(&mut self, value: f64) {
        self.call(Piece::Real, |inner| inner.real(value));
    }
    fn octets(&mut self, value: &[u8]) {
        self.call(Piece::Octets, |inner| inner.octets(value));
    }
    fn string_store(&mut self, value: &[u8]) {
        self.string(value, |inner| inner.string_store_replacing(value));
    }
    fn null(&mut self) {
        self.call(Piece::Bool, |inner| inner.null());
    }
    fn utf8_str(&mut self, value: &[u8]) {
        self.string(value, |inner| inner.utf8_str_replacing(value));
    }
    fn delayed(&mut self, bytes: &[u8]) {
        self.call(Piece::Raw, |inner| inner.delayed(bytes));
    }
    fn unread_delayed(&mut self) {
        self.call(Piece::Raw, |inner| inner.unread_delayed());
    }
    /// NCBI's exception of the read (`DelayedReadError`) as it leaves the writer's frames open
    /// (the member's own frame is not, `Encoder::delayed_read_failed`): a `CSerialException`
    /// with the name of each (`frame_names`, as `NullElement::exception`), another exception
    /// rethrown at each (`rethrown`); the bytes held are those NCBI's stream wrote before, and the
    /// calls after it write nothing (`call`, `string`). Where LOSAT does not place the frames (the
    /// XML stream and an exception other than a `CSerialException`), the member is not written
    /// (`unread_delayed`).
    fn delayed_read_failed(&mut self, error: &DelayedReadError) {
        if self.thrown.is_some() {
            return;
        }
        let stream = E::STREAM.filter(|stream| *stream != Stream::Json);
        let exception = match (stream, error) {
            (Some(_), DelayedReadError::Serial(exception)) => {
                let mut exception = exception.clone();
                let names = self.frame_names();
                let names: Vec<&str> = names.iter().rev().map(String::as_str).collect();
                exception.add_frame_infos(&names);
                Some(exception)
            }
            (Some(stream), DelayedReadError::Other(exception)) => {
                self.rethrown(stream, exception.clone())
            }
            (None, _) => None,
        };
        let Some(exception) = exception else {
            self.unread_delayed();
            return;
        };
        self.inner.delayed_read_failed(error);
        self.thrown = Some(Box::new(Thrown::DelayedRead(exception)));
    }
    /// NCBI's `WriteContainer` at an element that holds no object: the element's frame (pushed
    /// before NCBI's loop over the elements; none in NCBI's XML stream for an element of a named
    /// type, `name_element`), then NCBI's throw (`NullElement`) where the stream stands; the
    /// bytes held are those NCBI flushed, and the calls after it write nothing (`call`, `string`).
    /// An encoder that `Traced` does not follow (`Encoder::STREAM`), or the JSON stream, whose
    /// documents have no NULL element (`json::JsonEncoder::null_element`), keeps its error alone.
    fn null_element(&mut self, element_type: &str) {
        if self.thrown.is_some() {
            return;
        }
        self.end_element();
        self.frames.push(Frame::Element { shown: true });
        self.name_element(element_type);
        let element_frame = matches!(self.frames.last(), Some(Frame::Element { shown: true }));
        let position = self.inner.position();
        self.inner.null_element(element_type);
        let (Some(stream), Some(position)) =
            (E::STREAM.filter(|stream| *stream != Stream::Json), position)
        else {
            return;
        };
        let position = match position {
            StreamPosition::Line(line) => StreamPosition::Line(self.flushed_lines + line),
            StreamPosition::Byte(offset) => {
                StreamPosition::Byte(self.flushed_bytes + offset - self.base)
            }
        };
        self.thrown = Some(Box::new(Thrown::NullElement(NullElement {
            stream,
            element_frame,
            trace: self.stack_trace(),
            position,
            frames: self.frame_names(),
        })));
    }
    /// The posts of the read, made before the member starts (`Encoder::delayed_posts`), after
    /// what NCBI's buffer had written then.
    fn delayed_posts(&mut self, posts: &[ReplacedChar]) {
        if self.thrown.is_some() {
            return;
        }
        let written = self.buffer.as_ref().map_or(0, OutBuffer::written);
        self.posted.extend(posts.iter().map(|replaced| Posted {
            written,
            replaced: replaced.clone(),
        }));
    }
    fn buffered(&self) -> usize {
        self.inner.buffered()
    }
    /// Hands on what NCBI's buffer has written to its output so far, and with `all` (the
    /// document is complete) everything; with a post function first each byte that waits, after
    /// exactly the bytes that NCBI had written before posting it, the sink flushed.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        let mut counted = Counted {
            sink,
            bytes: 0,
            lines: 0,
        };
        let handed = self.hand_on_posting(&mut counted, all);
        self.flushed_bytes += counted.bytes;
        self.flushed_lines += counted.lines;
        handed
    }
    fn error(&self) -> Option<SerialError> {
        self.inner.error()
    }
}

impl<E: Encoder, P: FnMut(&ReplacedChar) -> std::io::Result<()>> Traced<E, P> {
    /// `flush_to` to the counting sink.
    fn hand_on_posting(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        if self.post.is_some() {
            for posted in std::mem::take(&mut self.posted) {
                self.hand_to(sink, posted.written)?;
                sink.flush()?;
                if let Some(post) = self.post.as_mut() {
                    post(&posted.replaced)?;
                }
            }
        }
        if all {
            if let Some(written) = self.buffer.as_ref().map(OutBuffer::written) {
                self.hand_to(sink, written)?;
            }
            self.foreign = 0;
            return self.inner.flush_to(sink, true);
        }
        match self.buffer.as_ref().map(OutBuffer::written) {
            Some(written) => self.hand_to(sink, written),
            None => self.inner.flush_to(sink, false),
        }
    }
}

#[cfg(test)]
#[path = "encoder_tests.rs"]
mod tests;
