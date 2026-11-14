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
//! ```c++
//! void CObjectOStream::Write(TConstObjectPtr object, TTypeInfo typeInfo)
//! {
//!     // root writer
//!     BEGIN_OBJECT_FRAME2(eFrameNamed, typeInfo);
//!
//!     WriteFileHeader(typeInfo);
//!
//!     WriteObject(object, typeInfo);
//!
//!     EndOfWrite();
//!
//!     END_OBJECT_FRAME();
//!
//!     if ( GetAutoSeparator() )
//!         Separator(*this);
//! }
//! ```
//!
//! The calls below are those pairs, one method each; a type writer calls them in the member
//! order of the ASN.1 module and skips an absent OPTIONAL member. `tag` is the member's or
//! variant's 0-based position in its ASN.1 SEQUENCE or CHOICE (absent OPTIONAL members still
//! count), which only the binary encoder writes (its context tag `[tag]`).

pub use super::out_buffer::Stream;
use super::out_buffer::{OutBuffer, Piece};

/// An error that NCBI's object stream raises while it writes a value, which LOSAT reports
/// instead of writing the document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerialError {
    /// `CObjectOStreamJson::WriteDouble2` on a NaN or an infinity
    /// (`real::JsonDoubleError`).
    JsonDouble(super::real::JsonDoubleError),
    /// A VisibleString byte outside 0x20..0x7E in the ASN.1 text or binary stream (in XML a
    /// control character other than TAB, LF and CR), which NCBI replaces by `#` with an error
    /// message (`eFNP_ReplaceAndWarn`, `ReplacedChar`): `visible_str` reports it here; a
    /// caller that posts NCBI's messages writes the document again through `Traced`
    /// (`visible_str_replacing`), the other callers refuse such strings before the search.
    NonVisibleChar(u8),
    /// A NaN or an infinity in the ASN.1 binary stream, which NCBI marks with its own form
    /// byte and LOSAT does not write (no BLAST report has one).
    NonFiniteReal,
    /// A byte >= 0x80 in a string of the XML streams of outfmt 5 / 14 / 16, whose US-ASCII
    /// encoding cannot hold it: NCBI's `CUtf8::SymbolToChar` throws (rc 255, or a document cut
    /// after the failed element); the callers refuse such strings first, so this does not
    /// occur. The archive's stream has no encoding and writes the byte.
    NonAsciiChar(u8),
    /// A call that the XML stream does not write outside the archive's object mode (a BOOLEAN,
    /// an ENUMERATED value or a CHOICE of outfmt 5).
    XmlUnsupported(&'static str),
    /// An OCTET STRING in a stream that does not write one (the archive of outfmt 11 is the
    /// only document with one; JSON and the XML of outfmt 5 / 14 / 16 do not write it).
    OctetString,
}

impl std::fmt::Display for SerialError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SerialError::JsonDouble(error) => f.write_str(error.message()),
            SerialError::NonFiniteReal => f.write_str(
                "a NaN or an infinite real number in an ASN.1 binary stream is not supported by LOSAT",
            ),
            SerialError::NonVisibleChar(byte) => write!(
                f,
                "a character 0x{byte:02x} outside the visible ASCII range in an ASN.1 VisibleString is not supported by LOSAT"
            ),
            SerialError::NonAsciiChar(byte) => write!(
                f,
                "a character 0x{byte:02x} outside the US-ASCII range in an XML string is not supported by LOSAT"
            ),
            SerialError::XmlUnsupported(what) => {
                write!(f, "an XML {what} is not supported by LOSAT")
            }
            SerialError::OctetString => {
                f.write_str("an OCTET STRING in this serialized format is not supported by LOSAT")
            }
        }
    }
}

impl std::error::Error for SerialError {}

/// Where NCBI's object stream stands when it replaces a byte (`GetPosition`): the line of the
/// ASN.1 text or XML stream (`m_Output.GetLine()`, 1 at the stream's start and one more for each
/// end of line written), the offset of the binary stream (`GetStreamPos()`, 0 at its start).
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:100-103
/// ```c++
/// string CObjectOStreamAsn::GetPosition(void) const
/// {
///     return "line "+NStr::SizetToString(m_Output.GetLine());
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:510-513
/// ```c++
/// string CObjectOStream::GetPosition(void) const
/// {
///     return "byte "+NStr::Int8ToString(NcbiStreamposToInt8(GetStreamPos()));
/// }
/// ```
/// (`CObjectOStreamXml::GetPosition`, objostrxml.cpp:192-195, is the text stream's.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StreamPosition {
    Line(u64),
    Byte(u64),
}

/// A byte of a VisibleString that NCBI's object stream writes as `#` and reports with
/// `ERR_POST_X(7, message)` of `Serial_IStream` (error 803.7: the function is in objistr.cpp).
///
/// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897
/// ```c++
/// char ReplaceVisibleChar(char c, EFixNonPrint fix_method,
///     const CObjectStack* io, const CTempString& str, char subst)
/// {
///     ...
///     if ( fix_method != eFNP_Replace ) {
///         string message;
///         if (io != NULL) {
///             message += io->GetStackTrace() + "\n";
///         }
///         message += "Bad char [0x" +
///                    NStr::NumericToString((unsigned char)c,0,16)+
///                    "] in string";
///         if (io != NULL) {
///             message += " at " + io->GetPosition();
///         }
///         if (!str.empty()) {
///             message += "\n" + str;
///         }
///         switch (fix_method) {
///         case eFNP_ReplaceAndWarn:
///             ...
///             ERR_POST_X(7, message);
///             ...
///     }
///     return subst;
/// }
/// ```
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
    /// (the hex digits in upper case without padding, `NumericToString(c, 0, 16)`) and the
    /// string when there is one.
    pub fn message(&self) -> Vec<u8> {
        let position = match self.position {
            StreamPosition::Line(line) => format!("line {line}"),
            StreamPosition::Byte(byte) => format!("byte {byte}"),
        };
        let mut message = format!(
            "{}\nBad char [0x{:X}] in string at {position}",
            self.trace, self.byte
        )
        .into_bytes();
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
    /// `WriteFileHeader(type)`: the start of a document whose top-level object is of the
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
    /// default does nothing.
    fn begin_named(&mut self, type_name: &str) {
        let _ = type_name;
    }
    /// `EndNamedType`.
    fn end_named(&mut self) {}
    /// `BeginClassMember(id)`: member `name` (position `tag`) of the current SEQUENCE.
    fn member(&mut self, name: &str, tag: u8);
    /// `EndClassMember`.
    fn end_member(&mut self);
    /// `BeginContainer`: a SEQUENCE OF (`set` false) or SET OF (`set` true) value starts.
    fn begin_list(&mut self, set: bool);
    /// `BeginContainerElement`: the next element of the current container.
    fn element(&mut self);
    /// `EndContainer`.
    fn end_list(&mut self);
    /// `BeginChoice` and `BeginChoiceVariant(id)`: a CHOICE value of the ASN.1 type
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
    /// A member or a choice variant: `.name`.
    Named(String),
    /// A container: nothing.
    List,
    /// A container element: `.E`, unless `shown` is false (an element of a named type in
    /// NCBI's XML stream, `Traced::name_element`).
    Element { shown: bool },
}

/// The encoder `E` (ASN.1 text, ASN.1 binary, or XML in the archive's object mode, whose
/// `variant_string_list` is the trait's sequence of calls; not the standard-XML mode of outfmt
/// 14 / 16, which writes it otherwise) with NCBI's object stack, for the posts of the bytes that
/// it replaces (`ReplacedChar`). A caller that posts NCBI's messages writes a document through
/// it when one of its strings has such a byte (the archive and the search strategy decide that
/// before they write, `report::archive`; the mask writers write a document again when the first
/// writing reported `SerialError::NonVisibleChar`); it writes the strings with
/// `visible_str_replacing`. The encoders by themselves keep no stack and replace nothing, so that
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
        }
    }
}

impl<E: Encoder, P> Traced<E, P> {
    /// `CObjectStack::GetStackTraceASN`: the document's type, then `.name` for each member or
    /// variant frame and `.E` for each container element.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-120
    /// ```c++
    ///     string stack = FetchFrameFromBottom(0).HasTypeInfo() ?
    ///         FetchFrameFromBottom(0).m_TypeInfo->GetName() : "?";
    ///     for ( size_t i = 1; i < GetStackDepth(); ++i ) {
    ///         const TFrame& frame = FetchFrameFromBottom(i);
    ///         switch ( frame.m_FrameType ) {
    ///         case TFrame::eFrameClassMember:
    ///         case TFrame::eFrameChoiceVariant:
    ///             {
    ///                 ...
    ///                     const CMemberId& id = *frame.m_MemberId;
    ///                     stack += '.';
    ///                     if ( !id.GetName().empty() ) {
    ///                         stack += id.GetName();
    ///                     }
    ///                     ...
    ///             }
    ///             break;
    ///         case TFrame::eFrameArrayElement:
    ///             stack += ".E";
    ///             break;
    ///         default:
    ///             break;
    ///         }
    ///     }
    ///     return stack;
    /// ```
    /// (A member without a name, `.[tag]`, is not in LOSAT's types.)
    pub fn stack_trace(&self) -> String {
        let mut stack = self.root.clone();
        for frame in &self.frames {
            match frame {
                Frame::Named(name) => {
                    stack.push('.');
                    stack.push_str(name);
                }
                Frame::Element { shown: true } => stack.push_str(".E"),
                Frame::Element { shown: false } | Frame::List => {}
            }
        }
        stack
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
    /// ```c++
    /// bool CObjectOStreamXml::WillHaveName(TTypeInfo elementType)
    /// {
    ///     while ( elementType->GetName().empty() ) {
    ///         if ( elementType->GetTypeFamily() != eTypeFamilyPointer )
    ///             return false;
    ///         elementType = CTypeConverter<CPointerTypeInfo>::SafeCast(
    ///             elementType)->GetPointedType();
    ///     }
    ///     // found named type
    ///     return true;
    /// }
    ///
    /// void CObjectOStreamXml::BeginContainerElement(TTypeInfo elementType)
    /// {
    ///     if ( !WillHaveName(elementType) ) {
    ///         BeginArrayElement(elementType);
    ///     }
    /// }
    /// ```
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
    /// `write` (one call of `inner`, which writes `piece`), followed in NCBI's buffer.
    fn call(&mut self, piece: Piece, write: impl FnOnce(&mut E)) {
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
        written
    }

    /// A string written by `write`: the bytes that `inner` replaced in it, with the stack trace,
    /// one shared copy of the string and what NCBI had written when it posted each.
    fn string(&mut self, string: &[u8], write: impl FnOnce(&mut E)) {
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
    fn hand_to(&mut self, sink: &mut dyn std::io::Write, target: u64) -> std::io::Result<()> {
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
    fn begin_document(&mut self, type_name: &str) {
        self.root = type_name.to_string();
        self.call(Piece::Header, |inner| inner.begin_document(type_name));
    }
    fn end_document(&mut self) {
        self.call(Piece::End, |inner| inner.end_document());
    }
    fn begin_sequence(&mut self, type_name: &str) {
        self.name_element(type_name);
        self.call(Piece::Markup, |inner| inner.begin_sequence(type_name));
    }
    fn end_sequence(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_sequence());
    }
    fn begin_named(&mut self, type_name: &str) {
        self.name_element(type_name);
        self.call(Piece::Markup, |inner| inner.begin_named(type_name));
    }
    fn end_named(&mut self) {
        self.call(Piece::Markup, |inner| inner.end_named());
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
        self.frames.push(Frame::Named(name.to_string()));
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
            self.foreign = 0;
            return self.inner.flush_to(sink, true);
        }
        match self.buffer.as_ref().map(OutBuffer::written) {
            Some(written) => self.hand_to(sink, written),
            None => self.inner.flush_to(sink, false),
        }
    }
}
