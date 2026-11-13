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

/// An error that NCBI's object stream raises while it writes a value, which LOSAT reports
/// instead of writing the document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerialError {
    /// `CObjectOStreamJson::WriteDouble2` on a NaN or an infinity
    /// (`real::JsonDoubleError`).
    JsonDouble(super::real::JsonDoubleError),
    /// A VisibleString byte outside 0x20..0x7E in the ASN.1 text or binary stream (in XML a
    /// control character other than TAB, LF and CR), which NCBI replaces by `#` with an error
    /// message (`eFNP_ReplaceAndWarn`, `ReplacedChar`): an encoder without `Frames` reports it
    /// here, and a caller that posts NCBI's messages writes the document again with `Frames`
    /// (`ReplacedChar`); the other callers refuse such strings before the search.
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
    /// `GetStackTrace()` (`Frames::stack_trace`).
    pub trace: String,
    pub byte: u8,
    pub position: StreamPosition,
    /// The string that holds the byte (empty in XML).
    pub string: Vec<u8>,
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
}

/// What an encoder keeps for the posts of `ReplacedChar`. `NoTrace`, the default, keeps nothing
/// and costs nothing: its encoder reports a replaced byte as `SerialError::NonVisibleChar`.
/// `Frames` follows NCBI's object stack and the stream's position and records each replaced
/// byte instead; a caller that posts NCBI's messages writes a document with `Frames` only when
/// its first writing reported the error.
pub trait Trace: Default {
    /// Whether the encoder records `ReplacedChar`s (`Frames`).
    const FRAMES: bool;
    /// `WriteFileHeader(type)`: the bottom frame, the document's type.
    fn document(&mut self, _type_name: &str) {}
    /// `BeginClassMember` / `BeginChoiceVariant`: a frame named `name`.
    fn push(&mut self, _name: &str) {}
    /// `EndClassMember` / `EndChoiceVariant`.
    fn pop(&mut self) {}
    /// `BeginContainer`.
    fn begin_list(&mut self) {}
    /// `BeginContainerElement` (and the `EndContainerElement` of the element before).
    fn element(&mut self) {}
    /// `EndContainer` (and the `EndContainerElement` of the last element).
    fn end_list(&mut self) {}
    /// The bytes of the document that the encoder handed to its sink and dropped.
    fn flushed(&mut self, _bytes: &[u8]) {}
    /// A replaced byte of `string` at `position` (`Frames` only).
    fn replaced(&mut self, _byte: u8, _position: StreamPosition, _string: &[u8]) {}
    /// The line of the stream after the document's bytes `buffered` that the encoder still
    /// holds (`Frames` only).
    fn line(&mut self, _buffered: &[u8]) -> u64 {
        0
    }
    /// The offset in the stream after the document's bytes `buffered` that the encoder still
    /// holds (`Frames` only).
    fn offset(&self, _buffered: &[u8]) -> u64 {
        0
    }
    /// The replaced bytes recorded so far, which the trace forgets.
    fn take(&mut self) -> Vec<ReplacedChar> {
        Vec::new()
    }
}

/// The encoders' default `Trace`: nothing kept.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoTrace;

impl Trace for NoTrace {
    const FRAMES: bool = false;
}

/// A frame of NCBI's object stack that `GetStackTraceASN` reads.
#[derive(Clone, Debug)]
enum Frame {
    /// A member or a choice variant: `.name`.
    Named(String),
    /// A container: nothing.
    List,
    /// A container element: `.E`.
    Element,
}

/// NCBI's object stack and the stream's position, for the posts of `ReplacedChar`.
#[derive(Clone, Debug, Default)]
pub struct Frames {
    root: String,
    frames: Vec<Frame>,
    /// The document's bytes that went to the sink, and the ends of line among them.
    flushed_bytes: u64,
    flushed_lines: u64,
    /// The ends of line before `counted` in the bytes still held: the bytes before the last end
    /// of line are final (the string wrap inserts its end of line into the last line only).
    counted: usize,
    counted_lines: u64,
    replaced: Vec<ReplacedChar>,
}

impl Frames {
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
                Frame::Element => stack.push_str(".E"),
                Frame::List => {}
            }
        }
        stack
    }
}

impl Trace for Frames {
    const FRAMES: bool = true;

    fn document(&mut self, type_name: &str) {
        self.root = type_name.to_string();
    }

    fn push(&mut self, name: &str) {
        self.frames.push(Frame::Named(name.to_string()));
    }

    fn pop(&mut self) {
        self.frames.pop();
    }

    fn begin_list(&mut self) {
        self.frames.push(Frame::List);
    }

    fn element(&mut self) {
        if matches!(self.frames.last(), Some(Frame::Element)) {
            self.frames.pop();
        }
        self.frames.push(Frame::Element);
    }

    fn end_list(&mut self) {
        if matches!(self.frames.last(), Some(Frame::Element)) {
            self.frames.pop();
        }
        self.frames.pop();
    }

    fn flushed(&mut self, bytes: &[u8]) {
        self.flushed_bytes += bytes.len() as u64;
        self.flushed_lines += bytes.iter().filter(|&&byte| byte == b'\n').count() as u64;
        self.counted = 0;
        self.counted_lines = 0;
    }

    fn replaced(&mut self, byte: u8, position: StreamPosition, string: &[u8]) {
        self.replaced.push(ReplacedChar {
            trace: self.stack_trace(),
            byte,
            position,
            string: string.to_vec(),
        });
    }

    fn line(&mut self, buffered: &[u8]) -> u64 {
        let end = buffered
            .iter()
            .rposition(|&byte| byte == b'\n')
            .map_or(0, |at| at + 1);
        if end > self.counted {
            self.counted_lines += buffered[self.counted..end]
                .iter()
                .filter(|&&byte| byte == b'\n')
                .count() as u64;
            self.counted = end;
        }
        1 + self.flushed_lines + self.counted_lines
    }

    fn offset(&self, buffered: &[u8]) -> u64 {
        self.flushed_bytes + buffered.len() as u64
    }

    fn take(&mut self) -> Vec<ReplacedChar> {
        std::mem::take(&mut self.replaced)
    }
}

/// The calls of NCBI's `CObjectOStream` that the serialized BLAST reports make.
pub trait Encoder {
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
    /// The bytes replaced so far by an encoder with `Frames` (`ReplacedChar`), which it forgets;
    /// none for the other encoders.
    fn take_replaced(&mut self) -> Vec<ReplacedChar> {
        Vec::new()
    }
}
