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

use super::object_stack::ObjectStack;
use super::visible::{BadChars, Position};

/// An error that NCBI's object stream raises while it writes a value, which LOSAT reports
/// instead of writing the document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerialError {
    /// `CObjectOStreamJson::WriteDouble2` on a NaN or an infinity
    /// (`real::JsonDoubleError`).
    JsonDouble(super::real::JsonDoubleError),
    /// A VisibleString byte outside 0x20..0x7E in the ASN.1 text or binary stream (in XML a
    /// control character other than TAB, LF and CR), which NCBI writes as `#` with an error
    /// message (`eFNP_ReplaceAndWarn`): a plain encoder reports it here; a writer that posts
    /// NCBI's messages writes the document again through `Traced` (`visible_str_replacing`, the
    /// messages from `visible::BadChars`), the others refuse it.
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

/// A byte of a VisibleString that the NCBI form of an encoder (`asn_text::AsnTextEncoder::ncbi`,
/// `asn_ber::AsnBinaryEncoder::ncbi`, the XML object stream) wrote as `#`
/// (`visible_str_replacing`), as the encoder knows it: `Traced` adds the stack trace and the
/// string (`visible::BadChars`, the one producer of NCBI's message).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Replaced {
    /// The byte's place in the string.
    pub index: usize,
    /// `GetPosition()` of the stream at the byte, counted from the document's start: the line of
    /// the ASN.1 text and XML streams, the bytes written of the binary stream.
    pub position: Position,
    /// The bytes of the document that NCBI's buffer (`stream_buffer::StreamBuffer`) has handed
    /// to its output when NCBI posts the message; 0 for an encoder without the model (the XML
    /// stream).
    pub after: usize,
    /// Whether NCBI's message quotes the string (the text and binary streams, not XML).
    pub with_string: bool,
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
    fn begin_document(&mut self, type_name: &'static str);
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
    fn member(&mut self, name: &'static str, tag: u8);
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
    fn variant(&mut self, type_name: &'static str, name: &'static str, tag: u8);
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
    /// `flush_to(sink, false)` that also keeps the last `keep` bytes, for a document that may be
    /// written again (`Traced`) from where the bytes handed out end. The default is `flush_to`.
    fn flush_keeping(&mut self, sink: &mut dyn std::io::Write, keep: usize) -> std::io::Result<()> {
        let _ = keep;
        self.flush_to(sink, false)
    }
    /// Hands the document's bytes before its byte `through` to `sink`, for the NCBI form of an
    /// encoder, which models where NCBI's buffer has handed its bytes to the output (a message of a
    /// replaced byte falls there among them, `Replaced::after`). The default hands nothing.
    fn flush_through(
        &mut self,
        sink: &mut dyn std::io::Write,
        through: usize,
    ) -> std::io::Result<()> {
        let _ = (sink, through);
        Ok(())
    }
    /// A CHOICE value of the type `type_name` whose variant `name` (position `tag`) is a
    /// SEQUENCE OF VisibleString: the strings of `items` in order. The default is `variant`,
    /// a list of string elements and `end_variant`; the XML stream in its standard-XML mode
    /// writes the type's tag around one member tag for each string, with no list around
    /// them (`<Target><subjects>a</subjects><subjects>b</subjects></Target>`).
    fn variant_string_list<I>(
        &mut self,
        type_name: &'static str,
        name: &'static str,
        tag: u8,
        items: I,
    ) where
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
    /// (`Replaced`), which the encoder forgets. `Traced` takes them after each string and keeps
    /// them as `visible::BadChars` (`Traced::take_bad_chars`).
    fn take_replaced(&mut self) -> Vec<Replaced> {
        Vec::new()
    }
    /// The strings written so far with bytes written as `#`, with NCBI's stack trace
    /// (`visible::BadChars`), which the encoder forgets: `Traced`'s; none for the others.
    fn take_bad_chars(&mut self) -> Vec<BadChars> {
        Vec::new()
    }
}

/// The encoder `E` with NCBI's object stack (`object_stack::ObjectStack`), for the messages of
/// the bytes that it writes as `#` (`visible::BadChars`): the NCBI form of the ASN.1 text or
/// binary encoder (which also models NCBI's buffer, so that each message knows the bytes NCBI
/// had handed out before it), or the XML encoder in the archive's object mode (whose
/// `variant_string_list` is the trait's sequence of calls). A writer writes a document with a
/// plain encoder, which keeps no stack and reports such a byte as `SerialError::NonVisibleChar`,
/// and writes it again through `Traced` only then (a cold path), so that the documents without
/// such a byte cost nothing more.
///
/// (NCBI's XML stream pushes no container element frame for an element of a named type
/// (`WillHaveName`, objostrxml.cpp:1186-1209), so its trace has no `.E` there; no XML document
/// that LOSAT traces has such a list with a string that can hold a control character.)
pub struct Traced<E> {
    inner: E,
    stack: ObjectStack,
    bad_chars: Vec<BadChars>,
}

impl<E: Encoder> Traced<E> {
    /// `inner` with an empty stack.
    pub fn new(inner: E) -> Self {
        Self {
            inner,
            stack: ObjectStack::default(),
            bad_chars: Vec::new(),
        }
    }

    /// Takes the bytes that `inner` replaced in `string`, which it just wrote, with the stack
    /// trace (`GetStackTraceASN`) and the string, as one `BadChars`.
    fn collect(&mut self, string: &[u8]) {
        let replaced = self.inner.take_replaced();
        let Some(first) = replaced.first() else {
            return;
        };
        let mut bad = BadChars::written(self.stack.trace(), string, first.with_string);
        for byte in replaced {
            bad.note_after(byte.index, byte.position, byte.after);
        }
        self.bad_chars.push(bad);
    }
}

impl<E: Encoder> Encoder for Traced<E> {
    fn writes_defaults(&self) -> bool {
        self.inner.writes_defaults()
    }
    fn begin_document(&mut self, type_name: &'static str) {
        self.stack.set_root(type_name);
        self.inner.begin_document(type_name);
    }
    fn end_document(&mut self) {
        self.inner.end_document();
    }
    fn begin_sequence(&mut self, type_name: &str) {
        self.inner.begin_sequence(type_name);
    }
    fn end_sequence(&mut self) {
        self.inner.end_sequence();
    }
    fn begin_named(&mut self, type_name: &str) {
        self.inner.begin_named(type_name);
    }
    fn end_named(&mut self) {
        self.inner.end_named();
    }
    // The frames of `CObjectStack` that `GetStackTraceASN` reads: one per member and variant,
    // one for the elements of a list (`E` from its first element; `""` adds nothing before it).
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:810-819
    // ```c++
    // void CObjectOStream::WriteContainer(const CContainerTypeInfo* cType,
    //                                     TConstObjectPtr containerPtr)
    // {
    //     BEGIN_OBJECT_FRAME2(eFrameArray, cType);
    //     BeginContainer(cType);
    //
    //     CContainerTypeInfo::CConstIterator i;
    //     if ( cType->InitIterator(i, containerPtr) ) {
    //         TTypeInfo elementType = cType->GetElementType();
    //         BEGIN_OBJECT_FRAME2(eFrameArrayElement, elementType);
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:906-911
    // ```c++
    // void CObjectOStream::WriteClassMember(const CMemberId& memberId,
    //                                       TTypeInfo memberType,
    //                                       TConstObjectPtr memberPtr)
    // {
    //     BEGIN_OBJECT_FRAME2(eFrameClassMember, memberId);
    //     BeginClassMember(memberId);
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/choice.cpp:451-453
    // ```c++
    //         variantInfo = choiceType->GetVariantInfo(index);
    //         BEGIN_OBJECT_FRAME_OF2(out, eFrameChoiceVariant, variantInfo->GetId());
    //         out.BeginChoiceVariant(choiceType, variantInfo->GetId());
    // ```
    fn member(&mut self, name: &'static str, tag: u8) {
        self.stack.push(name);
        self.inner.member(name, tag);
    }
    fn end_member(&mut self) {
        self.inner.end_member();
        self.stack.pop();
    }
    fn begin_list(&mut self, set: bool) {
        self.stack.push("");
        self.inner.begin_list(set);
    }
    fn element(&mut self) {
        self.stack.set_top("E");
        self.inner.element();
    }
    fn end_list(&mut self) {
        self.stack.pop();
        self.inner.end_list();
    }
    fn variant(&mut self, type_name: &'static str, name: &'static str, tag: u8) {
        self.stack.push(name);
        self.inner.variant(type_name, name, tag);
    }
    fn end_variant(&mut self) {
        self.inner.end_variant();
        self.stack.pop();
    }
    fn int(&mut self, value: i32) {
        self.inner.int(value);
    }
    fn big_int(&mut self, value: i64) {
        self.inner.big_int(value);
    }
    fn enumerated(&mut self, name: &str, value: i32) {
        self.inner.enumerated(name, value);
    }
    fn named_int(&mut self, name: &str, value: i32) {
        self.inner.named_int(name, value);
    }
    fn boolean(&mut self, value: bool) {
        self.inner.boolean(value);
    }
    fn visible_str(&mut self, value: &[u8]) {
        self.inner.visible_str_replacing(value);
        self.collect(value);
    }
    fn real(&mut self, value: f64) {
        self.inner.real(value);
    }
    fn octets(&mut self, value: &[u8]) {
        self.inner.octets(value);
    }
    fn string_store(&mut self, value: &[u8]) {
        self.inner.string_store_replacing(value);
        self.collect(value);
    }
    fn buffered(&self) -> usize {
        self.inner.buffered()
    }
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        self.inner.flush_to(sink, all)
    }
    fn flush_through(
        &mut self,
        sink: &mut dyn std::io::Write,
        through: usize,
    ) -> std::io::Result<()> {
        self.inner.flush_through(sink, through)
    }
    fn error(&self) -> Option<SerialError> {
        self.inner.error()
    }
    fn take_bad_chars(&mut self) -> Vec<BadChars> {
        std::mem::take(&mut self.bad_chars)
    }
}

/// A document that its writer writes through any encoder (`write_asn_text`, `write_asn_binary`).
pub(crate) trait SerialDocument {
    fn write<E: Encoder>(&self, enc: &mut E);
}

/// `document` in ASN.1 text appended to `out` by the plain encoder, and when it reported
/// `NonVisibleChar`, again by NCBI's form through `Traced` (`write_traced`): the encoder's error and
/// the strings with bytes written as `#`, each with its messages and the bytes of the document
/// NCBI had handed out before each (`visible::BadChars`, `Replaced::after`).
pub(crate) fn write_asn_text<D: SerialDocument>(
    document: &D,
    out: &mut Vec<u8>,
) -> (Option<SerialError>, Vec<BadChars>) {
    let start = out.len();
    let mut enc = super::asn_text::AsnTextEncoder::new(out);
    document.write(&mut enc);
    match enc.error() {
        Some(SerialError::NonVisibleChar(_)) => {
            out.truncate(start);
            write_traced(document, super::asn_text::AsnTextEncoder::ncbi(out))
        }
        error => (error, Vec::new()),
    }
}

/// `write_asn_text` in ASN.1 binary.
pub(crate) fn write_asn_binary<D: SerialDocument>(
    document: &D,
    out: &mut Vec<u8>,
) -> (Option<SerialError>, Vec<BadChars>) {
    let start = out.len();
    let mut enc = super::asn_ber::AsnBinaryEncoder::new(out);
    document.write(&mut enc);
    match enc.error() {
        Some(SerialError::NonVisibleChar(_)) => {
            out.truncate(start);
            write_traced(document, super::asn_ber::AsnBinaryEncoder::ncbi(out))
        }
        error => (error, Vec::new()),
    }
}

/// `document` through `Traced` of `enc` (NCBI's form of an encoder), the cold path of a document
/// with a byte that NCBI writes as `#`.
#[cold]
#[inline(never)]
pub(crate) fn write_traced<D: SerialDocument, E: Encoder>(
    document: &D,
    enc: E,
) -> (Option<SerialError>, Vec<BadChars>) {
    let mut enc = Traced::new(enc);
    document.write(&mut enc);
    (enc.error(), enc.take_bad_chars())
}
