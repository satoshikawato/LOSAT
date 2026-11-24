//! The reader under the type table: NCBI's object stack (`CObjectIStream`: its frames,
//! `ReadClassSequential`, `ReadChoiceSimple`, `ReadContainer`, `ReadPointer`, `HandleEOF`) for any
//! type of `table`, over the primitives of either encoding (`report::serial::asn_input::Input`):
//! ASN.1 text (`read_document`, `CObjectIStreamAsn`) and ASN.1 binary (`read_ber_document`,
//! `CObjectIStreamAsnBinary`), one object or objects one after another from one stream
//! (`ObjectStream`). The object stack of the table-driven readers: the generic `Value` of any
//! type (here, read without recursion to any depth), and the typed readers of the callers that
//! build their own model as they read (`super::typed`: makeblastdb's mask files, the Seq-id model
//! of `super::seq_id`) run on the same frames, members, pointers and errors. (rpsbproc's reader of
//! the BLAST archive, `report::serial::asn_text_reader` with `seqid::asn_text` for its Seq-ids,
//! and the database deflines' Seq-ids, `seqid::ber::decode_seq_id`, are readers of their own.)
//!
//! It reads as NCBI's reader reads:
//! - the members of a SEQUENCE in module order (named or by their tag `[3]` in text, by their
//!   context tag `a3` in binary: NCBI's automatic tags, which `scripts/asn_spec_gen.py` checks);
//!   an unknown or repeated member, an unknown variant and a missing mandatory member are NCBI's
//!   `CSerialException`; a member with a set flag (every member not held by `CRef<>`) adds
//!   `error while reading <member>` to a `CSerialException` that leaves its value
//!   (`ReadWithSetFlagMember`);
//! - NCBI's pointer forms where NCBI reads a `CRef<>` that is empty: an OPTIONAL member held by
//!   `CRef<>` (`REF`), each element of a list of `CRef<>`s, and every such member in a skip:
//!   `NULL` leaves it unset, `@n` is NCBI's error (no object is collected), `:Type value` naming
//!   the declared type is the value, in a frame of its own; a value of another class is NCBI's
//!   `CSerialException` where the program reports any with its own message (psiblast and
//!   tblastn), and refused where what NCBI does with it depends on the types its process has met
//!   before (`TypedReader::other_type`); an element of a list stated as `NULL` is held as the
//!   empty reference NCBI's list holds (`Value::Null`);
//! - INTEGER by the C++ storage of its member (`ReadInt4`, `ReadUint4`, `ReadInt8`, `ReadInt2`);
//!   ENUMERATED by name (text) or by a number that has a name, INTEGER with named values by name or
//!   any number;
//! - a VisibleString byte outside 0x20..0x7E as `#`, with NCBI's post 803.7
//!   (`Post::BadChar`), a UTF8String as it is; a StringStore as a VisibleString in text and as it
//!   is in binary; a UTF8String tag for a VisibleString in binary with NCBI's warning
//!   (`Post::Utf8Tag`, once per process);
//! - the delayed member (`Auth-list.names`) with the skip functions, kept as its bytes
//!   (`Value::Delayed`);
//! - the hooks of NCBI's hand-written classes that change what is read (`super::hooks`:
//!   `CSeq_descr::PostRead`, `CVariation_ref::PostRead`, the reserve hooks), and the Seq-id into
//!   LOSAT's Seq-id model (`super::seq_id`).
//!
//! Its errors are NCBI's classes (`ReadError`): a `CSerialException` (and an end of the input
//! inside an object, which NCBI's stack turns into one, `HandleEOF`), with the frames it leaves
//! (`AddFrameInfo`); a `CUtilException` of a number or a `CStringException` of a tag (text), with
//! the line NCBI's stack adds for each frame it leaves (`NCBI_RETHROW_SAME`, the report of
//! `CATCH_ALL`); an end of the input directly in the root object, a `CEofException` that NCBI's
//! stack passes on; the `std::exception`s of a reserve; and LOSAT's refusals of what its readers do
//! not take.

use super::read_xml::{XmlInput, XmlStyle, MISSING_VALUE};
use super::Encoding;
use super::{resolved, EnumValue, Kind, Member, Storage, Value, DELAYED};
use crate::ncbi_exception::{NcbiException, Severity, ThrowSite};
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{
    throw_error, utf8_tag_message, EnumItems, Fail, Input, Items, Names, StringKind, StringPost,
    TextInput, Utf8TagOnce, R,
};
use crate::report::serial::asn_text_reader::{AsnSource, PointerType, SliceSource};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::report::serial::object_frame::{
    frame_info, frame_name, FrameKind, CLASS_NOT_FOUND, EXPECTED_MEMBER, GET_REGISTERED_OBJECT,
    HANDLE_EOF, READ_INT2, READ_POINTER_FRAME_END, READ_POINTER_INCOMPATIBLE, SKIP_FILE_HEADER,
    SKIP_POINTER_FRAME_END,
};

/// The text of the byte block's post ([`Post::UnendedByteBlock`]).
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1554-1565
const UNENDED_BYTE_BLOCK: &str = "unended byte block";

/// A message that NCBI's reader posts while it reads, in the order posted.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Post {
    /// `ReplaceVisibleChar`'s post (803.7) of a byte read as `#`.
    BadChar(ReplacedChar),
    /// The error that a class's `PostRead` posts with its text (`CVariation_ref::PostRead`).
    Error(&'static str),
    /// The warning (803.10, once per run: `asn_input::Utf8TagOnce`) of a UTF8String tag for the
    /// VisibleString member at the stack trace (binary; `asn_input::utf8_tag_message`).
    Utf8Tag(String),
    /// The error post 803.4 of a byte block that a `std::exception` left with the stream in a
    /// good state (`TypedReader::octets`).
    UnendedByteBlock,
}

impl Post {
    /// The line that an application of NCBI posts on its standard error, with its diagnostic
    /// prefix `program` (the program's name in lower case, if it has one): the one owner of each
    /// post's text (`ReplacedChar::post`, `ncbi_exception::error_post_line`).
    pub(crate) fn line(&self, program: Option<&str>) -> Vec<u8> {
        match self {
            Post::BadChar(replaced) => replaced.post(program),
            Post::Error(text) => {
                crate::ncbi_exception::error_post_line(None, program, text.as_bytes())
            }
            Post::Utf8Tag(trace) => crate::ncbi_exception::post_line(
                Severity::Warning,
                Some("803.10"),
                program,
                utf8_tag_message(trace).as_bytes(),
            ),
            Post::UnendedByteBlock => crate::ncbi_exception::post_line(
                Severity::Error,
                Some("803.4"),
                program,
                UNENDED_BYTE_BLOCK.as_bytes(),
            ),
        }
    }

    /// The post's severity: a warning for the UTF8String tag, an error for the others.
    pub(crate) fn severity(&self) -> Severity {
        match self {
            Post::Utf8Tag(_) => Severity::Warning,
            Post::BadChar(_) | Post::Error(_) | Post::UnendedByteBlock => Severity::Error,
        }
    }
}

/// Why NCBI's reader stops: the class of what it throws.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ReadError {
    /// A `CSerialException`, with the frames it left (`AddFrameInfo`). NCBI's callers of the
    /// PSSM reader print their own line for it.
    Serial(NcbiException),
    /// Another `CException` that leaves the reader (`CUtilException`, `CStringException`), with
    /// the lines of NCBI's object stack.
    Exception(NcbiException),
    /// A `CEofException` that leaves the reader: the input ends directly in the root object.
    Eof,
    /// A `std::exception` that leaves the reader (its `what()`): `std::length_error` of a reserve
    /// hook (`hooks::pre_read_member`) or of the reserve of a string or byte block.
    Std(&'static str),
    /// `std::bad_alloc` of a string or byte block whose length passes every allocation
    /// (`asn_input::reserve`).
    OutOfMemory,
    /// What LOSAT's readers do not take (LOSAT's own limit, not an NCBI error): the caller refuses
    /// the input with this text.
    Unsupported(String),
}

impl From<Fail> for ReadError {
    fn from(fail: Fail) -> Self {
        match fail {
            Fail::Eof => ReadError::Eof,
            Fail::Serial(e) => ReadError::Serial(*e),
            Fail::Exception(e) => ReadError::Exception(*e),
            Fail::OutOfMemory => ReadError::OutOfMemory,
            Fail::Std(what) => ReadError::Std(what),
            Fail::Unsupported(what) => ReadError::Unsupported(what.0),
        }
    }
}

/// What the reader read, and the messages it posted before it stopped or finished.
pub(crate) type Read<T> = (std::result::Result<T, ReadError>, Vec<Post>);

impl Items for [Member] {
    fn count(&self) -> usize {
        self.len()
    }

    fn name(&self, index: usize) -> &'static str {
        self[index].name()
    }

    #[inline(always)]
    fn is_named(&self, index: usize, id: &[u8]) -> bool {
        self[index].name.is(id)
    }
}

/// The named values of an ENUMERATED (`integer` false) or an INTEGER with named values.
struct EnumType {
    integer: bool,
    values: &'static [EnumValue],
}

impl EnumItems for EnumType {
    fn integer(&self) -> bool {
        self.integer
    }

    fn value_of(&self, name: &[u8]) -> Option<i32> {
        self.values
            .iter()
            .find(|value| value.is_named(name))
            .map(|value| value.value)
    }

    fn has_value(&self, value: i32) -> bool {
        self.values.iter().any(|named| named.value == value)
    }

    fn type_name(&self) -> &'static str {
        super::read_xml::enum_type_name(self.values)
    }
}

/// The type of a frame that names none (a member, variant or element frame, and an array frame,
/// whose container type is anonymous: `Frame::type_name`).
pub(super) const NO_TYPE: u16 = u16::MAX;

/// A frame of NCBI's object stack (`CObjectStackFrame`), as the reader keeps it on its stack.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Frame {
    kind: FrameKind,
    /// Whether the frame belongs to a skip function (`SkipClassSequential`, ...).
    skip: bool,
    /// The named frame of `ReadPointer` (or `SkipPointer`) around a value of the type that
    /// `:Type` names, whose `END_OBJECT_FRAME` is that function's.
    pointer: bool,
    /// The type (the table's index) of a named, class or choice frame: the name of its
    /// `Object type= `, and the object whose `HasTypeInfo` NCBI's `GetParentObjectPtr` asks
    /// (`hooks::pre_read_variant`). [`NO_TYPE`] for the other frames.
    ty: u16,
    /// The member or variant (`SetTopMemberId`) of a member or variant frame.
    member: Option<&'static Member>,
    /// A `Seq-table` class frame's `num-rows` once read (the object's `GetNum_rows()`, an `int`:
    /// the member's Int4 storage).
    num_rows: Option<i32>,
}

impl Frame {
    /// Whether the frame is a class or a choice frame (`eFrameClass`, `eFrameChoice`).
    pub(super) fn is_class_or_choice(&self) -> bool {
        matches!(self.kind, FrameKind::Class | FrameKind::Choice)
    }

    /// The type of a class or choice frame.
    pub(super) fn ty(&self) -> Option<u16> {
        match self.kind {
            FrameKind::Class | FrameKind::Choice => Some(self.ty),
            _ => None,
        }
    }

    /// A `Seq-table` class frame's `num-rows`, once read.
    pub(super) fn num_rows(&self) -> Option<i32> {
        self.num_rows
    }

    /// The type the frame names (`Object type= `): the named, class and choice frames' type, the
    /// anonymous container type of an array frame (an implicit class such as `Seq-descr` has no
    /// frame of its own, `ReadImplicitMember`, classinfo.cpp:289-306: `Object type= ` names
    /// nothing, c07 of SURVEY_D3).
    fn type_name(&self) -> Option<&'static str> {
        match self.kind {
            FrameKind::Array => Some(""),
            FrameKind::Named | FrameKind::NamedType | FrameKind::Class | FrameKind::Choice => {
                Some(super::type_info(self.ty).name())
            }
            FrameKind::ClassMember | FrameKind::ChoiceVariant | FrameKind::ArrayElement => None,
        }
    }

    /// The member's or variant's name of a member or variant frame.
    fn member_name(&self) -> Option<&'static str> {
        self.member.map(Member::name)
    }

    /// `GetFrameInfo()` (`object_frame::frame_info`).
    fn info(&self) -> String {
        frame_info(self.kind, self.type_name(), self.member_name())
    }

    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves this frame in
    /// the stream `I` (`Input::frame_site`: `object_frame::text_frame_site` for the text and
    /// binary streams; `READ_POINTER_FRAME_END`, `SKIP_POINTER_FRAME_END`). The binary stream's
    /// primitives throw no `CException` but a `CSerialException` (and the end of the input), so
    /// only the text and XML streams' frames are ever named by their lines.
    fn site<I: Input>(&self) -> (&'static str, u32, &'static str) {
        match (self.pointer, self.skip) {
            (true, false) => READ_POINTER_FRAME_END,
            (true, true) => SKIP_POINTER_FRAME_END,
            (false, skip) => I::frame_site(self.kind, skip),
        }
    }
}

/// Reads the ASN.1 text document at the start of `text` as the root type `root` (an object of
/// `CObjectIStream::Read`: its file header, then the value; what follows it is not read), with
/// the messages posted.
pub(crate) fn read_document(text: &[u8], root: u16, names: Names) -> Read<Value> {
    ObjectStream::text(SliceSource::new(text), names).read_object(root)
}

/// [`read_document`] of the ASN.1 binary object at the start of `bytes` (the stream read as `MSerial_AsnBinary`
/// into the object: no file header).
pub(crate) fn read_ber_document(bytes: &[u8], root: u16, names: Names) -> Read<Value> {
    ObjectStream::binary(SliceSource::new(bytes), names).read_object(root)
}

/// [`read_document`] of the ASN.1 binary value of the type `ty` at the start of `bytes` (no file
/// header), as the tests make the value of a delayed member read from binary.
#[cfg(test)]
pub(super) fn read_ber_value(bytes: &[u8], ty: u16) -> Read<Value> {
    ObjectStream::binary(SliceSource::new(bytes), Names::DELAYED).read_object(ty)
}

// NCBI reference (598d8ae6): c++/src/serial/delaybuf.cpp:63-84
// NCBI reference (598d8ae6): c++/src/serial/member.cpp:489-499
/// The value of the delayed member `member` of the SEQUENCE `class` from the bytes that the
/// reader kept of it in `encoding` (`Value::Delayed`), as NCBI's `CDelayBuffer::Update` reads
/// them before a stream of another format writes the member (`write`): a new stream of that
/// format over the bytes (its line or byte count starts at them), under a class frame of `class`
/// and a member frame of `member` (`UpdateDelayedBuffer`), with what it posts.
pub(super) fn read_delayed(
    encoding: Encoding,
    bytes: &[u8],
    class: u16,
    member: &'static Member,
) -> Read<Value> {
    let src = SliceSource::new(bytes);
    match encoding {
        Encoding::Text => update(TextInput::new(src, Names::DELAYED), class, member),
        Encoding::Binary => update(BerReader::new(src, Names::DELAYED), class, member),
        Encoding::Xml => update(
            XmlInput::new(src, XmlStyle::Asn, Names::DELAYED),
            class,
            member,
        ),
    }
}

/// [`read_delayed`] of the stream `inp`.
fn update<I: Input>(inp: I, class: u16, member: &'static Member) -> Read<Value> {
    let mut reader = TypedReader::new(inp, Names::DELAYED);
    // The read that kept the bytes met their UTF8String tags (and posted for the first).
    reader.utf8_once = Utf8TagOnce::POSTED;
    reader.start_object();
    let value = reader.delayed_member(class, member);
    reader.inp.drain_posts(&mut reader.posts);
    (value.map_err(ReadError::from), reader.posts)
}

// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:122-145
/// The deepest nesting of object frames that a typed read follows (`ObjectStream::read_typed`):
/// a typed reader recurses on LOSAT's stack as the type nests (`super::typed`; a `Seq-loc` of
/// `mix` and `equiv` to any depth), so a frame past it is LOSAT's refusal
/// (`TypedReader::too_deep`). NCBI's object stream has no such limit (its frames are bounded by
/// memory: `PushFrameLong` doubles the stack), and the generic `Value` read, which keeps its
/// values on the heap, has none.
const TYPED_MAX_DEPTH: usize = 400;

/// Objects one after another from one input (`ObjectStream::read_object`), each read by a new
/// state of NCBI's stream at the byte after the last one the object before it took (NCBI gives the
/// bytes its stream did not use back when an object ends, `CIStreamBuffer::Close`,
/// strbuffer.cpp:125-139): its line count or its byte count and data limits start again, and its
/// frames and posts are its own.
pub(crate) struct ObjectStream<I: Input> {
    reader: TypedReader<I>,
    /// A read failed: NCBI's stream has its fail flags set, and its `EndOfData` is true from
    /// then on.
    failed: bool,
}

impl<S: AsnSource> ObjectStream<TextInput<S>> {
    /// Objects in ASN.1 text from `src`, which LOSAT's messages call by `names`.
    pub(crate) fn text(src: S, names: Names) -> Self {
        ObjectStream {
            reader: TypedReader::new(TextInput::new(src, names), names),
            failed: false,
        }
    }
}

impl<S: AsnSource> ObjectStream<BerReader<S>> {
    /// Objects in ASN.1 binary from `src`, which LOSAT's messages call by `names`.
    pub(crate) fn binary(src: S, names: Names) -> Self {
        ObjectStream {
            reader: TypedReader::new(BerReader::new(src, names), names),
            failed: false,
        }
    }
}

impl<S: AsnSource> ObjectStream<XmlInput<S>> {
    /// Objects in NCBI's XML of the `style` from `src`, which LOSAT's messages call by `names`.
    pub(crate) fn xml(src: S, style: XmlStyle, names: Names) -> Self {
        ObjectStream {
            reader: TypedReader::new(XmlInput::new(src, style, names), names),
            failed: false,
        }
    }
}

impl<I: Input> ObjectStream<I> {
    /// The next object, as the root type `root`, with the messages posted while it was read.
    pub(crate) fn read_object(&mut self, root: u16) -> Read<Value> {
        let reader = &mut self.reader;
        reader.start_object();
        let value = reader.document(root, |r| r.read_value(root, false));
        reader.inp.drain_posts(&mut reader.posts);
        self.failed |= value.is_err();
        (
            value.map_err(ReadError::from),
            std::mem::take(&mut reader.posts),
        )
    }

    /// The next object, as the root type `root`, read by a typed reader (`read`, on the reader's
    /// frames: `super::typed`), with the messages posted while it was read.
    pub(crate) fn read_typed<T>(
        &mut self,
        root: u16,
        read: impl FnOnce(&mut TypedReader<I>) -> R<T>,
    ) -> Read<T> {
        let reader = &mut self.reader;
        reader.start_object();
        reader.max_depth = TYPED_MAX_DEPTH;
        let value = reader.document(root, read).map_err(ReadError::from);
        reader.max_depth = usize::MAX;
        reader.inp.drain_posts(&mut reader.posts);
        self.failed |= value.is_err();
        (value, std::mem::take(&mut reader.posts))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:588-597
    /// `EndOfData()`: whether no object follows (`Input::end_of_data`), which a reader of objects
    /// one after another asks before each (a binary object that is not a SEQUENCE ends inside its
    /// frames, where the end of the input is NCBI's `CSerialException`). True after a failed
    /// read, whose fail flags NCBI's stream keeps.
    #[cfg(test)]
    pub(crate) fn end_of_data(&mut self) -> bool {
        self.failed || self.reader.inp.end_of_data()
    }

    /// The primitives of the stream after the objects read (`stream`: the text after them goes
    /// back to the caller).
    pub(crate) fn into_input(self) -> I {
        self.reader.inp
    }

    /// The primitives of the stream (where it stands: `Input::position`, the bytes it holds).
    pub(crate) fn input(&self) -> &I {
        &self.reader.inp
    }

    /// The primitives of the stream, to change how its source is read (tests).
    #[cfg(test)]
    pub(crate) fn input_mut(&mut self) -> &mut I {
        &mut self.reader.inp
    }
}

/// The object stack of one `CObjectIStream` over the primitives `I`: its frames, what it posted,
/// and the values the generic reader is reading.
pub(crate) struct TypedReader<I: Input> {
    pub(super) inp: I,
    frames: Vec<Frame>,
    posts: Vec<Post>,
    /// The values being read, outermost first (`Open`): NCBI's recursion of `ReadObject` as a list
    /// on the heap, so that a value nested to any depth is read in LOSAT's own stack.
    open: Vec<Open>,
    /// Whether the warning of a UTF8String tag was posted (`Utf8TagOnce`; one stream is all a
    /// run of LOSAT reads with it).
    utf8_once: Utf8TagOnce,
    /// What LOSAT's messages call the document, the reader and NCBI's part that allocates.
    names: Names,
    /// The deepest nesting of frames the reader follows (`TYPED_MAX_DEPTH` in a typed read).
    max_depth: usize,
    /// Whether a value of another class than the declared one is being read (`other_type`).
    reading_other: bool,
}

/// A value being read that waits for a value inside it (`TypedReader::read_value`).
enum Open {
    /// The end of `ReadObject` of the type `target`: the hooks of its class (not after a skip).
    Value { target: u16, skip: bool },
    /// A SEQUENCE reading the member at `reading` (`members`), whose value's frames are those
    /// above `depth`.
    Sequence {
        ty: u16,
        members: &'static [Member],
        skip: bool,
        values: Vec<(u16, Value)>,
        reading: usize,
        depth: usize,
    },
    /// A CHOICE reading its variant.
    Choice { index: u16 },
    /// A SEQUENCE OF or SET OF reading an element.
    List {
        element: u16,
        refs: bool,
        skip: bool,
        values: Vec<Value>,
    },
    /// A delayed member (`delayed`) from `start`.
    Delayed { start: usize },
    /// The value of `:Type` in `ReadPointer`'s named frame: the frame's end and
    /// `ReadOtherPointerEnd` follow it.
    Pointer,
    /// The value in the named frames of `count` named types (`ReadNamedType`, the XML stream's:
    /// `Input::NAMED_TYPES`), whose ends follow it.
    Named { count: u8 },
}

/// What the reader does next: read a value of a type (`skip`), or give a value read to the one
/// that holds it.
enum Step {
    Read(u16, bool),
    Done(Value),
}

/// Where `TypedReader::begin_pointer` found the value of an empty `CRef<>`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Pointee {
    /// `NULL`: there is none.
    Null,
    /// `eThisPointer`: it follows.
    This,
    /// `eOtherPointer` naming the declared type: it follows inside `ReadPointer`'s named frame.
    Other,
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:712-756
// (`eNullValue` and `eMissingValue` are not thrown for these types: `ExpectedMember` throws
// `eFormatError` with the default `eSerialVerifyData_Yes`. A member without a set flag, one held
// by `CRef<>`, is read by `ReadSimpleMember`, a delayed one by `ReadLongMember`, a skipped one by
// `SkipSimpleMember`: none of them adds a line.)
const READ_WITH_SET_FLAG_MEMBER: (&str, u32, &str) = (
    "src/serial/member.cpp",
    752,
    "ncbi::CMemberInfoFunctions::ReadWithSetFlagMember",
);

/// `ReadWithSetFlagMember`'s rethrow of a missing value (`eMissingValue`, the XML stream's
/// ENUMERATED without its attribute) of a mandatory member.
const MISSING_VALUE_MEMBER: ThrowSite = ThrowSite {
    file: "src/serial/member.cpp",
    line: 748,
    kind: "CSerialException::eFormatError",
    function: Some("ncbi::CMemberInfoFunctions::ReadWithSetFlagMember"),
};

/// What reading the value of the SEQUENCE member `member` makes of `fail` that leaves it: a
/// member with a set flag (one not held by `CRef<>`, `Member::is_ref`, and not delayed) is read
/// in `ReadWithSetFlagMember`, which rethrows a `CSerialException` with `error while reading
/// <member>` (a missing value of a mandatory member: `missing value <member>`,
/// `TypedReader::missing_value`); the others pass `fail` on. The one owner of that rule and its
/// rethrow, for the generic read (`TypedReader::unwind_values`) and the typed reads (`typed`); a
/// skip does not call it.
#[cold]
#[inline(never)]
pub(super) fn set_flag_error(fail: Fail, member: &Member) -> Fail {
    match fail {
        Fail::Serial(e) if !member.is_ref() && member.flags & DELAYED == 0 => {
            if e.kind() == MISSING_VALUE {
                return Fail::serial(e.rethrow(
                    MISSING_VALUE_MEMBER,
                    format!("missing value {}", member.name()),
                ));
            }
            let (file, line, function) = READ_WITH_SET_FLAG_MEMBER;
            let message = format!("error while reading {}", member.name());
            Fail::serial(e.rethrow_same(file, line, function, message))
        }
        other => other,
    }
}

impl<I: Input> TypedReader<I> {
    fn new(inp: I, names: Names) -> Self {
        TypedReader {
            inp,
            frames: Vec::new(),
            posts: Vec::new(),
            open: Vec::new(),
            utf8_once: Utf8TagOnce::default(),
            names,
            max_depth: usize::MAX,
            reading_other: false,
        }
    }

    /// The state of a new `CObjectIStream` at the next byte of the input: no frames, no values
    /// open (`Input::new_object`).
    fn start_object(&mut self) {
        self.inp.new_object();
        self.frames.clear();
        self.open.clear();
    }

    /// The refusal of a frame past `max_depth`.
    #[cold]
    #[inline(never)]
    fn too_deep(&self) -> Fail {
        Fail::unsupported(format!(
            "objects nested more than {} levels deep",
            self.max_depth
        ))
    }

    /// `BEGIN_OBJECT_FRAME` of a frame of `kind` (of the type `ty`, or [`NO_TYPE`]).
    #[inline(always)]
    pub(super) fn push(&mut self, kind: FrameKind, ty: u16, skip: bool) -> R<()> {
        if self.frames.len() >= self.max_depth {
            return Err(self.too_deep());
        }
        self.frames.push(Frame {
            kind,
            skip,
            pointer: false,
            ty,
            member: None,
            num_rows: None,
        });
        self.inp.frame_pushed(kind, ty, self.frames.len());
        Ok(())
    }

    /// `END_OBJECT_FRAME` of the top frame after what it holds was read.
    #[inline(always)]
    fn pop_frame(&mut self) -> Frame {
        let frame = self.frames.pop().expect("a frame to leave");
        self.inp.frame_popped(self.frames.len());
        frame
    }

    /// `END_OBJECT_FRAME` of the top frame without an error: the frame is dropped unread.
    #[inline(always)]
    fn drop_frame(&mut self) {
        let left = self.frames.pop();
        debug_assert!(left.is_some(), "a frame to leave");
        self.inp.frame_popped(self.frames.len());
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:189-214 (cited at
    // `report::serial::object_frame::text_frame_site`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608 (`HandleEOF`, cited at
    // `object_frame::HANDLE_EOF`)
    /// `END_OBJECT_FRAME`: the frame is left, and an error that leaves it gets NCBI's handling. A
    /// `CSerialException` takes the frame's name (`AddFrameInfo` of `GetFrameName()`); a
    /// `CEofException` is `HandleEOF`, which passes it on when fewer than 2 frames stay and throws
    /// its own `CSerialException` (with the frame's information, and without the frame's name)
    /// otherwise; another `CException` takes the frame's line (`NCBI_RETHROW_SAME` with
    /// `GetFrameInfo()`).
    #[inline(always)]
    pub(super) fn pop<T>(&mut self, mut result: R<T>) -> R<T> {
        match &mut result {
            Ok(_) => self.drop_frame(),
            Err(fail) => self.leave_top(fail),
        }
        result
    }

    /// `pop`'s error: the top frame is left by `fail`, which takes what the frame adds.
    #[cold]
    #[inline(never)]
    fn leave_top(&mut self, fail: &mut Fail) {
        let frame = self.pop_frame();
        let left = std::mem::replace(fail, Fail::Eof);
        *fail = self.leave(frame, left);
    }

    /// What leaving `frame` (the stack is below it now) adds to `fail`.
    #[cold]
    #[inline(never)]
    fn leave(&self, frame: Frame, fail: Fail) -> Fail {
        match fail {
            Fail::Serial(mut e) => {
                e.add_frame_info(frame_name(
                    frame.kind,
                    frame.type_name(),
                    frame.member_name(),
                ));
                Fail::Serial(e)
            }
            Fail::Eof if self.frames.len() >= 2 => Fail::serial(NcbiException::new(
                HANDLE_EOF,
                format!("{}: {}", self.inp.position(), frame.info()),
            )),
            Fail::Exception(e) => {
                let (file, line, function) = frame.site::<I>();
                Fail::exception(e.rethrow_same(file, line, function, frame.info()))
            }
            other => other,
        }
    }

    /// The frames above `depth` left by `fail`, from the top down.
    #[cold]
    #[inline(never)]
    fn unwind(&mut self, depth: usize, mut fail: Fail) -> Fail {
        while self.frames.len() > depth {
            let frame = self.pop_frame();
            fail = self.leave(frame, fail);
        }
        fail
    }

    /// `SetTopMemberId`: the member or variant of the top frame.
    #[inline(always)]
    pub(super) fn set_member(&mut self, member: &'static Member) {
        if let Some(top) = self.frames.last_mut() {
            top.member = Some(member);
        }
        self.inp.frame_member(member);
    }

    /// `GetPosition()` of the stream.
    pub(crate) fn position(&self) -> StreamPosition {
        self.inp.position()
    }

    /// `GetStackTrace()`: the path of the object being read (`PssmWithParameters.pssm.query.seq.
    /// descr.E.comment`), by the one writer of NCBI's object stack path.
    pub(crate) fn stack_trace(&self) -> String {
        object_stack_trace(
            self.frames
                .first()
                .and_then(Frame::type_name)
                .unwrap_or("?"),
            self.frames
                .iter()
                .skip(1)
                .filter_map(|frame| match frame.kind {
                    FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                        frame.member_name().map(StackStep::Member)
                    }
                    FrameKind::ArrayElement => Some(StackStep::Element),
                    _ => None,
                }),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:882-901
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:959-969
    /// The delayed member `member` of the SEQUENCE `class` in its frames (`read_delayed`).
    fn delayed_member(&mut self, class: u16, member: &'static Member) -> R<Value> {
        self.push(FrameKind::Class, class, false)?;
        let value = self.delayed_member_frame(member);
        self.pop(value)
    }

    fn delayed_member_frame(&mut self, member: &'static Member) -> R<Value> {
        self.push(FrameKind::ClassMember, NO_TYPE, false)?;
        self.set_member(member);
        let value = self.read_value(member.ty, false);
        self.pop(value)
    }

    /// `stream >> object` of the root type `root`: `SkipFileHeader`, then `Read` (here `read`),
    /// each in a named frame of its own.
    fn document<T>(&mut self, root: u16, read: impl FnOnce(&mut Self) -> R<T>) -> R<T> {
        self.push(FrameKind::Named, root, false)?;
        let header = self.file_header(super::type_info(root).name());
        self.pop(header)?;
        self.push(FrameKind::Named, root, false)?;
        let value = read(self);
        self.pop(value)
    }

    /// `SkipFileHeader` of the type: the type name of a text document must be the root's.
    fn file_header(&mut self, name: &'static str) -> R<()> {
        match self.inp.read_file_header()? {
            Some(id) if !id.is_empty() && id != name.as_bytes() => Err(throw_error(
                SKIP_FILE_HEADER,
                self.inp.position(),
                [
                    b"incompatible type ".as_slice(),
                    &id,
                    b"<>",
                    name.as_bytes(),
                ]
                .concat(),
            )),
            _ => Ok(()),
        }
    }

    /// `BeginClassMember` of the members `members` from `pos` on (`ReadClassSequential`'s loop):
    /// the next member, named in the member frame (`SetTopMemberId`), with the mandatory members
    /// before it checked (`ReadMissingMember`); `None` at the end of the class.
    #[inline(always)]
    pub(super) fn member_begin(
        &mut self,
        members: &'static [Member],
        pos: usize,
    ) -> R<Option<usize>> {
        let Some(index) = self.inp.begin_member(members, pos)? else {
            return Ok(None);
        };
        self.set_member(&members[index]);
        for missing in &members[pos..index] {
            self.expected_member(missing)?;
        }
        Ok(Some(index))
    }

    /// `ReadMissingSimpleMember` -> `ExpectedMember` of a mandatory member: "member X expected"
    /// (`CSerialException`, NCBI's default verification).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
    #[inline(always)]
    pub(super) fn expected_member(&self, member: &Member) -> R<()> {
        if member.optional() {
            return Ok(());
        }
        Err(self.member_expected(member))
    }

    #[cold]
    #[inline(never)]
    fn member_expected(&self, member: &Member) -> Fail {
        throw_error(
            EXPECTED_MEMBER,
            self.inp.position(),
            format!("member {} expected", member.name()),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/ptrinfo.cpp:210-227
    // NCBI reference (598d8ae6): c++/src/serial/ptrinfo.cpp:249-256
    /// Where the value of `member` is: a member held by `CRef<>` that is empty when NCBI reads it
    /// (OPTIONAL, or any one in a skip: `CPointerTypeInfo::ReadPointer`, `SkipPointer`) is read with
    /// the pointer forms (`begin_pointer`); the generated constructor allocates a mandatory one,
    /// read in place, as every other member is.
    #[inline(always)]
    pub(super) fn member_pointer(&mut self, member: &Member, skip: bool) -> R<Pointee> {
        if !member.is_ref() || !(skip || member.optional()) {
            return Ok(Pointee::This);
        }
        self.begin_pointer(member.ty, skip)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1133-1213
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1234-1272
    // NCBI reference (598d8ae6): c++/src/serial/classinfob.cpp:244-252
    // (`MapType` is `GetClassInfoByName`. The classes it knows are those whose type information
    // the process has built so far: the declared type is one of them, another type may or may not
    // be (makeblastdb 2.17.0: "class not found: Seq-interval" for `:Seq-interval` as the first
    // mask; a type it knows is read whole and then fails with "incompatible member type", or with
    // the first error inside it when the value is not of that type; psiblast 2.17.0 reads
    // `:PssmParameters` where params belongs and stops at `:Pssm` and `:Object-id`, probe
    // D3c-3b).)
    /// `ReadPointer` (`SkipPointer` in a skip) of an empty `CRef<>` to a value of the type
    /// `declared`: `NULL` is none; `@n` is NCBI's error (makeblastdb's and psiblast's streams
    /// collect no object); `:Type` naming `declared` is the value in a named frame of that type,
    /// pushed here (`end_pointer` leaves it); a value of another type is NCBI's error where what
    /// its process has read before does not show, and refused elsewhere (`other_type`).
    #[inline(always)]
    pub(super) fn begin_pointer(&mut self, declared: u16, skip: bool) -> R<Pointee> {
        match self.inp.read_pointer_type()? {
            PointerType::Null => Ok(Pointee::Null),
            PointerType::This => Ok(Pointee::This),
            other => self.begin_other_pointer(other, declared, skip),
        }
    }

    /// [`begin_pointer`](Self::begin_pointer) of an object index (`@n`) or a value of a named type
    /// (`:Type`).
    #[cold]
    #[inline(never)]
    fn begin_other_pointer(&mut self, kind: PointerType, declared: u16, skip: bool) -> R<Pointee> {
        match kind {
            PointerType::Null => Ok(Pointee::Null),
            PointerType::This => Ok(Pointee::This),
            PointerType::Object => {
                self.inp.read_object_pointer()?;
                Err(throw_error(
                    GET_REGISTERED_OBJECT,
                    self.inp.position(),
                    "invalid object index: NO_COLLECT defined",
                ))
            }
            PointerType::Other => {
                let name = self.inp.read_other_pointer()?;
                let declared_name = super::type_info(declared).name();
                if name != declared_name.as_bytes() {
                    return Err(self.other_type(&name, declared_name, skip));
                }
                self.push(FrameKind::Named, declared, skip)?;
                if let Some(top) = self.frames.last_mut() {
                    top.pointer = true;
                }
                Ok(Pointee::Other)
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1173-1219
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1341-1361
    /// `:name`, a value of the class `name`, where a `declared` belongs. NCBI's `MapType` finds
    /// the class when its process has built its type information, and throws "class not found"
    /// otherwise; found, the value is read as that class, and the read throws "incompatible member
    /// type" after it, or its first error inside it. Each is a `CSerialException`, which
    /// psiblast's and tblastn's PSSM reader report with a message of their own
    /// (`Names::own_serial_message`; NCBI 2.17.0 rc 1 for `:Object-id`, `:Pssm`, `:No-such-type`
    /// where params belongs and `:Bioseq` where the query belongs, text and binary, probe
    /// D3c-3b): there the value is read as the table's class of that name, if any, and the read
    /// ends with the error. A class outside the table is one that the PSSM's read never builds,
    /// and nothing before the read builds it: "class not found" (NCBI 2.17.0 with a number that is
    /// not one inside `:Blast4-cutoff`, `:Seq-submit`, `:Hit`; probe D3c-3b probe9). What NCBI
    /// does only when it knows the class (a message posted while the value is read, an error of
    /// another kind inside it: `:Object-id id x` is its exception report, exit status 255, in
    /// probe9), a value of another class inside it, a skip (`SkipPointer` skips a class it knows
    /// and stops at another) and the other streams (makeblastdb prints the exception's own text)
    /// are LOSAT's refusal.
    #[cold]
    #[inline(never)]
    fn other_type(&mut self, name: &[u8], declared: &str, skip: bool) -> Fail {
        if skip || !self.names.own_serial_message || self.reading_other {
            return self.other_type_refusal(name, declared);
        }
        let Some(class) = super::class_named(name) else {
            return Fail::serial(NcbiException::new(
                CLASS_NOT_FOUND,
                format!("class not found: {}", String::from_utf8_lossy(name)),
            ));
        };
        let posts = self.posts.len();
        self.reading_other = true;
        let read = self.read_other(class);
        self.reading_other = false;
        match read {
            _ if self.posts.len() > posts => {
                // NCBI's posts there depend on its registry: the refusal is LOSAT's alone.
                self.posts.truncate(posts);
                self.other_type_refusal(name, declared)
            }
            Ok(()) => throw_error(
                READ_POINTER_INCOMPATIBLE,
                self.inp.position(),
                "incompatible member type",
            ),
            Err(fail @ Fail::Serial(_)) => fail,
            Err(_) => self.other_type_refusal(name, declared),
        }
    }

    /// The value of the class `class` after `:name` (`other_type`): read in a named frame of the
    /// class, as `ReadPointer` reads it, to `ReadOtherPointerEnd`.
    fn read_other(&mut self, class: u16) -> R<()> {
        self.push(FrameKind::Named, class, false)?;
        if let Some(top) = self.frames.last_mut() {
            top.pointer = true;
        }
        let value = self.read_value(class, false);
        self.end_pointer(value).map(drop)
    }

    /// The refusal of `:name` (a value of the type `name`) where a `declared` belongs.
    #[cold]
    #[inline(never)]
    fn other_type_refusal(&self, name: &[u8], declared: &str) -> Fail {
        let name = String::from_utf8_lossy(name);
        Fail::unsupported(format!(
            "':{name}', a value of the type {name} where a {declared} belongs ({}): whether \
             {} reads it depends on the types its process has met before",
            self.inp.position(),
            self.names.allocator
        ))
    }

    /// The end of a value that `begin_pointer` found at `:Type` (`Pointee::Other`), read with
    /// `result`: the `END_OBJECT_FRAME` of the named frame and `ReadOtherPointerEnd`.
    #[inline(always)]
    pub(super) fn end_pointer<T>(&mut self, result: R<T>) -> R<T> {
        let value = self.pop(result)?;
        self.inp.read_other_pointer_end()?;
        Ok(value)
    }

    /// One value of the type `ty` (`ReadObject`), or its skip (`SkipObject`, `skip`), with the hooks
    /// of NCBI's classes that a read runs. The values inside it are read in turn from the list of
    /// the values being read (`open`), as NCBI's functions recurse; an error leaves NCBI's frames
    /// from the innermost out (`unwind`), as it passes `END_OBJECT_FRAME` of each, and the
    /// `ReadWithSetFlagMember` of each member it leaves, which takes a missing value
    /// (`missing_value`).
    pub(crate) fn read_value(&mut self, ty: u16, skip: bool) -> R<Value> {
        let (open, frames) = (self.open.len(), self.frames.len());
        let mut step = self.start(ty, skip);
        loop {
            step = match step {
                Ok(Step::Read(ty, skip)) => self.start(ty, skip),
                Ok(Step::Done(value)) if self.open.len() == open => return Ok(value),
                Ok(Step::Done(value)) => {
                    let held = self.open.pop().expect("a value being read");
                    self.resume(held, value)
                }
                Err(error) if I::MISSING_VALUES && is_missing_value(&error) => {
                    match self.missing_value_member(open) {
                        Some(at) => self.missing_value(at),
                        None => return Err(self.unwind_values(open, frames, error)),
                    }
                }
                Err(error) => return Err(self.unwind_values(open, frames, error)),
            };
        }
    }

    /// The position in `open` (from `floor`) of the member that takes a missing value: the
    /// innermost SEQUENCE member being read with a set flag (`ReadWithSetFlagMember`: read, not
    /// skipped or delayed, not held by `CRef<>`; `SkipSimpleMember`, `ReadLongMember` and
    /// `ReadSimpleMember` let it pass), when it is OPTIONAL or DEFAULT. A mandatory one rethrows
    /// it (`set_flag_error`) as any error leaves it.
    #[cold]
    #[inline(never)]
    fn missing_value_member(&self, floor: usize) -> Option<usize> {
        let at = self.open[floor..].iter().rposition(|held| match held {
            Open::Sequence {
                members,
                skip: false,
                reading,
                ..
            } => {
                let member = &members[*reading];
                !member.is_ref() && member.flags & DELAYED == 0
            }
            _ => false,
        })?;
        let Open::Sequence {
            members, reading, ..
        } = &self.open[floor + at]
        else {
            unreachable!("a SEQUENCE member")
        };
        members[*reading].optional().then_some(floor + at)
    }

    /// `ReadWithSetFlagMember`'s missing value in the OPTIONAL or DEFAULT member being read at
    /// `at` (`missing_value_member`): the values and frames inside the member are left (an error
    /// that passes no `CException` handler: the frames add nothing), the member is not set, and
    /// the class reads on (`EndClassMember`, the next member). What the stream had read stays
    /// read: the ENUMERATED's tag ended as self-closed (`ReadEnum`), the tags around it open.
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:714-756
    #[cold]
    #[inline(never)]
    fn missing_value(&mut self, at: usize) -> R<Step> {
        self.open.truncate(at + 1);
        let Some(Open::Sequence {
            ty,
            members,
            skip,
            values,
            reading,
            depth,
        }) = self.open.pop()
        else {
            unreachable!("the member's SEQUENCE")
        };
        while self.frames.len() > depth {
            self.drop_frame();
        }
        self.inp.end_member()?;
        self.next_member(ty, members, skip, values, reading + 1)
    }

    /// The values above `open` and the frames above `frames` left by `fail`: each SEQUENCE member
    /// being read with a set flag adds its line once the frames of its value are left.
    #[cold]
    #[inline(never)]
    fn unwind_values(&mut self, open: usize, frames: usize, mut fail: Fail) -> Fail {
        while self.open.len() > open {
            if let Some(Open::Sequence {
                members,
                reading,
                skip,
                depth,
                ..
            }) = self.open.pop()
            {
                fail = self.unwind(depth, fail);
                let member = &members[reading];
                if !skip {
                    fail = set_flag_error(fail, member);
                }
            }
        }
        self.unwind(frames, fail)
    }

    /// The start of a value of the type `ty` (`ReadObject`'s type function): a value without
    /// values inside it is read here, a SEQUENCE, CHOICE or SEQUENCE OF opens (`Open`). A Seq-id
    /// that is read (not skipped) is read into LOSAT's Seq-id model (`super::seq_id`).
    fn start(&mut self, ty: u16, skip: bool) -> R<Step> {
        let info = resolved(ty);
        let target = alias_target(ty);
        if I::NAMED_TYPES {
            let count = self.begin_names(ty, skip)?;
            if count > 0 {
                self.open.push(Open::Named { count });
            }
        }
        self.open.push(Open::Value { target, skip });
        let value = match info.kind {
            Kind::Choice(_) if target == super::table::named::SEQ_ID && !skip => {
                Value::SeqId(Box::new(super::seq_id::read_seq_id(self)?))
            }
            Kind::Sequence(members) => {
                self.push(FrameKind::Class, target, skip)?;
                self.inp.begin_class()?;
                self.push(FrameKind::ClassMember, NO_TYPE, skip)?;
                return self.next_member(target, members.members(), skip, Vec::new(), 0);
            }
            Kind::Choice(variants) => return self.start_choice(target, variants.members(), skip),
            Kind::List { element, set, refs } => return self.start_list(element, set, refs, skip),
            Kind::Alias(_) => unreachable!("resolved"),
            Kind::Integer(storage) => self.integer_value(storage, skip)?,
            Kind::Enumerated { integer, values } => {
                Value::Enum(self.enum_number(values.values(), integer, skip)?)
            }
            // `SkipBool` reads what `ReadBool` reads (and `SkipNull` what `ReadNull` reads, but
            // in XML: `Input::skip_null`).
            Kind::Boolean => Value::Bool(self.inp.read_bool()?),
            Kind::Null => {
                if skip {
                    self.inp.skip_null()?;
                } else {
                    self.inp.read_null()?;
                }
                Value::Null
            }
            Kind::Real => self.real(skip)?,
            Kind::String(kind) => Value::Str(self.string(kind, skip)?),
            Kind::Octets => self.octets(skip)?,
        };
        Ok(Step::Done(value))
    }

    /// A value read is given to the value being read that holds it (`held`).
    fn resume(&mut self, held: Open, value: Value) -> R<Step> {
        match held {
            Open::Value { target, skip } => {
                if skip {
                    return Ok(Step::Done(value));
                }
                // What the stream posted comes before what the hook posts.
                self.inp.drain_posts(&mut self.posts);
                super::hooks::post_read(target, value, &mut self.posts).map(Step::Done)
            }
            Open::Sequence {
                ty,
                members,
                skip,
                mut values,
                reading,
                ..
            } => {
                // `EndClassMember`, in the member's frame.
                self.inp.end_member()?;
                if !skip {
                    // The `num-rows` of a Seq-table, which the reserve hook of its columns' data
                    // reads from the object being read (`hooks::pre_read_variant`).
                    if ty == super::table::named::SEQ_TABLE && members[reading].name() == "num-rows"
                    {
                        let at = self.frames.len() - 2;
                        if let Value::Int(rows) = &value {
                            self.frames[at].num_rows = Some(*rows as i32);
                        }
                    }
                    values.push((reading as u16, value));
                }
                self.next_member(ty, members, skip, values, reading + 1)
            }
            Open::Choice { index } => {
                self.inp.end_variant()?;
                self.drop_frame();
                self.inp.end_choice()?;
                self.drop_frame();
                Ok(Step::Done(Value::Choice(index, Box::new(value))))
            }
            Open::List {
                element,
                refs,
                skip,
                mut values,
            } => {
                self.inp.end_element()?;
                if !skip {
                    values.push(value);
                }
                self.next_element(element, refs, skip, values)
            }
            Open::Delayed { start } => Ok(Step::Done(Value::Delayed {
                encoding: I::ENCODING,
                bytes: self.inp.delayed(start),
            })),
            Open::Pointer => {
                self.end_pointer(Ok(()))?;
                Ok(Step::Done(value))
            }
            Open::Named { count } => self.end_names(count, Ok(value)).map(Step::Done),
        }
    }

    /// `ReadNamedType` / `SkipNamedType` of the stream `I` (`Input::NAMED_TYPES`): a frame of its
    /// own (`FrameKind::NamedType`) and `BeginNamedType` for each named type of the alias chain
    /// from `ty` that is neither a SEQUENCE, a CHOICE nor an ENUMERATED (those have tags of their
    /// own: `BeginClass`, `BeginChoice`, `ReadEnum`), outermost first; the number of them. An
    /// error leaves the frames pushed here.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2041-2052
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1302-1313
    // (An alias (`CAliasTypeInfo`) reads its aliased type through `ReadNamedType`, a named
    // implicit class (`Seq-descr`) its list through `ReadImplicitMember`, classinfo.cpp:289-306;
    // a named string, number or OCTET STRING is an alias of the anonymous type: each is one such
    // frame.)
    pub(super) fn begin_names(&mut self, mut ty: u16, skip: bool) -> R<u8> {
        let mut count = 0;
        loop {
            let info = super::type_info(ty);
            if !info.name().is_empty()
                && !matches!(
                    info.kind,
                    Kind::Sequence(_) | Kind::Choice(_) | Kind::Enumerated { .. }
                )
            {
                if let Err(fail) = self.push(FrameKind::NamedType, ty, skip) {
                    return Err(self.unwind(self.frames.len() - usize::from(count), fail));
                }
                count += 1;
                if let Err(fail) = self.inp.begin_named() {
                    return Err(self.unwind(self.frames.len() - usize::from(count), fail));
                }
            }
            match info.kind {
                Kind::Alias(next) => ty = next,
                _ => return Ok(count),
            }
        }
    }

    /// `EndNamedType` and `END_OBJECT_FRAME` of the `count` named frames of `begin_names`, the
    /// innermost first, after the value inside them was read (`result`).
    pub(super) fn end_names<T>(&mut self, count: u8, mut result: R<T>) -> R<T> {
        for _ in 0..count {
            if result.is_ok() {
                if let Err(fail) = self.inp.end_named() {
                    result = Err(fail);
                }
            }
            result = self.pop(result);
        }
        result
    }

    /// `ReadClassSequential` / `SkipClassSequential` (`start`: the class frame, `BeginClass`, the
    /// member frame): the next member in module order (`member_begin` of the class type and `pos`)
    /// and a read's (not a skip's) member hooks before its value (`hooks::pre_read_member`), its
    /// pointer forms (`member_pointer`: a `NULL` leaves it unset); after the last, the member
    /// frame left, the mandatory members after it checked, `EndClass`, the class frame left.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:84-109
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397 (the binary stream's,
    // objistrasnb.cpp:1210-1233, is the same with its own `BeginClass`)
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:847-866 (`ReadLongMember`: a delayed
    // member's bytes start before its pointer forms, skipped by `SkipPointer`)
    fn next_member(
        &mut self,
        ty: u16,
        members: &'static [Member],
        skip: bool,
        mut values: Vec<(u16, Value)>,
        mut pos: usize,
    ) -> R<Step> {
        while let Some(index) = self.member_begin(members, pos)? {
            let member = &members[index];
            if !skip {
                super::hooks::pre_read_member(ty, index, &values)?;
            }
            let depth = self.frames.len();
            let delayed = member.flags & DELAYED != 0 && !skip;
            let start = if delayed {
                Some(self.inp.delay_start()?)
            } else {
                None
            };
            let pointee = self.member_pointer(member, skip || delayed)?;
            if pointee == Pointee::Null {
                // A delayed `NULL` is kept as the bytes read, as NCBI's buffer keeps it.
                if let Some(start) = start {
                    values.push((
                        index as u16,
                        Value::Delayed {
                            encoding: I::ENCODING,
                            bytes: self.inp.delayed(start),
                        },
                    ));
                }
                self.inp.end_member()?;
                pos = index + 1;
                continue;
            }
            self.open.push(Open::Sequence {
                ty,
                members,
                skip,
                values,
                reading: index,
                depth,
            });
            if let Some(start) = start {
                self.open.push(Open::Delayed { start });
            }
            if pointee == Pointee::Other {
                self.open.push(Open::Pointer);
            }
            return Ok(Step::Read(member.ty, skip || delayed));
        }
        self.drop_frame();
        for member in &members[pos..] {
            self.expected_member(member)?;
        }
        self.inp.end_class()?;
        self.drop_frame();
        Ok(Step::Done(Value::Sequence(values)))
    }

    /// `ReadChoiceSimple` / `SkipChoiceSimple`: the variant of the CHOICE `ty`; a read (not a skip)
    /// runs the class's variant hooks before its value (`hooks::pre_read_variant`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1501
    fn start_choice(&mut self, ty: u16, variants: &'static [Member], skip: bool) -> R<Step> {
        self.push(FrameKind::Choice, ty, skip)?;
        self.inp.begin_choice()?;
        self.push(FrameKind::ChoiceVariant, NO_TYPE, skip)?;
        let index = self.inp.begin_variant(variants)?;
        let variant = &variants[index];
        self.set_member(variant);
        if !skip {
            super::hooks::pre_read_variant(ty, index, &self.frames)?;
        }
        self.open.push(Open::Choice {
            index: index as u16,
        });
        Ok(Step::Read(variant.ty, skip))
    }

    /// `ReadContainer` / `SkipContainer` of a SEQUENCE OF or (`set`) SET OF.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1860-1875 (the XML stream's
    // `ReadContainer`: no array frame for a container without a name, `Input::ARRAY_FRAMES`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1921-1979 (its
    // `ReadContainerContents`: an element frame only for an element type that will have no name,
    // `Input::element_frame`)
    fn start_list(&mut self, element: u16, set: bool, refs: bool, skip: bool) -> R<Step> {
        if I::ARRAY_FRAMES {
            self.push(FrameKind::Array, NO_TYPE, skip)?;
        }
        self.inp.begin_container(set)?;
        if I::element_frame(element) {
            self.push(FrameKind::ArrayElement, NO_TYPE, skip)?;
        }
        let numbers = match resolved(element).kind {
            // The elements of an INTEGER or REAL list are read as the numbers themselves (no
            // `Value` each: a PSSM holds thousands); they have no hook.
            Kind::Integer(storage @ (Storage::Int2 | Storage::Int4)) => {
                // `integer` checked the storage's range.
                Value::Int4s(self.integers(element, storage, skip, |value| value as i32)?)
            }
            Kind::Integer(storage) => {
                Value::Ints(self.integers(element, storage, skip, |value| value)?)
            }
            Kind::Real => Value::Reals(self.reals(element, skip)?),
            _ => return self.next_element(element, refs, skip, Vec::new()),
        };
        self.end_list(element)?;
        Ok(Step::Done(numbers))
    }

    /// The element frame left, the container's end, the array frame left (the frames of
    /// `start_list`).
    #[inline(always)]
    fn end_list(&mut self, element: u16) -> R<()> {
        if I::element_frame(element) {
            self.drop_frame();
        }
        self.inp.end_container()?;
        if I::ARRAY_FRAMES {
            self.drop_frame();
        }
        Ok(())
    }

    /// The next element of a SEQUENCE OF or SET OF of values (of `CRef<>`s, `refs`: each with the
    /// pointer forms; a `NULL` one is held as `Value::Null`, the empty `CRef<>` that NCBI's list
    /// holds: `Seq-descr { NULL }` is not empty, and the program that uses the list dereferences
    /// it or not, `Value::Null`); after the last, the element frame left, the container's end, the
    /// array frame left.
    fn next_element(
        &mut self,
        element: u16,
        refs: bool,
        skip: bool,
        mut values: Vec<Value>,
    ) -> R<Step> {
        while self.inp.begin_element()? {
            let pointee = if refs {
                self.begin_pointer(element, skip)?
            } else {
                Pointee::This
            };
            if pointee == Pointee::Null {
                if !skip {
                    values.push(Value::Null);
                }
                self.inp.end_element()?;
                continue;
            }
            self.open.push(Open::List {
                element,
                refs,
                skip,
                values,
            });
            if pointee == Pointee::Other {
                self.open.push(Open::Pointer);
            }
            return Ok(Step::Read(element, skip));
        }
        self.end_list(element)?;
        Ok(Step::Done(Value::List(values)))
    }

    /// The elements of an INTEGER list of the type `element`, each `store`d as the list keeps
    /// it. Not inlined, as `reals`: the loop over a PSSM's scores keeps its own registers, which
    /// the walk of `start` around it would share (line RG-F3, PR 2: callgrind of `-in_pssm` text).
    #[inline(never)]
    fn integers<T>(
        &mut self,
        element: u16,
        storage: Storage,
        skip: bool,
        store: impl Fn(i64) -> T,
    ) -> R<Vec<T>> {
        let mut values = Vec::new();
        while self.inp.begin_element()? {
            let number = if I::NAMED_TYPES {
                let count = self.begin_names(element, skip)?;
                let number = self.integer_number(storage, skip);
                self.end_names(count, number)?
            } else {
                self.integer_number(storage, skip)?
            };
            if let Some(value) = number {
                values.push(store(value));
            }
            self.inp.end_element()?;
        }
        Ok(values)
    }

    /// The elements of a REAL list of the type `element` (not inlined, as `integers`).
    #[inline(never)]
    fn reals(&mut self, element: u16, skip: bool) -> R<Vec<f64>> {
        let mut values = Vec::new();
        while self.inp.begin_element()? {
            let number = if I::NAMED_TYPES {
                let count = self.begin_names(element, skip)?;
                let number = self.real_number(skip);
                self.end_names(count, number)?
            } else {
                self.real_number(skip)?
            };
            if let Some(value) = number {
                values.push(value);
            }
            self.inp.end_element()?;
        }
        Ok(values)
    }

    /// An INTEGER of storage `storage` (`ReadStd`), or its skip (`SkipSNumber`, `SkipUNumber`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1780-1808 (the skips of every size:
    // `SkipSNumber` for the signed ones, `SkipUNumber` for the unsigned ones)
    fn integer_value(&mut self, storage: Storage, skip: bool) -> R<Value> {
        Ok(self
            .integer_number(storage, skip)?
            .map_or(Value::Null, Value::Int))
    }

    /// `integer`'s number (`None` when skipped), for the lists that keep the numbers (no `Value`
    /// to make and drop for each element).
    #[inline(always)]
    fn integer_number(&mut self, storage: Storage, skip: bool) -> R<Option<i64>> {
        if skip {
            self.inp.skip_number(storage != Storage::Uint4)?;
            return Ok(None);
        }
        self.read_integer(storage).map(Some)
    }

    /// `ReadStd` of an INTEGER of storage `storage`.
    #[inline(always)]
    pub(super) fn read_integer(&mut self, storage: Storage) -> R<i64> {
        let value = match storage {
            Storage::Int4 => i64::from(self.inp.read_i32()?),
            Storage::Uint4 => i64::from(self.inp.read_u32()?),
            Storage::Int8 => self.inp.read_i64()?,
            Storage::Int2 => {
                let value = self.inp.read_i32()?;
                if i16::try_from(value).is_err() {
                    return Err(throw_error(
                        READ_INT2,
                        self.inp.position(),
                        "integer overflow",
                    ));
                }
                i64::from(value)
            }
        };
        Ok(value)
    }

    /// `CEnumeratedTypeInfo::ReadEnum` / `SkipEnum` over `ReadEnum` of the stream: an exception
    /// gets the line of the type's `NCBI_RETHROW_SAME` (an end of the input stays one).
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // (`SkipEnum`, 468-484, reads the value the same way and rethrows at 478.)
    pub(super) fn enum_number(
        &mut self,
        values: &'static [EnumValue],
        integer: bool,
        skip: bool,
    ) -> R<i32> {
        let items = EnumType { integer, values };
        let result = self.inp.read_enum(&items);
        result.map_err(|fail| enum_error(fail, skip))
    }

    /// `ReadDouble` / `SkipFNumber`.
    fn real(&mut self, skip: bool) -> R<Value> {
        Ok(self.real_number(skip)?.map_or(Value::Null, Value::Real))
    }

    /// `real`'s number (`None` when skipped), for the lists that keep the numbers.
    fn real_number(&mut self, skip: bool) -> R<Option<f64>> {
        if skip {
            self.inp.skip_real()?;
            return Ok(None);
        }
        Ok(Some(self.inp.read_real()?))
    }

    /// `ReadString` (`ReadStringStore`) / `SkipString` of a string of `kind`, with NCBI's posts.
    pub(super) fn string(&mut self, kind: StringKind, skip: bool) -> R<Vec<u8>> {
        let mut out = Vec::new();
        let mut posts = Vec::new();
        let result = if skip {
            self.inp.skip_string(kind, &mut posts)
        } else {
            self.inp.read_string(kind, &mut out, &mut posts)
        };
        for post in posts {
            self.post(post);
        }
        result?;
        Ok(out)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (cited at
    // `report::serial::encoder::ReplacedChar`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-310 (cited at
    // `asn_input::utf8_tag_message`)
    /// What reading a string posted, with the stack trace (`GetStackTrace()`) that NCBI's message
    /// names.
    #[cold]
    #[inline(never)]
    fn post(&mut self, post: StringPost) {
        match post {
            StringPost::BadChar {
                position,
                byte,
                string,
            } => self.posts.push(Post::BadChar(ReplacedChar {
                trace: self.stack_trace(),
                byte,
                position,
                string: string.into(),
            })),
            StringPost::Utf8Tag => {
                if self.utf8_once.first() {
                    self.posts.push(Post::Utf8Tag(self.stack_trace()));
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:610-614
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1554-1565
    // (A `std::exception` changes no state of the stream; the end of the data inside the bytes
    // and the stream's own errors leave it failed, and the block then posts nothing.)
    /// `ReadBytes` / `SkipByteBlock` of an OCTET STRING. A `std::exception` that leaves the byte
    /// block (its reserve, `asn_input::reserve`) passes its destructor with the stream in a good
    /// state, which posts that the block was not ended before the exception goes on.
    fn octets(&mut self, skip: bool) -> R<Value> {
        if skip {
            self.inp.skip_octets()?;
            return Ok(Value::Null);
        }
        let value = self.inp.read_octets();
        if matches!(value, Err(Fail::OutOfMemory | Fail::Std(_))) {
            self.inp.drain_posts(&mut self.posts);
            self.posts.push(Post::UnendedByteBlock);
        }
        Ok(Value::Octets(value?))
    }
}

/// What `CEnumeratedTypeInfo::ReadEnum` (`SkipEnum`, `skip`) makes of `fail`: its
/// `NCBI_RETHROW_SAME` line.
#[cold]
#[inline(never)]
fn enum_error(fail: Fail, skip: bool) -> Fail {
    let (line, function) = if skip {
        (478, "ncbi::CEnumeratedTypeInfo::SkipEnum")
    } else {
        (428, "ncbi::CEnumeratedTypeInfo::ReadEnum")
    };
    let rethrow = |e: Box<NcbiException>| {
        e.rethrow_same(
            "src/serial/enumerated.cpp",
            line,
            function,
            "invalid enum value",
        )
    };
    match fail {
        Fail::Serial(e) => Fail::serial(rethrow(e)),
        Fail::Exception(e) => Fail::exception(rethrow(e)),
        other => other,
    }
}

/// Whether `fail` is the `eMissingValue` of the XML stream (`read_xml::MISSING_VALUE`).
pub(super) fn is_missing_value(fail: &Fail) -> bool {
    matches!(fail, Fail::Serial(e) if e.kind() == MISSING_VALUE)
}

/// The type `ty` refers to after its aliases (the index of `resolved(ty)`): the type of its frames
/// and of the hooks.
pub(super) fn alias_target(mut ty: u16) -> u16 {
    while let Kind::Alias(target) = super::type_info(ty).kind {
        ty = target;
    }
    ty
}
