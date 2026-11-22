//! The generic reader of the type table: NCBI's `CObjectIStream` with its object stack, for any
//! type of `table` (`read_object`, `read_document`), over the steps of one encoding (`Front`):
//! ASN.1 text (`read_text`, `CObjectIStreamAsn`), binary ASN.1 (`read_ber`,
//! `CObjectIStreamAsnBinary`) and XML (`read_xml`, `CObjectIStreamXml`). The walk of a value
//! (the frames of NCBI's stack, the members in module order, the mandatory members, the lists,
//! the hooks of NCBI's classes) is this module's; how each step reads its bytes is the front
//! end's.
//!
//! It reads as NCBI's readers read:
//! - the members of a SEQUENCE in module order (by name, by tag in brackets or by context tag);
//!   an unknown or repeated member, an unknown variant and a missing mandatory member are NCBI's
//!   `CSerialException`;
//! - INTEGER by the C++ storage of its member (`ReadInt4`, `ReadUint4`, `ReadInt8`, `ReadInt2`);
//!   ENUMERATED by a name or a number that has a name, INTEGER with named values by name or any
//!   number;
//! - a VisibleString (and, in text, StringStore) byte outside 0x20..0x7E as `#`, with NCBI's
//!   `ERR_POST_X` with subcode 7 (`Post::BadChar`), a UTF8String as it is;
//! - the delayed member (`Auth-list.names`) with the skip functions;
//! - the hooks of NCBI's hand-written classes that change what is read (`super::hooks` of this
//!   file: `CSeq_descr::PostRead`, `CVariation_ref::PostRead`, the reserve hooks), and the Seq-id
//!   into LOSAT's Seq-id model.
//!
//! Its errors are NCBI's classes (`ReadError`): a `CSerialException` (and an end of the data
//! inside an object, which NCBI's stack turns into one), which its callers report as theirs; a
//! `CUtilException` of a number (`BadNumber`, `NumberOverflow`) or a `CStringException` of a tag,
//! with the line NCBI's stack adds for each frame it leaves (`NCBI_RETHROW_SAME`, the report of
//! `CATCH_ALL`); and an end of the data directly in the root object, a `CEofException` that NCBI's
//! stack passes on.

use super::{resolved, EnumValue, Kind, Member, Storage, StringKind, Value};
use crate::ncbi_exception::NcbiException;
use crate::report::serial::asn_text_reader::{AsnTextError, AsnTextErrorKind, Reader, SliceSource};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::report::serial::object_frame::{frame_info, FrameKind};

use super::read_text::Text;

/// A message that NCBI's reader posts while it reads, in the order posted.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Post {
    /// `ReplaceVisibleChar`'s `ERR_POST_X` of subcode 7 with a message, for a byte read as `#`.
    BadChar(ReplacedChar),
    /// `ERR_POST` of an `Error` with a text, from a class's `PostRead`
    /// (`CVariation_ref::PostRead`).
    Error(&'static str),
    /// The binary reader's warning of subcode 10 (`ERR_POST_X_ONCE`: once in the process) of a
    /// UTF8String read for a VisibleString, with the path of the member (`GetStackTraceASN`).
    Utf8ForVisible(String),
}

impl Post {
    /// The line that a BLAST application of NCBI posts on its standard error, with its
    /// diagnostic prefix (`SetDiagPostPrefix`: the program's name in lower case for the search
    /// programs; none for blast_formatter): the one owner of each post's text
    /// (`ReplacedChar::post`, `ncbi_exception::error_post_line`).
    pub(crate) fn line(&self, prefix: Option<&str>) -> Vec<u8> {
        match self {
            Post::BadChar(replaced) => replaced.post(prefix),
            Post::Error(text) => {
                crate::ncbi_exception::error_post_line(None, prefix, text.as_bytes())
            }
            Post::Utf8ForVisible(trace) => crate::ncbi_exception::post_line(
                crate::ncbi_exception::Severity::Warning,
                Some("803.10"),
                prefix,
                format!(
                    "CObjectIStreamAsnBinary: UTF8String data for VisibleString member {trace}, \
                     ASN.1 specification may need an update"
                )
                .as_bytes(),
            ),
        }
    }
}

/// Why NCBI's reader stops: the class of what it throws.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum ReadError {
    /// A `CSerialException` (with its message, for LOSAT's tests; NCBI's callers print their own
    /// line for it).
    Serial(String),
    /// Another `CException` that leaves the reader (`CUtilException`, `CStringException`), with
    /// the lines of NCBI's object stack.
    Exception(NcbiException),
    /// A `CEofException` that leaves the reader: the data ends directly in the root object.
    Eof,
    /// A `std::exception` that leaves the reader (its `what()`): `std::length_error` of a reserve
    /// hook (`hooks::pre_read_member`).
    Std(&'static str),
    /// A form of the data that LOSAT's reader does not take (`AsnTextErrorKind::Unsupported`,
    /// LOSAT's own limit, not an NCBI error): the caller refuses the input with this text.
    Unsupported(String),
}

/// What the reader read, and the messages it posted before it stopped or finished.
pub(crate) type Read<T> = (std::result::Result<T, ReadError>, Vec<Post>);

pub(super) type R<T> = std::result::Result<T, ReadError>;

impl From<AsnTextError> for ReadError {
    fn from(error: AsnTextError) -> Self {
        match error.kind {
            AsnTextErrorKind::Eof => ReadError::Eof,
            AsnTextErrorKind::Format => ReadError::Serial(error.message),
            AsnTextErrorKind::Unsupported => ReadError::Unsupported(error.message),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Frame {
    kind: FrameKind,
    /// The type the frame names (`Object type= `): the named, class, choice and array frames.
    type_name: Option<&'static str>,
    /// The member or variant (`SetTopMemberId`) of a member or variant frame.
    member: Option<&'static str>,
    /// Whether the frame belongs to a skip function (`SkipClassSequential`, ...).
    skip: bool,
    /// The type (the table's index) of a class or choice frame, the object whose
    /// `HasTypeInfo` NCBI's `GetParentObjectPtr` asks (`hooks::pre_read_variant`).
    ty: Option<u16>,
    /// A `Seq-table` class frame's `num-rows` once read (the object's `GetNum_rows()`).
    num_rows: Option<i64>,
}

impl Frame {
    /// Whether the frame is a class or a choice frame (`eFrameClass`, `eFrameChoice`).
    pub(super) fn is_class_or_choice(&self) -> bool {
        matches!(self.kind, FrameKind::Class | FrameKind::Choice)
    }

    /// The type of a class or choice frame.
    pub(super) fn ty(&self) -> Option<u16> {
        self.ty
    }

    /// A `Seq-table` class frame's `num-rows`, once read.
    pub(super) fn num_rows(&self) -> Option<i64> {
        self.num_rows
    }

    /// The frame's type (`GetFrameType`).
    pub(super) fn kind(&self) -> FrameKind {
        self.kind
    }

    /// The name of the type the frame names (`GetTypeInfo()->GetName()`); empty for none.
    pub(super) fn type_name(&self) -> &'static str {
        self.type_name.unwrap_or("")
    }

    /// The member or variant of a member or variant frame (`GetMemberId().GetName()`).
    pub(super) fn member(&self) -> &'static str {
        self.member.unwrap_or("")
    }
}

impl Frame {
    /// `GetFrameInfo()` of the frame (`object_frame::frame_info`).
    fn info(&self) -> String {
        frame_info(self.kind, self.type_name, self.member)
    }
}

/// `GetStackTrace()`: the path of the object being read (`PssmWithParameters.pssm.query.seq.
/// descr.E.comment`), by the one writer of NCBI's object stack path.
fn stack_trace(frames: &[Frame]) -> String {
    object_stack_trace(
        frames
            .first()
            .and_then(|frame| frame.type_name)
            .unwrap_or("?"),
        frames.iter().skip(1).filter_map(|frame| match frame.kind {
            FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                frame.member.map(StackStep::Member)
            }
            FrameKind::ArrayElement => Some(StackStep::Element),
            _ => None,
        }),
    )
}

/// What a step of a front end sees of the reader: NCBI's object stack (the frames, the innermost
/// last) and the messages posted so far.
pub(super) struct Cx<'a> {
    pub(super) frames: &'a [Frame],
    posts: &'a mut Vec<Post>,
}

impl Cx<'_> {
    /// `GetStackTrace()` of the frames.
    pub(super) fn trace(&self) -> String {
        stack_trace(self.frames)
    }

    /// `ReplaceVisibleChar`'s post of `byte` replaced at `position`, quoting `string`.
    pub(super) fn post_bad_char(&mut self, byte: u8, string: &[u8], position: StreamPosition) {
        let trace = self.trace();
        self.posts.push(Post::BadChar(ReplacedChar {
            trace,
            byte,
            position,
            string: string.into(),
        }));
    }

    /// A post other than a replaced byte.
    pub(super) fn post(&mut self, post: Post) {
        self.posts.push(post);
    }
}

/// The steps of one encoding of NCBI's object stream that the walk of a value calls, each as
/// that stream's function of the same name does it: where it ends, what it throws.
pub(super) trait Front {
    /// Whether the stream opens a tag of its own for a named type that is neither a SEQUENCE nor a
    /// CHOICE (`BeginNamedType`, `EndNamedType` in a frame of `FrameKind::NamedType`).
    const NAMED_TYPES: bool = false;
    /// Whether a SEQUENCE OF has a frame of its own (`eFrameArray`): `CObjectIStreamXml::
    /// ReadContainer` of a container without a name pushes none.
    const ARRAY_FRAMES: bool = true;
    /// Whether the elements of a list of `element`s have the frame of an element
    /// (`eFrameArrayElement`).
    fn element_frame(_element: u16) -> bool {
        true
    }
    /// `BeginNamedType` of the named type `ty`, after its frame was pushed.
    fn begin_named(&mut self, _ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }
    /// `EndNamedType`, before its frame is left.
    fn end_named(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }
    /// `BeginChoice` of the CHOICE `ty`, after its frame was pushed.
    fn begin_choice(&mut self, _ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }
    /// `EndChoice`, after the variant's frame was left.
    fn end_choice(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }
    /// The end of a list's element (`EndArrayElement`), with the element's frame still on the
    /// stack (when it has one).
    fn end_element(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }
    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves a frame of
    /// `kind` (of a skip function: `skip`): the place of the line that `NCBI_RETHROW_SAME` adds to
    /// an exception passing through it.
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str);
    /// Before an object: the bytes before it are not read again (one object in memory).
    fn start_object(&mut self);
    /// `ReadFileHeader`: the type name in the header (empty when the encoding has none).
    fn file_header(&mut self, cx: &mut Cx<'_>) -> R<Vec<u8>>;
    /// `BeginClass` of the SEQUENCE `ty`.
    fn begin_class(&mut self, ty: u16, cx: &mut Cx<'_>) -> R<()>;
    /// `BeginClassMember`: the position of the next member (one of `members` from `pos`), or
    /// `None` at the end of the class. `block` is the class's state (`true` before its first
    /// member).
    fn next_member(
        &mut self,
        members: &'static [Member],
        pos: usize,
        block: &mut bool,
        cx: &mut Cx<'_>,
    ) -> R<Option<usize>>;
    /// `EndClassMember` after a member's value.
    fn end_member(&mut self, cx: &mut Cx<'_>) -> R<()>;
    /// `EndClass`, after the member frame was left and the mandatory members checked.
    fn end_class(&mut self, cx: &mut Cx<'_>) -> R<()>;
    /// `BeginChoiceVariant` of the CHOICE `ty`: the position of the variant.
    fn begin_variant(&mut self, ty: u16, variants: &'static [Member], cx: &mut Cx<'_>) -> R<usize>;
    /// `EndChoiceVariant` after the variant's value.
    fn end_variant(&mut self, cx: &mut Cx<'_>) -> R<()>;
    /// `BeginContainer` of a SEQUENCE OF (`set` false) or SET OF.
    fn begin_list(&mut self, set: bool, cx: &mut Cx<'_>) -> R<()>;
    /// `BeginContainerElement`: whether an element of `element` follows. `block` is the list's
    /// state (`true` before its first element).
    fn next_element(&mut self, element: u16, block: &mut bool, cx: &mut Cx<'_>) -> R<bool>;
    /// `EndContainer`, after the element frame was left.
    fn end_list(&mut self, cx: &mut Cx<'_>) -> R<()>;
    /// An INTEGER of the type `ty` and storage `storage` (`ReadStd`), or its skip (`None`).
    fn integer(&mut self, ty: u16, storage: Storage, skip: bool, cx: &mut Cx<'_>)
        -> R<Option<i64>>;
    /// `ReadEnum` of the ENUMERATED or INTEGER with named values `ty` (the skip reads the same).
    fn enumerated(
        &mut self,
        ty: u16,
        values: &'static [EnumValue],
        integer: bool,
        skip: bool,
        cx: &mut Cx<'_>,
    ) -> R<i32>;
    /// `ReadBool` / `SkipBool`.
    fn boolean(&mut self, ty: u16, skip: bool, cx: &mut Cx<'_>) -> R<bool>;
    /// `ReadNull` / `SkipNull`.
    fn null(&mut self, ty: u16, skip: bool, cx: &mut Cx<'_>) -> R<()>;
    /// `ReadDouble` / `SkipFNumber` (`None`).
    fn real(&mut self, ty: u16, skip: bool, cx: &mut Cx<'_>) -> R<Option<f64>>;
    /// `ReadString` / `SkipString` of a string of `kind` (an empty one when skipped), with the
    /// posts of the bytes replaced.
    fn string(&mut self, ty: u16, kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>>;
    /// `ReadBytes` / `SkipByteBlock` of an OCTET STRING (`None` when skipped).
    fn octets(&mut self, ty: u16, skip: bool, cx: &mut Cx<'_>) -> R<Option<Vec<u8>>>;
    /// `StartDelayBuffer` before a delayed member's value is skipped: the position its buffer
    /// starts at.
    fn delay_start(&mut self, cx: &mut Cx<'_>) -> R<usize>;
    /// `EndDelayBuffer` after the delayed member of the type `ty` was skipped from `start`: the
    /// member's value.
    fn delay_end(&mut self, start: usize, ty: u16, cx: &mut Cx<'_>) -> R<Value>;
}

/// Reads the ASN.1 text document at the start of `text` as the root type `root` (an object of
/// `CObjectIStream::Read`: its file header, then the value; what follows it is not read), with
/// the messages posted.
pub(crate) fn read_document(text: &[u8], root: u16) -> Read<Value> {
    read_object(
        &mut Text::new(Reader::new(SliceSource::new(text), "PSSM", "PSSM reader")),
        root,
    )
}

/// Reads one object of the root type `root` from where the front end `front` stands (`in >>
/// object`: `CObjectIStream::Read`, its file header, then the value), with the messages posted.
/// The front end stands after the object's last byte when it is read (a stream of objects reads
/// the next one from there, `stream::ObjectStream`).
pub(super) fn read_object<F: Front>(front: &mut F, root: u16) -> Read<Value> {
    front.start_object();
    let mut reader = TypedReader {
        r: front,
        frames: Vec::new(),
        posts: Vec::new(),
        open: Vec::new(),
    };
    let value = reader.document(root);
    (value, reader.posts)
}

/// Reads one value of the type `ty` from where `front` stands, as the value of a member that
/// NCBI keeps in a delay buffer and parses when it is first used (`CDelayBuffer::Update`: a new
/// object stream over the buffer's bytes, no frame of the stream that read the buffer); its
/// posts go to those of `cx`.
pub(super) fn read_delayed<F: Front>(front: &mut F, ty: u16, cx: &mut Cx<'_>) -> R<Value> {
    let mut reader = TypedReader {
        r: front,
        frames: Vec::new(),
        posts: Vec::new(),
        open: Vec::new(),
    };
    let name = super::type_info(ty).name();
    reader.push(FrameKind::Named, Some(name), false);
    let value = reader.value(ty, false);
    let value = reader.pop(value);
    cx.posts.append(&mut reader.posts);
    value
}

/// `read_delayed` for a stream whose tags are named by the frames of the object that holds the
/// member (`SkipStackTagName`): the new stream's stack starts with `outer`, the frames of the
/// class and of the member that was read into the buffer (`CDelayBuffer::Update`, delaybuf.cpp:
/// 63-84, `BEGIN_OBJECT_FRAME_OF2` of the class and of the member).
pub(super) fn read_delayed_in<F: Front>(
    front: &mut F,
    outer: &[Frame],
    ty: u16,
    cx: &mut Cx<'_>,
) -> R<Value> {
    let mut reader = TypedReader {
        r: front,
        frames: outer.to_vec(),
        posts: Vec::new(),
        open: Vec::new(),
    };
    let value = reader.value(ty, false);
    cx.posts.append(&mut reader.posts);
    value
}

struct TypedReader<'r, F: Front> {
    r: &'r mut F,
    frames: Vec<Frame>,
    posts: Vec<Post>,
    /// The values being read, outermost first (`Open`): NCBI's recursion of `ReadObject` as a list
    /// on the heap, so that a value nested to any depth is read in LOSAT's own stack.
    open: Vec<Open>,
}

/// A value being read that waits for a value inside it (`TypedReader::value`).
enum Open {
    /// The end of `ReadObject` of the type `target`: the hooks of its class (not after a skip).
    Value { target: u16, skip: bool },
    /// A SEQUENCE reading the member at `reading` (`members`).
    Sequence {
        ty: u16,
        members: &'static [Member],
        skip: bool,
        values: Vec<(u16, Value)>,
        reading: usize,
        block_start: bool,
    },
    /// A CHOICE reading its variant.
    Choice { index: u16 },
    /// A SEQUENCE OF or SET OF reading an element.
    List {
        element: u16,
        skip: bool,
        values: Vec<Value>,
        block_start: bool,
    },
    /// A delayed member of the type `ty` (`delayed`) from `start`.
    Delayed { start: usize, ty: u16 },
    /// The named types around a value (`ReadNamedType`): the end of `count` of them, the
    /// innermost first.
    Named { count: usize },
}

/// What the reader does next: read a value of a type (`skip`), or give a value read to the one
/// that holds it.
enum Step {
    Read(u16, bool),
    Done(Value),
}

/// The context of a front end's step: the reader's frames and posts.
macro_rules! cx {
    ($reader:expr) => {
        &mut Cx {
            frames: &$reader.frames,
            posts: &mut $reader.posts,
        }
    };
}

impl<F: Front> TypedReader<'_, F> {
    fn push(&mut self, kind: FrameKind, type_name: Option<&'static str>, skip: bool) {
        self.frames.push(Frame {
            kind,
            type_name,
            member: None,
            skip,
            ty: None,
            num_rows: None,
        });
    }

    /// `push` of a class or choice frame of the type `ty`.
    fn push_object(&mut self, kind: FrameKind, ty: u16, type_name: &'static str, skip: bool) {
        self.push(kind, Some(type_name), skip);
        if let Some(frame) = self.frames.last_mut() {
            frame.ty = Some(ty);
        }
    }

    /// `END_OBJECT_FRAME`: the frame is left, and an error that leaves it gets NCBI's handling
    /// (`Front::frame_site`): a `CException` other than a `CSerialException` takes the frame's
    /// line; an end of the data (`CEofException`, `HandleEOF`, objistr.cpp:599-608) becomes a
    /// `CSerialException` unless fewer than 2 frames stay.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608
    fn pop<T>(&mut self, result: R<T>) -> R<T> {
        let frame = self.frames.pop().expect("a frame to leave");
        result.map_err(|error| match error {
            ReadError::Exception(e) => {
                let (file, line, function) = F::frame_site(frame.kind, frame.skip);
                ReadError::Exception(e.rethrow_same(file, line, function, frame.info()))
            }
            ReadError::Eof if self.frames.len() >= 2 => {
                ReadError::Serial(format!("unexpected end of the data:{}", frame.info()))
            }
            other => other,
        })
    }

    fn set_member(&mut self, name: &'static str) {
        if let Some(top) = self.frames.last_mut() {
            top.member = Some(name);
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-900 and 959-969
    fn document(&mut self, root: u16) -> R<Value> {
        let name = super::type_info(root).name();
        self.push(FrameKind::Named, Some(name), false);
        let header = self.r.file_header(cx!(self)).and_then(|id| {
            if !id.is_empty() && id != name.as_bytes() {
                return Err(ReadError::Serial(format!(
                    "incompatible type {}<>{name}",
                    String::from_utf8_lossy(&id)
                )));
            }
            Ok(())
        });
        self.pop(header)?;
        self.push(FrameKind::Named, Some(name), false);
        let value = self.value(root, false);
        self.pop(value)
    }

    /// One value of the type `ty` (`ReadObject`), or its skip (`SkipObject`, `skip`), with the
    /// hooks of NCBI's classes that a read runs. The values inside it are read in turn from the
    /// list of the values being read (`open`), as NCBI's functions recurse; an error leaves
    /// NCBI's frames from the innermost out (`pop`), as it passes `END_OBJECT_FRAME` of each.
    fn value(&mut self, ty: u16, skip: bool) -> R<Value> {
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
                Err(mut error) => {
                    self.open.truncate(open);
                    while self.frames.len() > frames {
                        error = match self.pop::<()>(Err(error)) {
                            Err(error) => error,
                            Ok(()) => unreachable!("an error stays one"),
                        };
                    }
                    return Err(error);
                }
            };
        }
    }

    /// The start of a value of the type `ty` (`ReadObject`'s type function): a value without
    /// values inside it is read here, a SEQUENCE, CHOICE or SEQUENCE OF opens (`Open`).
    fn start(&mut self, ty: u16, skip: bool) -> R<Step> {
        let info = resolved(ty);
        let target = alias_target(ty);
        let names = self.begin_names(ty, skip)?;
        if names > 0 {
            self.open.push(Open::Named { count: names });
        }
        self.open.push(Open::Value { target, skip });
        let value = match info.kind {
            Kind::Sequence(members) => {
                self.push_object(FrameKind::Class, target, info.name(), skip);
                self.r.begin_class(target, cx!(self))?;
                self.push(FrameKind::ClassMember, None, skip);
                return self.next_member(target, members.members(), skip, Vec::new(), 0, true);
            }
            Kind::Choice(variants) => {
                return self.start_choice(target, info.name(), variants.members(), skip)
            }
            Kind::List { element, set } => return self.start_list(element, set, skip),
            Kind::Alias(_) => unreachable!("resolved"),
            Kind::Integer(storage) => self
                .r
                .integer(target, storage, skip, cx!(self))?
                .map_or(Value::Null, Value::Int),
            Kind::Enumerated { integer, values } => {
                self.enumerated(target, values.values(), integer, skip)?
            }
            Kind::Boolean => Value::Bool(self.r.boolean(target, skip, cx!(self))?),
            Kind::Null => {
                self.r.null(target, skip, cx!(self))?;
                Value::Null
            }
            Kind::Real => self
                .r
                .real(target, skip, cx!(self))?
                .map_or(Value::Null, Value::Real),
            Kind::String(kind) => Value::Str(self.r.string(target, kind, skip, cx!(self))?),
            Kind::Octets => self
                .r
                .octets(target, skip, cx!(self))?
                .map_or(Value::Null, Value::Octets),
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
                super::hooks::post_read(target, value, &mut self.posts).map(Step::Done)
            }
            Open::Sequence {
                ty,
                members,
                skip,
                mut values,
                reading,
                block_start,
            } => {
                self.r.end_member(cx!(self))?;
                if !skip {
                    // The `num-rows` of a Seq-table, which the reserve hook of its columns' data
                    // reads from the object being read (`hooks::pre_read_variant`).
                    if ty == super::table::named::SEQ_TABLE && members[reading].name() == "num-rows"
                    {
                        let at = self.frames.len() - 2;
                        if let Value::Int(rows) = &value {
                            self.frames[at].num_rows = Some(*rows);
                        }
                    }
                    values.push((reading as u16, value));
                }
                self.next_member(ty, members, skip, values, reading + 1, block_start)
            }
            Open::Choice { index } => {
                self.r.end_variant(cx!(self))?;
                self.frames.pop();
                self.r.end_choice(cx!(self))?;
                self.frames.pop();
                Ok(Step::Done(Value::Choice(index, Box::new(value))))
            }
            Open::List {
                element,
                skip,
                mut values,
                block_start,
            } => {
                self.r.end_element(cx!(self))?;
                if !skip {
                    values.push(value);
                }
                self.next_element(element, skip, values, block_start)
            }
            Open::Delayed { start, ty } => self.r.delay_end(start, ty, cx!(self)).map(Step::Done),
            Open::Named { count } => {
                self.end_names(count)?;
                Ok(Step::Done(value))
            }
        }
    }

    /// `ReadNamedType` / `SkipNamedType` of each named type of the alias chain from `ty` that a
    /// stream opens a tag for (`Front::NAMED_TYPES`), outermost first: its frame, then
    /// `BeginNamedType`. The number of them. A named SEQUENCE or CHOICE has its own tag (its
    /// class or choice calls) and so has a named ENUMERATED (`ReadEnum` of a type with a name).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2039-2052 and objistr.cpp:1302-1313
    // (An alias (`CAliasTypeInfo`) reads through `ReadNamedType` of its aliased type, a named
    // implicit class (`Seq-descr`) through `ReadImplicitMember`, a named string, number or
    // OCTET STRING as an alias of the anonymous type: each is one such frame.)
    fn begin_names(&mut self, mut ty: u16, skip: bool) -> R<usize> {
        let mut count = 0;
        if !F::NAMED_TYPES {
            return Ok(0);
        }
        loop {
            let info = super::type_info(ty);
            if !info.name().is_empty()
                && !matches!(
                    info.kind,
                    Kind::Sequence(_) | Kind::Choice(_) | Kind::Enumerated { .. }
                )
            {
                self.push(FrameKind::NamedType, Some(info.name()), skip);
                count += 1;
                self.r.begin_named(ty, cx!(self))?;
            }
            match info.kind {
                Kind::Alias(next) => ty = next,
                _ => return Ok(count),
            }
        }
    }

    /// `EndNamedType` and the frame left, for each of the `count` named types, the innermost first.
    fn end_names(&mut self, count: usize) -> R<()> {
        for _ in 0..count {
            self.r.end_named(cx!(self))?;
            self.frames.pop();
        }
        Ok(())
    }

    /// `ReadClassSequential` / `SkipClassSequential` (`start`: the class frame, `BeginClass`, the
    /// member frame): the next member in module order (`BeginClassMember` from `pos`) with the
    /// mandatory members before it checked, and a read's (not a skip's) member hooks before
    /// its value (`hooks::pre_read_member`); after the last, the member frame left, the mandatory
    /// members after it checked, `EndClass`, the class frame left.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:86-111
    fn next_member(
        &mut self,
        ty: u16,
        members: &'static [Member],
        skip: bool,
        values: Vec<(u16, Value)>,
        pos: usize,
        mut block_start: bool,
    ) -> R<Step> {
        if let Some(index) = self
            .r
            .next_member(members, pos, &mut block_start, cx!(self))?
        {
            let member = &members[index];
            self.set_member(member.name());
            for missing in &members[pos..index] {
                expected_member(missing)?;
            }
            if !skip {
                super::hooks::pre_read_member(ty, index, &values)?;
            }
            let delayed = member.flags & super::DELAYED != 0 && !skip;
            let start = if delayed {
                Some(self.r.delay_start(cx!(self))?)
            } else {
                None
            };
            self.open.push(Open::Sequence {
                ty,
                members,
                skip,
                values,
                reading: index,
                block_start,
            });
            if let Some(start) = start {
                self.open.push(Open::Delayed {
                    start,
                    ty: member.ty,
                });
                return Ok(Step::Read(member.ty, true));
            }
            return Ok(Step::Read(member.ty, skip));
        }
        self.frames.pop();
        for member in &members[pos..] {
            expected_member(member)?;
        }
        self.r.end_class(cx!(self))?;
        self.frames.pop();
        Ok(Step::Done(Value::Sequence(values)))
    }

    /// `ReadChoiceSimple` / `SkipChoiceSimple`: the variant of the CHOICE `ty`; a read (not a skip)
    /// runs the class's variant hooks before its value (`hooks::pre_read_variant`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1501
    fn start_choice(
        &mut self,
        ty: u16,
        name: &'static str,
        variants: &'static [Member],
        skip: bool,
    ) -> R<Step> {
        self.push_object(FrameKind::Choice, ty, name, skip);
        self.r.begin_choice(ty, cx!(self))?;
        self.push(FrameKind::ChoiceVariant, None, skip);
        let index = self.r.begin_variant(ty, variants, cx!(self))?;
        let variant = &variants[index];
        self.set_member(variant.name());
        if !skip {
            super::hooks::pre_read_variant(ty, index, &self.frames)?;
        }
        self.open.push(Open::Choice {
            index: index as u16,
        });
        Ok(Step::Read(variant.ty, skip))
    }

    /// `ReadContainer` / `SkipContainer`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274
    // The container's type information is anonymous (an implicit class such as `Seq-descr` has
    // no frame of its own, `ReadImplicitMember`, classinfo.cpp:289-306): `Object type= ` names
    // nothing (c07 of SURVEY_D3).
    fn start_list(&mut self, element: u16, set: bool, skip: bool) -> R<Step> {
        if F::ARRAY_FRAMES {
            self.push(FrameKind::Array, Some(""), skip);
        }
        self.r.begin_list(set, cx!(self))?;
        if F::element_frame(element) {
            self.push(FrameKind::ArrayElement, None, skip);
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
            _ => return self.next_element(element, skip, Vec::new(), true),
        };
        self.end_list_frames(element)?;
        Ok(Step::Done(numbers))
    }

    /// The next element of a SEQUENCE OF or SET OF of values; after the last, the element frame
    /// left, `EndContainer`, the array frame left.
    fn next_element(
        &mut self,
        element: u16,
        skip: bool,
        values: Vec<Value>,
        mut block_start: bool,
    ) -> R<Step> {
        if self.r.next_element(element, &mut block_start, cx!(self))? {
            self.open.push(Open::List {
                element,
                skip,
                values,
                block_start,
            });
            return Ok(Step::Read(element, skip));
        }
        self.end_list_frames(element)?;
        Ok(Step::Done(Value::List(values)))
    }

    /// The element frame left, `EndContainer`, the array frame left (the frames a list has).
    fn end_list_frames(&mut self, element: u16) -> R<()> {
        if F::element_frame(element) {
            self.frames.pop();
        }
        self.r.end_list(cx!(self))?;
        if F::ARRAY_FRAMES {
            self.frames.pop();
        }
        Ok(())
    }

    /// The elements of an INTEGER list of the type `element`, each `store`d as the list keeps it.
    fn integers<T>(
        &mut self,
        element: u16,
        storage: Storage,
        skip: bool,
        store: impl Fn(i64) -> T,
    ) -> R<Vec<T>> {
        let target = alias_target(element);
        let mut values = Vec::new();
        let mut block_start = true;
        while self.r.next_element(element, &mut block_start, cx!(self))? {
            let names = self.begin_names(element, skip)?;
            if let Some(value) = self.r.integer(target, storage, skip, cx!(self))? {
                values.push(store(value));
            }
            self.end_names(names)?;
            self.r.end_element(cx!(self))?;
        }
        Ok(values)
    }

    /// The elements of a REAL list of the type `element`.
    fn reals(&mut self, element: u16, skip: bool) -> R<Vec<f64>> {
        let target = alias_target(element);
        let mut values = Vec::new();
        let mut block_start = true;
        while self.r.next_element(element, &mut block_start, cx!(self))? {
            let names = self.begin_names(element, skip)?;
            if let Some(value) = self.r.real(target, skip, cx!(self))? {
                values.push(value);
            }
            self.end_names(names)?;
            self.r.end_element(cx!(self))?;
        }
        Ok(values)
    }

    /// `CEnumeratedTypeInfo::ReadEnum` / `SkipEnum` over `ReadEnum` of the stream; an exception
    /// that is not a `CSerialException` gets the line of the type's `NCBI_RETHROW_SAME`.
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // (`SkipEnum`, 468-484, reads the value the same way and rethrows at 478.)
    fn enumerated(
        &mut self,
        ty: u16,
        values: &'static [EnumValue],
        integer: bool,
        skip: bool,
    ) -> R<Value> {
        let result = self.r.enumerated(ty, values, integer, skip, cx!(self));
        result.map(Value::Enum).map_err(|error| match error {
            ReadError::Exception(e) => {
                let (line, function) = if skip {
                    (478, "ncbi::CEnumeratedTypeInfo::SkipEnum")
                } else {
                    (428, "ncbi::CEnumeratedTypeInfo::ReadEnum")
                };
                ReadError::Exception(e.rethrow_same(
                    "src/serial/enumerated.cpp",
                    line,
                    function,
                    "invalid enum value",
                ))
            }
            other => other,
        })
    }
}

/// `ReadMissingSimpleMember` -> `ExpectedMember` of a mandatory member: "member X expected"
/// (`CSerialException`, NCBI's default verification).
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
fn expected_member(member: &Member) -> R<()> {
    if member.optional() {
        return Ok(());
    }
    Err(ReadError::Serial(format!(
        "member {} expected",
        member.name()
    )))
}

/// The type `ty` refers to after its aliases (the index of `resolved(ty)`), for the hooks.
pub(super) fn alias_target(mut ty: u16) -> u16 {
    while let Kind::Alias(target) = super::type_info(ty).kind {
        ty = target;
    }
    ty
}
