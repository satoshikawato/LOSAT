//! The generic reader under the type table: NCBI's object stack (`CObjectIStream`: its frames,
//! `ReadClassSequential`, `ReadChoiceSimple`, `ReadContainer`, `HandleEOF`) for any type of
//! `table`, over the primitives of either encoding (`report::serial::asn_input::Input`): ASN.1
//! text (`read_document`, `CObjectIStreamAsn`) and ASN.1 binary (`read_ber_document`,
//! `CObjectIStreamAsnBinary`), one object or objects one after another from one stream
//! (`ObjectStream`).
//!
//! It reads as NCBI's reader reads:
//! - the members of a SEQUENCE in module order (named or by their tag `[3]` in text, by their
//!   context tag `a3` in binary: NCBI's automatic tags, which `scripts/asn_spec_gen.py` checks);
//!   an unknown or repeated member, an unknown variant and a missing mandatory member are NCBI's
//!   `CSerialException`;
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
//!   LOSAT's Seq-id model.
//!
//! Its errors are NCBI's classes (`ReadError`): a `CSerialException` (and an end of the input
//! inside an object, which NCBI's stack turns into one, `HandleEOF`), with the frames it leaves
//! (`AddFrameInfo`); a `CUtilException` of a number or a `CStringException` of a tag (text), with
//! the line NCBI's stack adds for each frame it leaves (`NCBI_RETHROW_SAME`, the report of
//! `CATCH_ALL`); an end of the input directly in the root object, a `CEofException` that NCBI's
//! stack passes on; the `std::exception`s of a reserve; and LOSAT's refusals of what its readers do
//! not take.

use super::{resolved, EnumValue, Kind, Member, Storage, Value};
use crate::ncbi_exception::{NcbiException, Severity};
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{
    read_with_set_flag_member, throw_error, utf8_tag_message, EnumItems, Fail, Input, Items, Names,
    StringKind, StringPost, TextInput, Utf8TagOnce, R,
};
use crate::report::serial::asn_text_reader::{AsnSource, SliceSource};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep};
use crate::report::serial::object_frame::{
    frame_info, frame_name, text_frame_site, FrameKind, EXPECTED_MEMBER, HANDLE_EOF, READ_INT2,
    SKIP_FILE_HEADER,
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
    /// The line that a BLAST application of NCBI posts on its standard error, with its
    /// diagnostic prefix `program` (the program's name in lower case, if it has one): the one owner of each
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
            .find(|value| value.name().as_bytes() == name)
            .map(|value| value.value)
    }

    fn has_value(&self, value: i32) -> bool {
        self.values.iter().any(|named| named.value == value)
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

    /// `GetFrameInfo()` (`object_frame::frame_info`).
    fn info(&self) -> String {
        frame_info(self.kind, self.type_name, self.member)
    }

    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves this frame
    /// (`object_frame::text_frame_site`; a binary stream's frames see no other `CException`).
    fn site(&self) -> (&'static str, u32, &'static str) {
        text_frame_site(self.kind, self.skip)
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

// NCBI reference (598d8ae6): c++/src/serial/delaybuf.cpp:63-84
/// The value of the type `ty` that a delayed member read from ASN.1 binary holds
/// (`CDelayBuffer::Update`: a new binary stream over the bytes kept, `Value::Delayed`), which a
/// writer reads before it writes it (`write`).
pub(super) fn read_ber_value(bytes: &[u8], ty: u16) -> Read<Value> {
    let mut reader = TypedReader::new(BerReader::new(SliceSource::new(bytes), Names::DELAYED));
    // The read that kept the bytes met their UTF8String tags (and posted for the first).
    reader.utf8_once = Utf8TagOnce::POSTED;
    let value = reader.value(ty, false).map_err(ReadError::from);
    (value, reader.posts)
}

/// Objects one after another from one input (`ObjectStream::read_object`), each read by a new
/// state of NCBI's stream at the byte after the last one the object before it took (NCBI gives the
/// bytes its stream did not use back when an object ends, `CIStreamBuffer::Close`,
/// strbuffer.cpp:125-139): its line count or its byte count and data limits start again, and its
/// frames and posts are its own.
pub(crate) struct ObjectStream<I: Input> {
    reader: TypedReader<I>,
}

impl<S: AsnSource> ObjectStream<TextInput<S>> {
    /// Objects in ASN.1 text from `src`, which LOSAT's messages call by `names`.
    pub(crate) fn text(src: S, names: Names) -> Self {
        ObjectStream {
            reader: TypedReader::new(TextInput::new(src, names)),
        }
    }
}

impl<S: AsnSource> ObjectStream<BerReader<S>> {
    /// Objects in ASN.1 binary from `src`, which LOSAT's messages call by `names`.
    pub(crate) fn binary(src: S, names: Names) -> Self {
        ObjectStream {
            reader: TypedReader::new(BerReader::new(src, names)),
        }
    }
}

impl<I: Input> ObjectStream<I> {
    /// The next object, as the root type `root`, with the messages posted while it was read.
    pub(crate) fn read_object(&mut self, root: u16) -> Read<Value> {
        let reader = &mut self.reader;
        reader.inp.new_object();
        reader.frames.clear();
        reader.open.clear();
        let value = reader.document(root).map_err(ReadError::from);
        (value, std::mem::take(&mut reader.posts))
    }

    /// One value of the type `ty` at the stream's position, inside an object that a hand grammar
    /// reads over the same input (rpsbproc's archive reader, `rpsbproc::archive`), with the
    /// messages posted. NCBI reads the whole object with its typed reader; the value is read as
    /// that reader reads it, below the frames of `root` (the object's named frame) and of `steps`
    /// (a member or variant frame for each name, an element frame for each `E`), which its posts'
    /// stack traces name. The frames of the classes and containers between them (which show in
    /// no trace) are not there: a failure that leaves the value ends the object, whose reader
    /// only asks whether it is LOSAT's refusal (`ReadError::Unsupported`).
    pub(crate) fn read_value(
        &mut self,
        ty: u16,
        root: &'static str,
        steps: &[StackStep<'static>],
    ) -> Read<Value> {
        let reader = &mut self.reader;
        reader.frames.clear();
        reader.open.clear();
        reader.push(FrameKind::Named, Some(root), false);
        for step in steps {
            match *step {
                StackStep::Member(name) => {
                    reader.push(FrameKind::ClassMember, None, false);
                    reader.set_member(name);
                }
                StackStep::Element => reader.push(FrameKind::ArrayElement, None, false),
            }
        }
        let value = reader.value(ty, false).map_err(ReadError::from);
        (value, std::mem::take(&mut reader.posts))
    }

    /// The primitives of the stream, for the hand grammar that reads around `read_value`.
    pub(crate) fn input_mut(&mut self) -> &mut I {
        &mut self.reader.inp
    }

    /// `EndOfData()`: whether no object follows (`Input::end_of_data`), which a reader of objects
    /// one after another asks before each (a binary object that is not a SEQUENCE ends inside its
    /// frames, where the end of the input is NCBI's `CSerialException`).
    #[cfg(test)]
    pub(crate) fn end_of_data(&mut self) -> bool {
        self.reader.inp.end_of_data()
    }

    /// The primitives of the stream (where it stands: `Input::position`).
    #[cfg(test)]
    pub(crate) fn input(&self) -> &I {
        &self.reader.inp
    }
}

struct TypedReader<I: Input> {
    inp: I,
    frames: Vec<Frame>,
    posts: Vec<Post>,
    /// The values being read, outermost first (`Open`): NCBI's recursion of `ReadObject` as a list
    /// on the heap, so that a value nested to any depth is read in LOSAT's own stack.
    open: Vec<Open>,
    /// Whether the warning of a UTF8String tag was posted (`Utf8TagOnce`; one stream is all a
    /// run of LOSAT reads with it).
    utf8_once: Utf8TagOnce,
}

/// A value being read that waits for a value inside it (`TypedReader::value`).
enum Open {
    /// The end of `ReadObject` of the type `target`: the hooks of its class (not after a skip).
    Value { target: u16, skip: bool },
    /// A SEQUENCE reading the member at `reading` (`members`), inside the member's frame, the
    /// `frames`-th.
    Sequence {
        ty: u16,
        members: &'static [Member],
        skip: bool,
        values: Vec<(u16, Value)>,
        reading: usize,
        frames: usize,
    },
    /// A CHOICE reading its variant.
    Choice { index: u16 },
    /// A SEQUENCE OF or SET OF reading an element.
    List {
        element: u16,
        skip: bool,
        values: Vec<Value>,
    },
    /// A delayed member (`delayed`) from `start`.
    Delayed { start: usize },
}

/// What the reader does next: read a value of a type (`skip`), or give a value read to the one
/// that holds it.
enum Step {
    Read(u16, bool),
    Done(Value),
}

impl<I: Input> TypedReader<I> {
    fn new(inp: I) -> Self {
        TypedReader {
            inp,
            frames: Vec::new(),
            posts: Vec::new(),
            open: Vec::new(),
            utf8_once: Utf8TagOnce::default(),
        }
    }

    /// `BEGIN_OBJECT_FRAME` of a frame of `kind`.
    fn push(&mut self, kind: FrameKind, type_name: Option<&'static str>, skip: bool) {
        self.frames.push(Frame {
            kind,
            type_name,
            member: None,
            skip,
            ty: None,
            num_rows: None,
        });
        self.inp.set_stack_depth(self.frames.len());
    }

    /// `push` of a class or choice frame of the type `ty`.
    fn push_object(&mut self, kind: FrameKind, ty: u16, type_name: &'static str, skip: bool) {
        self.push(kind, Some(type_name), skip);
        if let Some(frame) = self.frames.last_mut() {
            frame.ty = Some(ty);
        }
    }

    /// `END_OBJECT_FRAME` of the top frame after what it holds was read.
    fn pop_frame(&mut self) -> Frame {
        let frame = self.frames.pop().expect("a frame to leave");
        self.inp.set_stack_depth(self.frames.len());
        frame
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
    fn pop<T>(&mut self, result: R<T>) -> R<T> {
        let frame = self.pop_frame();
        result.map_err(|fail| self.leave(frame, fail))
    }

    /// What leaving `frame` (the stack is below it now) adds to `fail`.
    #[cold]
    #[inline(never)]
    fn leave(&self, frame: Frame, fail: Fail) -> Fail {
        match fail {
            Fail::Serial(mut e) => {
                e.add_frame_info(frame_name(frame.kind, frame.type_name, frame.member));
                Fail::Serial(e)
            }
            Fail::Eof if self.frames.len() >= 2 => Fail::serial(NcbiException::new(
                HANDLE_EOF,
                format!("{}: {}", self.inp.position(), frame.info()),
            )),
            Fail::Exception(e) => {
                let (file, line, function) = frame.site();
                Fail::exception(e.rethrow_same(file, line, function, frame.info()))
            }
            other => other,
        }
    }

    /// The frames above the `depth`-th left (`leave`) with `error`.
    fn unwind(&mut self, depth: usize, mut error: Fail) -> Fail {
        while self.frames.len() > depth {
            let frame = self.pop_frame();
            error = self.leave(frame, error);
        }
        error
    }

    fn set_member(&mut self, name: &'static str) {
        if let Some(top) = self.frames.last_mut() {
            top.member = Some(name);
        }
    }

    /// `GetStackTrace()`: the path of the object being read (`PssmWithParameters.pssm.query.seq.
    /// descr.E.comment`), by the one writer of NCBI's object stack path.
    fn stack_trace(&self) -> String {
        object_stack_trace(
            self.frames
                .first()
                .and_then(|frame| frame.type_name)
                .unwrap_or("?"),
            self.frames
                .iter()
                .skip(1)
                .filter_map(|frame| match frame.kind {
                    FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                        frame.member.map(StackStep::Member)
                    }
                    FrameKind::ArrayElement => Some(StackStep::Element),
                    _ => None,
                }),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:882-901
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:959-969
    /// `stream >> object` of the root type `root`: `SkipFileHeader`, then `Read`, each in a named
    /// frame of its own.
    fn document(&mut self, root: u16) -> R<Value> {
        let name = super::type_info(root).name();
        self.push(FrameKind::Named, Some(name), false);
        let header = self.file_header(name);
        self.pop(header)?;
        self.push(FrameKind::Named, Some(name), false);
        let value = self.value(root, false);
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

    /// One value of the type `ty` (`ReadObject`), or its skip (`SkipObject`, `skip`), with the hooks
    /// of NCBI's classes that a read runs. The values inside it are read in turn from the list of
    /// the values being read (`open`), as NCBI's functions recurse; an error leaves NCBI's frames
    /// from the innermost out (`pop`), as it passes `END_OBJECT_FRAME` of each.
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
                    // From the innermost value out: the frames inside each member being read, then
                    // the member's `ReadWithSetFlagMember` (`member_error`).
                    while self.open.len() > open {
                        if let Some(Open::Sequence {
                            members,
                            skip,
                            reading,
                            frames: depth,
                            ..
                        }) = self.open.pop()
                        {
                            error = self.unwind(depth, error);
                            if !skip {
                                error = member_error(error, &members[reading]);
                            }
                        }
                    }
                    error = self.unwind(frames, error);
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
        self.open.push(Open::Value { target, skip });
        let value = match info.kind {
            Kind::Sequence(members) => {
                self.push_object(FrameKind::Class, target, info.name(), skip);
                self.inp.begin_class()?;
                self.push(FrameKind::ClassMember, None, skip);
                return self.next_member(target, members.members(), skip, Vec::new(), 0);
            }
            Kind::Choice(variants) => {
                return self.start_choice(target, info.name(), variants.members(), skip)
            }
            Kind::List { element, set } => return self.start_list(element, set, skip),
            Kind::Alias(_) => unreachable!("resolved"),
            Kind::Integer(storage) => self.integer(storage, skip)?,
            Kind::Enumerated { integer, values } => {
                self.enumerated(values.values(), integer, skip)?
            }
            // `SkipBool` and `SkipNull` read what `ReadBool` and `ReadNull` read.
            Kind::Boolean => Value::Bool(self.inp.read_bool()?),
            Kind::Null => {
                self.inp.read_null()?;
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
                super::hooks::post_read(target, value, &mut self.posts).map(Step::Done)
            }
            Open::Sequence {
                ty,
                members,
                skip,
                mut values,
                reading,
                frames: _,
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
                            self.frames[at].num_rows = Some(*rows);
                        }
                    }
                    values.push((reading as u16, value));
                }
                self.next_member(ty, members, skip, values, reading + 1)
            }
            Open::Choice { index } => {
                self.inp.end_variant()?;
                self.pop_frame();
                self.inp.end_choice()?;
                self.pop_frame();
                Ok(Step::Done(Value::Choice(index, Box::new(value))))
            }
            Open::List {
                element,
                skip,
                mut values,
            } => {
                if !skip {
                    values.push(value);
                }
                self.next_element(element, skip, values)
            }
            Open::Delayed { start } => Ok(Step::Done(Value::Delayed {
                encoding: I::ENCODING,
                bytes: self.inp.delayed(start),
            })),
        }
    }

    /// `ReadClassSequential` / `SkipClassSequential` (`start`: the class frame, `BeginClass`, the
    /// member frame): the next member in module order (`BeginClassMember` of the class type and `pos`) with
    /// the mandatory members before it checked, and a read's (not a skip's) member hooks before
    /// its value (`hooks::pre_read_member`); after the last, the member frame left, the mandatory
    /// members after it checked, `EndClass`, the class frame left.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:84-109
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397 (the binary stream's,
    // objistrasnb.cpp:1210-1233, is the same with its own `BeginClass`)
    fn next_member(
        &mut self,
        ty: u16,
        members: &'static [Member],
        skip: bool,
        values: Vec<(u16, Value)>,
        pos: usize,
    ) -> R<Step> {
        if let Some(index) = self.inp.begin_member(members, pos)? {
            let member = &members[index];
            self.set_member(member.name());
            for missing in &members[pos..index] {
                self.expected_member(missing)?;
            }
            if !skip {
                super::hooks::pre_read_member(ty, index, &values)?;
            }
            let delayed = member.flags & super::DELAYED != 0 && !skip;
            let start = if delayed {
                Some(self.inp.delay_start()?)
            } else {
                None
            };
            self.open.push(Open::Sequence {
                ty,
                members,
                skip,
                values,
                reading: index,
                frames: self.frames.len(),
            });
            if let Some(start) = start {
                self.open.push(Open::Delayed { start });
                return Ok(Step::Read(member.ty, true));
            }
            return Ok(Step::Read(member.ty, skip));
        }
        self.pop_frame();
        for member in &members[pos..] {
            self.expected_member(member)?;
        }
        self.inp.end_class()?;
        self.pop_frame();
        Ok(Step::Done(Value::Sequence(values)))
    }

    /// `ReadMissingSimpleMember` -> `ExpectedMember` of a mandatory member: "member X expected"
    /// (`CSerialException`, NCBI's default verification).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
    fn expected_member(&self, member: &Member) -> R<()> {
        if member.optional() {
            return Ok(());
        }
        Err(throw_error(
            EXPECTED_MEMBER,
            self.inp.position(),
            format!("member {} expected", member.name()),
        ))
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
        self.inp.begin_choice()?;
        self.push(FrameKind::ChoiceVariant, None, skip);
        let index = self.inp.begin_variant(variants)?;
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

    /// `ReadContainer` / `SkipContainer` of a SEQUENCE OF or (`set`) SET OF.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274
    // The container's type information is anonymous (an implicit class such as `Seq-descr` has
    // no frame of its own, `ReadImplicitMember`, classinfo.cpp:289-306): `Object type= ` names
    // nothing (c07 of SURVEY_D3).
    fn start_list(&mut self, element: u16, set: bool, skip: bool) -> R<Step> {
        self.push(FrameKind::Array, Some(""), skip);
        self.inp.begin_container(set)?;
        self.push(FrameKind::ArrayElement, None, skip);
        let numbers = match resolved(element).kind {
            // The elements of an INTEGER or REAL list are read as the numbers themselves (no
            // `Value` each: a PSSM holds thousands); they have no hook.
            Kind::Integer(storage @ (Storage::Int2 | Storage::Int4)) => {
                // `integer` checked the storage's range.
                Value::Int4s(self.integers(storage, skip, |value| value as i32)?)
            }
            Kind::Integer(storage) => Value::Ints(self.integers(storage, skip, |value| value)?),
            Kind::Real => Value::Reals(self.reals(skip)?),
            _ => return self.next_element(element, skip, Vec::new()),
        };
        self.pop_frame();
        self.inp.end_container()?;
        self.pop_frame();
        Ok(Step::Done(numbers))
    }

    /// The next element of a SEQUENCE OF or SET OF of values; after the last, the element frame
    /// left, the container's end, the array frame left.
    fn next_element(&mut self, element: u16, skip: bool, values: Vec<Value>) -> R<Step> {
        if self.inp.begin_element()? {
            self.open.push(Open::List {
                element,
                skip,
                values,
            });
            return Ok(Step::Read(element, skip));
        }
        self.pop_frame();
        self.inp.end_container()?;
        self.pop_frame();
        Ok(Step::Done(Value::List(values)))
    }

    /// The elements of an INTEGER list, each `store`d as the list keeps it.
    fn integers<T>(&mut self, storage: Storage, skip: bool, store: impl Fn(i64) -> T) -> R<Vec<T>> {
        let mut values = Vec::new();
        while self.inp.begin_element()? {
            if let Some(value) = self.integer_number(storage, skip)? {
                values.push(store(value));
            }
        }
        Ok(values)
    }

    /// The elements of a REAL list.
    fn reals(&mut self, skip: bool) -> R<Vec<f64>> {
        let mut values = Vec::new();
        while self.inp.begin_element()? {
            if let Some(value) = self.real_number(skip)? {
                values.push(value);
            }
        }
        Ok(values)
    }

    /// An INTEGER of storage `storage` (`ReadStd`), or its skip (`SkipSNumber`, `SkipUNumber`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1780-1808 (the skips of every size:
    // `SkipSNumber` for the signed ones, `SkipUNumber` for the unsigned ones)
    fn integer(&mut self, storage: Storage, skip: bool) -> R<Value> {
        Ok(self
            .integer_number(storage, skip)?
            .map_or(Value::Null, Value::Int))
    }

    /// `integer`'s number (`None` when skipped), for the lists that keep the numbers (no `Value`
    /// to make and drop for each element).
    fn integer_number(&mut self, storage: Storage, skip: bool) -> R<Option<i64>> {
        if skip {
            self.inp.skip_number(storage != Storage::Uint4)?;
            return Ok(None);
        }
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
        Ok(Some(value))
    }

    /// `CEnumeratedTypeInfo::ReadEnum` / `SkipEnum` over `ReadEnum` of the stream: an exception
    /// gets the line of the type's `NCBI_RETHROW_SAME` (an end of the input stays one).
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // (`SkipEnum`, 468-484, reads the value the same way and rethrows at 478.)
    fn enumerated(&mut self, values: &'static [EnumValue], integer: bool, skip: bool) -> R<Value> {
        let items = EnumType { integer, values };
        let result = self.inp.read_enum(&items);
        result.map(Value::Enum).map_err(|fail| {
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
        })
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
    fn string(&mut self, kind: StringKind, skip: bool) -> R<Vec<u8>> {
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
            self.posts.push(Post::UnendedByteBlock);
        }
        Ok(Value::Octets(value?))
    }
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:715-755
/// What reading `member` makes of `error`: `ReadWithSetFlagMember` of a member with a set flag
/// (`Member::has_set_flag`; not a delayed one, whose value is read later, `ReadLongMember`, and
/// not a skipped one, `SkipSimpleMember`: neither catches), `asn_input::read_with_set_flag_member`.
fn member_error(error: Fail, member: &Member) -> Fail {
    let set_flag = member.has_set_flag() && member.flags & super::DELAYED == 0;
    read_with_set_flag_member(error, set_flag, member.name())
}

/// The type `ty` refers to after its aliases (the index of `resolved(ty)`), for the hooks.
fn alias_target(mut ty: u16) -> u16 {
    while let Kind::Alias(target) = super::type_info(ty).kind {
        ty = target;
    }
    ty
}
