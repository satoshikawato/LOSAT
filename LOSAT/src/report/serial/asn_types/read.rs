//! The generic reader of ASN.1 text under the type table: NCBI's `CObjectIStreamAsn` with its
//! object stack, for any type of `table` (`read_document`), over the lexical layer of
//! `report::serial::asn_text_reader`.
//!
//! It reads as NCBI's reader reads:
//! - the members of a SEQUENCE in module order, named, or by their tag in brackets (`[3]`); an
//!   unknown or repeated member, an unknown variant and a missing mandatory member are NCBI's
//!   `CSerialException`;
//! - INTEGER by the C++ storage of its member (`GetInt4`, `GetUint4` without `-`, `GetInt8`,
//!   `ReadInt2`); ENUMERATED by name or by a number that has a name, INTEGER with named values by
//!   name or any number;
//! - a VisibleString (and StringStore) byte outside 0x20..0x7E as `#`, with NCBI's
//!   `ERR_POST_X` with subcode 7 (`Post::BadChar`), a UTF8String as it is;
//! - the delayed member (`Auth-list.names`) with the skip functions, kept as its text;
//! - the hooks of NCBI's hand-written classes that change what is read (`super::hooks` of this
//!   file: `CSeq_descr::PostRead`, `CVariation_ref::PostRead`), and the Seq-id into LOSAT's
//!   Seq-id model.
//!
//! Its errors are NCBI's classes (`ReadError`): a `CSerialException` (and an end of the text
//! inside an object, which NCBI's stack turns into one), which its callers report as theirs; a
//! `CUtilException` of a number (`BadNumber`, `NumberOverflow`) or a `CStringException` of a tag,
//! with the line NCBI's stack adds for each frame it leaves (`NCBI_RETHROW_SAME`, the report of
//! `CATCH_ALL`); and an end of the text directly in the root object, a `CEofException` that NCBI's
//! stack passes on.

use super::{resolved, EnumValue, Kind, Member, Storage, StringKind, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_text_reader::{
    real_of_parts, AsnSource, AsnTextError, AsnTextErrorKind, MemberId, Reader, SliceSource,
    BAD_NUMBER, NUMBER_OVERFLOW,
};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::report::serial::object_frame::{frame_info, text_frame_site, FrameKind};

/// A message that NCBI's reader posts while it reads, in the order posted.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Post {
    /// `ReplaceVisibleChar`'s `ERR_POST_X` of subcode 7 with a message, for a byte read as `#`.
    BadChar(ReplacedChar),
    /// `ERR_POST` of an `Error` with a text, from a class's `PostRead`
    /// (`CVariation_ref::PostRead`).
    Error(&'static str),
}

impl Post {
    /// The line that a BLAST application of NCBI posts on its standard error, with its
    /// diagnostic prefix `program` (the program's name in lower case): the one owner of each
    /// post's text (`ReplacedChar::post`, `ncbi_exception::error_post_line`).
    pub(crate) fn line(&self, program: &str) -> Vec<u8> {
        match self {
            Post::BadChar(replaced) => replaced.post(Some(program)),
            Post::Error(text) => {
                crate::ncbi_exception::error_post_line(None, Some(program), text.as_bytes())
            }
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
    /// A `CEofException` that leaves the reader: the text ends directly in the root object.
    Eof,
    /// A `std::exception` that leaves the reader (its `what()`): `std::length_error` of a reserve
    /// hook (`hooks::pre_read_member`).
    Std(&'static str),
    /// A form of the text that LOSAT's lexical reader does not take (`AsnTextErrorKind::Unsupported`,
    /// LOSAT's own limit, not an NCBI error): the caller refuses the input with this text.
    Unsupported(String),
}

/// What the reader read, and the messages it posted before it stopped or finished.
pub(crate) type Read<T> = (std::result::Result<T, ReadError>, Vec<Post>);

type R<T> = std::result::Result<T, ReadError>;

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
}

impl Frame {
    /// `GetFrameInfo()` (`object_frame::frame_info`).
    fn info(&self) -> String {
        frame_info(self.kind, self.type_name, self.member)
    }

    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves this frame
    /// (`object_frame::text_frame_site`).
    fn site(&self) -> (&'static str, u32, &'static str) {
        text_frame_site(self.kind, self.skip)
    }
}

/// Reads the ASN.1 text document at the start of `text` as the root type `root` (an object of
/// `CObjectIStream::Read`: its file header, then the value; what follows it is not read), with
/// the messages posted.
pub(crate) fn read_document(text: &[u8], root: u16) -> Read<Value> {
    let mut reader = TypedReader {
        r: Reader::new(SliceSource::new(text), "PSSM", "PSSM reader"),
        frames: Vec::new(),
        posts: Vec::new(),
        open: Vec::new(),
    };
    let value = reader.document(root);
    (value, reader.posts)
}

struct TypedReader<'a> {
    r: Reader<SliceSource<'a>>,
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
    /// A delayed member (`delayed`) from `start`.
    Delayed { start: usize },
}

/// What the reader does next: read a value of a type (`skip`), or give a value read to the one
/// that holds it.
enum Step {
    Read(u16, bool),
    Done(Value),
}

impl TypedReader<'_> {
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
    /// (`Frame::site`): a `CException` other than a `CSerialException` takes the frame's line; an
    /// end of the text (`CEofException`, `HandleEOF`, objistr.cpp:599-608) becomes a
    /// `CSerialException` unless fewer than 2 frames stay.
    fn pop<T>(&mut self, result: R<T>) -> R<T> {
        let frame = self.frames.pop().expect("a frame to leave");
        result.map_err(|error| match error {
            ReadError::Exception(e) => {
                let (file, line, function) = frame.site();
                ReadError::Exception(e.rethrow_same(file, line, function, frame.info()))
            }
            ReadError::Eof if self.frames.len() >= 2 => {
                ReadError::Serial(format!("unexpected end of the text:{}", frame.info()))
            }
            other => other,
        })
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

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-959 and 959-969
    fn document(&mut self, root: u16) -> R<Value> {
        let name = super::type_info(root).name();
        self.push(FrameKind::Named, Some(name), false);
        let header = self
            .r
            .file_header()
            .map_err(ReadError::from)
            .and_then(|id| {
                let id = self.r.bytes(id);
                if !id.is_empty() && id != name.as_bytes() {
                    return Err(ReadError::Serial(format!(
                        "incompatible type {}<>{name}",
                        String::from_utf8_lossy(id)
                    )));
                }
                Ok(())
            });
        self.pop(header)?;
        self.push(FrameKind::Named, Some(name), false);
        let value = self.value(root, false);
        self.pop(value)
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
        self.open.push(Open::Value { target, skip });
        let value = match info.kind {
            Kind::Sequence(members) => {
                self.push_object(FrameKind::Class, target, info.name(), skip);
                self.r.start_block()?;
                self.push(FrameKind::ClassMember, None, skip);
                return self.next_member(target, members.members(), skip, Vec::new(), 0, true);
            }
            Kind::Choice(variants) => {
                return self.start_choice(target, info.name(), variants.members(), skip)
            }
            Kind::List { element, .. } => return self.start_list(element, skip),
            Kind::Alias(_) => unreachable!("resolved"),
            Kind::Integer(storage) => self.integer(storage, skip)?,
            Kind::Enumerated { integer, values } => {
                self.enumerated(values.values(), integer, skip)?
            }
            Kind::Boolean => Value::Bool(self.r.boolean()?),
            Kind::Null => {
                self.r.null()?;
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
                block_start,
            } => {
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
                self.frames.pop();
                self.frames.pop();
                Ok(Step::Done(Value::Choice(index, Box::new(value))))
            }
            Open::List {
                element,
                skip,
                mut values,
                block_start,
            } => {
                if !skip {
                    values.push(value);
                }
                self.next_element(element, skip, values, block_start)
            }
            Open::Delayed { start } => Ok(Step::Done(Value::Delayed(
                self.r.bytes(start..self.r.pos).to_vec(),
            ))),
        }
    }

    /// `ReadClassSequential` / `SkipClassSequential` (`start`: the class frame, `{`, the member
    /// frame): the next member in module order (`BeginClassMember` of `classType` at `pos`) with
    /// the mandatory members before it checked, and a read's (not a skip's) member hooks before its
    /// value (`hooks::pre_read_member`); after the last, the member frame left, the mandatory
    /// members after it checked, `}`, the class frame left.
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
        if self.r.next_element(&mut block_start)? {
            let index = self.item_index(members, pos, true)?;
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
                Some(self.delay_start()?)
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
                self.open.push(Open::Delayed { start });
                return Ok(Step::Read(member.ty, true));
            }
            return Ok(Step::Read(member.ty, skip));
        }
        self.frames.pop();
        for member in &members[pos..] {
            expected_member(member)?;
        }
        self.r.end_block()?;
        self.frames.pop();
        Ok(Step::Done(Value::Sequence(values)))
    }

    /// `ReadMemberId` after `SkipWhiteSpace` and `GetMemberIndex` of `classType`, `id` and `pos`
    /// (from `pos`, `in_class`) or `GetChoiceIndex` of `choiceType` and `id`: the member or variant
    /// named or tagged.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:407-424
    // `GetAltItemIndex` (366-387) tries the name with its first letter in upper case and takes it
    // only for an item without a prefix, which no member of an ASN.1 module is (its names start
    // with a lower-case letter): an upper-case first letter is never a member here (c21 of
    // SURVEY_D3: `Biomol`).
    fn item_index(&mut self, items: &'static [Member], pos: usize, in_class: bool) -> R<usize> {
        let found = match self.r.member_or_tag()? {
            MemberId::Name(range) => {
                let id = self.r.bytes(range);
                if id.is_empty() && !in_class {
                    return Err(ReadError::Serial("choice variant id expected".into()));
                }
                items[pos..]
                    .iter()
                    .position(|item| item.name().as_bytes() == id)
                    .map(|at| at + pos)
            }
            MemberId::Tag(range) => {
                // Parse the tag with `NStr::StringToInt` using default flags: its `CStringException` leaves
                // the reader (NCBI 2.17.0: `[3x]`, `[99999999999]`, probe e03, e07).
                let tag = crate::ncbi_str::string_to_int(self.r.bytes(range), 0)
                    .map_err(ReadError::Exception)?;
                items[pos..]
                    .iter()
                    .position(|item| i32::from(item.tag) == tag)
                    .map(|at| at + pos)
            }
        };
        found.ok_or_else(|| ReadError::Serial("unexpected member or variant".into()))
    }

    /// The delayed member: `StartDelayBuffer` (white space skipped), then the value skipped
    /// (`Open::Delayed`), the text kept (`Value::Delayed`).
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1648-1652
    // (psiblast and tblastn set no hook on their stream: `ShouldParseDelayBuffer` is false,
    // objistr.cpp:745-760.)
    fn delay_start(&mut self) -> R<usize> {
        self.r.skip_white_space()?;
        Ok(self.r.pos)
    }

    /// `ReadChoiceSimple` / `SkipChoiceSimple`: the variant of the CHOICE `ty`; a read (not a skip)
    /// runs the class's variant hooks before its value (`hooks::pre_read_variant`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1501
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1475 (`BeginChoiceVariant`,
    // cited at `Reader::choice`: `NextElement` at the choice's start, then the variant's id)
    fn start_choice(
        &mut self,
        ty: u16,
        name: &'static str,
        variants: &'static [Member],
        skip: bool,
    ) -> R<Step> {
        self.push_object(FrameKind::Choice, ty, name, skip);
        self.push(FrameKind::ChoiceVariant, None, skip);
        let mut block_start = true;
        if !self.r.next_element(&mut block_start)? {
            return Err(ReadError::Serial("choice variant id expected".into()));
        }
        let index = self.item_index(variants, 0, false)?;
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
    fn start_list(&mut self, element: u16, skip: bool) -> R<Step> {
        self.push(FrameKind::Array, Some(""), skip);
        self.r.start_block()?;
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
            _ => return self.next_element(element, skip, Vec::new(), true),
        };
        self.frames.pop();
        self.r.end_block()?;
        self.frames.pop();
        Ok(Step::Done(numbers))
    }

    /// The next element of a SEQUENCE OF or SET OF of values; after the last, the element frame
    /// left, `}`, the array frame left.
    fn next_element(
        &mut self,
        element: u16,
        skip: bool,
        values: Vec<Value>,
        mut block_start: bool,
    ) -> R<Step> {
        if self.r.next_element(&mut block_start)? {
            self.open.push(Open::List {
                element,
                skip,
                values,
                block_start,
            });
            return Ok(Step::Read(element, skip));
        }
        self.frames.pop();
        self.r.end_block()?;
        self.frames.pop();
        Ok(Step::Done(Value::List(values)))
    }

    /// `elements` of an INTEGER list, each `store`d as the list keeps it.
    fn integers<T>(&mut self, storage: Storage, skip: bool, store: impl Fn(i64) -> T) -> R<Vec<T>> {
        let mut values = Vec::new();
        let mut block_start = true;
        while self.r.next_element(&mut block_start)? {
            if let Some(value) = self.integer_number(storage, skip)? {
                values.push(store(value));
            }
        }
        Ok(values)
    }

    /// `elements` of a REAL list.
    fn reals(&mut self, skip: bool) -> R<Vec<f64>> {
        let mut values = Vec::new();
        let mut block_start = true;
        while self.r.next_element(&mut block_start)? {
            if let Some(value) = self.real_number(skip)? {
                values.push(value);
            }
        }
        Ok(values)
    }

    fn util_error(&self, site: ThrowSite, what: &str) -> ReadError {
        ReadError::Exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.r.line),
        ))
    }

    /// `GetInt4`, `GetUint4` (`signed` false: no `-`) and `GetInt8` after `SkipWhiteSpace`, the
    /// digits checked against `max` as NCBI checks them.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-668
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let r = &mut self.r;
        let (pos, c) = r.white_space_at(r.pos)?;
        // `SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
        let (pos, mut c) = r.ws_at(pos, c)?;
        let mut pos = pos + 1;
        let mut negative = false;
        if signed && c == b'-' {
            negative = true;
            c = r.char_at(pos)?;
            pos += 1;
        } else if c == b'+' {
            c = r.char_at(pos)?;
            pos += 1;
        }
        if !c.is_ascii_digit() {
            r.pos = pos;
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        let mut overflow = false;
        let (end, _) = r.src.scan(pos, |d| {
            if !d.is_ascii_digit() {
                return false;
            }
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                overflow = true;
                return false;
            }
            n = n * 10 + d;
            true
        });
        if overflow {
            // `GetInt4` takes the digit before it checks the overflow.
            r.pos = end + 1;
            return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
        }
        r.pos = end;
        // `-Int4(n)` of n = 2^31 is the least Int4 (and likewise in 64 bits).
        Ok(if negative {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    fn int4(&mut self) -> R<i32> {
        Ok(self.get_number(i32::MAX as u64, true)? as i32)
    }

    fn uint4(&mut self) -> R<u32> {
        Ok(self.get_number(u64::from(u32::MAX), false)? as u32)
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
            self.skip_number(storage != Storage::Uint4)?;
            return Ok(None);
        }
        let value = match storage {
            Storage::Int4 => i64::from(self.int4()?),
            Storage::Uint4 => i64::from(self.uint4()?),
            Storage::Int8 => self.get_number(i64::MAX as u64, true)?,
            Storage::Int2 => {
                let value = self.int4()?;
                if i16::try_from(value).is_err() {
                    return Err(ReadError::Serial("integer overflow".into()));
                }
                i64::from(value)
            }
        };
        Ok(Some(value))
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1048
    // (`SkipUNumber`, 1050-1074, takes only `+` and peeks past the digits with `PeekCharNoEOF`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let r = &mut self.r;
        let (start, mut c) = r.white_space_at(r.pos)?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = match r.src.byte_at(start + 1) {
                Some(c) => c,
                None => return Err(Self::end_of_number(r, start)),
            };
            i = 2;
        }
        if !c.is_ascii_digit() {
            r.pos = start;
            return Err(ReadError::Serial(format!(
                "bad {} integer in line {}",
                if signed { "signed" } else { "unsigned" },
                r.line
            )));
        }
        let (end, c) = r.src.scan(start + i, |c| c.is_ascii_digit());
        // `SkipSNumber` peeks past the digits with `PeekChar`, which fails at the end of the
        // text; `SkipUNumber` with `PeekCharNoEOF`.
        if signed && c.is_none() {
            return Err(Self::end_of_number(r, start));
        }
        r.pos = end;
        Ok(())
    }

    /// The end of the text inside a number that starts at `start`: `PeekChar`'s, at `start`.
    #[cold]
    fn end_of_number(r: &mut Reader<SliceSource<'_>>, start: usize) -> ReadError {
        r.pos = start;
        r.end_of_text().into()
    }

    /// `CEnumeratedTypeInfo::ReadEnum` / `SkipEnum` over `ReadEnum` of the stream: a lower-case
    /// name of the list, or a number (one of the list for an ENUMERATED); an exception that is not
    /// a `CSerialException` gets the line of the type's `NCBI_RETHROW_SAME`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // (`SkipEnum`, 468-484, reads the value the same way and rethrows at 478.)
    fn enumerated(&mut self, values: &'static [EnumValue], integer: bool, skip: bool) -> R<Value> {
        let result = self.enum_value(values, integer);
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

    fn enum_value(&mut self, values: &'static [EnumValue], integer: bool) -> R<i32> {
        let c = self.r.skip_white_space()?;
        let id = self.r.scan_end_of_id(c.is_ascii_lowercase())?;
        if !id.is_empty() {
            let id = self.r.bytes(id);
            return values
                .iter()
                .find(|value| value.name().as_bytes() == id)
                .map(|value| value.value)
                .ok_or_else(|| ReadError::Serial("invalid value of enumerated type".into()));
        }
        let value = self.int4()?;
        if !integer && !values.iter().any(|named| named.value == value) {
            return Err(ReadError::Serial("invalid value of enumerated type".into()));
        }
        Ok(value)
    }

    /// `ReadDouble` / `SkipFNumber`: a word, or `{ mantissa, base, exponent }` whose base is read
    /// with `ReadUint4` (a `CUtilException` for a `-`) and exponent with `ReadInt4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (cited at `Reader::real`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
    fn real(&mut self, skip: bool) -> R<Value> {
        Ok(self.real_number(skip)?.map_or(Value::Null, Value::Real))
    }

    /// `real`'s number (`None` when skipped), for the lists that keep the numbers.
    fn real_number(&mut self, skip: bool) -> R<Option<f64>> {
        if skip {
            if self.r.skip_white_space()? != b'{' {
                self.r.scan_end_of_id(true)?;
                return Ok(None);
            }
            self.r.pos += 1;
            self.skip_number(true)?;
            self.r.expect(b',')?;
            let base = self.uint4()?;
            self.r.expect(b',')?;
            self.skip_number(true)?;
            self.r.expect(b'}')?;
            if base != 2 && base != 10 {
                return Err(ReadError::Serial(
                    "illegal REAL base (must be 2 or 10)".into(),
                ));
            }
            return Ok(None);
        }
        if let Some(word) = self.r.real_word()? {
            return Ok(Some(word));
        }
        let (mantissa, is_negative) = self.r.real_mantissa()?;
        self.r.expect(b',')?;
        let base = self.uint4()?;
        self.r.expect(b',')?;
        let exp = self.int4()?;
        self.r.expect(b'}')?;
        real_of_parts(mantissa, is_negative, base, exp)
            .map(Some)
            .ok_or_else(|| ReadError::Serial("illegal REAL base (must be 2 or 10)".into()))
    }

    /// `ReadStringValue` with the fix of `kind` / `SkipString` (`Reader::visible_string`): for a
    /// VisibleString each byte outside 0x20..0x7E read as `#`, with NCBI's message.
    fn string(&mut self, kind: StringKind, skip: bool) -> R<Vec<u8>> {
        let mut out = Vec::new();
        let mut replaced = Vec::new();
        let result = self.r.visible_string(
            &mut out,
            kind != StringKind::Utf8,
            skip,
            &mut |line, byte, piece| replaced.push((line, byte, piece.to_vec())),
        );
        for (line, byte, piece) in replaced {
            self.post_bad_char(line, byte, &piece);
        }
        result?;
        Ok(out)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (cited at
    // `report::serial::encoder::ReplacedChar`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:807-810
    // (`CObjectIStreamAsn` keeps `GetPosition` of its base; `CIStreamBuffer::GetLine` is the line
    // of the reader.)
    fn post_bad_char(&mut self, line: usize, byte: u8, string: &[u8]) {
        self.posts.push(Post::BadChar(ReplacedChar {
            trace: self.stack_trace(),
            byte,
            position: StreamPosition::Line(line as u64),
            string: string.into(),
        }));
    }

    /// `ReadBytes` / `SkipByteBlock` (the skip takes `H` or `B` after white space).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
    fn octets(&mut self, skip: bool) -> R<Value> {
        if !skip {
            return Ok(Value::Octets(self.r.octets()?));
        }
        self.r.expect(b'\'')?;
        loop {
            let c = self.r.peek_char(0)?;
            self.r.pos += 1;
            match c {
                b'0'..=b'9' | b'A'..=b'F' | b'a'..=b'f' => {}
                b'\'' => break,
                b'\r' | b'\n' => self.r.skip_end_of_line(c),
                _ => {
                    self.r.pos -= 1;
                    return Err(ReadError::Serial(format!(
                        "bad char in octet string: #{}",
                        c as i8
                    )));
                }
            }
        }
        match self.r.skip_white_space()? {
            b'H' | b'B' => {
                self.r.pos += 1;
                Ok(Value::Null)
            }
            _ => Err(ReadError::Serial("'H' or 'B' expected".into())),
        }
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
fn alias_target(mut ty: u16) -> u16 {
    while let Kind::Alias(target) = super::type_info(ty).kind {
        ty = target;
    }
    ty
}
