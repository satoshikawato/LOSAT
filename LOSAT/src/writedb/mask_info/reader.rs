//! The reader of makeblastdb's mask files (`makeblastdb.cpp` `s_ReadObject` and
//! `x_ProcessMaskData`): the format guess, then `Blast-db-mask-info` objects one after another
//! in ASN.1 text or binary.
//!
//! NCBI reads an object with its serial library (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`,
//! driven by the type information of the generated classes). The two encodings share their
//! structure: this module is the part of `CObjectIStream` that is common to both
//! (`ReadClassSequential`, `ReadChoiceSimple`, `ReadContainer`, the object frames); the two
//! encodings' primitives are `report::serial::asn_input::TextInput` and
//! `report::serial::asn_ber_reader::BerReader`, which the reader of the type table
//! (`report::serial::asn_types::read`) shares.
//!
//! What NCBI's reader throws for a file it cannot read is what LOSAT reports (`Fail`, then
//! `MaskFileError::Ncbi`): the exception of the primitive at its place in NCBI's source, with
//! what each object frame it leaves adds (`END_OBJECT_FRAME`: the frame's name for a
//! `CSerialException`, its line for another `CException`; `ReadWithSetFlagMember` and
//! `CEnumeratedTypeInfo::ReadEnum` add theirs), as makeblastdb's `CATCH_ALL` prints it. What it
//! posts while it reads (a byte of a VisibleString read as `#`, a UTF8String tag) is handed out
//! in order (`MaskFileReader::take_posts`).
//!
//! The end of the file is NCBI's `CEofException` (`Fail::Eof`). Every object frame turns it into
//! a `CSerialException` (`HandleEOF`) unless fewer than two frames remain after it, in which case
//! the `CEofException` goes on to `x_ProcessMaskData`, which treats it as the end of the file and
//! drops the object being read. That happens at the type name of a text object, at its opening
//! brace, and at the tag, the length and the closing bytes of the outermost binary SEQUENCE.

use std::fmt;
use std::io::Read;
use std::path::{Path, PathBuf};

use super::guess::guess_format;
use super::schema::{
    ChoiceDef, ClassDef, EnumDef, Member, Ty, BLAST_DB_MASK_INFO, BLAST_FILTER_PROGRAM,
    BLAST_MASK_LIST, DATE, DATE_STD, DBTAG, GIIMPORT_ID, ID_PAT, ID_PAT_ID, NA_STRAND,
    OBJECT_ID_ID, PACKED_SEQPNT, PATENT_SEQ_ID, PDB_SEQ_ID, SEQ_BOND, SEQ_ID, SEQ_INTERVAL,
    SEQ_INTERVAL_FROM, SEQ_INTERVAL_ID, SEQ_INTERVAL_STRAND, SEQ_INTERVAL_TO, SEQ_LOC,
    SEQ_LOC_BOND, SEQ_LOC_EMPTY, SEQ_LOC_EQUIV, SEQ_LOC_INT, SEQ_LOC_MIX, SEQ_LOC_NULL,
    SEQ_LOC_PACKED_INT, SEQ_LOC_PACKED_PNT, SEQ_LOC_PNT, SEQ_LOC_WHOLE, SEQ_POINT, TEXTSEQ_ID,
    T_OBJECT_ID,
};
use super::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskPart, NullElement};
use crate::ncbi_exception::{NcbiException, Severity};
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{
    throw_error, utf8_tag_message, Fail, Input, Names, StringKind, StringPost, TextInput, R,
};
use crate::report::serial::asn_text_reader::{PointerType, StreamSource};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep};
use crate::report::serial::object_frame::{
    frame_info, frame_name, text_frame_site, FrameKind, EXPECTED_MEMBER, GET_REGISTERED_OBJECT,
    HANDLE_EOF, READ_POINTER_FRAME_END, SKIP_FILE_HEADER,
};
use crate::seqid::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};
use crate::writedb::WriteDbError;

/// The encoding of a mask file, as `CFormatGuess` finds it. Anything else is an error when the
/// file is opened ("Unknown encoding").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MaskFileFormat {
    /// `eBinaryASN`: `MSerial_AsnBinary`.
    BinaryAsn,
    /// `eTextASN`: `MSerial_AsnText`.
    TextAsn,
}

/// Why a mask file could not be read.
#[derive(Debug)]
pub(crate) enum MaskFileError {
    /// What makeblastdb ends with, as NCBI's code writes it: the "Unknown encoding ..."
    /// `CInvalidDataException` of the format guess, and what NCBI's serial reader throws for a
    /// file it cannot read: its `CSerialException` or `CUtilException` (`Exception`, "Error: "
    /// and the report, status 255), `std::bad_alloc` of a string's length (`OutOfMemory`,
    /// status 4) or `std::length_error` (`Other` with its `what()`, status 255).
    Ncbi(WriteDbError),
    /// NCBI reads it; LOSAT does not.
    Unsupported(String),
    /// The file ends where `next_list` expected the continuation that `more TRUE` announces.
    /// `x_ProcessMaskData` has no handler for the `CEofException` that this raises outside
    /// of its first `try` block: `CATCH_ALL` ignores a `CIOException` that is not a flush
    /// failure, so makeblastdb leaves with status 0 and no database, after the
    /// "No sequences matched any of the masks provided." message of its build object.
    EofInContinuation { path: PathBuf },
}

impl fmt::Display for MaskFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MaskFileError::Ncbi(e) => write!(f, "{e}"),
            MaskFileError::Unsupported(what) => {
                write!(f, "mask data that LOSAT does not read: {what}")
            }
            MaskFileError::EofInContinuation { path } => write!(
                f,
                "mask data in '{}': the file ends where the continuation of 'more TRUE' was expected",
                path.display()
            ),
        }
    }
}

impl std::error::Error for MaskFileError {}

/// A message that NCBI's reader posts while it reads a mask file.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum MaskPost {
    /// The post 803.7 of a byte read as `#`.
    BadChar(ReplacedChar),
    /// The warning 803.10 of a UTF8String tag for the VisibleString member at the stack trace
    /// (makeblastdb writes the first of its run, `asn_input::Utf8TagOnce`).
    Utf8Tag(String),
}

impl MaskPost {
    /// The line makeblastdb writes on its standard error (`ReplacedChar::post`,
    /// `ncbi_exception::post_line`).
    pub(crate) fn line(&self) -> Vec<u8> {
        match self {
            MaskPost::BadChar(replaced) => replaced.post(Some("makeblastdb")),
            MaskPost::Utf8Tag(trace) => crate::ncbi_exception::post_line(
                Severity::Warning,
                Some("803.10"),
                Some("makeblastdb"),
                utf8_tag_message(trace).as_bytes(),
            ),
        }
    }

    /// Sets the stack trace that the post names (`Rd::read_str`).
    fn set_trace(&mut self, trace: String) {
        match self {
            MaskPost::BadChar(replaced) => replaced.trace = trace,
            MaskPost::Utf8Tag(text) => *text = trace,
        }
    }
}

/// A frame of NCBI's object stack (`object_frame`). The reader keeps only the depth of its stack
/// (`Rd::depth`): each frame is a value of the function that reads what it holds (`Rd::frame`),
/// which hands it to `Rd::exit` when it is left, for the error that leaves it (`Rd::leave`) and
/// the stack traces that wait for it (`Trace`).
#[derive(Clone, Copy, Debug)]
struct Frame {
    kind: FrameKind,
    /// The named frame of `ReadPointer`, whose `END_OBJECT_FRAME` is not the text reader's for
    /// its kind (`end`).
    read_pointer: bool,
    /// The type of a named, class, choice or container frame (`Object type= `).
    type_name: Option<&'static str>,
    /// The member or variant of a member or variant frame (`SetTopMemberId`), as the place of its
    /// name in the schema (`Member::name`, `Variant::name`): one word to store for each member.
    member: Option<&'static &'static str>,
}

impl Frame {
    /// A frame of `kind` of the type `type_name`.
    fn new(kind: FrameKind, type_name: Option<&'static str>) -> Frame {
        Frame {
            kind,
            read_pointer: false,
            type_name,
            member: None,
        }
    }

    /// The member or variant of a member or variant frame.
    fn member(self) -> Option<&'static str> {
        self.member.copied()
    }

    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves the frame: the
    /// text reader's for its kind (`text_frame_site`), `READ_POINTER_FRAME_END` for the named
    /// frame of `ReadPointer`.
    fn end(self) -> (&'static str, u32, &'static str) {
        if self.read_pointer {
            READ_POINTER_FRAME_END
        } else {
            text_frame_site(self.kind, false)
        }
    }
}

/// Where [`Rd::begin_pointer`] found the value of a `CRef<>`.
#[derive(Clone, Copy)]
enum Pointee {
    /// `NULL`: there is none.
    Null,
    /// `eThisPointer`: it follows.
    This,
    /// `eOtherPointer`: it follows inside `ReadPointer`'s named frame of the declared type,
    /// entered at this depth.
    Other(usize, &'static str),
}

/// The `GetStackTrace()` of the posts of one VisibleString (`Rd::read_str`), built as the frames
/// that were on the stack when it was read give their steps: a member or variant frame when its
/// member changes or it is left, any other frame when it is left; the frame at the bottom names
/// the type, and fills the trace into the posts.
struct Trace {
    /// The posts in `Rd::posts`.
    posts: std::ops::Range<usize>,
    /// The frames below this depth have not given their steps yet.
    level: usize,
    /// The steps given so far, from the top of the stack down.
    steps: Vec<StackStep<'static>>,
}

/// How far `CSeq_loc::CheckId` walks one Seq-loc (`x_CheckId`, cited at `read_seq_loc_parts`):
/// each Seq-loc starts a walk of its own (`my_id`, NULL), takes the Seq-ids of its parts in
/// order, and hands its own to the Seq-loc around it (`x_UpdateId` of `id` and `my_id`) once its
/// walk is done. makeblastdb's `GetId` of a mask is this walk of the mask's Seq-loc.
#[derive(Debug)]
enum IdWalk {
    /// Every Seq-id met so far is `Equals` to this one (`None`: none yet).
    Same(Option<ParsedSeqId>),
    /// `x_CheckId` returned false: two Seq-ids differ, or a `feat` came. Nothing after it is
    /// walked, and `GetId` returns NULL.
    Stopped,
    /// The walk reached an element that is an empty `CRef<>`: `CCoreException` `eNullPtr`.
    Null(NullElement),
}

impl IdWalk {
    /// `x_UpdateId` of `id` and the address of `seq_id`, for a part.
    fn update(&mut self, seq_id: ParsedSeqId) {
        match self {
            IdWalk::Same(None) => *self = IdWalk::Same(Some(seq_id)),
            // `Equals` compares the values: the typed Seq-ids are equal exactly when it does.
            IdWalk::Same(Some(first)) if *first != seq_id => *self = IdWalk::Stopped,
            _ => {}
        }
    }

    /// `CheckId` of `id` and `may_throw`, called on an element (`*li`) of a `mix` or `equiv`:
    /// the element's own walk `inner`, if this walk gets that far.
    fn nested(&mut self, inner: IdWalk) {
        if !matches!(self, IdWalk::Same(_)) {
            return;
        }
        match inner {
            IdWalk::Same(Some(seq_id)) => self.update(seq_id),
            // `x_UpdateId` of `id` and `NULL` keeps `id`.
            IdWalk::Same(None) => {}
            stopped_or_null => *self = stopped_or_null,
        }
    }

    /// An element that is an empty `CRef<>`, dereferenced if this walk gets that far.
    fn null(&mut self, null: NullElement) {
        if matches!(self, IdWalk::Same(_)) {
            *self = IdWalk::Null(null);
        }
    }
}

/// The lists of `CRef`s that a mask file holds, as `NullElement::list` names them.
const LIST_MASKS: &str = "the masks of a Blast-mask-list";
const LIST_MIX: &str = "a Seq-loc-mix (mix)";
const LIST_EQUIV: &str = "a Seq-loc-equiv (equiv)";
const LIST_PACKED_INT: &str = "a Packed-seqint (packed-int)";

/// A frame that the Seq-loc reader (`Rd::read_seq_loc_parts`) entered and has not left. The
/// Seq-loc is the one type of a mask file that holds values of its own type (`mix`, `equiv`): the
/// frames of a nest of them are kept on a list on the heap (`Rd::loc_frames`), two bytes each,
/// not on the call stack, so that a Seq-loc nested to any depth is read without recursion.
#[derive(Clone, Copy, Debug)]
enum LocFrame {
    /// The choice frame of a Seq-loc (`read_choice`),
    Choice,
    /// its variant frame, with the variant once `begin_variant` has read it,
    Variant(Option<u8>),
    /// the container frame of a `mix` or `equiv` (`read_container`),
    Array,
    /// its element frame,
    Element,
    /// and the named frame of `ReadPointer` around an element stated as `:Seq-loc`
    /// (`begin_pointer`).
    Pointer,
}

impl LocFrame {
    /// The frame, as `read_choice`, `read_container` and `end_pointer` make it.
    fn frame(self) -> Frame {
        match self {
            LocFrame::Choice => Frame::new(FrameKind::Choice, Some(SEQ_LOC.name)),
            LocFrame::Variant(index) => Frame {
                member: index.map(|i| &SEQ_LOC.variants[usize::from(i)].name),
                ..Frame::new(FrameKind::ChoiceVariant, None)
            },
            LocFrame::Array => Frame::new(FrameKind::Array, Some("")),
            LocFrame::Element => Frame::new(FrameKind::ArrayElement, None),
            LocFrame::Pointer => Frame {
                read_pointer: true,
                ..Frame::new(FrameKind::Named, Some(SEQ_LOC.name))
            },
        }
    }
}

/// A `mix` or `equiv` that the Seq-loc reader is reading the elements of: the walk of its Seq-loc
/// (`None` while it is `IdWalk::Same(None)`, as it is until an element gives it a Seq-id, so that a
/// level of a deep nest holds one word) and which of the two it is.
struct LocLevel {
    walk: Option<Box<IdWalk>>,
    equiv: bool,
}

impl LocLevel {
    fn walk(&mut self) -> &mut IdWalk {
        self.walk
            .get_or_insert_with(|| Box::new(IdWalk::Same(None)))
    }

    /// [`IdWalk::nested`] of the walk of an element: one that met no Seq-id changes nothing.
    fn nested(&mut self, inner: IdWalk) {
        if !matches!(inner, IdWalk::Same(None)) {
            self.walk().nested(inner);
        }
    }

    /// The list, as [`NullElement::list`] names it.
    fn list(&self) -> &'static str {
        if self.equiv {
            LIST_EQUIV
        } else {
            LIST_MIX
        }
    }
}

/// What the Seq-loc reader does next: read a Seq-loc, the next element of the innermost `mix` or
/// `equiv` (or its end), or the end of a Seq-loc whose value it has read, with its walk.
enum LocStep {
    Begin,
    Element,
    End(IdWalk),
}

const IN_A_LEVEL: &str = "an element frame belongs to a mix or equiv";

/// The state of one `CObjectIStream`: the primitives, the object stack and what was posted.
pub(super) struct Rd<I: Input> {
    pub(super) inp: I,
    /// `GetStackDepth()`: the frames entered and not left.
    depth: usize,
    /// The stack traces that wait for frames still on the stack (empty but for a VisibleString
    /// that posted).
    traces: Vec<Trace>,
    /// The messages posted since the caller last took them (`MaskFileReader::take_posts`).
    pub(super) posts: Vec<MaskPost>,
    /// The lists of the Seq-loc reader (`read_seq_loc_parts`), kept empty between its calls, so
    /// that the masks of a file reuse their room.
    loc_frames: Vec<LocFrame>,
    loc_levels: Vec<LocLevel>,
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:712-756
// (`eNullValue` and `eMissingValue` are not thrown for these types: `ExpectedMember` throws
// `eFormatError` with the default `eSerialVerifyData_Yes`.)
const READ_WITH_SET_FLAG_MEMBER: (&str, u32, &str) = (
    "src/serial/member.cpp",
    752,
    "ncbi::CMemberInfoFunctions::ReadWithSetFlagMember",
);

/// What reading `member` makes of `fail`: `ReadWithSetFlagMember` rethrows a `CSerialException`
/// of a member with a set flag (`Member::set_flag`).
#[cold]
#[inline(never)]
fn member_error(fail: Fail, member: &Member) -> Fail {
    match fail {
        Fail::Serial(e) if member.set_flag => {
            let (file, line, function) = READ_WITH_SET_FLAG_MEMBER;
            let message = format!("error while reading {}", member.name);
            Fail::serial(e.exception().rethrow_same(file, line, function, message))
        }
        other => other,
    }
}

impl<I: Input> Rd<I> {
    pub(super) fn new(inp: I) -> Self {
        Rd {
            inp,
            depth: 0,
            traces: Vec::new(),
            posts: Vec::new(),
            loc_frames: Vec::new(),
            loc_levels: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:189-214 (cited at
    // `report::serial::object_frame::text_frame_site`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608 (`HandleEOF`, cited at
    // `HANDLE_EOF`)
    /// One `BEGIN_OBJECT_FRAME` ... `END_OBJECT_FRAME` of a frame of `kind` (of the type
    /// `type_name`): an error that leaves it gets NCBI's handling. A `CSerialException` takes the
    /// frame's name (`AddFrameInfo` of `GetFrameName()`); a `CEofException` is `HandleEOF`, which
    /// passes it on when fewer than 2 frames stay and throws its own `CSerialException` (with the
    /// frame's information, and without the frame's name) otherwise; another `CException` takes
    /// the frame's line (`NCBI_RETHROW_SAME` with `GetFrameInfo()`, at the `END_OBJECT_FRAME` of
    /// the text reader: only its numbers and tags throw one).
    fn frame<T>(
        &mut self,
        kind: FrameKind,
        type_name: Option<&'static str>,
        f: impl FnOnce(&mut Self) -> R<T>,
    ) -> R<T> {
        self.frame_with(Frame::new(kind, type_name), |s, _| f(s))
    }

    /// [`frame`](Self::frame) of `frame`, which `f` may change (`set_member`).
    #[inline(always)]
    fn frame_with<T>(
        &mut self,
        mut frame: Frame,
        f: impl FnOnce(&mut Self, &mut Frame) -> R<T>,
    ) -> R<T> {
        let depth = self.enter();
        let result = f(self, &mut frame);
        self.exit(depth, frame, result)
    }

    /// `BEGIN_OBJECT_FRAME`: the depth below the new frame, which [`exit`](Self::exit) takes
    /// once what the frame holds has been read. Every `enter` is followed by its `exit`, with no
    /// return between them. The stack has no bound but memory: the frames of the one type that
    /// nests in itself, the Seq-loc, are kept on the heap ([`read_seq_loc_parts`](Self::read_seq_loc_parts)).
    #[inline(always)]
    fn enter(&mut self) -> usize {
        let depth = self.depth;
        self.depth = depth + 1;
        self.inp.set_stack_depth(depth + 1);
        depth
    }

    /// `END_OBJECT_FRAME` of `frame`, which [`enter`](Self::enter) entered at `depth`, after
    /// what it holds was read with `result`: the traces that wait for it take its step, and an
    /// error gets [`leave`](Self::leave) (`exit_traced`, the one place that needs the frame).
    #[inline(always)]
    fn exit<T>(&mut self, depth: usize, frame: Frame, result: R<T>) -> R<T> {
        self.depth = depth;
        self.inp.set_stack_depth(depth);
        match result {
            Ok(value) if self.traces.is_empty() => Ok(value),
            result => self.exit_traced(depth, frame, result),
        }
    }

    /// [`exit`](Self::exit) of a frame that traces wait for or that an error leaves.
    #[cold]
    #[inline(never)]
    fn exit_traced<T>(&mut self, depth: usize, frame: Frame, result: R<T>) -> R<T> {
        if !self.traces.is_empty() {
            self.trace_frame(depth, frame);
        }
        result.map_err(|fail| self.leave(frame, fail))
    }

    /// What leaving `frame` (the stack is below it now) adds to `fail`.
    #[cold]
    #[inline(never)]
    fn leave(&mut self, frame: Frame, fail: Fail) -> Fail {
        match fail {
            Fail::Serial(mut e) => {
                e.add_frame(frame_name(frame.kind, frame.type_name, frame.member()));
                Fail::Serial(e)
            }
            Fail::Eof if self.depth >= 2 => {
                let info = frame_info(frame.kind, frame.type_name, frame.member());
                Fail::serial(NcbiException::new(
                    HANDLE_EOF,
                    format!("{}: {info}", self.inp.position()),
                ))
            }
            Fail::Exception(e) => {
                let (file, line, function) = frame.end();
                let info = frame_info(frame.kind, frame.type_name, frame.member());
                Fail::exception(e.rethrow_same(file, line, function, info))
            }
            other => other,
        }
    }

    /// `SetTopMemberId`: the member or variant of the top frame, `frame`. A trace that waits for
    /// the frame takes the member it had.
    #[inline(always)]
    fn set_member(&mut self, frame: &mut Frame, name: &'static &'static str) {
        if !self.traces.is_empty() {
            self.trace_frame(self.depth - 1, *frame);
        }
        frame.member = Some(name);
    }

    /// The step of `frame`, at `depth`, for the traces that wait for it (`GetStackTrace()`: the
    /// type of the bottom frame, then the member of each member or variant frame and `E` of each
    /// element frame above it). The bottom frame completes them.
    #[cold]
    #[inline(never)]
    fn trace_frame(&mut self, depth: usize, frame: Frame) {
        let mut i = 0;
        while i < self.traces.len() {
            if self.traces[i].level != depth + 1 {
                i += 1;
                continue;
            }
            if depth == 0 {
                let trace = self.traces.remove(i);
                self.complete_trace(trace, frame.type_name);
                continue;
            }
            let trace = &mut self.traces[i];
            trace.level = depth;
            match frame.kind {
                FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                    trace.steps.extend(frame.member().map(StackStep::Member));
                }
                FrameKind::ArrayElement => trace.steps.push(StackStep::Element),
                _ => {}
            }
            i += 1;
        }
    }

    /// Fills `trace`, whose bottom frame is of the type `root`, into its posts.
    fn complete_trace(&mut self, trace: Trace, root: Option<&'static str>) {
        let text = object_stack_trace(root.unwrap_or("?"), trace.steps.into_iter().rev());
        for post in &mut self.posts[trace.posts] {
            post.set_trace(text.clone());
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1038
    /// `ReadMissingSimpleMember` of the members `missing` (`ExpectedMember` of the first that is
    /// neither `OPTIONAL` nor has a `DEFAULT`).
    #[inline]
    fn missing_members(&self, missing: &[Member]) -> R<()> {
        match missing.iter().find(|m| !m.optional) {
            Some(m) => Err(self.expected_member(m)),
            None => Ok(()),
        }
    }

    /// `ExpectedMember` of the member `m`.
    #[cold]
    #[inline(never)]
    fn expected_member(&self, m: &Member) -> Fail {
        throw_error(
            EXPECTED_MEMBER,
            self.inp.position(),
            format!("member {} expected", m.name),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397 (`ReadClassSequential`;
    // the binary stream's, objistrasnb.cpp:1210-1233, is the same with its own `BeginClass`)
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:84-109
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:871-876
    /// A SEQUENCE read member by member, in declaration order (`ReadClassSequential`): a member
    /// that is not in the file is an error unless it is `OPTIONAL` or has a `DEFAULT`. `read`
    /// reads the value of member `index`; a member with a set flag reads it in
    /// `ReadWithSetFlagMember` (cited at [`READ_WITH_SET_FLAG_MEMBER`]). An `OPTIONAL` member held
    /// by `CRef<>` (`Member::pointer`) is read by `read_pointer`: `read` is not called for a
    /// `NULL`, which leaves the member unset.
    pub(super) fn read_class(
        &mut self,
        def: &'static ClassDef,
        mut read: impl FnMut(&mut Self, usize) -> R<()>,
    ) -> R<()> {
        self.frame(FrameKind::Class, Some(def.name), |s| {
            s.inp.begin_class()?;
            let mut pos = 0usize;
            s.frame_with(Frame::new(FrameKind::ClassMember, None), |s, frame| {
                while let Some(index) = s.inp.begin_member(def, pos)? {
                    let member = &def.members[index];
                    s.set_member(frame, &member.name);
                    if index > pos {
                        s.missing_members(&def.members[pos..index])?;
                    }
                    // `ReadPointer` of a member held by `CRef<>` around `read`, which is called at
                    // this one place. Such a member has no set flag (`Member::set_flag`): the
                    // errors of `begin_pointer` are not `ReadWithSetFlagMember`'s.
                    let pointee = if member.pointer {
                        s.begin_pointer(member.ty.type_name())?
                    } else {
                        Pointee::This
                    };
                    if !matches!(pointee, Pointee::Null) {
                        let result = read(s, index);
                        if let Err(fail) = s.end_pointer(pointee, result) {
                            return Err(member_error(fail, member));
                        }
                    }
                    pos = index + 1;
                    s.inp.end_member()?;
                }
                Ok(())
            })?;
            s.missing_members(&def.members[pos..])?;
            s.inp.end_class()
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1500
    /// A CHOICE: `read` reads the value of the variant.
    pub(super) fn read_choice<T>(
        &mut self,
        def: &'static ChoiceDef,
        read: impl FnOnce(&mut Self, usize) -> R<T>,
    ) -> R<T> {
        self.frame(FrameKind::Choice, Some(def.name), |s| {
            s.inp.begin_choice()?;
            let value = s.frame_with(Frame::new(FrameKind::ChoiceVariant, None), |s, frame| {
                let index = s.inp.begin_variant(def)?;
                s.set_member(frame, &def.variants[index].name);
                let value = read(s, index)?;
                s.inp.end_variant()?;
                Ok(value)
            })?;
            s.inp.end_choice()?;
            Ok(value)
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1318-1346
    // (The containers of these types are the anonymous `SEQUENCE OF` and `SET OF` of a member or
    // of an implicit class such as `Packed-seqint`, which NCBI reads as its member without a
    // class frame: the container frame's type has no name, its element frame names none,
    // checked with makeblastdb 2.17.0's traces.)
    /// A `SEQUENCE OF` or `SET OF`: `read` reads one element.
    pub(super) fn read_container(
        &mut self,
        set: bool,
        mut read: impl FnMut(&mut Self) -> R<()>,
    ) -> R<()> {
        self.frame(FrameKind::Array, Some(""), |s| {
            s.inp.begin_container(set)?;
            s.frame(FrameKind::ArrayElement, None, |s| {
                while s.inp.begin_element()? {
                    read(s)?;
                }
                Ok(())
            })?;
            s.inp.end_container()
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/ptrinfo.cpp:210-227
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1133-1213
    // NCBI reference (598d8ae6): c++/src/serial/classinfob.cpp:244-252
    // (`MapType` is `GetClassInfoByName`. The classes it knows are those whose type information
    // the process has built so far: the declared type is one of them, another type may or may not
    // be (makeblastdb 2.17.0: "class not found: Seq-interval" for `:Seq-interval` as the first
    // mask; a type it knows is read whole and then fails with "incompatible member type", or with
    // the first error inside it when the value is not of that type).)
    /// `CPointerTypeInfo::ReadPointer` of a `CRef<>` to a value of the type named `declared` that
    /// is empty when NCBI reads it: an `OPTIONAL` member held by `CRef<>` (the generated classes
    /// preallocate the mandatory ones) and each element of a list of `CRef`s (`AddElementIn`
    /// pushes an empty `CRef`, stltypes.hpp:246-263). `read` reads the value and takes what it
    /// needs of it (`This`, and `Other` naming `declared`, inside `ReadPointer`'s named frame);
    /// `false` is a `NULL`, which it does not read. `@n` is NCBI's error; a value of another type
    /// is refused (`Fail::Unsupported`), as what NCBI does with it depends on what its process
    /// has read before.
    fn read_pointer(
        &mut self,
        declared: &'static str,
        read: impl FnOnce(&mut Self) -> R<()>,
    ) -> R<bool> {
        let pointee = self.begin_pointer(declared)?;
        if let Pointee::Null = pointee {
            return Ok(false);
        }
        let result = read(self);
        self.end_pointer(pointee, result).map(|()| true)
    }

    /// [`read_pointer`](Self::read_pointer) up to the value: where it is ([`Pointee`]), the named
    /// frame of `Other` entered. [`end_pointer`](Self::end_pointer) follows the value. A mask
    /// file's pointers are `NULL` or the value itself; the others are read by
    /// [`begin_other_pointer`](Self::begin_other_pointer).
    #[inline(always)]
    fn begin_pointer(&mut self, declared: &'static str) -> R<Pointee> {
        match self.inp.read_pointer_type()? {
            PointerType::Null => Ok(Pointee::Null),
            PointerType::This => Ok(Pointee::This),
            other => self.begin_other_pointer(other, declared),
        }
    }

    /// [`begin_pointer`](Self::begin_pointer) of an object index (`@n`) or a value of a named type
    /// (`:Type`).
    #[cold]
    #[inline(never)]
    fn begin_other_pointer(&mut self, kind: PointerType, declared: &'static str) -> R<Pointee> {
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
                if name != declared.as_bytes() {
                    return Err(self.other_type(&name, declared));
                }
                Ok(Pointee::Other(self.enter(), declared))
            }
        }
    }

    /// The refusal of `:name` (a value of the type `name`) where a `declared` belongs.
    #[cold]
    #[inline(never)]
    fn other_type(&self, name: &[u8], declared: &'static str) -> Fail {
        let name = String::from_utf8_lossy(name);
        Fail::unsupported(format!(
            "':{name}', a value of the type {name} where a {declared} belongs ({}): whether \
             makeblastdb reads it depends on the types its process has met before",
            self.inp.position()
        ))
    }

    /// [`read_pointer`](Self::read_pointer) after the value of `pointee`, read with `result`: the
    /// `END_OBJECT_FRAME` of the named frame of `Other` and `ReadOtherPointerEnd`.
    #[inline(always)]
    fn end_pointer<T>(&mut self, pointee: Pointee, result: R<T>) -> R<T> {
        match pointee {
            Pointee::Other(depth, declared) => {
                let frame = Frame {
                    read_pointer: true,
                    ..Frame::new(FrameKind::Named, Some(declared))
                };
                let value = self.exit(depth, frame, result)?;
                self.inp.read_other_pointer_end()?;
                Ok(value)
            }
            Pointee::Null | Pointee::This => result,
        }
    }

    /// The element stated as `NULL` that `read_pointer` just read, of the list `list`.
    fn null_element(&self, list: &'static str) -> NullElement {
        NullElement {
            list,
            position: self.inp.position().to_string(),
        }
    }

    /// A `VisibleString`: each byte outside `' '..='~'` read as `#`, with NCBI's posts
    /// (`StringPost`, here with the stack trace).
    fn read_str(&mut self) -> R<Vec<u8>> {
        let mut posts = Vec::new();
        let mut out = Vec::new();
        let result = self
            .inp
            .read_string(StringKind::Visible, &mut out, &mut posts)
            .map(|()| out);
        if !posts.is_empty() {
            let start = self.posts.len();
            for post in posts {
                self.posts.push(match post {
                    StringPost::Utf8Tag => MaskPost::Utf8Tag(String::new()),
                    StringPost::BadChar {
                        position,
                        byte,
                        string,
                    } => MaskPost::BadChar(ReplacedChar {
                        trace: String::new(),
                        byte,
                        position,
                        string: string.into(),
                    }),
                });
            }
            // `GetStackTrace()` now: the frames on the stack give it as they are left.
            let trace = Trace {
                posts: start..self.posts.len(),
                level: self.depth,
                steps: Vec::new(),
            };
            if self.depth == 0 {
                self.complete_trace(trace, None);
            } else {
                self.traces.push(trace);
            }
        }
        result
    }

    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // (A `CEofException` stays one through `NCBI_RETHROW_SAME`: `HandleEOF` does not report it.)
    /// An ENUMERATED or an INTEGER with named values (`CEnumeratedTypeInfo::ReadEnum`).
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        const READ_ENUM: (&str, u32, &str) = (
            "src/serial/enumerated.cpp",
            428,
            "ncbi::CEnumeratedTypeInfo::ReadEnum",
        );
        let (file, line, function) = READ_ENUM;
        self.inp.read_enum(def).map_err(|fail| match fail {
            Fail::Serial(e) => {
                Fail::serial(
                    e.exception()
                        .rethrow_same(file, line, function, "invalid enum value"),
                )
            }
            Fail::Exception(e) => {
                Fail::exception(e.rethrow_same(file, line, function, "invalid enum value"))
            }
            other => other,
        })
    }

    /// Reads a value of type `ty` and appends a canonical encoding of it to `out`: two values
    /// give the same bytes exactly when they are equal as `CSerialObject::Equals` compares
    /// them, apart from members that are left out in one and stated with their `DEFAULT` in the
    /// other (those differ here). The encoding is self-delimiting (every item states its kind,
    /// a string its length, a list its count).
    pub(super) fn read_value(&mut self, ty: &'static Ty, out: &mut Vec<u8>) -> R<()> {
        match ty {
            Ty::Int32 => {
                let v = self.inp.read_i32()?;
                out.push(b'i');
                out.extend_from_slice(&i64::from(v).to_le_bytes());
            }
            Ty::Int64 => {
                let v = self.inp.read_i64()?;
                out.push(b'i');
                out.extend_from_slice(&v.to_le_bytes());
            }
            Ty::Uint32 => {
                let v = self.inp.read_u32()?;
                out.push(b'i');
                out.extend_from_slice(&i64::from(v).to_le_bytes());
            }
            Ty::Bool => {
                let v = self.inp.read_bool()?;
                out.push(b'b');
                out.push(u8::from(v));
            }
            Ty::Str => {
                let v = self.read_str()?;
                out.push(b's');
                out.extend_from_slice(&(v.len() as u64).to_le_bytes());
                out.extend_from_slice(&v);
            }
            Ty::Null => {
                self.inp.read_null()?;
                out.push(b'n');
            }
            Ty::Enum(def) => {
                let v = self.read_enum(def)?;
                out.push(b'e');
                out.extend_from_slice(&v.to_le_bytes());
            }
            Ty::Class(def) => {
                out.push(b'S');
                let mut next = 0usize;
                self.read_class(def, |s, index| {
                    for _ in next..index {
                        out.push(0);
                    }
                    out.push(1);
                    next = index + 1;
                    s.read_value(def.members[index].ty, out)
                })?;
                for _ in next..def.members.len() {
                    out.push(0);
                }
            }
            Ty::Choice(def) => {
                out.push(b'C');
                self.read_choice(def, |s, index| {
                    out.push(index as u8);
                    s.read_value(def.variants[index].ty, out)
                })?;
            }
            Ty::SeqOf(element) | Ty::SetOf(element) => {
                // The lists read here hold numbers (`Int-fuzz.alt`); the lists of `CRef`s (`masks`,
                // `mix`, `equiv`, `packed-int`) are read by their own readers, through
                // `read_pointer`.
                out.push(b'L');
                let count_at = out.len();
                out.extend_from_slice(&[0u8; 8]);
                let mut count = 0u64;
                self.read_container(matches!(ty, Ty::SetOf(_)), |s| {
                    count += 1;
                    s.read_value(element, out)
                })?;
                out[count_at..count_at + 8].copy_from_slice(&count.to_le_bytes());
            }
        }
        Ok(())
    }

    /// A `Seq-id` of any alternative, as the value that NCBI's `CSeq_id` holds (`ParsedSeqId`: the
    /// one Seq-id type of LOSAT; its members at a DEFAULT or OPTIONAL are `None` when the file does
    /// not state them, as NCBI keeps their set flags). Text and binary files take the same path:
    /// the frames and members are the reader's, the primitives the encoding's.
    fn read_seq_id(&mut self) -> R<ParsedSeqId> {
        self.read_choice(&SEQ_ID, |s, index| {
            Ok(match index {
                0 => ParsedSeqId::Local(s.read_object_id()?),
                1 => ParsedSeqId::Gibbsq(s.inp.read_i32()?),
                2 => ParsedSeqId::Gibbmt(s.inp.read_i32()?),
                3 => ParsedSeqId::Giim(s.read_giimport_id()?),
                4 => ParsedSeqId::Genbank(s.read_textseq_id()?),
                5 => ParsedSeqId::Embl(s.read_textseq_id()?),
                6 => ParsedSeqId::Pir(s.read_textseq_id()?),
                7 => ParsedSeqId::Swissprot(s.read_textseq_id()?),
                8 => ParsedSeqId::Patent(s.read_patent_seq_id()?),
                9 => ParsedSeqId::Other(s.read_textseq_id()?),
                10 => ParsedSeqId::General(s.read_dbtag()?),
                11 => ParsedSeqId::Gi(s.inp.read_i64()?),
                12 => ParsedSeqId::Ddbj(s.read_textseq_id()?),
                13 => ParsedSeqId::Prf(s.read_textseq_id()?),
                14 => ParsedSeqId::Pdb(s.read_pdb_seq_id()?),
                15 => ParsedSeqId::Tpg(s.read_textseq_id()?),
                16 => ParsedSeqId::Tpe(s.read_textseq_id()?),
                17 => ParsedSeqId::Tpd(s.read_textseq_id()?),
                18 => ParsedSeqId::Gpipe(s.read_textseq_id()?),
                _ => ParsedSeqId::NamedAnnotTrack(s.read_textseq_id()?),
            })
        })
    }

    /// An `Object-id` (`T_OBJECT_ID`).
    fn read_object_id(&mut self) -> R<ObjectId> {
        let Ty::Choice(def) = &T_OBJECT_ID else {
            unreachable!("Object-id is a CHOICE")
        };
        self.read_choice(def, |s, index| {
            Ok(if index == OBJECT_ID_ID {
                ObjectId::Id(s.inp.read_i32()?)
            } else {
                ObjectId::Str(s.read_str()?)
            })
        })
    }

    fn read_dbtag(&mut self) -> R<Dbtag> {
        let (mut db, mut tag) = (None, None);
        self.read_class(&DBTAG, |s, index| {
            match index {
                0 => db = Some(s.read_str()?),
                _ => tag = Some(s.read_object_id()?),
            }
            Ok(())
        })?;
        Ok(Dbtag {
            db: mandatory(db),
            tag: mandatory(tag),
        })
    }

    fn read_textseq_id(&mut self) -> R<TextseqId> {
        let mut id = TextseqId::default();
        self.read_class(&TEXTSEQ_ID, |s, index| {
            match index {
                0 => id.name = Some(s.read_str()?),
                1 => id.accession = Some(s.read_str()?),
                2 => id.release = Some(s.read_str()?),
                _ => id.version = Some(s.inp.read_i32()?),
            }
            Ok(())
        })?;
        Ok(id)
    }

    fn read_giimport_id(&mut self) -> R<GiimportId> {
        let (mut id, mut db, mut release) = (None, None, None);
        self.read_class(&GIIMPORT_ID, |s, index| {
            match index {
                0 => id = Some(s.inp.read_i32()?),
                1 => db = Some(s.read_str()?),
                _ => release = Some(s.read_str()?),
            }
            Ok(())
        })?;
        Ok(GiimportId {
            id: mandatory(id),
            db,
            release,
        })
    }

    fn read_patent_seq_id(&mut self) -> R<PatentSeqId> {
        let (mut seqid, mut cit) = (None, None);
        self.read_class(&PATENT_SEQ_ID, |s, index| {
            match index {
                0 => seqid = Some(s.inp.read_i32()?),
                _ => cit = Some(s.read_id_pat()?),
            }
            Ok(())
        })?;
        Ok(PatentSeqId {
            seqid: mandatory(seqid),
            cit: mandatory(cit),
        })
    }

    fn read_id_pat(&mut self) -> R<IdPat> {
        let (mut country, mut id, mut doc_type) = (None, None, None);
        self.read_class(&ID_PAT, |s, index| {
            match index {
                0 => country = Some(s.read_str()?),
                1 => {
                    id = Some(s.read_choice(&ID_PAT_ID, |s, variant| {
                        Ok(if variant == 0 {
                            PatentNumber::Number(s.read_str()?)
                        } else {
                            PatentNumber::AppNumber(s.read_str()?)
                        })
                    })?);
                }
                _ => doc_type = Some(s.read_str()?),
            }
            Ok(())
        })?;
        Ok(IdPat {
            country: mandatory(country),
            id: mandatory(id),
            doc_type,
        })
    }

    fn read_pdb_seq_id(&mut self) -> R<PdbSeqId> {
        let (mut mol, mut chain, mut rel, mut chain_id) = (None, None, None, None);
        self.read_class(&PDB_SEQ_ID, |s, index| {
            match index {
                0 => mol = Some(s.read_str()?),
                1 => chain = Some(s.inp.read_i32()?),
                2 => rel = Some(s.read_date()?),
                _ => chain_id = Some(s.read_str()?),
            }
            Ok(())
        })?;
        Ok(PdbSeqId {
            mol: mandatory(mol),
            chain,
            rel,
            chain_id,
        })
    }

    fn read_date(&mut self) -> R<Date> {
        self.read_choice(&DATE, |s, variant| {
            if variant == 0 {
                return Ok(Date::Str(s.read_str()?));
            }
            let mut std = DateStd {
                year: 0,
                month: None,
                day: None,
                season: None,
                hour: None,
                minute: None,
                second: None,
            };
            let mut year = None;
            s.read_class(&DATE_STD, |s, index| {
                match index {
                    0 => year = Some(s.inp.read_i32()?),
                    1 => std.month = Some(s.inp.read_i32()?),
                    2 => std.day = Some(s.inp.read_i32()?),
                    3 => std.season = Some(s.read_str()?),
                    4 => std.hour = Some(s.inp.read_i32()?),
                    5 => std.minute = Some(s.inp.read_i32()?),
                    _ => std.second = Some(s.inp.read_i32()?),
                }
                Ok(())
            })?;
            std.year = mandatory(year);
            Ok(Date::Std(std))
        })
    }

    /// A `Seq-interval`, and its Seq-id.
    fn read_seq_interval(&mut self) -> R<(MaskInterval, ParsedSeqId)> {
        let (mut from, mut to, mut strand, mut id) = (0u32, 0u32, None, None);
        let mut scratch = Vec::new();
        self.read_class(&SEQ_INTERVAL, |s, index| {
            match index {
                SEQ_INTERVAL_FROM => from = s.inp.read_u32()?,
                SEQ_INTERVAL_TO => to = s.inp.read_u32()?,
                SEQ_INTERVAL_STRAND => {
                    let v = s.read_enum(&NA_STRAND)?;
                    // The values of Na-strand are 0..=4 and 255.
                    strand = Some(v as u8);
                }
                SEQ_INTERVAL_ID => id = Some(s.read_seq_id()?),
                _ => {
                    // fuzz-from, fuzz-to: read, not used.
                    scratch.clear();
                    s.read_value(SEQ_INTERVAL.members[index].ty, &mut scratch)?;
                }
            }
            Ok(())
        })?;
        Ok((MaskInterval { from, to, strand }, mandatory(id)))
    }

    /// A `Seq-point` (also the members `a` and `b` of a `Seq-bond`), and its Seq-id.
    fn read_seq_point(&mut self) -> R<(MaskInterval, ParsedSeqId)> {
        let (mut point, mut strand, mut id) = (0u32, None, None);
        let mut scratch = Vec::new();
        self.read_class(&SEQ_POINT, |s, index| {
            match index {
                0 => point = s.inp.read_u32()?,
                1 => strand = Some(s.read_enum(&NA_STRAND)? as u8),
                2 => id = Some(s.read_seq_id()?),
                _ => {
                    // fuzz: read, not used.
                    scratch.clear();
                    s.read_value(SEQ_POINT.members[index].ty, &mut scratch)?;
                }
            }
            Ok(())
        })?;
        let id = mandatory(id);
        let interval = MaskInterval {
            from: point,
            to: point,
            strand,
        };
        Ok((interval, id))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3528-3620
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:403-422
    // NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_loc.hpp:926-945
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1555,1559-1600
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1500 (`ReadChoiceSimple`, cited
    // at `read_choice`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274 (`ReadContainer`, cited at
    // `read_container`)
    /// One `Seq-loc` (any variant, nested `mix` and `equiv` included): appends the ranges
    /// `CSeq_loc_CI` visits to `parts`, in its order (a `whole` Seq-loc is the interval 0 to
    /// 4294967294, whose `CRange` is `TRange::GetWhole()`; a point the interval of one position),
    /// and returns its walk ([`IdWalk`]) as `CSeq_loc::CheckId` walks it. An element of
    /// `packed-int`, `mix` or `equiv` stated as `NULL` (an empty `CRef<>`) is where the walk
    /// dereferences it, if it gets that far (`**ii`, `(*li)->CheckId`).
    ///
    /// NCBI reads a `mix` or `equiv` by recursion (`ReadChoiceSimple`, `ReadContainer` and
    /// `ReadPointer` of each element), so a Seq-loc nested deep enough overflows its stack
    /// (PD-LOSAT-NCBI-DEFECTS exception 25). This reader keeps the frames of the Seq-locs, mixes and
    /// equivs open (`Rd::loc_frames`) and the walks of the mixes and equivs (`Rd::loc_levels`) on
    /// the heap, and reads every depth with no recursion, as far as memory allows: it enters and
    /// leaves the frames in NCBI's order, the steps of `read_choice`, `read_container` and
    /// `read_pointer` for a Seq-loc, and an error leaves every frame still open, innermost first,
    /// as each of NCBI's functions passes it on. A variant that holds no Seq-loc is read by
    /// [`read_seq_loc_variant`](Self::read_seq_loc_variant).
    fn read_seq_loc_parts(&mut self, parts: &mut Vec<MaskPart>) -> R<IdWalk> {
        let mut frames = std::mem::take(&mut self.loc_frames);
        let mut levels = std::mem::take(&mut self.loc_levels);
        let result = match self.read_seq_loc_steps(parts, &mut frames, &mut levels) {
            Ok(walk) => Ok(walk),
            Err(mut fail) => {
                while let Some(frame) = frames.pop() {
                    let depth = self.depth - 1;
                    fail = match self.exit(depth, frame.frame(), Err::<(), Fail>(fail)) {
                        Err(fail) => fail,
                        Ok(()) => unreachable!("an error leaves its frames as an error"),
                    };
                }
                Err(fail)
            }
        };
        debug_assert!(frames.is_empty());
        levels.clear();
        self.loc_frames = frames;
        self.loc_levels = levels;
        result
    }

    /// [`read_seq_loc_parts`](Self::read_seq_loc_parts) up to its end or its error: `frames` holds
    /// the frames it entered and has not left, `levels` the mixes and equivs it reads.
    fn read_seq_loc_steps(
        &mut self,
        parts: &mut Vec<MaskPart>,
        frames: &mut Vec<LocFrame>,
        levels: &mut Vec<LocLevel>,
    ) -> R<IdWalk> {
        let mut step = LocStep::Begin;
        loop {
            step = match step {
                LocStep::Begin => {
                    // `read_choice` of a Seq-loc
                    self.enter_loc(frames, LocFrame::Choice);
                    self.inp.begin_choice()?;
                    self.enter_loc(frames, LocFrame::Variant(None));
                    let index = self.inp.begin_variant(&SEQ_LOC)?;
                    self.set_loc_variant(frames, index);
                    if index == SEQ_LOC_MIX || index == SEQ_LOC_EQUIV {
                        // `read_container` of its elements, each a Seq-loc by `read_pointer`
                        let equiv = index == SEQ_LOC_EQUIV;
                        self.enter_loc(frames, LocFrame::Array);
                        self.inp.begin_container(equiv)?;
                        self.enter_loc(frames, LocFrame::Element);
                        levels.push(LocLevel { walk: None, equiv });
                        LocStep::Element
                    } else {
                        let mut walk = IdWalk::Same(None);
                        self.read_seq_loc_variant(index, parts, &mut walk)?;
                        LocStep::End(walk)
                    }
                }
                LocStep::Element => {
                    if self.inp.begin_element()? {
                        match self.begin_pointer(SEQ_LOC.name)? {
                            Pointee::Null => {
                                let level = levels.last_mut().expect(IN_A_LEVEL);
                                let null = self.null_element(level.list());
                                level.walk().null(null);
                                LocStep::Element
                            }
                            Pointee::This => LocStep::Begin,
                            Pointee::Other(..) => {
                                // `begin_pointer` entered the named frame of `ReadPointer`
                                frames.push(LocFrame::Pointer);
                                LocStep::Begin
                            }
                        }
                    } else {
                        self.exit_loc(frames);
                        self.inp.end_container()?;
                        self.exit_loc(frames);
                        let level = levels.pop().expect(IN_A_LEVEL);
                        LocStep::End(level.walk.map_or(IdWalk::Same(None), |walk| *walk))
                    }
                }
                LocStep::End(walk) => {
                    // the end of `read_choice`
                    self.inp.end_variant()?;
                    self.exit_loc(frames);
                    self.inp.end_choice()?;
                    self.exit_loc(frames);
                    let Some(level) = levels.last_mut() else {
                        return Ok(walk);
                    };
                    // the end of `read_pointer` of an element of `level`
                    level.nested(walk);
                    if matches!(frames.last(), Some(LocFrame::Pointer)) {
                        self.exit_loc(frames);
                        self.inp.read_other_pointer_end()?;
                    }
                    LocStep::Element
                }
            };
        }
    }

    /// `BEGIN_OBJECT_FRAME` of `frame`, a frame of the Seq-loc reader, kept on `frames`.
    #[inline(always)]
    fn enter_loc(&mut self, frames: &mut Vec<LocFrame>, frame: LocFrame) {
        self.enter();
        frames.push(frame);
    }

    /// `END_OBJECT_FRAME` of the top frame of `frames` once what it holds has been read (`exit`).
    #[inline(always)]
    fn exit_loc(&mut self, frames: &mut Vec<LocFrame>) {
        let frame = frames.pop().expect("the Seq-loc reader's frame is open");
        let depth = self.depth - 1;
        let left: R<()> = self.exit(depth, frame.frame(), Ok(()));
        debug_assert!(left.is_ok());
    }

    /// `SetTopMemberId` of the Seq-loc's variant `index` on its variant frame, the top of `frames`.
    #[inline(always)]
    fn set_loc_variant(&mut self, frames: &mut [LocFrame], index: usize) {
        let top = frames.last_mut().expect("the variant frame is open");
        let mut frame = top.frame();
        self.set_member(&mut frame, &SEQ_LOC.variants[index].name);
        // Seq-loc has 11 variants.
        *top = LocFrame::Variant(Some(index as u8));
    }

    /// The value of the variant `index` of a Seq-loc that holds no Seq-loc (every variant but
    /// `mix` and `equiv`, which [`read_seq_loc_steps`](Self::read_seq_loc_steps) reads): its
    /// ranges appended to `parts`, its Seq-ids walked in `walk`.
    fn read_seq_loc_variant(
        &mut self,
        index: usize,
        parts: &mut Vec<MaskPart>,
        walk: &mut IdWalk,
    ) -> R<()> {
        match index {
            SEQ_LOC_NULL => {
                self.inp.read_null()?;
                parts.push(MaskPart::Null);
            }
            SEQ_LOC_EMPTY => {
                let seq_id = self.read_seq_id()?;
                walk.update(seq_id);
                parts.push(MaskPart::Empty);
            }
            SEQ_LOC_WHOLE => {
                let seq_id = self.read_seq_id()?;
                walk.update(seq_id);
                parts.push(MaskPart::Interval(MaskInterval {
                    from: 0,
                    to: u32::MAX - 1,
                    strand: None,
                }));
            }
            SEQ_LOC_INT => {
                let (interval, seq_id) = self.read_seq_interval()?;
                walk.update(seq_id);
                parts.push(MaskPart::Interval(interval));
            }
            SEQ_LOC_PACKED_INT => {
                self.read_container(false, |s| {
                    let read = s.read_pointer(SEQ_INTERVAL.name, |s| {
                        let (interval, seq_id) = s.read_seq_interval()?;
                        walk.update(seq_id);
                        parts.push(MaskPart::Interval(interval));
                        Ok(())
                    })?;
                    if !read {
                        walk.null(s.null_element(LIST_PACKED_INT));
                    }
                    Ok(())
                })?;
            }
            SEQ_LOC_PNT => {
                let (interval, seq_id) = self.read_seq_point()?;
                walk.update(seq_id);
                parts.push(MaskPart::Interval(interval));
            }
            SEQ_LOC_PACKED_PNT => {
                let (mut strand, mut seq_id, mut points) = (None, None, Vec::new());
                let mut scratch = Vec::new();
                self.read_class(&PACKED_SEQPNT, |s, member| {
                    match member {
                        0 => strand = Some(s.read_enum(&NA_STRAND)? as u8),
                        1 => seq_id = Some(s.read_seq_id()?),
                        3 => s.read_container(false, |s| {
                            points.push(s.inp.read_u32()?);
                            Ok(())
                        })?,
                        _ => {
                            // fuzz: read, not used.
                            scratch.clear();
                            s.read_value(PACKED_SEQPNT.members[member].ty, &mut scratch)?;
                        }
                    }
                    Ok(())
                })?;
                walk.update(mandatory(seq_id));
                parts.extend(points.into_iter().map(|p| {
                    MaskPart::Interval(MaskInterval {
                        from: p,
                        to: p,
                        strand,
                    })
                }));
            }
            SEQ_LOC_MIX | SEQ_LOC_EQUIV => {
                unreachable!("read_seq_loc_steps reads the Seq-locs of a mix or equiv")
            }
            SEQ_LOC_BOND => {
                // `b` is an OPTIONAL `CRef<>` (`read_class`): a `NULL` leaves it unset, which
                // `x_CheckId` and `CSeq_loc_CI` skip (`CanGetB`, `IsSetB`).
                self.read_class(&SEQ_BOND, |s, _| {
                    let (interval, seq_id) = s.read_seq_point()?;
                    walk.update(seq_id);
                    parts.push(MaskPart::Interval(interval));
                    Ok(())
                })?;
            }
            _ => {
                // `feat`: read; `x_CheckId` fails, so the Seq-loc has no Seq-id.
                let mut scratch = Vec::new();
                self.read_value(SEQ_LOC.variants[index].ty, &mut scratch)?;
                *walk = IdWalk::Stopped;
            }
        }
        Ok(())
    }

    /// One `Seq-loc` of `masks`: the ranges `CSeq_loc_CI` visits and `CSeq_loc::GetId()`, or the
    /// element stated as `NULL` that `GetId` dereferences.
    fn read_seq_loc(&mut self) -> R<Result<MaskLoc, NullElement>> {
        let mut parts = Vec::new();
        let id = match self.read_seq_loc_parts(&mut parts)? {
            IdWalk::Same(id) => id,
            IdWalk::Stopped => None,
            IdWalk::Null(null) => return Ok(Err(null)),
        };
        Ok(Ok(MaskLoc { id, parts }))
    }

    /// A `Blast-mask-list`; the caller checks the end of the class frame. Its masks are kept up to
    /// the first one whose `GetId` dereferences an element stated as `NULL` (the mask itself, or
    /// one inside it), which is `MaskList::null_element`; the rest is read and dropped.
    fn read_mask_list_members(&mut self) -> R<MaskList> {
        let mut masks = Vec::new();
        let mut null_element = None;
        let mut more = None;
        self.read_class(&BLAST_MASK_LIST, |s, index| {
            if index == 0 {
                s.read_container(false, |s| {
                    let read = s.read_pointer(SEQ_LOC.name, |s| {
                        let loc = s.read_seq_loc()?;
                        if null_element.is_none() {
                            match loc {
                                Ok(loc) => masks.push(loc),
                                Err(null) => null_element = Some(null),
                            }
                        }
                        Ok(())
                    })?;
                    if !read && null_element.is_none() {
                        null_element = Some(s.null_element(LIST_MASKS));
                    }
                    Ok(())
                })
            } else {
                more = Some(s.inp.read_bool()?);
                Ok(())
            }
        })?;
        Ok(MaskList {
            masks,
            more: mandatory(more),
            null_element,
        })
    }

    /// The read of `file` into the object through `MSerial_Asn...`, for a `Blast-db-mask-info`.
    fn read_info(&mut self) -> R<MaskInfo> {
        let mut algo_id = 0;
        let mut algo_program = 0;
        let mut algo_options = Vec::new();
        let mut masks = None;
        self.read_top(&BLAST_DB_MASK_INFO, |s| {
            s.read_class(&BLAST_DB_MASK_INFO, |s, index| {
                match index {
                    0 => algo_id = s.inp.read_i32()?,
                    1 => algo_program = s.read_enum(&BLAST_FILTER_PROGRAM)?,
                    2 => algo_options = s.read_str()?,
                    _ => masks = Some(s.read_mask_list_members()?),
                }
                Ok(())
            })
        })?;
        Ok(MaskInfo {
            algo_id,
            algo_program,
            algo_options,
            masks: mandatory(masks),
        })
    }

    /// The read of `file` into the object through `MSerial_Asn...`, for a `Blast-mask-list` (the
    /// continuation).
    fn read_list(&mut self) -> R<MaskList> {
        let mut list = None;
        self.read_top(&BLAST_MASK_LIST, |s| {
            list = Some(s.read_mask_list_members()?);
            Ok(())
        })?;
        Ok(mandatory(list))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:882-901
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:971-975
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:959-969
    /// The read of `file` into the object through `MSerial_Asn...`, of the type `def`:
    /// `SkipFileHeader`, then `Read` with `read`, each in a frame of its own.
    fn read_top(&mut self, def: &'static ClassDef, read: impl FnOnce(&mut Self) -> R<()>) -> R<()> {
        // Every frame entered was left (`exit`), and with it every stack trace completed.
        debug_assert!(self.depth == 0 && self.traces.is_empty());
        self.depth = 0;
        self.frame(FrameKind::Named, Some(def.name), |s| {
            match s.inp.read_file_header()? {
                Some(name) if !name.is_empty() && name != def.name.as_bytes() => Err(throw_error(
                    SKIP_FILE_HEADER,
                    s.inp.position(),
                    [
                        b"incompatible type ".as_slice(),
                        &name,
                        b"<>",
                        def.name.as_bytes(),
                    ]
                    .concat(),
                )),
                _ => Ok(()),
            }
        })?;
        self.frame(FrameKind::Named, Some(def.name), read)
    }
}

/// The value of a member that `ReadClassSequential` read: a member that is neither `OPTIONAL` nor
/// has a `DEFAULT` was read when `Rd::read_class` returns (`missing_members`).
fn mandatory<T>(value: Option<T>) -> T {
    value.expect("ReadClassSequential has read every mandatory member")
}

/// One `Seq-id` of an ASN.1 binary value, read by the mask reader's builder, for the test that
/// compares it with `seqid::ber::decode_seq_id` (DECISIONS 35).
#[cfg(test)]
pub(super) fn seq_id_from_binary(data: &[u8]) -> R<ParsedSeqId> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let mut rd = Rd::new(BerReader::new(SliceSource::new(data), Names::MASK));
    rd.frame(FrameKind::Named, Some("Seq-id"), |s| s.read_seq_id())
}

/// [`seq_id_from_binary`] for ASN.1 text, for the test that compares the builder with
/// `seqid::asn_text`.
#[cfg(test)]
pub(super) fn seq_id_from_text(data: &[u8]) -> R<ParsedSeqId> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let mut rd = Rd::new(TextInput::new(SliceSource::new(data), Names::MASK));
    rd.frame(FrameKind::Named, Some("Seq-id"), |s| s.read_seq_id())
}

/// Where the reader is in its file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Open,
    /// The clean end of the file.
    Ended,
    /// An error: every later call is an error again, never a clean end.
    Failed,
    /// The file ended where a continuation was expected.
    FailedEof,
}

fn unknown_encoding(name: &Path) -> MaskFileError {
    MaskFileError::Ncbi(WriteDbError::Options(format!(
        "Unknown encoding for mask data in '{}'",
        name.display()
    )))
}

/// The bytes of a mask file, read in blocks as its objects are read.
type FileSource = StreamSource<Box<dyn Read>>;

/// Where the objects of a mask file come from: one reader for all objects of the file (only the
/// bytes of the object being read are kept), each object read by a new state of NCBI's stream
/// (`Input::new_object`).
enum Source {
    /// ASN.1 text.
    Text(Box<Rd<TextInput<FileSource>>>),
    /// ASN.1 binary.
    Binary(Box<Rd<BerReader<FileSource>>>),
}

/// A mask file read from the disk.
pub(crate) struct MaskFileReader {
    path: PathBuf,
    source: Source,
    /// Set by the first error or end of file (NCBI stops reading the file there).
    state: State,
    /// What NCBI's reader posted while it read the objects, not yet taken (`take_posts`).
    posts: Vec<MaskPost>,
}

impl fmt::Debug for MaskFileReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut debug = f.debug_struct("MaskFileReader");
        debug.field("path", &self.path);
        match &self.source {
            Source::Text(_) => debug.field("format", &MaskFileFormat::TextAsn),
            Source::Binary(rd) => debug
                .field("format", &MaskFileFormat::BinaryAsn)
                .field("offset", &rd.inp.offset()),
        };
        debug.field("state", &self.state).finish()
    }
}

impl MaskFileReader {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:933-947
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:160-181
    /// [`open_named`](Self::open_named) for a file that messages call by its path.
    #[cfg(test)]
    pub(crate) fn open(path: &Path) -> Result<Self, MaskFileError> {
        Self::open_named(path, path)
    }

    /// Opens the file and guesses its format from its first bytes; messages call the file
    /// `name` (the name on makeblastdb's command line). The text of the error for a file in
    /// neither format is NCBI's. A file that cannot be opened or read is a failed stream for the
    /// format guess (`x_TestInput`), which answers `eUnknown`: the same error. The file is then
    /// read in blocks as its objects are read (a read error ends the file, as the bad bit of
    /// NCBI's stream does).
    pub(crate) fn open_named(path: &Path, name: &Path) -> Result<Self, MaskFileError> {
        match std::fs::File::open(path) {
            Ok(file) => Self::from_reader(name, Box::new(file)),
            Err(_) => Err(unknown_encoding(name)),
        }
    }

    /// The reader of a file whose bytes are `data`; `path` is only used in messages.
    #[cfg(test)]
    pub(crate) fn from_bytes(path: &Path, data: Vec<u8>) -> Result<Self, MaskFileError> {
        Self::from_reader(path, Box::new(std::io::Cursor::new(data)))
    }

    /// The reader of the bytes `input` gives, read one byte per fill after the format guess (the
    /// bytes it keeps end at every byte); `path` is only used in messages.
    #[cfg(test)]
    pub(crate) fn from_read(path: &Path, input: Box<dyn Read>) -> Result<Self, MaskFileError> {
        let mut reader = Self::from_reader(path, input)?;
        match &mut reader.source {
            Source::Text(rd) => rd.inp.r.src.fill_one_byte_at_a_time(),
            Source::Binary(rd) => rd.inp.src.fill_one_byte_at_a_time(),
        }
        Ok(reader)
    }

    /// The reader of the file `input`; `path` is only used in messages.
    fn from_reader(path: &Path, mut input: Box<dyn Read>) -> Result<Self, MaskFileError> {
        let Some((format, sample)) = guess_format(&mut input) else {
            return Err(unknown_encoding(path));
        };
        let src = StreamSource::new(input, &sample);
        let source = match format {
            MaskFileFormat::TextAsn => {
                Source::Text(Box::new(Rd::new(TextInput::new(src, Names::MASK))))
            }
            MaskFileFormat::BinaryAsn => {
                Source::Binary(Box::new(Rd::new(BerReader::new(src, Names::MASK))))
            }
        };
        Ok(MaskFileReader {
            path: path.to_path_buf(),
            source,
            state: State::Open,
            posts: Vec::new(),
        })
    }

    /// The bytes of the file that the reader keeps (the stream's buffer).
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        match &self.source {
            Source::Text(rd) => rd.inp.buffered(),
            Source::Binary(rd) => rd.inp.src.buffered(),
        }
    }

    #[cfg(test)]
    pub(crate) fn format(&self) -> MaskFileFormat {
        match self.source {
            Source::Text(_) => MaskFileFormat::TextAsn,
            Source::Binary(_) => MaskFileFormat::BinaryAsn,
        }
    }

    fn run<T>(
        &mut self,
        text: impl FnOnce(&mut Rd<TextInput<FileSource>>) -> R<T>,
        binary: impl FnOnce(&mut Rd<BerReader<FileSource>>) -> R<T>,
    ) -> R<T> {
        match &mut self.source {
            Source::Text(rd) => {
                rd.inp.new_object();
                let result = text(rd);
                self.posts.append(&mut rd.posts);
                result
            }
            Source::Binary(rd) => {
                rd.inp.new_object();
                let result = binary(rd);
                self.posts.append(&mut rd.posts);
                result
            }
        }
    }

    /// The messages that NCBI's reader posted while it read the objects returned so far (and the
    /// one that failed), in the order posted: makeblastdb writes them as they come, before what
    /// it then does with the object or the error.
    pub(crate) fn take_posts(&mut self) -> Vec<MaskPost> {
        std::mem::take(&mut self.posts)
    }

    /// A call after the reader stopped (makeblastdb does not read on after an error).
    fn stopped(&self) -> MaskFileError {
        if self.state == State::FailedEof {
            return MaskFileError::EofInContinuation {
                path: self.path.clone(),
            };
        }
        MaskFileError::Unsupported(format!(
            "mask data in '{}' after the error that stopped its reader",
            self.path.display()
        ))
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-263 (`CATCH_ALL`)
    /// The error of the reader's `fail`: what makeblastdb's `CATCH_ALL` reports for what NCBI's
    /// reader threw.
    fn fail(&mut self, fail: Fail) -> MaskFileError {
        self.state = if matches!(fail, Fail::Eof) {
            State::FailedEof
        } else {
            State::Failed
        };
        match fail {
            Fail::Serial(e) => MaskFileError::Ncbi(WriteDbError::Exception(e.exception())),
            Fail::Exception(e) => MaskFileError::Ncbi(WriteDbError::Exception(*e)),
            Fail::OutOfMemory => MaskFileError::Ncbi(WriteDbError::OutOfMemory),
            Fail::Std(what) => MaskFileError::Ncbi(WriteDbError::Other(what.to_string())),
            Fail::Unsupported(what) => MaskFileError::Unsupported(what.0),
            Fail::Eof => MaskFileError::EofInContinuation {
                path: self.path.clone(),
            },
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:947-960
    /// The next `Blast-db-mask-info`; `None` at the end of the file. A `CEofException` of the
    /// serial reader (nothing left but white space and comments, or a file that ends inside the
    /// type name or the first bytes of an object) is the end, as in `x_ProcessMaskData`; a
    /// partly read object is dropped.
    pub(crate) fn next_info(&mut self) -> Result<Option<MaskInfo>, MaskFileError> {
        match self.state {
            State::Open => {}
            State::Ended => return Ok(None),
            State::Failed | State::FailedEof => return Err(self.stopped()),
        }
        match self.run(|rd| rd.read_info(), |rd| rd.read_info()) {
            Ok(info) => Ok(Some(info)),
            Err(Fail::Eof) => {
                self.state = State::Ended;
                Ok(None)
            }
            Err(fail) => Err(self.fail(fail)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:994-998
    /// The `Blast-mask-list` that follows an object with `more TRUE`. A file that ends there is
    /// [`MaskFileError::EofInContinuation`]. After an error every call is an error again.
    pub(crate) fn next_list(&mut self) -> Result<MaskList, MaskFileError> {
        match self.state {
            State::Open => {}
            State::Ended => {
                return Err(MaskFileError::EofInContinuation {
                    path: self.path.clone(),
                })
            }
            State::Failed | State::FailedEof => return Err(self.stopped()),
        }
        match self.run(|rd| rd.read_list(), |rd| rd.read_list()) {
            Ok(list) => Ok(list),
            Err(fail) => Err(self.fail(fail)),
        }
    }
}
