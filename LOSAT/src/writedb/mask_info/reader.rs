//! The reader of makeblastdb's mask files (`makeblastdb.cpp` `s_ReadObject` and
//! `x_ProcessMaskData`): the format guess, then `Blast-db-mask-info` objects one after another
//! in ASN.1 text or binary.
//!
//! NCBI reads an object with its serial library (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`,
//! driven by the type information of the generated classes). The two encodings share their
//! structure: this module is the part of `CObjectIStream` that is common to both
//! (`ReadClassSequential`, `ReadChoiceSimple`, `ReadContainer`, the object frames); `text.rs`
//! and `ber.rs` are the two encodings' primitives.
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

use super::ber::BerIn;
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
use super::text::{TextIn, TextSource};
use super::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskPart, NullElement};
use crate::ncbi_exception::{NcbiException, Severity, ThrowSite};
use crate::report::serial::asn_text_reader::{PointerType, StreamSource};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::report::serial::object_frame::{
    frame_info, frame_name, text_frame_site, FrameKind, READ_POINTER_FRAME_END,
};
use crate::rpsbproc::{HOST_DEPENDENT_BYTES, USER_SPACE_BYTES};
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

/// What a primitive or a frame reports: the class of what NCBI's reader throws. The exceptions
/// are boxed (`Fail::serial`, `Fail::exception`), which keeps every primitive's result (`R`) as
/// small as a `String`: no exception is built on the paths that read a valid file.
#[derive(Debug)]
pub(super) enum Fail {
    /// NCBI's `CEofException` (see the module comment).
    Eof,
    /// A `CSerialException`, with the frame stack that `END_OBJECT_FRAME` adds to it.
    Serial(Box<NcbiException>),
    /// Another `CException` (the `CUtilException` of a number in text, the `CStringException` of
    /// a member's tag), with a line for each frame it leaves (`NCBI_RETHROW_SAME`).
    Exception(Box<NcbiException>),
    /// `std::bad_alloc` of the string that a length announces (`GetChars`' `reserve`).
    OutOfMemory,
    /// `std::length_error` of that `reserve`: `what()` is [`STRING_CREATE_LENGTH_ERROR`].
    LengthError,
    /// NCBI reads it; LOSAT does not.
    Unsupported(String),
}

pub(super) type R<T> = Result<T, Fail>;

impl Fail {
    /// [`Fail::Serial`] of `e`.
    #[cold]
    pub(super) fn serial(e: NcbiException) -> Fail {
        Fail::Serial(Box::new(e))
    }

    /// [`Fail::Exception`] of `e`.
    #[cold]
    pub(super) fn exception(e: NcbiException) -> Fail {
        Fail::Exception(Box::new(e))
    }
}

/// `ThrowError(flag, message)` of NCBI's object stream at `site` (`ThrowError1`,
/// objistr.cpp:828-853): a `CSerialException` whose message is the stream's position
/// (`GetPosition()`: `line N` in text, `byte N` in binary) and `message`. It returns the
/// [`Fail`], not an `R`: the path that calls it ends there (`return Err(throw_error(..))`), so
/// the reading path around it keeps nothing for its return.
#[cold]
#[inline(never)]
pub(super) fn throw_error(site: ThrowSite, position: &str, message: impl fmt::Display) -> Fail {
    Fail::serial(NcbiException::new(site, format!("{position}: {message}")))
}

/// A `ThrowError` of `CObjectIStream` (objistr.cpp) at `line` of `function`.
const fn istr_site(line: u32, kind: &'static str, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/serial/objistr.cpp",
        line,
        kind,
        function: Some(function),
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608
// ```c++
// void CObjectIStream::HandleEOF(CEofException& expt)
// {
//     string msg(TopFrame().GetFrameInfo());
//     PopFrame();
//     if (GetStackDepth() < 2) {
//         NCBI_RETHROW_SAME(expt,msg);
//     } else {
//         ThrowError(fEOF, msg);
//     }
// }
// ```
const HANDLE_EOF: ThrowSite = istr_site(
    606,
    "CSerialException::eEOF",
    "ncbi::CObjectIStream::HandleEOF",
);
const SKIP_FILE_HEADER: ThrowSite = istr_site(
    896,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::SkipFileHeader",
);
const EXPECTED_MEMBER: ThrowSite = istr_site(
    1037,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::ExpectedMember",
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:873-880
// ```c++
// const CReadObjectInfo&
// CObjectIStream::GetRegisteredObject(CReadObjectInfo::TObjectIndex index)
// {
//     if ( !m_Objects ) {
//         ThrowError(fFormatError,"invalid object index: NO_COLLECT defined");
//     }
//     return m_Objects->GetRegisteredObject(index);
// }
// ```
// (makeblastdb's streams collect no objects: `m_Objects` is null, checked with makeblastdb
// 2.17.0 on `@0` and `@1`.)
const GET_REGISTERED_OBJECT: ThrowSite = istr_site(
    877,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::GetRegisteredObject",
);

/// `what()` of the `std::length_error` that libstdc++'s `_M_create` throws for a capacity over
/// `max_size()` (`rpsbproc::prosite::CxxString`).
pub(super) const STRING_CREATE_LENGTH_ERROR: &str = "basic_string::_M_create";

/// The bytes of the file that one string of `length` bytes asks NCBI for (`GetChars`:
/// `str.reserve(count)` before the bytes are read, when they are not all in the input buffer):
/// a length over `max_size()` (2^63 - 1, `rpsbproc::prosite::CxxString`) is `std::length_error`;
/// an allocation of [`USER_SPACE_BYTES`] or more `std::bad_alloc` on every host; one of
/// [`HOST_DEPENDENT_BYTES`] or more succeeds or fails by the host's memory, which LOSAT refuses
/// (DECISIONS 40 of line PS: memory that NCBI builds from a length field of the input, as
/// `losat rpsbproc` does). Below that the allocation succeeds and the read goes on to the end of
/// the file.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:496-523
// ```c++
// void CIStreamBuffer::GetChars(string& str, size_t count)
//     THROWS1((CIOException))
// {
//     // cache pos
//     const char* pos = m_CurrentPos;
//     size_t in_buffer = m_DataEndPos - pos;
//     if ( in_buffer >= count ) {
//         // simplest case - plain copy
//         str.assign(pos, count);
//         m_CurrentPos = pos + count;
//         return;
//     }
//     str.reserve(count);
//     str.assign(pos, in_buffer);
//     for ( ;; ) {
//         count -= in_buffer;
//         m_CurrentPos = pos += in_buffer;
//         pos = FillBuffer(pos);
// ```
// (checked with makeblastdb 2.17.0 on lengths 2^62 - 1, 2^62, 2^63 - 2 and 2^63 - 1: "BLAST ran
// out of memory", status 4; 2^63: "Error: basic_string::_M_create", status 255.)
pub(super) fn reserve_string(length: u64) -> R<()> {
    if length > (1 << 63) - 1 {
        return Err(Fail::LengthError);
    }
    let bytes = length + 1;
    if bytes >= USER_SPACE_BYTES {
        return Err(Fail::OutOfMemory);
    }
    if bytes >= HOST_DEPENDENT_BYTES {
        return Err(host_dependent_string(length));
    }
    Ok(())
}

/// The refusal of a string of `length` bytes that `reserve_string` leaves to the host.
#[cold]
#[inline(never)]
fn host_dependent_string(length: u64) -> Fail {
    let bytes = length + 1;
    Fail::Unsupported(format!(
        "a string whose length ({length} bytes) passes the end of the file: makeblastdb \
         allocates {bytes} bytes for it, which succeeds or fails by the host's memory"
    ))
}

/// What reading a VisibleString posts, in the order NCBI posts it.
#[derive(Debug)]
pub(super) enum StringPost {
    /// The warning of a UTF8String tag for a VisibleString (binary,
    /// `ExpectStringTag`, objistrasnb.cpp:296-305).
    Utf8Tag,
    /// `ReplaceVisibleChar`'s post of a byte read as `#`: where, the byte, and the string (or
    /// piece of it) that the message quotes.
    BadChar {
        position: StreamPosition,
        byte: u8,
        string: Vec<u8>,
    },
}

/// A message that NCBI's reader posts while it reads a mask file.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum MaskPost {
    /// `ERR_POST_X(7)` of a byte read as `#`.
    BadChar(ReplacedChar),
    /// `ERR_POST_X_ONCE(10, Warning << ...)` of a UTF8String tag for the VisibleString member at
    /// the stack trace: once per process.
    Utf8Tag(String),
}

impl MaskPost {
    /// The line makeblastdb writes on its standard error (`ReplacedChar::post`,
    /// `ncbi_exception::post_line`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-310
    // ```c++
    //         static CSafeStatic<NCBI_PARAM_TYPE(SERIAL, READ_ANY_VISIBLESTRING_TAG)> sx_ReadAny;
    //         if ( sx_ReadAny->Get() ) {
    //             // may be eUTF8String
    //             TByte alt_tag = MakeTagByte(eUniversal,ePrimitive,eUTF8String);
    //             if ( PeekTagByte() == alt_tag ) {
    //                 // optionally issue a warning
    //                 if ( sx_ReadAny->Get() == 1 ) {
    //                     ERR_POST_X_ONCE(10, Warning<<
    //                                     "CObjectIStreamAsnBinary: UTF8String data for VisibleString member "<<GetStackTraceASN()<<", ASN.1 specification may need an update");
    //                 }
    // ```
    pub(crate) fn line(&self) -> Vec<u8> {
        match self {
            MaskPost::BadChar(replaced) => replaced.post(Some("makeblastdb")),
            MaskPost::Utf8Tag(trace) => crate::ncbi_exception::post_line(
                Severity::Warning,
                Some("803.10"),
                Some("makeblastdb"),
                format!(
                    "CObjectIStreamAsnBinary: UTF8String data for VisibleString member {trace}, \
                     ASN.1 specification may need an update"
                )
                .as_bytes(),
            ),
        }
    }
}

/// The primitives of one encoding (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`).
pub(super) trait Input {
    /// `GetPosition()`: `line N` in text, `byte N` (from the start of the object) in binary.
    fn position(&self) -> String;
    /// The depth of NCBI's object stack (`GetStackDepth()`), which the text reader's
    /// `SkipWhiteSpace` asks at the end of the text.
    #[inline]
    fn set_stack_depth(&mut self, _depth: usize) {}
    /// `ReadFileHeader`: the type name that precedes an object in text; binary files have none.
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>>;
    fn begin_class(&mut self) -> R<()>;
    /// `BeginClassMember(classType, pos)`: the next member, searched from `pos` on, or `None`
    /// at the end of the class.
    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>>;
    fn end_member(&mut self) -> R<()>;
    fn end_class(&mut self) -> R<()>;
    fn begin_choice(&mut self) -> R<()>;
    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize>;
    fn end_variant(&mut self) -> R<()>;
    fn end_choice(&mut self) -> R<()>;
    fn begin_container(&mut self, set: bool) -> R<()>;
    fn begin_element(&mut self) -> R<bool>;
    fn end_container(&mut self) -> R<()>;
    fn read_i32(&mut self) -> R<i32>;
    fn read_i64(&mut self) -> R<i64>;
    fn read_u32(&mut self) -> R<u32>;
    fn read_bool(&mut self) -> R<bool>;
    fn read_null(&mut self) -> R<()>;
    /// A `VisibleString` with each byte outside `' '..='~'` read as `#` (`eFNP_ReplaceAndWarn`,
    /// the default of `SERIAL_WRONG_CHARS_READ`), and what NCBI posts while it reads it, in
    /// `posts` (also when it then fails).
    fn read_string(&mut self, posts: &mut Vec<StringPost>) -> R<Vec<u8>>;
    /// `ReadEnum(values)`: the value (`CEnumeratedTypeInfo::ReadEnum` adds its line to an error).
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32>;
    /// `ReadPointerType()`, where NCBI reads a `CRef<>` that is empty at that moment
    /// (`Rd::read_pointer`): a `NULL` is taken, the `@` or `:` that starts the other two forms too.
    fn read_pointer_type(&mut self) -> R<PointerType>;
    /// `ReadObjectPointer()`: the object index after `@` (text: `ReadInt8`, `TObjectIndex` being
    /// 64 bits) or of the object reference tag (binary: `ReadStdSigned`). Its value does not matter:
    /// `GetRegisteredObject` throws for every index.
    fn read_object_pointer(&mut self) -> R<()>;
    /// `ReadOtherPointer()`: the name of the type whose value follows.
    fn read_other_pointer(&mut self) -> R<Vec<u8>>;
    /// `ReadOtherPointerEnd()`, after that value.
    fn read_other_pointer_end(&mut self) -> R<()>;
}

/// The deepest nesting of object frames that LOSAT follows (NCBI's stack is only limited by memory).
const MAX_DEPTH: usize = 400;

/// The frames that a reader reserves room for at its start: a mask file's objects nest about 20
/// deep (a `Dbtag` Seq-id of an interval of a `packed-int` mask), so the stack is not grown while
/// they are read.
const FRAMES_RESERVED: usize = 32;

/// The refusal of a frame past [`MAX_DEPTH`].
#[cold]
#[inline(never)]
fn too_deep() -> Fail {
    Fail::Unsupported(format!("objects nested more than {MAX_DEPTH} levels deep"))
}

/// A frame of NCBI's object stack as the mask reader keeps it (`object_frame`).
#[derive(Clone, Copy, Debug)]
struct Frame {
    kind: FrameKind,
    /// The named frame of `ReadPointer`, whose `END_OBJECT_FRAME` is not the text reader's for
    /// its kind (`end`).
    read_pointer: bool,
    /// The type of a named, class, choice or container frame (`Object type= `).
    type_name: Option<&'static str>,
    /// The member or variant of a member or variant frame (`SetTopMemberId`).
    member: Option<&'static str>,
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
    /// `eOtherPointer`: it follows inside `ReadPointer`'s named frame, entered at this depth.
    Other(usize),
}

/// How far `CSeq_loc::CheckId` walks one Seq-loc (`x_CheckId`, quoted at `read_seq_loc_parts`):
/// each Seq-loc starts a walk of its own (`my_id`, NULL), takes the Seq-ids of its parts in
/// order, and hands its own to the Seq-loc around it (`x_UpdateId(id, my_id)`) once its walk is
/// done. makeblastdb's `GetId` of a mask is this walk of the mask's Seq-loc.
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
    /// `x_UpdateId(id, &seq_id)` of a part.
    fn update(&mut self, seq_id: ParsedSeqId) {
        match self {
            IdWalk::Same(None) => *self = IdWalk::Same(Some(seq_id)),
            // `Equals` compares the values: the typed Seq-ids are equal exactly when it does.
            IdWalk::Same(Some(first)) if *first != seq_id => *self = IdWalk::Stopped,
            _ => {}
        }
    }

    /// `(*li)->CheckId(id, may_throw)` of an element of a `mix` or `equiv`: the element's own walk
    /// `inner`, if this walk gets that far.
    fn nested(&mut self, inner: IdWalk) {
        if !matches!(self, IdWalk::Same(_)) {
            return;
        }
        match inner {
            IdWalk::Same(Some(seq_id)) => self.update(seq_id),
            // `x_UpdateId(id, NULL)` keeps `id`.
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

/// The state of one `CObjectIStream`: the primitives, the object stack and what was posted.
pub(super) struct Rd<I: Input> {
    pub(super) inp: I,
    frames: Vec<Frame>,
    /// The messages posted since the caller last took them (`MaskFileReader::take_posts`).
    pub(super) posts: Vec<MaskPost>,
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:712-756
// ```c++
// void CMemberInfoFunctions::ReadWithSetFlagMember(CObjectIStream& in,
//                                                  const CMemberInfo* memberInfo,
//                                                  TObjectPtr classPtr)
// {
//     _ASSERT(!memberInfo->CanBeDelayed());
//     _ASSERT(memberInfo->HaveSetFlag());
//     memberInfo->UpdateSetFlagYes(classPtr);
//     try {
//         in.ReadObject(memberInfo->GetItemPtr(classPtr),
//                       memberInfo->GetTypeInfo());
//     }
//     catch (CSerialException& e) {
//         if (e.GetErrCode() == CSerialException::eNullValue) {
//             ...
//         } else if (e.GetErrCode() == CSerialException::eMissingValue) {
//             ...
//         } else {
//             NCBI_RETHROW_SAME(e,
//                 "error while reading " + memberInfo->GetId().GetName());
//         }
//     }
// }
// ```
// (`eNullValue` and `eMissingValue` are not thrown for these types: `ExpectedMember` throws
// `eFormatError` with the default `eSerialVerifyData_Yes`.)
const READ_WITH_SET_FLAG_MEMBER: (&str, u32, &str) = (
    "src/serial/member.cpp",
    752,
    "ncbi::CMemberInfoFunctions::ReadWithSetFlagMember",
);

impl<I: Input> Rd<I> {
    pub(super) fn new(inp: I) -> Self {
        Rd {
            inp,
            frames: Vec::with_capacity(FRAMES_RESERVED),
            posts: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:189-214 (quoted at
    // `report::serial::object_frame::text_frame_site`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608 (`HandleEOF`, quoted at
    // `HANDLE_EOF`)
    /// One `BEGIN_OBJECT_FRAME` ... `END_OBJECT_FRAME` of a frame of `kind` (of the type
    /// `type_name`): an error that leaves it gets NCBI's handling. A `CSerialException` takes the
    /// frame's name (`AddFrameInfo(GetFrameName())`); a `CEofException` is `HandleEOF`, which
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
        let depth = self.enter(Frame::new(kind, type_name))?;
        let result = f(self);
        self.exit(depth, result)
    }

    /// `BEGIN_OBJECT_FRAME` of `frame`: the depth below it, which [`exit`](Self::exit) takes
    /// once what the frame holds has been read. Every `enter` is followed by its `exit`, with no
    /// return between them.
    #[inline(always)]
    fn enter(&mut self, frame: Frame) -> R<usize> {
        let depth = self.frames.len();
        if depth >= MAX_DEPTH {
            return Err(too_deep());
        }
        self.frames.push(frame);
        self.inp.set_stack_depth(depth + 1);
        Ok(depth)
    }

    /// `END_OBJECT_FRAME` of the frame that [`enter`](Self::enter) entered at `depth`, after
    /// what it holds was read with `result` (an error gets [`leave`](Self::leave)).
    #[inline(always)]
    fn exit<T>(&mut self, depth: usize, result: R<T>) -> R<T> {
        match result {
            Ok(value) => {
                self.frames.truncate(depth);
                self.inp.set_stack_depth(depth);
                Ok(value)
            }
            Err(fail) => Err(self.leave(fail)),
        }
    }

    /// `END_OBJECT_FRAME` of the top frame for `fail`: the frame is left, then `fail` gets what
    /// leaving it adds.
    #[cold]
    #[inline(never)]
    fn leave(&mut self, fail: Fail) -> Fail {
        let frame = self.frames.pop().expect("a frame to leave");
        self.inp.set_stack_depth(self.frames.len());
        match fail {
            Fail::Serial(mut e) => {
                e.add_frame_info(frame_name(frame.kind, frame.type_name, frame.member));
                Fail::Serial(e)
            }
            Fail::Eof if self.frames.len() >= 2 => {
                let info = frame_info(frame.kind, frame.type_name, frame.member);
                Fail::serial(NcbiException::new(
                    HANDLE_EOF,
                    format!("{}: {info}", self.inp.position()),
                ))
            }
            Fail::Exception(e) => {
                let (file, line, function) = frame.end();
                let info = frame_info(frame.kind, frame.type_name, frame.member);
                Fail::exception(e.rethrow_same(file, line, function, info))
            }
            other => other,
        }
    }

    /// `SetTopMemberId`: the member or variant of the top frame.
    fn set_member(&mut self, name: &'static str) {
        if let Some(top) = self.frames.last_mut() {
            top.member = Some(name);
        }
    }

    /// `GetStackTrace()`: the path of the value being read (`Blast-db-mask-info.algo-options`).
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

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1038
    // ```c++
    // bool CObjectIStream::ExpectedMember(const CMemberInfo* memberInfo)
    // {
    //     const CItemInfo* info = CItemsInfo::FindNextMandatory(memberInfo);
    //     if (info) {
    //         if (GetVerifyData() == eSerialVerifyData_Yes) {
    //             ThrowError(fFormatError,
    //                     "member "+info->GetId().ToString()+" expected");
    // ```
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
            &self.inp.position(),
            format_args!("member {} expected", m.name),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397 (`ReadClassSequential`;
    // the binary stream's, objistrasnb.cpp:1210-1233, is the same with its own `BeginClass`)
    // ```c++
    //     BEGIN_OBJECT_FRAME3(eFrameClass, classType, classPtr);
    //     StartBlock();
    //
    //     ReadClassSequentialContentsBegin(classType);
    //
    //     TMemberIndex index;
    //     while ( (index = BeginClassMember(classType,*pos)) != kInvalidMember ) {
    //         ReadClassSequentialContentsMember(classPtr);
    //     }
    //
    //     ReadClassSequentialContentsEnd(classPtr);
    //
    //     EndBlock();
    //     END_OBJECT_FRAME();
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:84-107
    // ```c++
    // #define ClassSequentialContentsBegin(classType) \
    // { \
    //     CClassTypeInfo::CIterator pos(classType); \
    //     BEGIN_OBJECT_FRAME(eFrameClassMember); \
    //     {
    //
    // #define ClassSequentialContentsMember(Func, Args) \
    //     { \
    //         const CMemberInfo* memberInfo = classType->GetMemberInfo(index); \
    //         SetTopMemberId(memberInfo->GetId()); \
    //         for ( TMemberIndex i = *pos; i < index; ++i ) { \
    //             classType->GetMemberInfo(i)->NCBI_NAME2(Func,MissingMember)Args; \
    //         } \
    //         { \
    //             memberInfo->NCBI_NAME2(Func,Member)Args; \
    //         } \
    //         pos.SetIndex(index + 1); \
    //     }
    //
    // #define ClassSequentialContentsEnd(Func, Args) \
    //     } \
    //     END_OBJECT_FRAME(); \
    //     for ( ; pos.Valid(); ++pos ) { \
    //         classType->GetMemberInfo(pos)->NCBI_NAME2(Func,MissingMember)Args; \
    //     } \
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:871-876
    // ```c++
    // void CMemberInfoFunctions::ReadMissingSimpleMember(CObjectIStream& in,
    //                                                    const CMemberInfo* memberInfo,
    //                                                    TObjectPtr classPtr)
    // {
    //     _ASSERT(!memberInfo->Optional());
    //     in.ExpectedMember(memberInfo);
    // ```
    /// A SEQUENCE read member by member, in declaration order (`ReadClassSequential`): a member
    /// that is not in the file is an error unless it is `OPTIONAL` or has a `DEFAULT`. `read`
    /// reads the value of member `index`; a member with a set flag reads it in
    /// `ReadWithSetFlagMember` (quoted at [`READ_WITH_SET_FLAG_MEMBER`]). An `OPTIONAL` member held
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
            s.frame(FrameKind::ClassMember, None, |s| {
                while let Some(index) = s.inp.begin_member(def, pos)? {
                    let member = &def.members[index];
                    s.set_member(member.name);
                    s.missing_members(&def.members[pos..index])?;
                    // `ReadPointer` of a member held by `CRef<>` around `read`, which is called at
                    // this one place. Such a member has no set flag (`Member::set_flag`): the
                    // errors of `begin_pointer` are not `ReadWithSetFlagMember`'s.
                    let pointee = if member.pointer {
                        s.begin_pointer(member.ty.type_name())?
                    } else {
                        Pointee::This
                    };
                    let result = match pointee {
                        Pointee::Null => Ok(()),
                        Pointee::This | Pointee::Other(_) => {
                            let result = read(s, index);
                            s.end_pointer(pointee, result)
                        }
                    };
                    if member.set_flag {
                        result.map_err(|fail| match fail {
                            Fail::Serial(e) => {
                                let (file, line, function) = READ_WITH_SET_FLAG_MEMBER;
                                let message = format!("error while reading {}", member.name);
                                Fail::serial(e.rethrow_same(file, line, function, message))
                            }
                            other => other,
                        })?;
                    } else {
                        result?;
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
    // ```c++
    // void CObjectIStream::ReadChoiceSimple(const CChoiceTypeInfo* choiceType,
    //                                       TObjectPtr choicePtr)
    // {
    //     BEGIN_OBJECT_FRAME3(eFrameChoice, choiceType, choicePtr);
    //     BeginChoice(choiceType);
    //     BEGIN_OBJECT_FRAME(eFrameChoiceVariant);
    //     TMemberIndex index = BeginChoiceVariant(choiceType);
    //     if ( index == kInvalidMember ) {
    //         ThrowError(fFormatError,"choice variant id expected");
    //     }
    //     const CVariantInfo* variantInfo = choiceType->GetVariantInfo(index);
    //     SetTopMemberId(variantInfo->GetId());
    //
    //     variantInfo->ReadVariant(*this, choicePtr);
    //
    //     EndChoiceVariant();
    //     END_OBJECT_FRAME();
    //     EndChoice();
    //     END_OBJECT_FRAME();
    // ```
    /// A CHOICE: `read` reads the value of the variant.
    pub(super) fn read_choice<T>(
        &mut self,
        def: &'static ChoiceDef,
        read: impl FnOnce(&mut Self, usize) -> R<T>,
    ) -> R<T> {
        self.frame(FrameKind::Choice, Some(def.name), |s| {
            s.inp.begin_choice()?;
            let value = s.frame(FrameKind::ChoiceVariant, None, |s| {
                let index = s.inp.begin_variant(def)?;
                s.set_member(def.variants[index].name);
                let value = read(s, index)?;
                s.inp.end_variant()?;
                Ok(value)
            })?;
            s.inp.end_choice()?;
            Ok(value)
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1318-1346
    // ```c++
    // void CObjectIStream::ReadContainer(const CContainerTypeInfo* containerType,
    //                                    TObjectPtr containerPtr)
    // {
    //     BEGIN_OBJECT_FRAME2(eFrameArray, containerType);
    //     BeginContainer(containerType);
    //
    //     TTypeInfo elementType = containerType->GetElementType();
    //     BEGIN_OBJECT_FRAME2(eFrameArrayElement, elementType);
    //     ...
    //     while ( BeginContainerElement(elementType) ) {
    //         ...
    //         EndContainerElement();
    //     }
    //
    //     END_OBJECT_FRAME();
    //
    //     EndContainer();
    //     END_OBJECT_FRAME();
    // }
    // ```
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
    // ```c++
    // void CPointerTypeInfo::ReadPointer(CObjectIStream& in,
    //                                    TTypeInfo objectType,
    //                                    TObjectPtr objectPtr)
    // {
    //     ...
    //     TObjectPtr pointedPtr = pointerType->GetObjectPointer(objectPtr);
    //     if ( pointedPtr ) {
    //         //pointedType->SetDefault(pointedPtr);
    //         in.ReadObject(pointedPtr, pointedType);
    //     }
    //     else {
    //         pointerType->SetObjectPointer(objectPtr,
    //                                       in.ReadPointer(pointedType).first);
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1133-1213
    // ```c++
    // pair<TObjectPtr, TTypeInfo> CObjectIStream::ReadPointer(TTypeInfo declaredType)
    // {
    //     ...
    //     switch ( ReadPointerType() ) {
    //     case eNullPointer:
    //         return pair<TObjectPtr, TTypeInfo>((TObjectPtr)0, declaredType);
    //     case eObjectPointer:
    //         {
    //             TObjectIndex index = ReadObjectPointer();
    //             const CReadObjectInfo& info = GetRegisteredObject(index);
    //             ...
    //         }
    //     case eThisPointer:
    //         {
    //             ...
    //             RegisterObject(objectPtr, declaredType);
    //             ReadObject(objectPtr, declaredType);
    //             ...
    //             return make_pair(objectPtr, declaredType);
    //         }
    //     case eOtherPointer:
    //         {
    //             string className = ReadOtherPointer();
    //             objectType = MapType(className);
    //
    //             BEGIN_OBJECT_FRAME2(eFrameNamed, objectType);
    //             ...
    //             RegisterObject(objectPtr, objectType);
    //             ReadObject(objectPtr, objectType);
    //             ...
    //             END_OBJECT_FRAME();
    //
    //             ReadOtherPointerEnd();
    //             break;
    //         }
    //     ...
    //     }
    //     while ( objectType != declaredType ) {
    //         // try to check parent class pointer
    //         if ( objectType->GetTypeFamily() != eTypeFamilyClass ) {
    //             ThrowError(fFormatError,"incompatible member type");
    //         }
    //         ...
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/classinfob.cpp:244-252
    // ```c++
    // TTypeInfo CClassTypeInfoBase::GetClassInfoByName(const string& name)
    // {
    //     TClassesByName& classes = ClassesByName();
    //     pair<TClassesByName::iterator, TClassesByName::iterator> i =
    //         classes.equal_range(name);
    //     if (  i.first == i.second ) {
    //         NCBI_THROW_FMT(CSerialException, eInvalidData,
    //                        "class not found: "<<name);
    //     }
    // ```
    // (`MapType` is `GetClassInfoByName`. The classes it knows are those whose type information
    // the process has built so far: the declared type is one of them, another type may or may not
    // be (makeblastdb 2.17.0: "class not found: Seq-interval" for `:Seq-interval` as the first
    // mask; a type it knows is read whole and then fails with "incompatible member type", or with
    // the first error inside it when the value is not of that type).)
    /// `CPointerTypeInfo::ReadPointer` of a `CRef<>` to a value of the type named `declared` that
    /// is empty when NCBI reads it: an `OPTIONAL` member held by `CRef<>` (the generated classes
    /// preallocate the mandatory ones) and each element of a list of `CRef`s (`AddElementIn`
    /// pushes an empty `CRef`, stltypes.hpp:246-263). `read` reads the value (`This`, and `Other`
    /// naming `declared`, inside `ReadPointer`'s named frame); `None` is a `NULL`. `@n` is NCBI's
    /// error; a value of another type is refused (`Fail::Unsupported`), as what NCBI does with
    /// it depends on what its process has read before.
    fn read_pointer<T>(
        &mut self,
        declared: &'static str,
        read: impl FnOnce(&mut Self) -> R<T>,
    ) -> R<Option<T>> {
        let pointee = self.begin_pointer(declared)?;
        if let Pointee::Null = pointee {
            return Ok(None);
        }
        let result = read(self);
        self.end_pointer(pointee, result).map(Some)
    }

    /// [`read_pointer`](Self::read_pointer) up to the value: where it is ([`Pointee`]), the named
    /// frame of `Other` entered. [`end_pointer`](Self::end_pointer) follows the value.
    fn begin_pointer(&mut self, declared: &'static str) -> R<Pointee> {
        match self.inp.read_pointer_type()? {
            PointerType::Null => Ok(Pointee::Null),
            PointerType::This => Ok(Pointee::This),
            PointerType::Object => {
                self.inp.read_object_pointer()?;
                Err(throw_error(
                    GET_REGISTERED_OBJECT,
                    &self.inp.position(),
                    "invalid object index: NO_COLLECT defined",
                ))
            }
            PointerType::Other => {
                let name = self.inp.read_other_pointer()?;
                if name != declared.as_bytes() {
                    return Err(self.other_type(&name, declared));
                }
                let frame = Frame {
                    read_pointer: true,
                    ..Frame::new(FrameKind::Named, Some(declared))
                };
                Ok(Pointee::Other(self.enter(frame)?))
            }
        }
    }

    /// The refusal of `:name` (a value of the type `name`) where a `declared` belongs.
    #[cold]
    #[inline(never)]
    fn other_type(&self, name: &[u8], declared: &'static str) -> Fail {
        let name = String::from_utf8_lossy(name);
        Fail::Unsupported(format!(
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
            Pointee::Other(depth) => {
                let value = self.exit(depth, result)?;
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
            position: self.inp.position(),
        }
    }

    /// A `VisibleString`: each byte outside `' '..='~'` read as `#`, with NCBI's posts
    /// (`StringPost`, here with the stack trace).
    fn read_str(&mut self) -> R<Vec<u8>> {
        let mut posts = Vec::new();
        let result = self.inp.read_string(&mut posts);
        if !posts.is_empty() {
            let trace = self.stack_trace();
            for post in posts {
                self.posts.push(match post {
                    StringPost::Utf8Tag => MaskPost::Utf8Tag(trace.clone()),
                    StringPost::BadChar {
                        position,
                        byte,
                        string,
                    } => MaskPost::BadChar(ReplacedChar {
                        trace: trace.clone(),
                        byte,
                        position,
                        string: string.into(),
                    }),
                });
            }
        }
        result
    }

    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:416-432
    // ```c++
    // void CEnumeratedTypeInfo::ReadEnum(CObjectIStream& in,
    //                                    TTypeInfo objectType,
    //                                    TObjectPtr objectPtr)
    // {
    //     const CEnumeratedTypeInfo* enumType =
    //         CTypeConverter<CEnumeratedTypeInfo>::SafeCast(objectType);
    //     try {
    //         enumType->m_ValueType->SetValueInt(objectPtr,
    //                                            in.ReadEnum(enumType->Values()));
    //     }
    //     catch ( CException& e ) {
    //         NCBI_RETHROW_SAME(e,"invalid enum value");
    //     }
    // ```
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
                Fail::serial(e.rethrow_same(file, line, function, "invalid enum value"))
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

    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3541-3620
    // ```c++
    // bool CSeq_loc::x_CheckId(const CSeq_id*& id, bool may_throw) const
    // {
    //     switch ( Which() ) {
    //     case e_not_set:
    //     case e_Null:
    //         {
    //             // no Seq-id
    //             return true;
    //         }
    //     case e_Empty:
    //         {
    //             return x_UpdateId(id, &GetEmpty(), may_throw);
    //         }
    //     case e_Whole:
    //         {
    //             return x_UpdateId(id, &GetWhole(), may_throw);
    //         }
    //     case e_Int: ... x_UpdateId(id, &loc.GetId(), may_throw);
    //     case e_Pnt: ... x_UpdateId(id, &pnt.GetId(), may_throw);
    //     case e_Packed_int:
    //         {
    //             // Check ID of each interval
    //             const CPacked_seqint& ints = GetPacked_int();
    //             ITERATE ( CPacked_seqint::Tdata, ii, ints.Get() ) {
    //                 const CSeq_interval& loc = **ii;
    //                 if ( !x_UpdateId(id, &loc.GetId(), may_throw) ) {
    //                     return false;
    //                 }
    //             }
    //             return true;
    //         }
    //     case e_Packed_pnt: ... x_UpdateId(id, &pnts.GetId(), may_throw);
    //     case e_Mix:
    //         {
    //             // Check ID of each sub-location.
    //             const CSeq_loc_mix& mix = GetMix();
    //             ITERATE( CSeq_loc_mix::Tdata, li, mix.Get() ) {
    //                 if ( !(*li)->CheckId(id, may_throw) ) {
    //                     return false;
    //                 }
    //             }
    //             return true;
    //         }
    //     case e_Bond: ... A, then B when it is set ...
    //     case e_Equiv: ... (*li)->CheckId(id, may_throw) for each, as e_Mix ...
    //     case e_Feat:
    //     default:
    //         if (may_throw) { NCBI_THROW_FMT(CSeqLocException, eUnsupported, ...); }
    //         return false;
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:403-422
    // ```c++
    // bool CSeq_loc::x_UpdateId(const CSeq_id*& total_id, const CSeq_id* id,
    //                           bool may_throw) const
    // {
    //     if ( total_id == id ) {
    //         return true;
    //     }
    //
    //     if ( !total_id ) {
    //         total_id = id;
    //     } else if ( (id  &&  !total_id->Equals(*id)) ) {
    //         if (may_throw) {
    //             ...
    //         } else {
    //             return false;
    //         }
    //     }
    //     return true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_loc.hpp:926-945
    // ```c++
    // inline
    // bool CSeq_loc::CheckId(const CSeq_id*& id, bool may_throw) const
    // {
    //     const CSeq_id* my_id = m_IdCache.load(memory_order_acquire);
    //     if ( my_id == NULL ) {
    //         if ( !x_CheckId(my_id, may_throw) ) {
    //             return false;
    //         }
    //         m_IdCache.store(my_id, memory_order_release);
    //     }
    //     return x_UpdateId(id, my_id, may_throw);
    // }
    //
    //
    // inline
    // const CSeq_id* CSeq_loc::GetId(void) const
    // {
    //     const CSeq_id* sip = NULL;
    //     return CheckId(sip, false) ? sip : NULL;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:1439-1555,1559-1600
    // ```c++
    // void CSeq_loc_CI_Impl::x_ProcessLocation(const CSeq_loc& loc)
    // {
    //     switch ( loc.Which() ) {
    //     case CSeq_loc::e_not_set:
    //     case CSeq_loc::e_Null:
    //     case CSeq_loc::e_Empty:
    //         ... (a range only with eEmpty_Allow: TRange::GetEmpty(), the Seq-id of `empty`) ...
    //     case CSeq_loc::e_Whole:
    //         ... info.m_Range = TRange::GetWhole(); ...
    //     case CSeq_loc::e_Int:  x_ProcessInterval(loc.GetInt(), loc);
    //     case CSeq_loc::e_Pnt:  x_ProcessPoint(loc.GetPnt(), loc);
    //     case CSeq_loc::e_Packed_int: ... x_ProcessInterval(**ii, loc) for each ...
    //     case CSeq_loc::e_Packed_pnt:
    //         ... strand of the pack; for each point: info.m_Range.Set(*pi, *pi) ...
    //     case CSeq_loc::e_Mix:   ... x_ProcessLocation(**li) for each ...
    //     case CSeq_loc::e_Equiv: ... x_ProcessLocation(**li) for each part ...
    //     case CSeq_loc::e_Bond:
    //         x_ProcessPoint(bond.GetA(), loc);
    //         if ( bond.IsSetB() ) { x_ProcessPoint(bond.GetB(), loc); }
    //     default:
    //         NCBI_THROW_FMT(CSeqLocException, eUnsupported, ...);
    //     }
    // }
    // ```
    /// One `Seq-loc` (any variant, nested `mix` and `equiv` included): appends the ranges
    /// `CSeq_loc_CI` visits to `parts`, in its order (a `whole` Seq-loc is the interval 0 to
    /// 4294967294, whose `CRange` is `TRange::GetWhole()`; a point the interval of one position),
    /// and returns how far `CSeq_loc::CheckId` walks it ([`IdWalk`]). An element of `packed-int`,
    /// `mix` or `equiv` stated as `NULL` (an empty `CRef<>`) is where the walk dereferences it,
    /// if it gets that far (`**ii`, `(*li)->CheckId`).
    fn read_seq_loc_parts(&mut self, parts: &mut Vec<MaskPart>) -> R<IdWalk> {
        let mut walk = IdWalk::Same(None);
        self.read_choice(&SEQ_LOC, |s, index| {
            match index {
                SEQ_LOC_NULL => {
                    s.inp.read_null()?;
                    parts.push(MaskPart::Null);
                }
                SEQ_LOC_EMPTY => {
                    let seq_id = s.read_seq_id()?;
                    walk.update(seq_id);
                    parts.push(MaskPart::Empty);
                }
                SEQ_LOC_WHOLE => {
                    let seq_id = s.read_seq_id()?;
                    walk.update(seq_id);
                    parts.push(MaskPart::Interval(MaskInterval {
                        from: 0,
                        to: u32::MAX - 1,
                        strand: None,
                    }));
                }
                SEQ_LOC_INT => {
                    let (interval, seq_id) = s.read_seq_interval()?;
                    walk.update(seq_id);
                    parts.push(MaskPart::Interval(interval));
                }
                SEQ_LOC_PACKED_INT => {
                    s.read_container(false, |s| {
                        match s.read_pointer(SEQ_INTERVAL.name, |s| s.read_seq_interval())? {
                            Some((interval, seq_id)) => {
                                walk.update(seq_id);
                                parts.push(MaskPart::Interval(interval));
                            }
                            None => walk.null(s.null_element(LIST_PACKED_INT)),
                        }
                        Ok(())
                    })?;
                }
                SEQ_LOC_PNT => {
                    let (interval, seq_id) = s.read_seq_point()?;
                    walk.update(seq_id);
                    parts.push(MaskPart::Interval(interval));
                }
                SEQ_LOC_PACKED_PNT => {
                    let (mut strand, mut seq_id, mut points) = (None, None, Vec::new());
                    let mut scratch = Vec::new();
                    s.read_class(&PACKED_SEQPNT, |s, member| {
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
                    let list = if index == SEQ_LOC_MIX {
                        LIST_MIX
                    } else {
                        LIST_EQUIV
                    };
                    s.read_container(index == SEQ_LOC_EQUIV, |s| {
                        match s.read_pointer(SEQ_LOC.name, |s| s.read_seq_loc_parts(parts))? {
                            Some(inner) => walk.nested(inner),
                            None => walk.null(s.null_element(list)),
                        }
                        Ok(())
                    })?;
                }
                SEQ_LOC_BOND => {
                    // `b` is an OPTIONAL `CRef<>` (`read_class`): a `NULL` leaves it unset, which
                    // `x_CheckId` and `CSeq_loc_CI` skip (`CanGetB`, `IsSetB`).
                    s.read_class(&SEQ_BOND, |s, _| {
                        let (interval, seq_id) = s.read_seq_point()?;
                        walk.update(seq_id);
                        parts.push(MaskPart::Interval(interval));
                        Ok(())
                    })?;
                }
                _ => {
                    // `feat`: read; `x_CheckId` fails, so the Seq-loc has no Seq-id.
                    let mut scratch = Vec::new();
                    s.read_value(SEQ_LOC.variants[index].ty, &mut scratch)?;
                    walk = IdWalk::Stopped;
                }
            }
            Ok(())
        })?;
        Ok(walk)
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
                    let loc = match s.read_pointer(SEQ_LOC.name, |s| s.read_seq_loc())? {
                        Some(loc) => loc,
                        None => Err(s.null_element(LIST_MASKS)),
                    };
                    if null_element.is_none() {
                        match loc {
                            Ok(loc) => masks.push(loc),
                            Err(null) => null_element = Some(null),
                        }
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

    /// `file >> MSerial_Asn... >> *obj` for a `Blast-db-mask-info`.
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

    /// `file >> MSerial_Asn... >> *obj` for a `Blast-mask-list` (the continuation).
    fn read_list(&mut self) -> R<MaskList> {
        let mut list = None;
        self.read_top(&BLAST_MASK_LIST, |s| {
            list = Some(s.read_mask_list_members()?);
            Ok(())
        })?;
        Ok(mandatory(list))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-901
    // ```c++
    // // root reader
    // void CObjectIStream::SkipFileHeader(TTypeInfo typeInfo)
    // {
    //     ResetState();
    //     ...
    //     BEGIN_OBJECT_FRAME2(eFrameNamed, typeInfo);
    //
    //     string name = ReadFileHeader();
    //     const string& tname = typeInfo->GetName();
    //     if ( !name.empty() && !tname.empty() && name != tname ) {
    //         ThrowError(fFormatError,
    //                    "incompatible type "+name+"<>"+typeInfo->GetName());
    //     }
    //
    //     END_OBJECT_FRAME();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:971-975
    // ```c++
    // void CObjectIStream::Read(TObjectPtr object, TTypeInfo typeInfo)
    // {
    //     // root object
    //     SkipFileHeader(typeInfo);
    //     Read(object, typeInfo, eNoFileHeader);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:959-969
    // ```c++
    // void CObjectIStream::Read(TObjectPtr object, TTypeInfo typeInfo, ENoFileHeader)
    // {
    //     // root object
    //     BEGIN_OBJECT_FRAME2(eFrameNamed, typeInfo);
    //
    //     ReadObject(object, typeInfo);
    //
    //     EndOfRead();
    //
    //     END_OBJECT_FRAME();
    // }
    // ```
    /// `file >> MSerial_Asn... >> *obj` of the type `def`: `SkipFileHeader`, then `Read` with
    /// `read`, each in a frame of its own.
    fn read_top(&mut self, def: &'static ClassDef, read: impl FnOnce(&mut Self) -> R<()>) -> R<()> {
        self.frames.clear();
        self.frame(FrameKind::Named, Some(def.name), |s| {
            match s.inp.read_file_header()? {
                Some(name) if !name.is_empty() && name != def.name.as_bytes() => Err(throw_error(
                    SKIP_FILE_HEADER,
                    &s.inp.position(),
                    format_args!(
                        "incompatible type {}<>{}",
                        String::from_utf8_lossy(&name),
                        def.name
                    ),
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
    let mut rd = Rd::new(BerIn::new(data));
    rd.frame(FrameKind::Named, Some("Seq-id"), |s| s.read_seq_id())
}

/// [`seq_id_from_binary`] for ASN.1 text, for the test that compares the builder with
/// `seqid::asn_text`.
#[cfg(test)]
pub(super) fn seq_id_from_text(data: &[u8]) -> R<ParsedSeqId> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let mut rd = Rd::new(TextIn::new(SliceSource::new(data)));
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

/// Where the objects of a mask file come from.
enum Source {
    /// ASN.1 text, read from the file in blocks by one reader for all objects (`text.rs`): only
    /// the bytes of the object being read are kept.
    Text(Box<Rd<TextIn<TextSource>>>),
    /// ASN.1 binary: the bytes of the file and the offset of the next object (`ber.rs` reads a
    /// slice).
    Binary { data: Vec<u8>, offset: usize },
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
            Source::Binary { data, offset } => debug
                .field("format", &MaskFileFormat::BinaryAsn)
                .field("bytes", &data.len())
                .field("offset", offset),
        };
        debug.field("state", &self.state).finish()
    }
}

impl MaskFileReader {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:933-943
    // ```c++
    //         CNcbiIfstream mask_file(file_list[i].c_str(), ios::binary);
    //         CFormatGuess::EFormat mask_file_format = CFormatGuess::eUnknown;
    //         {{
    //              CFormatGuess fg(mask_file);
    //              fg.GetFormatHints().AddPreferredFormat(CFormatGuess::eBinaryASN);
    //              fg.GetFormatHints().AddPreferredFormat(CFormatGuess::eTextASN);
    //              fg.GetFormatHints().DisableAllNonpreferred();
    //              mask_file_format = fg.GuessFormat();
    //          }}
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:160-176
    // ```c++
    //     switch (fmt) {
    //     ...
    //     default:
    //         NCBI_THROW(CInvalidDataException, eInvalidInput, string("Unknown encoding for ") + msg);
    // ```
    /// [`open_named`](Self::open_named) for a file that messages call by its path.
    #[cfg(test)]
    pub(crate) fn open(path: &Path) -> Result<Self, MaskFileError> {
        Self::open_named(path, path)
    }

    /// Opens the file and guesses its format from its first bytes; messages call the file
    /// `name` (the name on makeblastdb's command line). The text of the error for a file in
    /// neither format is NCBI's. A file that cannot be opened or read is a failed stream for the
    /// format guess (`x_TestInput`), which answers `eUnknown`: the same error. A text file is
    /// then read in blocks as its objects are read; a binary file is read whole (a read error
    /// there is the same error).
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

    /// The reader of the file `input`; `path` is only used in messages.
    fn from_reader(path: &Path, mut input: Box<dyn Read>) -> Result<Self, MaskFileError> {
        let Some((format, sample)) = guess_format(&mut input) else {
            return Err(unknown_encoding(path));
        };
        let source = match format {
            MaskFileFormat::TextAsn => Source::Text(Box::new(Rd::new(TextIn::new(
                StreamSource::new(input, &sample),
            )))),
            MaskFileFormat::BinaryAsn => {
                let mut data = sample;
                if input.read_to_end(&mut data).is_err() {
                    return Err(unknown_encoding(path));
                }
                Source::Binary { data, offset: 0 }
            }
        };
        Ok(MaskFileReader {
            path: path.to_path_buf(),
            source,
            state: State::Open,
            posts: Vec::new(),
        })
    }

    /// The bytes of the file that the reader keeps: the stream's buffer of a text file, the
    /// whole of a binary one.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        match &self.source {
            Source::Text(rd) => rd.inp.buffered(),
            Source::Binary { data, .. } => data.len(),
        }
    }

    #[cfg(test)]
    pub(crate) fn format(&self) -> MaskFileFormat {
        match self.source {
            Source::Text(_) => MaskFileFormat::TextAsn,
            Source::Binary { .. } => MaskFileFormat::BinaryAsn,
        }
    }

    fn run<T>(
        &mut self,
        text: impl FnOnce(&mut Rd<TextIn<TextSource>>) -> R<T>,
        binary: impl FnOnce(&mut Rd<BerIn<'_>>) -> R<T>,
    ) -> R<T> {
        match &mut self.source {
            Source::Text(rd) => {
                rd.inp.new_object();
                let result = text(rd);
                self.posts.append(&mut rd.posts);
                result
            }
            Source::Binary { data, offset } => {
                let start = (*offset).min(data.len());
                let mut rd = Rd::new(BerIn::new(&data[start..]));
                let result = binary(&mut rd);
                *offset = start + rd.inp.offset();
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
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    //     ...
    //     catch (const std::bad_alloc&) {                                         \
    //         LOG_POST(Error << "BLAST ran out of memory");                       \
    //         exit_code = BLAST_OUT_OF_MEMORY;                                    \
    //     }                                                                       \
    //     catch (const std::exception& e) {                                       \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    /// The error of the reader's `fail`: what makeblastdb's `CATCH_ALL` reports for what NCBI's
    /// reader threw.
    fn fail(&mut self, fail: Fail) -> MaskFileError {
        self.state = if matches!(fail, Fail::Eof) {
            State::FailedEof
        } else {
            State::Failed
        };
        match fail {
            Fail::Serial(e) | Fail::Exception(e) => {
                MaskFileError::Ncbi(WriteDbError::Exception(*e))
            }
            Fail::OutOfMemory => MaskFileError::Ncbi(WriteDbError::OutOfMemory),
            Fail::LengthError => {
                MaskFileError::Ncbi(WriteDbError::Other(STRING_CREATE_LENGTH_ERROR.to_string()))
            }
            Fail::Unsupported(what) => MaskFileError::Unsupported(what),
            Fail::Eof => MaskFileError::EofInContinuation {
                path: self.path.clone(),
            },
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:947-960
    // ```c++
    //         int algo_id = -1;
    //         while (true) {
    //             CRef<CBlast_db_mask_info> first_obj;
    //
    //             try {
    //                 s_ReadObject(mask_file, mask_file_format, first_obj, "mask data in '" + file_list[i] + "'");
    //             }
    //             catch (CEofException&) {
    //                 // must be end of file
    //                 break;
    //             }
    // ```
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

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:994-996
    // ```c++
    //                 if (! masks->GetMore())
    //                     break;
    //
    //                 s_ReadObject(mask_file, mask_file_format, masks, "mask data (continuation)");
    // ```
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
