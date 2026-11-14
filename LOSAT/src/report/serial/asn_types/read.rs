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
//!   `ERR_POST_X(7)` (`Post::BadChar`), a UTF8String as it is;
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

use super::{resolved, Kind, Member, Storage, StringKind, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_text_reader::{
    real_of_parts, AsnSource, AsnTextError, AsnTextErrorKind, MemberId, Reader, SliceSource,
};
use crate::report::serial::encoder::{ReplacedChar, StreamPosition};

/// A message that NCBI's reader posts while it reads, in the order posted.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Post {
    /// `ReplaceVisibleChar`'s `ERR_POST_X(7, message)` of a byte read as `#`.
    BadChar(ReplacedChar),
    /// `ERR_POST(Error << text)` of a class's `PostRead` (`CVariation_ref::PostRead`).
    Error(&'static str),
}

impl Post {
    /// The line that a BLAST application of NCBI posts on its standard error (its diagnostic
    /// prefix `program`, the program's name in lower case), with its line feed.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5672-5695 (quoted at
    // `ncbi_exception::NcbiException::report_default_ex`: `<severity>: (<code>.<subcode>) `)
    pub(crate) fn line(&self, program: &str) -> Vec<u8> {
        let mut line = b"Error: ".to_vec();
        match self {
            Post::BadChar(replaced) => {
                let (code, subcode) = ReplacedChar::ERR_CODE;
                line.extend_from_slice(format!("({code}.{subcode}) [{program}] ").as_bytes());
                line.extend_from_slice(&replaced.message());
            }
            Post::Error(text) => {
                line.extend_from_slice(format!("[{program}] {text}").as_bytes());
            }
        }
        line.push(b'\n');
        line
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
}

/// What the reader read, and the messages it posted before it stopped or finished.
pub(crate) type Read<T> = (std::result::Result<T, ReadError>, Vec<Post>);

type R<T> = std::result::Result<T, ReadError>;

impl From<AsnTextError> for ReadError {
    fn from(error: AsnTextError) -> Self {
        match error.kind {
            AsnTextErrorKind::Eof => ReadError::Eof,
            AsnTextErrorKind::Format | AsnTextErrorKind::Unsupported => {
                ReadError::Serial(error.message)
            }
        }
    }
}

/// A frame of NCBI's object stack (`CObjectStackFrame::EFrameType`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FrameKind {
    Named,
    Class,
    ClassMember,
    Choice,
    ChoiceVariant,
    Array,
    ArrayElement,
}

#[derive(Clone, Copy, Debug)]
struct Frame {
    kind: FrameKind,
    /// The type the frame names (`Object type= `): the named, class, choice and array frames.
    type_name: Option<&'static str>,
    /// The member or variant (`SetTopMemberId`) of a member or variant frame.
    member: Option<&'static str>,
    /// Whether the frame belongs to a skip function (`SkipClassSequential`, ...).
    skip: bool,
}

// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:366-381
// ```c++
// const char* CObjectStackFrame::GetFrameTypeName(void) const
// {
//     const char* s;
//     switch (GetFrameType())
//     {
//     default:                  s = "UNKNOWN";             break;
//     case eFrameOther:         s = "eFrameOther";         break;
//     case eFrameNamed:         s = "eFrameNamed";         break;
//     case eFrameArray:         s = "eFrameArray";         break;
//     case eFrameArrayElement:  s = "eFrameArrayElement";  break;
//     case eFrameClass:         s = "eFrameClass";         break;
//     case eFrameClassMember:   s = "eFrameClassMember";   break;
//     case eFrameChoice:        s = "eFrameChoice";        break;
//     case eFrameChoiceVariant: s = "eFrameChoiceVariant"; break;
//     }
//     return s;
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:399-410
// ```c++
// string CObjectStackFrame::GetFrameInfo(void) const
// {
//     string info(" Frame type= ");
//     info += GetFrameTypeName();
//     if (m_TypeInfo) {
//         info += ", Object type= " + m_TypeInfo->GetName();
//     }
//     if (m_MemberId) {
//         info += ", Member name= " + m_MemberId->GetName();
//     }
//     return info;
// }
// ```
impl Frame {
    fn info(&self) -> String {
        let kind = match self.kind {
            FrameKind::Named => "eFrameNamed",
            FrameKind::Class => "eFrameClass",
            FrameKind::ClassMember => "eFrameClassMember",
            FrameKind::Choice => "eFrameChoice",
            FrameKind::ChoiceVariant => "eFrameChoiceVariant",
            FrameKind::Array => "eFrameArray",
            FrameKind::ArrayElement => "eFrameArrayElement",
        };
        let mut info = format!(" Frame type= {kind}");
        if let Some(name) = self.type_name {
            info.push_str(", Object type= ");
            info.push_str(name);
        }
        if let Some(member) = self.member {
            info.push_str(", Member name= ");
            info.push_str(member);
        }
        info
    }

    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves this frame: the
    /// place of the line that `NCBI_RETHROW_SAME` adds to an exception passing through it.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
    // ```c++
    // #define END_OBJECT_FRAME_OF(Stream) \
    //     } catch (CSerialException& s_expt) { \
    //         std::string msg((Stream).TopFrame().GetFrameName()); \
    //         (Stream).PopFrame(); \
    //         s_expt.AddFrameInfo(msg); \
    //         throw; \
    //     } catch ( CEofException& e_expt ) { \
    //         (Stream).HandleEOF(e_expt); \
    //     } catch (CException& expt) { \
    //         std::string msg((Stream).TopFrame().GetFrameInfo()); \
    //         (Stream).PopFrame(); \
    //         NCBI_RETHROW_SAME(expt,msg); \
    //     } \
    //     (Stream).PopFrame()
    // ```
    // The lines are those of `END_OBJECT_FRAME()` in the functions that push the frames:
    // `CObjectIStream::Read` (objistr.cpp:959-969), `CObjectIStreamAsn::ReadClassSequential`
    // (objistrasn.cpp:1380-1397, the member frame ends in `ReadClassSequentialContentsEnd` at
    // 1393), `CObjectIStream::ReadChoiceSimple` (objistr.cpp:1482-1501),
    // `CObjectIStreamAsn::ReadContainer` (objistrasn.cpp:1245-1274), and the skip functions
    // `SkipClassSequential` (objistrasn.cpp:1417-1433), `SkipChoiceSimple` (objistr.cpp:1503-1522)
    // and `SkipContainer` (objistrasn.cpp:1276-1292).
    fn site(&self) -> (&'static str, u32, &'static str) {
        const ASN: &str = "src/serial/objistrasn.cpp";
        const ISTR: &str = "src/serial/objistr.cpp";
        match (self.kind, self.skip) {
            (FrameKind::Named, _) => (ISTR, 968, "ncbi::CObjectIStream::Read"),
            (FrameKind::Class, false) => {
                (ASN, 1396, "ncbi::CObjectIStreamAsn::ReadClassSequential")
            }
            (FrameKind::ClassMember, false) => {
                (ASN, 1393, "ncbi::CObjectIStreamAsn::ReadClassSequential")
            }
            (FrameKind::Choice, false) => (ISTR, 1500, "ncbi::CObjectIStream::ReadChoiceSimple"),
            (FrameKind::ChoiceVariant, false) => {
                (ISTR, 1498, "ncbi::CObjectIStream::ReadChoiceSimple")
            }
            (FrameKind::Array, false) => (ASN, 1273, "ncbi::CObjectIStreamAsn::ReadContainer"),
            (FrameKind::ArrayElement, false) => {
                (ASN, 1270, "ncbi::CObjectIStreamAsn::ReadContainer")
            }
            (FrameKind::Class, true) => (ASN, 1432, "ncbi::CObjectIStreamAsn::SkipClassSequential"),
            (FrameKind::ClassMember, true) => {
                (ASN, 1429, "ncbi::CObjectIStreamAsn::SkipClassSequential")
            }
            (FrameKind::Choice, true) => (ISTR, 1521, "ncbi::CObjectIStream::SkipChoiceSimple"),
            (FrameKind::ChoiceVariant, true) => {
                (ISTR, 1519, "ncbi::CObjectIStream::SkipChoiceSimple")
            }
            (FrameKind::Array, true) => (ASN, 1291, "ncbi::CObjectIStreamAsn::SkipContainer"),
            (FrameKind::ArrayElement, true) => {
                (ASN, 1288, "ncbi::CObjectIStreamAsn::SkipContainer")
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:570-584
// ```c++
// void CIStreamBuffer::BadNumber(void)
//     THROWS1((CUtilException))
// {
//     m_Error = "bad number";
//     NCBI_THROW_FMT(CUtilException, eWrongData,
//                    "bad number in line " << GetLine());
// }
//
//
// void CIStreamBuffer::NumberOverflow(void)
//     THROWS1((CUtilException))
// {
//     m_Error = "number overflow";
//     NCBI_THROW_FMT(CUtilException, eWrongData,
//                    "number overflow in line " << GetLine());
// }
// ```
const BAD_NUMBER: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 574,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::BadNumber"),
};
const NUMBER_OVERFLOW: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 583,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::NumberOverflow"),
};

/// Reads the ASN.1 text document at the start of `text` as the root type `root` (an object of
/// `CObjectIStream::Read`: its file header, then the value; what follows it is not read), with
/// the messages posted.
pub(crate) fn read_document(text: &[u8], root: u16) -> Read<Value> {
    let mut reader = TypedReader {
        r: Reader::new(SliceSource::new(text), "PSSM", "PSSM reader"),
        frames: Vec::new(),
        posts: Vec::new(),
    };
    let value = reader.document(root);
    (value, reader.posts)
}

struct TypedReader<'a> {
    r: Reader<SliceSource<'a>>,
    frames: Vec<Frame>,
    posts: Vec<Post>,
}

impl TypedReader<'_> {
    fn push(&mut self, kind: FrameKind, type_name: Option<&'static str>, skip: bool) {
        self.frames.push(Frame {
            kind,
            type_name,
            member: None,
            skip,
        });
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

    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-118
    // ```c++
    // string CObjectStack::GetStackTraceASN(void) const
    // {
    //     ...
    //     string stack = FetchFrameFromBottom(0).HasTypeInfo() ?
    //         FetchFrameFromBottom(0).m_TypeInfo->GetName() : "?";
    //     for ( size_t i = 1; i < GetStackDepth(); ++i ) {
    //         const TFrame& frame = FetchFrameFromBottom(i);
    //         switch ( frame.m_FrameType ) {
    //         case TFrame::eFrameClassMember:
    //         case TFrame::eFrameChoiceVariant:
    //             {
    //                 if ( !frame.m_MemberId ) {
    //                     _ASSERT(i == GetStackDepth() - 1);
    //                 }
    //                 else {
    //                     const CMemberId& id = *frame.m_MemberId;
    //                     stack += '.';
    //                     if ( !id.GetName().empty() ) {
    //                         stack += id.GetName();
    //                     }
    //                     ...
    //         case TFrame::eFrameArrayElement:
    //             stack += ".E";
    //             break;
    // ```
    /// `GetStackTrace()`: the path of the object being read (`PssmWithParameters.pssm.query.seq.
    /// descr.E.comment`).
    fn stack_trace(&self) -> String {
        let mut stack = self
            .frames
            .first()
            .and_then(|frame| frame.type_name)
            .unwrap_or("?")
            .to_string();
        for frame in self.frames.iter().skip(1) {
            match frame.kind {
                FrameKind::ClassMember | FrameKind::ChoiceVariant => {
                    if let Some(member) = frame.member {
                        stack.push('.');
                        stack.push_str(member);
                    }
                }
                FrameKind::ArrayElement => stack.push_str(".E"),
                _ => {}
            }
        }
        stack
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-900 and 959-969
    // ```c++
    // void CObjectIStream::SkipFileHeader(TTypeInfo typeInfo)
    // {
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
    // ...
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
    fn document(&mut self, root: u16) -> R<Value> {
        let name = super::type_info(root).name;
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

    /// One value of the type `ty` (`ReadObject`), or its skip (`SkipObject`, `skip`), then the
    /// hooks of NCBI's classes that a read runs.
    fn value(&mut self, ty: u16, skip: bool) -> R<Value> {
        let info = resolved(ty);
        let target = alias_target(ty);
        let value = match info.kind {
            Kind::Sequence(members) => self.sequence(info.name, members, skip)?,
            Kind::Choice(variants) => self.choice(info.name, variants, skip)?,
            Kind::List { element, .. } => self.list(element, skip)?,
            Kind::Alias(_) => unreachable!("resolved"),
            Kind::Integer(storage) => self.integer(storage, skip)?,
            Kind::Enumerated { integer, values } => self.enumerated(values, integer, skip)?,
            Kind::Boolean => Value::Bool(self.r.boolean()?),
            Kind::Null => {
                self.r.null()?;
                Value::Null
            }
            Kind::Real => self.real(skip)?,
            Kind::String(kind) => Value::Str(self.string(kind, skip)?),
            Kind::Octets => self.octets(skip)?,
        };
        if skip {
            return Ok(value);
        }
        super::hooks::post_read(target, value, &mut self.posts)
    }

    /// `ReadClassSequential` / `SkipClassSequential`: the class frame, `{`, the member frame, the
    /// members in module order (`BeginClassMember(classType, pos)`) with the mandatory members
    /// before each one checked, the member frame left, the mandatory members after the last one
    /// checked, `}`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:86-111
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
    fn sequence(&mut self, name: &'static str, members: &'static [Member], skip: bool) -> R<Value> {
        self.push(FrameKind::Class, Some(name), skip);
        let result = self.sequence_in_frame(members, skip);
        self.pop(result)
    }

    fn sequence_in_frame(&mut self, members: &'static [Member], skip: bool) -> R<Value> {
        self.r.start_block()?;
        self.push(FrameKind::ClassMember, None, skip);
        let result = self.members(members, skip);
        let (values, pos) = self.pop(result)?;
        for member in &members[pos..] {
            expected_member(member)?;
        }
        self.r.end_block()?;
        Ok(Value::Sequence(values))
    }

    fn members(&mut self, members: &'static [Member], skip: bool) -> R<(Vec<(u16, Value)>, usize)> {
        let mut values = Vec::new();
        let mut pos = 0;
        let mut block_start = true;
        while self.r.next_element(&mut block_start)? {
            let index = self.item_index(members, pos, true)?;
            let member = &members[index];
            self.set_member(member.name);
            for missing in &members[pos..index] {
                expected_member(missing)?;
            }
            let value = if member.flags & super::DELAYED != 0 && !skip {
                self.delayed(member.ty)?
            } else {
                self.value(member.ty, skip)?
            };
            if !skip {
                values.push((index as u16, value));
            }
            pos = index + 1;
        }
        Ok((values, pos))
    }

    /// `ReadMemberId(SkipWhiteSpace())` and `GetMemberIndex(classType, id, pos)` (from `pos`,
    /// `in_class`) or `GetChoiceIndex(choiceType, id)`: the member or variant named or tagged.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:407-424
    // ```c++
    //     TMemberIndex idx;
    //     if (!id.empty()  &&  isdigit((unsigned char) id[0])) {
    //         idx = classType->GetMembers().Find
    //             (CMemberId::TTag(NStr::StringToInt(id)), CAsnBinaryDefs::eContextSpecific, pos);
    //     }
    //     else {
    //         idx = classType->GetMembers().Find(id, pos);
    //         if (idx == kInvalidMember) {
    //             idx = GetAltItemIndex(classType,id,pos);
    //         }
    //     }
    //     return idx;
    // ```
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
                    .position(|item| item.name.as_bytes() == id)
                    .map(|at| at + pos)
            }
            MemberId::Tag(range) => {
                // NStr::StringToInt(id) with the default flags: its `CStringException` leaves
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

    /// The delayed member: `StartDelayBuffer` (white space skipped), the value skipped, the text
    /// kept (`Value::Delayed`).
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    // ```c++
    // void CMemberInfoFunctions::ReadLongMember(CObjectIStream& in,
    //                                           const CMemberInfo* memberInfo,
    //                                           TObjectPtr classPtr)
    // {
    //     if ( memberInfo->CanBeDelayed() ) {
    //         CDelayBuffer& buffer = memberInfo->GetDelayBuffer(classPtr);
    //         if ( !buffer ) {
    //             if (!in.ShouldParseDelayBuffer()) {
    //                 memberInfo->UpdateSetFlagYes(classPtr);
    //                 in.StartDelayBuffer();
    //                 memberInfo->GetTypeInfo()->SkipData(in);
    //                 in.EndDelayBuffer(buffer, memberInfo, classPtr);
    //                 return;
    //             }
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1648-1652
    // ```c++
    // void CObjectIStreamAsn::StartDelayBuffer(void)
    // {
    //     SkipWhiteSpace();
    //     CObjectIStream::StartDelayBuffer();
    // }
    // ```
    // (psiblast and tblastn set no hook on their stream: `ShouldParseDelayBuffer` is false,
    // objistr.cpp:745-760.)
    fn delayed(&mut self, ty: u16) -> R<Value> {
        self.r.skip_white_space()?;
        let start = self.r.pos;
        self.value(ty, true)?;
        Ok(Value::Delayed(self.r.bytes(start..self.r.pos).to_vec()))
    }

    /// `ReadChoiceSimple` / `SkipChoiceSimple`.
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1501
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
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1475 (`BeginChoiceVariant`,
    // quoted at `Reader::choice`: `NextElement` at the choice's start, then the variant's id)
    fn choice(&mut self, name: &'static str, variants: &'static [Member], skip: bool) -> R<Value> {
        self.push(FrameKind::Choice, Some(name), skip);
        self.push(FrameKind::ChoiceVariant, None, skip);
        let result = self.variant(variants, skip);
        let result = self.pop(result);
        self.pop(result)
    }

    fn variant(&mut self, variants: &'static [Member], skip: bool) -> R<Value> {
        let mut block_start = true;
        if !self.r.next_element(&mut block_start)? {
            return Err(ReadError::Serial("choice variant id expected".into()));
        }
        let index = self.item_index(variants, 0, false)?;
        let variant = &variants[index];
        self.set_member(variant.name);
        let value = self.value(variant.ty, skip)?;
        Ok(Value::Choice(index as u16, Box::new(value)))
    }

    /// `ReadContainer` / `SkipContainer`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274
    // ```c++
    // void CObjectIStreamAsn::ReadContainer(const CContainerTypeInfo* contType,
    //                                       TObjectPtr contPtr)
    // {
    //     BEGIN_OBJECT_FRAME2(eFrameArray, contType);
    //     StartBlock();
    //
    //     BEGIN_OBJECT_FRAME(eFrameArrayElement);
    //     ...
    //     while ( NextElement() ) {
    //         ...
    //             contType->AddElement(contPtr, *this);
    //     }
    //     ...
    //     END_OBJECT_FRAME();
    //
    //     EndBlock();
    //     END_OBJECT_FRAME();
    // }
    // ```
    // The container's type information is anonymous (an implicit class such as `Seq-descr` has
    // no frame of its own, `ReadImplicitMember`, classinfo.cpp:289-306): `Object type= ` names
    // nothing (c07 of SURVEY_D3).
    fn list(&mut self, element: u16, skip: bool) -> R<Value> {
        self.push(FrameKind::Array, Some(""), skip);
        let result = self.list_in_frame(element, skip);
        self.pop(result)
    }

    fn list_in_frame(&mut self, element: u16, skip: bool) -> R<Value> {
        self.r.start_block()?;
        self.push(FrameKind::ArrayElement, None, skip);
        let result = self.elements(element, skip);
        let values = self.pop(result)?;
        self.r.end_block()?;
        Ok(Value::List(values))
    }

    fn elements(&mut self, element: u16, skip: bool) -> R<Vec<Value>> {
        let mut values = Vec::new();
        let mut block_start = true;
        while self.r.next_element(&mut block_start)? {
            let value = self.value(element, skip)?;
            if !skip {
                values.push(value);
            }
        }
        Ok(values)
    }

    /// `CIStreamBuffer::SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
    fn get_after_ws(&mut self) -> R<u8> {
        loop {
            let c = self.r.peek_char(0)?;
            self.r.pos += 1;
            if !crate::ncbi_str::is_space(c) {
                return Ok(c);
            }
        }
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
    // ```c++
    // Int4 CIStreamBuffer::GetInt4(void)
    //     THROWS1((CIOException,CUtilException))
    // {
    //     bool sign;
    //     char c = SkipWs();
    //     switch ( c ) {
    //     case '-':
    //         sign = true;
    //         c = GetChar();
    //         break;
    //     case '+':
    //         sign = false;
    //         c = GetChar();
    //         break;
    //     default:
    //         sign = false;
    //         break;
    //     }
    //
    //     Uint4 n = c - '0';
    //     if ( n > 9 )
    //         BadNumber();
    //
    //     // overflow limits
    //     const Uint4 kMaxBeforeMul = kMax_I4/10;
    //     const Uint1 kMaxLimitAdd = Uint1(kMax_I4%10 + sign);
    //
    //     for ( ;; ) {
    //         Uint1 d = (Uint1)(PeekCharNoEOF() - '0');
    //         if ( d > 9 )
    //             break;
    //         SkipChar();
    //
    //         // check multiplication overflow
    //         if ( n > kMaxBeforeMul || (n == kMaxBeforeMul && d > kMaxLimitAdd) )
    //             NumberOverflow();
    //
    //         n = n * 10 + d;
    //     }
    //     if ( sign )
    //         return -Int4(n);
    //     else
    //         return n;
    // }
    // ```
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        self.r.skip_white_space()?;
        let mut c = self.get_after_ws()?;
        let mut negative = false;
        if signed && c == b'-' {
            negative = true;
            c = self.r.peek_char(0)?;
            self.r.pos += 1;
        } else if c == b'+' {
            c = self.r.peek_char(0)?;
            self.r.pos += 1;
        }
        if !c.is_ascii_digit() {
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        loop {
            let d = self.r.peek_at(0);
            if !d.is_ascii_digit() {
                break;
            }
            self.r.pos += 1;
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
            }
            n = n * 10 + d;
        }
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
    // ```c++
    // Int2 CObjectIStream::ReadInt2(void)
    // {
    //     Int4 data = ReadInt4();
    //     Int2 ret = Int2(data);
    //     if ( ret != data )
    //         ThrowError(fOverflow, "integer overflow");
    //     return ret;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1780-1808 (the skips of every size:
    // `SkipSNumber` for the signed ones, `SkipUNumber` for the unsigned ones)
    fn integer(&mut self, storage: Storage, skip: bool) -> R<Value> {
        if skip {
            self.skip_number(storage != Storage::Uint4)?;
            return Ok(Value::Null);
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
        Ok(Value::Int(value))
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1048
    // ```c++
    // void CObjectIStreamAsn::SkipSNumber(void)
    // {
    //     size_t i;
    //     char c = SkipWhiteSpace();
    //     switch ( c ) {
    //     case '-':
    //     case '+':
    //         c = m_Input.PeekChar(1);
    //         // next char
    //         i = 2;
    //         break;
    //     default:
    //         // next char
    //         i = 1;
    //         break;
    //     }
    //     if ( c < '0' || c > '9' ) {
    //         ThrowError(fFormatError, "bad signed integer in line "
    //             + NStr::SizetToString(m_Input.GetLine()));
    //     }
    //     while ( (c = m_Input.PeekChar(i)) >= '0' && c <= '9' ) {
    //         ++i;
    //     }
    //     m_Input.SkipChars(i);
    // }
    // ```
    // (`SkipUNumber`, 1050-1074, takes only `+` and peeks past the digits with `PeekCharNoEOF`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        let mut c = self.r.skip_white_space()?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = self.r.peek_char(1)?;
            i = 2;
        }
        if !c.is_ascii_digit() {
            return Err(ReadError::Serial(format!(
                "bad {} integer in line {}",
                if signed { "signed" } else { "unsigned" },
                self.r.line
            )));
        }
        loop {
            let c = if signed {
                self.r.peek_char(i)?
            } else {
                self.r.peek_at(i)
            };
            if !c.is_ascii_digit() {
                break;
            }
            i += 1;
        }
        self.r.pos += i;
        Ok(())
    }

    /// `CEnumeratedTypeInfo::ReadEnum` / `SkipEnum` over `ReadEnum` of the stream: a lower-case
    /// name of the list, or a number (one of the list for an ENUMERATED); an exception that is not
    /// a `CSerialException` gets the line of the type's `NCBI_RETHROW_SAME`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    // ```c++
    // TEnumValueType CObjectIStreamAsn::ReadEnum(const CEnumeratedTypeValues& values)
    // {
    //     // not integer
    //     CTempString id = ReadLCaseId(SkipWhiteSpace());
    //     if ( !id.empty() ) {
    //         // enum element by name
    //         return values.FindValue(id);
    //     }
    //     // enum element by value
    //     TEnumValueType value = m_Input.GetInt4();
    //     if ( !values.IsInteger() ) // check value
    //         values.FindName(value, false);
    //
    //     return value;
    // }
    // ```
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
    // (`SkipEnum`, 468-484, reads the value the same way and rethrows at 478.)
    fn enumerated(
        &mut self,
        values: &'static [(&'static str, i32)],
        integer: bool,
        skip: bool,
    ) -> R<Value> {
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

    fn enum_value(&mut self, values: &'static [(&'static str, i32)], integer: bool) -> R<i32> {
        let c = self.r.skip_white_space()?;
        let id = self.r.scan_end_of_id(c.is_ascii_lowercase())?;
        if !id.is_empty() {
            let id = self.r.bytes(id);
            return values
                .iter()
                .find(|(name, _)| name.as_bytes() == id)
                .map(|&(_, number)| number)
                .ok_or_else(|| ReadError::Serial("invalid value of enumerated type".into()));
        }
        let value = self.int4()?;
        if !integer && !values.iter().any(|&(_, number)| number == value) {
            return Err(ReadError::Serial("invalid value of enumerated type".into()));
        }
        Ok(value)
    }

    /// `ReadDouble` / `SkipFNumber`: a word, or `{ mantissa, base, exponent }` whose base is read
    /// with `ReadUint4` (a `CUtilException` for a `-`) and exponent with `ReadInt4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (quoted at `Reader::real`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
    // ```c++
    // void CObjectIStreamAsn::SkipFNumber(void)
    // {
    //     if (PeekChar(true) != '{') {
    //         ScanEndOfId(true);
    //         return;
    //     }
    //     Expect('{', true);
    //     SkipSNumber();
    //     Expect(',', true);
    //     unsigned base = ReadUint4();
    //     Expect(',', true);
    //     SkipSNumber();
    //     Expect('}', true);
    //     if ( base != 2 && base != 10 )
    //         ThrowError(fFormatError, "illegal REAL base (must be 2 or 10)");
    // }
    // ```
    fn real(&mut self, skip: bool) -> R<Value> {
        if skip {
            if self.r.skip_white_space()? != b'{' {
                self.r.scan_end_of_id(true)?;
                return Ok(Value::Null);
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
            return Ok(Value::Null);
        }
        if let Some(word) = self.r.real_word()? {
            return Ok(Value::Real(word));
        }
        let (mantissa, is_negative) = self.r.real_mantissa()?;
        self.r.expect(b',')?;
        let base = self.uint4()?;
        self.r.expect(b',')?;
        let exp = self.int4()?;
        self.r.expect(b'}')?;
        real_of_parts(mantissa, is_negative, base, exp)
            .map(Value::Real)
            .ok_or_else(|| ReadError::Serial("illegal REAL base (must be 2 or 10)".into()))
    }

    /// `ReadStringValue` with the fix of `kind` / `SkipString`: the string's lines joined, `""` as
    /// one quote, and for a VisibleString each byte outside 0x20..0x7E read as `#` with NCBI's
    /// message, which quotes the piece of the string being appended (`AppendStringData`: up to the
    /// end of the line, the quote, or 128 bytes; none when skipping).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:877-912 (`AppendStringData`)
    // ```c++
    //     const char* data = m_Input.GetCurrentPos();
    //     if ( fix_method != eFNP_Allow ) {
    //         size_t done = 0, valid = 0;
    //         for ( size_t i = 0; i < count; ++i ) {
    //             char c = data[i];
    //             if ( !GoodVisibleChar(c) ) {
    //                 ...
    //                 if (valid == 0) {
    //                     if ( i > done ) {
    //                         s.append(data + done, i - done);
    //                     }
    //                     c = ReplaceVisibleChar(c, fix_method, this, CTempString(data,count), x_FixCharsSubst());
    //                     if (c != 0) {
    //                         s += c;
    //                     }
    //                     done = i + 1;
    //                 }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-986 (`ReadStringValue`, quoted
    // at `Reader::string`: a piece is appended at an end of line, at the quote, and every 128
    // bytes; after `""` the next piece starts with the second quote)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1093-1155
    // ```c++
    // void CObjectIStreamAsn::SkipString(EStringType type)
    // {
    //     ...
    //             default:
    //                 if (type == eStringTypeVisible) {
    //                     if ( !GoodVisibleChar(c) ) {
    //                         ...
    //                         if (valid == 0) {
    //                             c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, kEmptyStr, x_FixCharsSubst());
    //                         }
    // ```
    // (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0, objistr.hpp:1163: `valid` stays 0;
    // the fix is `eFNP_ReplaceAndWarn`, the default of `SERIAL_WRONG_CHARS_READ`, objistr.cpp:302.)
    fn string(&mut self, kind: StringKind, skip: bool) -> R<Vec<u8>> {
        let visible = kind != StringKind::Utf8;
        self.r.expect(b'"')?;
        let start_line = self.r.line;
        let unclosed = |line: usize| {
            ReadError::Serial(format!(
                "unclosed string starts at line {start_line} ({line})"
            ))
        };
        let mut out = Vec::new();
        let mut i = 0usize;
        loop {
            let Some(c) = self.r.src.byte_at(self.r.pos + i) else {
                return Err(unclosed(self.r.line));
            };
            match c {
                b'\r' | b'\n' => {
                    self.append(&mut out, i, visible, skip);
                    self.r.pos += i + 1;
                    i = 0;
                    self.r.skip_end_of_line(c);
                }
                b'"' => {
                    self.append(&mut out, i, visible, skip);
                    self.r.pos += i + 1;
                    if self.r.peek_at(0) != b'"' {
                        return Ok(out);
                    }
                    // `""`: the second quote starts the next piece (`SkipString` takes both).
                    if skip {
                        self.r.pos += 1;
                        i = 0;
                    } else {
                        i = 1;
                    }
                }
                _ => {
                    if skip && visible && !(b' '..=b'~').contains(&c) {
                        self.post_bad_char(c, &[]);
                    }
                    i += 1;
                    if i == 128 {
                        self.append(&mut out, i, visible, skip);
                        self.r.pos += i;
                        i = 0;
                    }
                }
            }
        }
    }

    /// `AppendStringData` of the `count` bytes at the reader's position (not taken).
    fn append(&mut self, out: &mut Vec<u8>, count: usize, visible: bool, skip: bool) {
        if skip {
            return;
        }
        let data = self.r.bytes(self.r.pos..self.r.pos + count).to_vec();
        if !visible {
            out.extend_from_slice(&data);
            return;
        }
        let mut done = 0;
        for (i, &c) in data.iter().enumerate() {
            if !(b' '..=b'~').contains(&c) {
                out.extend_from_slice(&data[done..i]);
                self.post_bad_char(c, &data);
                out.push(b'#');
                done = i + 1;
            }
        }
        out.extend_from_slice(&data[done..]);
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (quoted at
    // `report::serial::encoder::ReplacedChar`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:807-810
    // ```c++
    // string CObjectIStream::GetPosition(void) const
    // {
    //     return "line "+NStr::SizetToString(m_Input.GetLine());
    // }
    // ```
    // (`CObjectIStreamAsn` keeps `GetPosition` of its base; `CIStreamBuffer::GetLine` is the line
    // of the reader.)
    fn post_bad_char(&mut self, byte: u8, string: &[u8]) {
        self.posts.push(Post::BadChar(ReplacedChar {
            trace: self.stack_trace(),
            byte,
            position: StreamPosition::Line(self.r.line as u64),
            string: string.into(),
        }));
    }

    /// `ReadBytes` / `SkipByteBlock` (the skip takes `H` or `B` after white space).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
    // ```c++
    // void CObjectIStreamAsn::SkipByteBlock(void)
    // {
    //     Expect('\'', true);
    //     for ( ;; ) {
    //         char c = GetChar();
    //         if ( c >= '0' && c <= '9' ) {
    //             continue;
    //         }
    //         ...
    //         else if ( c == '\'' ) {
    //             break;
    //         }
    //         else if ( c == '\r' || c == '\n' ) {
    //             SkipEndOfLine(c);
    //         }
    //         else {
    //             m_Input.UngetChar(c);
    //             ThrowError(fFormatError, "bad char in octet string: #"
    //                 + NStr::IntToString(c));
    //         }
    //     }
    //     Expect('H', 'B', true);
    // }
    // ```
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
// ```c++
// bool CObjectIStream::ExpectedMember(const CMemberInfo* memberInfo)
// {
//     const CItemInfo* info = CItemsInfo::FindNextMandatory(memberInfo);
//     if (info) {
//         if (GetVerifyData() == eSerialVerifyData_Yes) {
//             ThrowError(fFormatError,
//                     "member "+info->GetId().ToString()+" expected");
// ```
fn expected_member(member: &Member) -> R<()> {
    if member.optional() {
        return Ok(());
    }
    Err(ReadError::Serial(format!(
        "member {} expected",
        member.name
    )))
}

/// The type `ty` refers to after its aliases (the index of `resolved(ty)`), for the hooks.
fn alias_target(mut ty: u16) -> u16 {
    while let Kind::Alias(target) = super::type_info(ty).kind {
        ty = target;
    }
    ty
}
