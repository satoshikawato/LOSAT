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
//! Where an error is raised matters for one thing, the end of the file. NCBI's input buffer throws
//! a `CEofException` when a read passes the end. Every object frame (`END_OBJECT_FRAME`) turns it
//! into a `CSerialException` (an error) unless fewer than two frames remain after it, in which case
//! the `CEofException` goes on to `x_ProcessMaskData`, which treats it as the end of the file and drops the
//! object being read. That happens at the type name of a text object, at its opening brace, and at the tag, the
//! length and the closing bytes of the outermost binary SEQUENCE; the `Fail::Eof` below is that exception.

use std::fmt;
use std::io::Read;
use std::path::{Path, PathBuf};

use super::ber::BerIn;
use super::guess::guess_format;
use super::schema::{
    ChoiceDef, ClassDef, EnumDef, Ty, BLAST_DB_MASK_INFO, BLAST_FILTER_PROGRAM, BLAST_MASK_LIST,
    DATE, DATE_STD, DBTAG, GIIMPORT_ID, ID_PAT, ID_PAT_ID, NA_STRAND, OBJECT_ID_ID, PACKED_SEQPNT,
    PATENT_SEQ_ID, PDB_SEQ_ID, SEQ_BOND, SEQ_ID, SEQ_INTERVAL, SEQ_INTERVAL_FROM, SEQ_INTERVAL_ID,
    SEQ_INTERVAL_STRAND, SEQ_INTERVAL_TO, SEQ_LOC, SEQ_LOC_BOND, SEQ_LOC_EMPTY, SEQ_LOC_EQUIV,
    SEQ_LOC_INT, SEQ_LOC_MIX, SEQ_LOC_NULL, SEQ_LOC_PACKED_INT, SEQ_LOC_PACKED_PNT, SEQ_LOC_PNT,
    SEQ_LOC_WHOLE, SEQ_POINT, TEXTSEQ_ID, T_OBJECT_ID,
};
use super::text::{TextIn, TextSource};
use super::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskPart};
use crate::report::serial::asn_text_reader::StreamSource;
use crate::report::serial::visible::BadChar;
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
    /// An error whose text NCBI's own code writes (the "Unknown encoding ..."
    /// `CInvalidDataException`); carries the `WriteDbError` makeblastdb would end with.
    Ncbi(WriteDbError),
    /// The file is not a valid object for NCBI's serial reader; NCBI's serial-library text is
    /// not reproduced (makeblastdb ends with status 255 and a `CSerialException` or
    /// `CUtilException` text that names the NCBI build directory).
    Malformed { path: PathBuf, message: String },
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
            MaskFileError::Malformed { path, message } => {
                write!(f, "mask data in '{}': {message}", path.display())
            }
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

/// What a primitive or a frame reports.
#[derive(Debug)]
pub(super) enum Fail {
    /// NCBI's `CEofException` (see the module comment).
    Eof,
    /// The serial reader throws: `CSerialException`, `CUtilException` (bad number, overflow).
    Format(String),
    /// NCBI reads it; LOSAT does not.
    Unsupported(String),
}

pub(super) type R<T> = Result<T, Fail>;

pub(super) fn format_err<T>(message: impl Into<String>) -> R<T> {
    Err(Fail::Format(message.into()))
}

/// The primitives of one encoding (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`).
pub(super) trait Input {
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
    /// A `VisibleString` (line breaks of a text string gone) with each byte outside `' '..='~'`
    /// stored as `#`; `notes` gets what NCBI posts for it, in its order (also when the read then
    /// fails): a UTF8String tag of a binary file, then each replaced byte.
    fn read_string(&mut self, notes: &mut Vec<StringNote>) -> R<Vec<u8>>;
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32>;
}

/// What a string read posts in NCBI (`Input::read_string`).
#[derive(Debug)]
pub(super) enum StringNote {
    /// A binary file gave a VisibleString member the tag of a UTF8String.
    Utf8Tag,
    /// A byte stored as `#` (`report::serial::visible`).
    Bad(BadChar),
}

/// A diagnostic that NCBI's reader posts while it reads a mask file (with `ERR_POST`, standard
/// error, in the order read; the status stays 0). `process_mask_data` writes them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum MaskDiagnostic {
    /// `ReplaceVisibleChar`'s message for a byte stored as `#` in a VisibleString (`ERR_POST_X(7,
    /// ...)`, an error): the trace, "Bad char [0x..] in string at line N" (text) or "at byte N"
    /// (binary), the chunk.
    BadChar(Vec<u8>),
    /// A UTF8String tag for the VisibleString member whose trace this is (`ERR_POST_X_ONCE(10,
    /// Warning << ...)`, once per process); the string is still read as a VisibleString.
    Utf8Tag(String),
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-305
// ```c++
//             if ( PeekTagByte() == alt_tag ) {
//                 // optionally issue a warning
//                 if ( sx_ReadAny->Get() == 1 ) {
//                     ERR_POST_X_ONCE(10, Warning<<
//                                     "CObjectIStreamAsnBinary: UTF8String data for VisibleString member "<<GetStackTraceASN()<<", ASN.1 specification may need an update");
//                 }
// ```
impl MaskDiagnostic {
    /// The lines NCBI's application `program` posts (without the final newline).
    pub(crate) fn text(&self, program: &str) -> Vec<u8> {
        match self {
            MaskDiagnostic::BadChar(message) => {
                let mut text = format!("Error: (803.7) [{program}] ").into_bytes();
                text.extend_from_slice(message);
                text
            }
            MaskDiagnostic::Utf8Tag(trace) => format!(
                "Warning: (803.10) [{program}] CObjectIStreamAsnBinary: UTF8String data for VisibleString member {trace}, ASN.1 specification may need an update"
            )
            .into_bytes(),
        }
    }
}

/// The deepest nesting of object frames that LOSAT follows (NCBI's stack is only limited by memory).
const MAX_DEPTH: usize = 400;

/// The state of one `CObjectIStream`: the primitives, and the number of object frames.
pub(super) struct Rd<I: Input> {
    pub(super) inp: I,
    depth: usize,
    /// The names of `GetStackTraceASN`: the root type, then a name per class-member and
    /// choice-variant frame (the member or variant being read) and `E` per array-element frame.
    trace: Vec<&'static str>,
    /// What NCBI posted while the object was read, in order.
    pub(super) diagnostics: Vec<MaskDiagnostic>,
}

impl<I: Input> Rd<I> {
    pub(super) fn new(inp: I) -> Self {
        Rd {
            inp,
            depth: 0,
            trace: Vec::new(),
            diagnostics: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-119
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
    //                 ...
    //                     const CMemberId& id = *frame.m_MemberId;
    //                     stack += '.';
    //                     if ( !id.GetName().empty() ) {
    //                         stack += id.GetName();
    //                     ...
    //         case TFrame::eFrameArrayElement:
    //             stack += ".E";
    //             break;
    //         default:
    //             break;
    // ```
    /// `GetStackTraceASN`.
    fn stack_trace(&self) -> String {
        self.trace.join(".")
    }

    /// `read` with `name` on the trace (a member, a variant or `E`).
    fn named<T>(&mut self, name: &'static str, read: impl FnOnce(&mut Self) -> R<T>) -> R<T> {
        self.trace.push(name);
        let result = read(self);
        self.trace.pop();
        result
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:189-199
    // ```c++
    // #define BEGIN_OBJECT_FRAME_OFx(Stream, Args) \
    //     (Stream).PushFrame Args; \
    //     try {
    //
    // #define END_OBJECT_FRAME_OF(Stream) \
    //     } catch (CSerialException& s_expt) { \
    //     ...
    //     } catch ( CEofException& e_expt ) { \
    //         (Stream).HandleEOF(e_expt); \
    // ```
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
    /// One `BEGIN_OBJECT_FRAME` ... `END_OBJECT_FRAME`.
    pub(super) fn frame<T>(&mut self, f: impl FnOnce(&mut Self) -> R<T>) -> R<T> {
        if self.depth >= MAX_DEPTH {
            return Err(Fail::Unsupported(format!(
                "objects nested more than {MAX_DEPTH} levels deep"
            )));
        }
        self.depth += 1;
        let result = f(self);
        self.depth -= 1;
        match result {
            Err(Fail::Eof) if self.depth >= 2 => format_err("unexpected end of file"),
            other => other,
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1398-1425
    // ```c++
    // void CObjectIStream::ReadClassSequential(const CClassTypeInfo* classType,
    //                                          TObjectPtr classPtr)
    // {
    //     TMemberIndex prevIndex = kInvalidMember;
    //     BEGIN_OBJECT_FRAME3(eFrameClass, classType, classPtr);
    //     BeginClass(classType);
    //
    //     ReadClassSequentialContentsBegin(classType);
    //
    //     TMemberIndex index;
    //     while ( (index = BeginClassMember(classType, *pos)) != kInvalidMember ) {
    //         ...
    //         ReadClassSequentialContentsMember(classPtr);
    //
    //         EndClassMember();
    //     }
    //
    //     ReadClassSequentialContentsEnd(classPtr);
    //
    //     EndClass();
    //     END_OBJECT_FRAME();
    // }
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
    /// A SEQUENCE read member by member, in declaration order (`ReadClassSequential`): a member
    /// that is not in the file is an error unless it is `OPTIONAL` or has a `DEFAULT`. `read`
    /// reads the value of member `index`.
    pub(super) fn read_class(
        &mut self,
        def: &'static ClassDef,
        mut read: impl FnMut(&mut Self, usize) -> R<()>,
    ) -> R<()> {
        self.frame(|s| {
            s.inp.begin_class()?;
            let mut pos = 0usize;
            s.frame(|s| {
                while let Some(index) = s.inp.begin_member(def, pos)? {
                    for m in &def.members[pos..index] {
                        if !m.optional {
                            return format_err(format!("member {} expected", m.name));
                        }
                    }
                    s.named(def.members[index].name, |s| read(s, index))?;
                    pos = index + 1;
                    s.inp.end_member()?;
                }
                Ok(())
            })?;
            for m in &def.members[pos..] {
                if !m.optional {
                    return format_err(format!("member {} expected", m.name));
                }
            }
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
    //     ...
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
        self.frame(|s| {
            s.inp.begin_choice()?;
            let value = s.frame(|s| {
                let index = s.inp.begin_variant(def)?;
                let value = s.named(def.variants[index].name, |s| read(s, index))?;
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
    /// A `SEQUENCE OF` or `SET OF`: `read` reads one element.
    pub(super) fn read_container(
        &mut self,
        set: bool,
        mut read: impl FnMut(&mut Self) -> R<()>,
    ) -> R<()> {
        self.frame(|s| {
            s.inp.begin_container(set)?;
            s.frame(|s| {
                while s.inp.begin_element()? {
                    s.named("E", &mut read)?;
                }
                Ok(())
            })?;
            s.inp.end_container()
        })
    }

    /// A `VisibleString`, with each byte outside `' '..='~'` stored as `#`; what NCBI posts for
    /// it goes to `diagnostics` with the trace of the frames.
    fn read_str(&mut self) -> R<Vec<u8>> {
        let mut notes = Vec::new();
        let result = self.inp.read_string(&mut notes);
        if !notes.is_empty() {
            let trace = self.stack_trace();
            for note in notes {
                self.diagnostics.push(match note {
                    StringNote::Utf8Tag => MaskDiagnostic::Utf8Tag(trace.clone()),
                    StringNote::Bad(bad) => MaskDiagnostic::BadChar(bad.message(&trace)),
                });
            }
        }
        result
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
                let v = self.inp.read_enum(def)?;
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
        match (db, tag) {
            (Some(db), Some(tag)) => Ok(Dbtag { db, tag }),
            _ => format_err("member expected"),
        }
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
        match id {
            Some(id) => Ok(GiimportId { id, db, release }),
            None => format_err("member id expected"),
        }
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
        match (seqid, cit) {
            (Some(seqid), Some(cit)) => Ok(PatentSeqId { seqid, cit }),
            _ => format_err("member expected"),
        }
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
        match (country, id) {
            (Some(country), Some(id)) => Ok(IdPat {
                country,
                id,
                doc_type,
            }),
            _ => format_err("member expected"),
        }
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
        match mol {
            Some(mol) => Ok(PdbSeqId {
                mol,
                chain,
                rel,
                chain_id,
            }),
            None => format_err("member mol expected"),
        }
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
            match year {
                Some(year) => {
                    std.year = year;
                    Ok(Date::Std(std))
                }
                None => format_err("member year expected"),
            }
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
                    let v = s.inp.read_enum(&NA_STRAND)?;
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
        let Some(id) = id else {
            return format_err("member id expected");
        };
        Ok((MaskInterval { from, to, strand }, id))
    }

    /// A `Seq-point` (also the members `a` and `b` of a `Seq-bond`), and its Seq-id.
    fn read_seq_point(&mut self) -> R<(MaskInterval, ParsedSeqId)> {
        let (mut point, mut strand, mut id) = (0u32, None, None);
        let mut scratch = Vec::new();
        self.read_class(&SEQ_POINT, |s, index| {
            match index {
                0 => point = s.inp.read_u32()?,
                1 => strand = Some(s.inp.read_enum(&NA_STRAND)? as u8),
                2 => id = Some(s.read_seq_id()?),
                _ => {
                    // fuzz: read, not used.
                    scratch.clear();
                    s.read_value(SEQ_POINT.members[index].ty, &mut scratch)?;
                }
            }
            Ok(())
        })?;
        let Some(id) = id else {
            return format_err("member id expected");
        };
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
    //     case e_Packed_int: ... (each interval) ...
    //     case e_Packed_pnt: ... x_UpdateId(id, &pnts.GetId(), may_throw);
    //     case e_Mix: ... (*li)->CheckId(id, may_throw) for each sub-location ...
    //     case e_Bond: ... A, then B when it is set ...
    //     case e_Equiv: ... (*li)->CheckId(id, may_throw) for each ...
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
    // NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_loc.hpp:939-945
    // ```c++
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
    /// and follows `CSeq_loc::GetId()` in `id` (`x_CheckId`: the first Seq-id, kept while the others
    /// are `Equals` to it; `Err` once two differ or a `feat` comes, which has none).
    fn read_seq_loc_parts(
        &mut self,
        parts: &mut Vec<MaskPart>,
        id: &mut Result<Option<ParsedSeqId>, ()>,
    ) -> R<()> {
        fn update(total: &mut Result<Option<ParsedSeqId>, ()>, id: ParsedSeqId) {
            match total {
                Ok(None) => *total = Ok(Some(id)),
                // `Equals` compares the values: the typed Seq-ids are equal exactly when it does.
                Ok(Some(first)) if *first != id => *total = Err(()),
                _ => {}
            }
        }
        self.read_choice(&SEQ_LOC, |s, index| {
            match index {
                SEQ_LOC_NULL => {
                    s.inp.read_null()?;
                    parts.push(MaskPart::Null);
                }
                SEQ_LOC_EMPTY => {
                    let seq_id = s.read_seq_id()?;
                    update(id, seq_id);
                    parts.push(MaskPart::Empty);
                }
                SEQ_LOC_WHOLE => {
                    let seq_id = s.read_seq_id()?;
                    update(id, seq_id);
                    parts.push(MaskPart::Interval(MaskInterval {
                        from: 0,
                        to: u32::MAX - 1,
                        strand: None,
                    }));
                }
                SEQ_LOC_INT => {
                    let (interval, seq_id) = s.read_seq_interval()?;
                    update(id, seq_id);
                    parts.push(MaskPart::Interval(interval));
                }
                SEQ_LOC_PACKED_INT => {
                    s.read_container(false, |s| {
                        let (interval, seq_id) = s.read_seq_interval()?;
                        update(id, seq_id);
                        parts.push(MaskPart::Interval(interval));
                        Ok(())
                    })?;
                }
                SEQ_LOC_PNT => {
                    let (interval, seq_id) = s.read_seq_point()?;
                    update(id, seq_id);
                    parts.push(MaskPart::Interval(interval));
                }
                SEQ_LOC_PACKED_PNT => {
                    let (mut strand, mut seq_id, mut points) = (None, None, Vec::new());
                    let mut scratch = Vec::new();
                    s.read_class(&PACKED_SEQPNT, |s, member| {
                        match member {
                            0 => strand = Some(s.inp.read_enum(&NA_STRAND)? as u8),
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
                    let Some(seq_id) = seq_id else {
                        return format_err("member id expected");
                    };
                    update(id, seq_id);
                    parts.extend(points.into_iter().map(|p| {
                        MaskPart::Interval(MaskInterval {
                            from: p,
                            to: p,
                            strand,
                        })
                    }));
                }
                SEQ_LOC_MIX | SEQ_LOC_EQUIV => {
                    s.read_container(index == SEQ_LOC_EQUIV, |s| s.read_seq_loc_parts(parts, id))?;
                }
                SEQ_LOC_BOND => {
                    s.read_class(&SEQ_BOND, |s, _| {
                        let (interval, seq_id) = s.read_seq_point()?;
                        update(id, seq_id);
                        parts.push(MaskPart::Interval(interval));
                        Ok(())
                    })?;
                }
                _ => {
                    // `feat`: read; `x_CheckId` fails, so the Seq-loc has no Seq-id.
                    let mut scratch = Vec::new();
                    s.read_value(SEQ_LOC.variants[index].ty, &mut scratch)?;
                    *id = Err(());
                }
            }
            Ok(())
        })
    }

    /// One `Seq-loc` of `masks`: the ranges `CSeq_loc_CI` visits and `CSeq_loc::GetId()`.
    fn read_seq_loc(&mut self) -> R<MaskLoc> {
        let mut parts = Vec::new();
        let mut id = Ok(None);
        self.read_seq_loc_parts(&mut parts, &mut id)?;
        Ok(MaskLoc {
            id: id.ok().flatten(),
            parts,
        })
    }

    /// A `Blast-mask-list`; the caller checks the end of the class frame.
    fn read_mask_list_members(&mut self) -> R<MaskList> {
        let mut masks = Vec::new();
        let mut more = None;
        self.read_class(&BLAST_MASK_LIST, |s, index| {
            if index == 0 {
                s.read_container(false, |s| {
                    masks.push(s.read_seq_loc()?);
                    Ok(())
                })
            } else {
                more = Some(s.inp.read_bool()?);
                Ok(())
            }
        })?;
        let Some(more) = more else {
            return format_err("member more expected");
        };
        Ok(MaskList { masks, more })
    }

    /// `file >> MSerial_Asn... >> *obj` for a `Blast-db-mask-info`.
    fn read_info(&mut self) -> R<MaskInfo> {
        let mut algo_id = 0;
        let mut algo_program = 0;
        let mut algo_options = Vec::new();
        let mut masks = None;
        self.read_top(BLAST_DB_MASK_INFO.name, |s| {
            s.read_class(&BLAST_DB_MASK_INFO, |s, index| {
                match index {
                    0 => algo_id = s.inp.read_i32()?,
                    1 => algo_program = s.inp.read_enum(&BLAST_FILTER_PROGRAM)?,
                    2 => algo_options = s.read_str()?,
                    _ => masks = Some(s.read_mask_list_members()?),
                }
                Ok(())
            })
        })?;
        let Some(masks) = masks else {
            return format_err("member masks expected");
        };
        Ok(MaskInfo {
            algo_id,
            algo_program,
            algo_options,
            masks,
        })
    }

    /// `file >> MSerial_Asn... >> *obj` for a `Blast-mask-list` (the continuation).
    fn read_list(&mut self) -> R<MaskList> {
        let mut list = None;
        self.read_top(BLAST_MASK_LIST.name, |s| {
            list = Some(s.read_mask_list_members()?);
            Ok(())
        })?;
        list.ok_or_else(|| Fail::Format("no Blast-mask-list".into()))
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
    fn read_top(
        &mut self,
        type_name: &'static str,
        read: impl FnOnce(&mut Self) -> R<()>,
    ) -> R<()> {
        self.depth = 0;
        self.trace.clear();
        self.trace.push(type_name);
        if let Some(name) = self.frame(|s| s.inp.read_file_header())? {
            if !name.is_empty() && name != type_name.as_bytes() {
                return format_err(format!(
                    "incompatible type {}<>{type_name}",
                    String::from_utf8_lossy(&name)
                ));
            }
        }
        self.frame(read)
    }
}

/// One `Seq-id` of an ASN.1 binary value, read by the mask reader's builder, for the test that
/// compares it with `seqid::ber::decode_seq_id` (DECISIONS 35).
#[cfg(test)]
pub(super) fn seq_id_from_binary(data: &[u8]) -> R<ParsedSeqId> {
    let mut rd = Rd::new(BerIn::new(data));
    rd.frame(|s| s.read_seq_id())
}

/// [`seq_id_from_binary`] for ASN.1 text, for the test that compares the builder with
/// `seqid::asn_text`.
#[cfg(test)]
pub(super) fn seq_id_from_text(data: &[u8]) -> R<ParsedSeqId> {
    use crate::report::serial::asn_text_reader::SliceSource;
    let mut rd = Rd::new(TextIn::new(SliceSource::new(data)));
    rd.frame(|s| s.read_seq_id())
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
    /// What NCBI posted while it read the objects so far (`take_diagnostics`).
    diagnostics: Vec<MaskDiagnostic>,
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
    /// Opens the file and guesses its format from its first bytes. The text of the error for a
    /// file in neither format is NCBI's. A file that cannot be opened or read is a failed stream
    /// for the format guess (`x_TestInput`), which answers `eUnknown`: the same error. A text
    /// file is then read in blocks as its objects are read; a binary file is read whole (a read
    /// error there is the same error).
    #[cfg(test)]
    pub(crate) fn open(path: &Path) -> Result<Self, MaskFileError> {
        Self::open_named(path, path)
    }

    /// [`open`](Self::open) for a file that messages call `name` (the name on makeblastdb's
    /// command line).
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
            diagnostics: Vec::new(),
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
                self.diagnostics.append(&mut rd.diagnostics);
                result
            }
            Source::Binary { data, offset } => {
                let start = (*offset).min(data.len());
                let mut rd = Rd::new(BerIn::new(&data[start..]));
                let result = binary(&mut rd);
                *offset = start + rd.inp.offset();
                self.diagnostics.append(&mut rd.diagnostics);
                result
            }
        }
    }

    /// What NCBI's reader posted (standard error) while it read the objects returned since the
    /// last call, and the one that failed: `next_info` and `next_list` read on after them.
    pub(crate) fn take_diagnostics(&mut self) -> Vec<MaskDiagnostic> {
        std::mem::take(&mut self.diagnostics)
    }

    fn stopped(&self) -> MaskFileError {
        if self.state == State::FailedEof {
            return MaskFileError::EofInContinuation {
                path: self.path.clone(),
            };
        }
        MaskFileError::Malformed {
            path: self.path.clone(),
            message: "reader stopped after an earlier error".into(),
        }
    }

    fn fail(&mut self, fail: Fail) -> MaskFileError {
        self.state = if matches!(fail, Fail::Eof) {
            State::FailedEof
        } else {
            State::Failed
        };
        match fail {
            Fail::Format(message) => MaskFileError::Malformed {
                path: self.path.clone(),
                message,
            },
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
