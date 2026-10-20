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

use super::guess::guess_format;
use super::schema::{
    ChoiceDef, ClassDef, EnumDef, Ty, BLAST_DB_MASK_INFO, BLAST_FILTER_PROGRAM, BLAST_MASK_LIST,
    NA_STRAND, OBJECT_ID_ID, SEQ_ID, SEQ_ID_LOCAL, SEQ_INTERVAL, SEQ_INTERVAL_FROM,
    SEQ_INTERVAL_ID, SEQ_INTERVAL_STRAND, SEQ_INTERVAL_TO, SEQ_LOC, SEQ_LOC_INT,
    SEQ_LOC_PACKED_INT, T_OBJECT_ID,
};
use super::{MaskInfo, MaskInterval, MaskList, MaskLoc, MaskSeqId};
use crate::writedb::WriteDbError;

/// The encoding of a mask file, as `CFormatGuess` finds it. Anything else is an error when the
/// file is opened ("Unknown encoding").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskFileFormat {
    /// `eBinaryASN`: `MSerial_AsnBinary`.
    BinaryAsn,
    /// `eTextASN`: `MSerial_AsnText`.
    TextAsn,
}

/// Why a mask file could not be read.
#[derive(Debug)]
pub enum MaskFileError {
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
    /// A `VisibleString` as it is in the file (line breaks of a text string gone), and whether a
    /// binary file gave it the tag of a UTF8String.
    fn read_string(&mut self) -> R<(Vec<u8>, bool)>;
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32>;
    /// Bytes read so far (for the "file is not at its end" decisions of the caller).
    fn offset(&self) -> usize;
}

// NCBI reference (598d8ae6): c++/include/serial/objistr.hpp:1160-1163
// ```c++
// inline
// bool GoodVisibleChar(char c);
//
// #define SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING 0
// ```
// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:455-459
// ```c++
// inline
// bool GoodVisibleChar(char c)
// {
//     return c >= ' ' && c <= '~';
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:307-311
// ```c++
// NCBI_PARAM_ENUM_DECL(EFixNonPrint, SERIAL, WRONG_CHARS_READ);
// NCBI_PARAM_ENUM_DEF(EFixNonPrint, SERIAL, WRONG_CHARS_READ, eFNP_ReplaceAndWarn);
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1898
// ```c++
// char ReplaceVisibleChar(char c, EFixNonPrint fix_method,
//     const CObjectStack* io, const CTempString& str, char subst)
// {
//     ...
//     if ( fix_method != eFNP_Replace ) {
//         ...
//         case eFNP_ReplaceAndWarn:
//             ERR_POST_X(7, message);
//             break;
//     ...
//     return subst;
// }
// ```
/// A `VisibleString` is read with the default `eFNP_ReplaceAndWarn`: every byte outside
/// `' '..='~'` becomes `m_NonPrintSubst` (`'#'`, objistr.cpp:463), and NCBI posts a diagnostic
/// per byte (`Error: (803.7) [makeblastdb] Blast-db-mask-info.algo-options` / `Bad char [0x9]
/// in string at line 1` / the chunk); a UTF8String tag for a VisibleString posts a warning once
/// (`Warning: (803.10) ... UTF8String data for VisibleString member`). Both leave the status
/// at 0. LOSAT does not write those lines, so the reader refuses such a string
/// ([`Rd::read_str`]); this function is the value NCBI stores, which the refusal does not need
/// but the tests compare.
pub(super) fn fix_visible_chars(s: &mut [u8]) {
    for c in s {
        if !(b' '..=b'~').contains(c) {
            *c = b'#';
        }
    }
}

/// The deepest nesting of object frames that LOSAT follows (NCBI's stack is only limited by memory).
const MAX_DEPTH: usize = 400;

/// What `Seq-loc` reading found, besides the value.
struct Collected {
    /// The first variant of a `Seq-loc` that LOSAT does not model; it is raised after the whole
    /// object has been read, so that an error NCBI raises later in the object comes first.
    unsupported: Option<String>,
}

/// The state of one `CObjectIStream`: the primitives, and the number of object frames.
pub(super) struct Rd<I: Input> {
    pub(super) inp: I,
    depth: usize,
    collected: Collected,
}

impl<I: Input> Rd<I> {
    pub(super) fn new(inp: I) -> Self {
        Rd {
            inp,
            depth: 0,
            collected: Collected { unsupported: None },
        }
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
                    read(s, index)?;
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
                    read(s)?;
                }
                Ok(())
            })?;
            s.inp.end_container()
        })
    }

    /// A `VisibleString` member `name`. A byte outside `' '..='~'`, or the tag of a UTF8String
    /// in a binary file, is read (NCBI stores `'#'` for the byte), but the object is then refused
    /// as [`Fail::Unsupported`] once it has been read to its end, so that an error NCBI raises
    /// later in the object comes first: NCBI's diagnostics for it are not reproduced.
    fn read_str(&mut self, name: &'static str) -> R<Vec<u8>> {
        let (mut s, utf8_tag) = self.inp.read_string()?;
        if utf8_tag {
            self.collected.unsupported.get_or_insert_with(|| {
                format!("a UTF8String tag for the VisibleString member `{name}` (NCBI warns and reads it)")
            });
        }
        if let Some(&byte) = s.iter().find(|c| !(b' '..=b'~').contains(c)) {
            self.collected.unsupported.get_or_insert_with(|| {
                format!(
                    "byte 0x{byte:02x} in the VisibleString member `{name}` (NCBI replaces it with '#' and logs it)"
                )
            });
        }
        fix_visible_chars(&mut s);
        Ok(s)
    }

    /// Reads a value of type `ty` and appends a canonical encoding of it to `out`: two values
    /// give the same bytes exactly when they are equal as `CSerialObject::Equals` compares
    /// them, apart from members that are left out in one and stated with their `DEFAULT` in the
    /// other (those differ here). The encoding is self-delimiting (every item states its kind,
    /// a string its length, a list its count).
    pub(super) fn read_value(
        &mut self,
        ty: &'static Ty,
        name: &'static str,
        out: &mut Vec<u8>,
    ) -> R<()> {
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
                let v = self.read_str(name)?;
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
                    s.read_value(def.members[index].ty, def.members[index].name, out)
                })?;
                for _ in next..def.members.len() {
                    out.push(0);
                }
            }
            Ty::Choice(def) => {
                out.push(b'C');
                self.read_choice(def, |s, index| {
                    out.push(index as u8);
                    s.read_value(def.variants[index].ty, def.variants[index].name, out)
                })?;
            }
            Ty::SeqOf(element) | Ty::SetOf(element) => {
                out.push(b'L');
                let count_at = out.len();
                out.extend_from_slice(&[0u8; 8]);
                let mut count = 0u64;
                self.read_container(matches!(ty, Ty::SetOf(_)), |s| {
                    count += 1;
                    s.read_value(element, name, out)
                })?;
                out[count_at..count_at + 8].copy_from_slice(&count.to_le_bytes());
            }
        }
        Ok(())
    }

    /// A `Seq-id`: `local id N` is a [`MaskSeqId::LocalId`], any other value is read in full.
    fn read_seq_id(&mut self) -> R<MaskSeqId> {
        let mut canon = Vec::new();
        let mut local = None;
        self.read_choice(&SEQ_ID, |s, index| {
            canon.push(index as u8);
            if index == SEQ_ID_LOCAL {
                let Ty::Choice(object_id) = &T_OBJECT_ID else {
                    unreachable!("Object-id is a CHOICE")
                };
                s.read_choice(object_id, |s, v| {
                    if v == OBJECT_ID_ID {
                        local = Some(s.inp.read_i32()?);
                        Ok(())
                    } else {
                        canon.push(v as u8);
                        s.read_value(
                            object_id.variants[v].ty,
                            object_id.variants[v].name,
                            &mut canon,
                        )
                    }
                })
            } else {
                s.read_value(
                    SEQ_ID.variants[index].ty,
                    SEQ_ID.variants[index].name,
                    &mut canon,
                )
            }
        })?;
        Ok(match local {
            Some(n) => MaskSeqId::LocalId(n),
            None => MaskSeqId::Other(canon),
        })
    }

    /// A `Seq-interval`, and its Seq-id.
    fn read_seq_interval(&mut self) -> R<(MaskInterval, MaskSeqId)> {
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
                    s.read_value(
                        SEQ_INTERVAL.members[index].ty,
                        SEQ_INTERVAL.members[index].name,
                        &mut scratch,
                    )?;
                }
            }
            Ok(())
        })?;
        let Some(id) = id else {
            return format_err("member id expected");
        };
        Ok((MaskInterval { from, to, strand }, id))
    }

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
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3541-3590
    // ```c++
    // bool CSeq_loc::x_CheckId(const CSeq_id*& id, bool may_throw) const
    // {
    //     switch ( Which() ) {
    //     ...
    //     case e_Int:
    //         {
    //             const CSeq_interval& loc = GetInt();
    //             return x_UpdateId(id, &loc.GetId(), may_throw);
    //         }
    //     ...
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
    /// `CSeq_loc::GetId()` over the Seq-ids of the intervals of one Seq-loc.
    fn common_id(ids: Vec<MaskSeqId>) -> R<Option<MaskSeqId>> {
        let mut ids = ids.into_iter();
        let Some(first) = ids.next() else {
            // An empty packed-int has no Seq-id.
            return Ok(None);
        };
        for id in ids {
            if id != first {
                // `Equals` compares the values: two different `local id` or a `local id` and
                // another kind of Seq-id are different; two others may still be equal for
                // NCBI (a member at its DEFAULT), which LOSAT does not model.
                return match (&first, &id) {
                    (MaskSeqId::Other(_), MaskSeqId::Other(_)) => Err(Fail::Unsupported(
                        "a Seq-loc with two different Seq-ids that are not both `local id`".into(),
                    )),
                    _ => Ok(None),
                };
            }
        }
        Ok(Some(first))
    }

    /// One `Seq-loc` of `masks`: `int` and `packed-int` are read into a [`MaskLoc`]; every other
    /// variant is read to check it, then remembered as not supported.
    fn read_seq_loc(&mut self) -> R<Option<MaskLoc>> {
        let mut loc = None;
        self.read_choice(&SEQ_LOC, |s, index| {
            match index {
                SEQ_LOC_INT => {
                    let (interval, id) = s.read_seq_interval()?;
                    loc = Some(MaskLoc {
                        id: Self::common_id(vec![id])?,
                        intervals: vec![interval],
                    });
                }
                SEQ_LOC_PACKED_INT => {
                    let mut intervals = Vec::new();
                    let mut ids = Vec::new();
                    s.read_container(false, |s| {
                        let (interval, id) = s.read_seq_interval()?;
                        intervals.push(interval);
                        ids.push(id);
                        Ok(())
                    })?;
                    loc = Some(MaskLoc {
                        id: Self::common_id(ids)?,
                        intervals,
                    });
                }
                _ => {
                    let mut scratch = Vec::new();
                    s.read_value(
                        SEQ_LOC.variants[index].ty,
                        SEQ_LOC.variants[index].name,
                        &mut scratch,
                    )?;
                    s.collected.unsupported.get_or_insert_with(|| {
                        format!(
                            "Seq-loc variant `{}` in mask data",
                            SEQ_LOC.variants[index].name
                        )
                    });
                }
            }
            Ok(())
        })?;
        Ok(loc)
    }

    /// A `Blast-mask-list`; the caller checks the end of the class frame.
    fn read_mask_list_members(&mut self) -> R<MaskList> {
        let mut masks = Vec::new();
        let mut more = None;
        self.read_class(&BLAST_MASK_LIST, |s, index| {
            if index == 0 {
                s.read_container(false, |s| {
                    if let Some(loc) = s.read_seq_loc()? {
                        masks.push(loc);
                    }
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
                    2 => algo_options = s.read_str("algo-options")?,
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
    fn read_top(&mut self, type_name: &str, read: impl FnOnce(&mut Self) -> R<()>) -> R<()> {
        self.depth = 0;
        self.collected.unsupported = None;
        if let Some(name) = self.frame(|s| s.inp.read_file_header())? {
            if !name.is_empty() && name != type_name.as_bytes() {
                return format_err(format!(
                    "incompatible type {}<>{type_name}",
                    String::from_utf8_lossy(&name)
                ));
            }
        }
        self.frame(read)?;
        match self.collected.unsupported.take() {
            Some(what) => Err(Fail::Unsupported(what)),
            None => Ok(()),
        }
    }
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

/// A mask file read from the disk.
pub struct MaskFileReader {
    path: PathBuf,
    format: MaskFileFormat,
    data: Vec<u8>,
    /// The offset of the next object; `data.len()` once the file is over.
    offset: usize,
    /// Set by the first error or end of file (NCBI stops reading the file there).
    state: State,
}

impl fmt::Debug for MaskFileReader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MaskFileReader")
            .field("path", &self.path)
            .field("format", &self.format)
            .field("bytes", &self.data.len())
            .field("offset", &self.offset)
            .field("state", &self.state)
            .finish()
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
    /// Reads the file and guesses its format. The text of the error for a file in neither
    /// format is NCBI's. A file that cannot be opened or read is a failed stream for the format
    /// guess (`x_TestInput`), which answers `eUnknown`: the same error.
    pub fn open(path: &Path) -> Result<Self, MaskFileError> {
        Self::open_named(path, path)
    }

    /// [`open`](Self::open) for a file that messages call `name` (the name on makeblastdb's
    /// command line).
    pub fn open_named(path: &Path, name: &Path) -> Result<Self, MaskFileError> {
        let mut data = Vec::new();
        if std::fs::File::open(path)
            .and_then(|mut f| f.read_to_end(&mut data))
            .is_err()
        {
            // A read that fails half way leaves the bytes read; NCBI's sample would have them.
            // An empty or unreadable file has no sample at all.
            return Err(unknown_encoding(name));
        }
        Self::from_bytes(name, data)
    }

    /// The reader of a file whose bytes are `data`; `path` is only used in messages.
    pub(crate) fn from_bytes(path: &Path, data: Vec<u8>) -> Result<Self, MaskFileError> {
        let Some(format) = guess_format(&data) else {
            return Err(unknown_encoding(path));
        };
        Ok(MaskFileReader {
            path: path.to_path_buf(),
            format,
            data,
            offset: 0,
            state: State::Open,
        })
    }

    pub fn format(&self) -> MaskFileFormat {
        self.format
    }

    fn run<T>(
        &mut self,
        text: impl FnOnce(&mut Rd<super::text::TextIn<'_>>) -> R<T>,
        binary: impl FnOnce(&mut Rd<super::ber::BerIn<'_>>) -> R<T>,
    ) -> R<T> {
        let start = self.offset.min(self.data.len());
        let data = &self.data[start..];
        let (result, used) = match self.format {
            MaskFileFormat::TextAsn => {
                let mut rd = Rd::new(super::text::TextIn::new(data));
                let result = text(&mut rd);
                (result, rd.inp.offset())
            }
            MaskFileFormat::BinaryAsn => {
                let mut rd = Rd::new(super::ber::BerIn::new(data));
                let result = binary(&mut rd);
                (result, rd.inp.offset())
            }
        };
        self.offset = start + used;
        result
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
    pub fn next_info(&mut self) -> Result<Option<MaskInfo>, MaskFileError> {
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
    pub fn next_list(&mut self) -> Result<MaskList, MaskFileError> {
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
