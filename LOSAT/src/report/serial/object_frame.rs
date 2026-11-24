//! The frames of NCBI's object stack (`CObjectStackFrame`) as its readers report them: the one
//! owner of a frame's texts (`GetFrameInfo`, which `NCBI_RETHROW_SAME` adds to a `CException`
//! that leaves the frame, and `GetFrameName`, which `CSerialException::AddFrameInfo` adds to the
//! exception's frame stack) and of the places of the `END_OBJECT_FRAME` of NCBI's ASN.1 text
//! reader, and the throws of the object stack itself (with the getter's throw of an unset member,
//! which with `ExpectedMember`'s is a check of the data that NCBI's settings turn off:
//! `DataCheck`). The table reader
//! (`report::serial::asn_types::read`), which reads makeblastdb's mask files too
//! (`writedb::mask_info`), keeps its stack of these frames.

use crate::ncbi_exception::{NcbiException, ThrowSite};

/// A frame of NCBI's object stack (`CObjectStackFrame::EFrameType`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FrameKind {
    Named,
    Class,
    ClassMember,
    Choice,
    ChoiceVariant,
    Array,
    ArrayElement,
}

// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:366-381
// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:399-410
/// `GetFrameInfo()` of a frame of `kind` whose type (the named, class, choice and array frames)
/// is `type_name` (empty for an anonymous type) and whose member or variant (a member or variant
/// frame, once `SetTopMemberId` named it) is `member`.
pub(crate) fn frame_info(kind: FrameKind, type_name: Option<&str>, member: Option<&str>) -> String {
    let kind = match kind {
        FrameKind::Named => "eFrameNamed",
        FrameKind::Class => "eFrameClass",
        FrameKind::ClassMember => "eFrameClassMember",
        FrameKind::Choice => "eFrameChoice",
        FrameKind::ChoiceVariant => "eFrameChoiceVariant",
        FrameKind::Array => "eFrameArray",
        FrameKind::ArrayElement => "eFrameArrayElement",
    };
    let mut info = format!(" Frame type= {kind}");
    if let Some(name) = type_name {
        info.push_str(", Object type= ");
        info.push_str(name);
    }
    if let Some(member) = member {
        info.push_str(", Member name= ");
        info.push_str(member);
    }
    info
}

// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:412-452
/// `GetFrameName()`: what `CSerialException::AddFrameInfo` adds for the frame (the member's name,
/// `[]` for a container and its element, the type of a named frame, nothing for a class or a
/// choice). Every member of LOSAT's types has a name.
pub(crate) fn frame_name<'a>(
    kind: FrameKind,
    type_name: Option<&'a str>,
    member: Option<&'a str>,
) -> &'a str {
    match kind {
        FrameKind::ClassMember | FrameKind::ChoiceVariant => member.unwrap_or(""),
        FrameKind::Array | FrameKind::ArrayElement => "[]",
        FrameKind::Named => type_name.unwrap_or(""),
        FrameKind::Class | FrameKind::Choice => "",
    }
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

// The throws of the object stack itself, which both encodings share.
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608
pub(crate) const HANDLE_EOF: ThrowSite = istr_site(
    606,
    "CSerialException::eEOF",
    "ncbi::CObjectIStream::HandleEOF",
);
pub(crate) const SKIP_FILE_HEADER: ThrowSite = istr_site(
    896,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::SkipFileHeader",
);
/// `ExpectedMember`'s throw of a mandatory member that the input skips (the read check,
/// `DataCheck::Read`): every reader's (the table reader's, `asn_types::read`; the header records',
/// `seqid::ber`).
pub(crate) const EXPECTED_MEMBER: ThrowSite = istr_site(
    1037,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::ExpectedMember",
);

// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:300-307
// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:317-320
// NCBI reference (598d8ae6): c++/src/serial/datatool/classstr.cpp:849-855
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:805-810
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:836-842
// The generated getters call `ThrowUnassigned` with the index, so the place is the `NCBI_THROW` of
// serialobject.cpp:306 (`file_name` and `file_line` are 0). Its function is a `const` method,
// whose `__PRETTY_FUNCTION__` ends with `const`, not `)`: `ParseCurrFunctName` names no
// function then.
/// `ThrowUnassigned`, the throw of the getter of a member that is not set (the get check,
/// `DataCheck::Get`): every getter's (the PSSM query's feature table, `annot_index::table`;
/// rpsbproc's).
pub(crate) const UNASSIGNED_MEMBER: ThrowSite = ThrowSite {
    file: "src/serial/serialobject.cpp",
    line: 306,
    kind: "CUnassignedMember::eGet",
    function: None,
};

/// A check of the data of NCBI's serial library that a setting of the process turns off
/// (`blastinput::ncbi_environment::ReadVerification`): the read of a mandatory member
/// (`SERIAL_VERIFY_DATA_READ`) and the get of an unset member (`SERIAL_VERIFY_DATA_GET`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DataCheck {
    Read,
    Get,
}

impl DataCheck {
    /// The check that threw `exception` first (the oldest of its chain, which the frames it left
    /// rethrew: `NcbiException::first_site`), `None` for an exception of no check. The one test
    /// of which errors depend on the settings.
    pub(crate) fn of(exception: &NcbiException) -> Option<DataCheck> {
        match exception.first_site() {
            site if site == EXPECTED_MEMBER => Some(DataCheck::Read),
            site if site == UNASSIGNED_MEMBER => Some(DataCheck::Get),
            _ => None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1206-1219
/// `ReadPointer`'s check of a value of a named class that is not the declared type (a class of
/// the module has no parent class: the second throw).
pub(crate) const READ_POINTER_INCOMPATIBLE: ThrowSite = istr_site(
    1217,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::ReadPointer",
);
// NCBI reference (598d8ae6): c++/src/serial/classinfob.cpp:244-252
/// `GetClassInfoByName` of a name that the registry of classes does not hold.
pub(crate) const CLASS_NOT_FOUND: ThrowSite = ThrowSite {
    file: "src/serial/classinfob.cpp",
    line: 251,
    kind: "CSerialException::eInvalidData",
    function: Some("ncbi::CClassTypeInfoBase::GetClassInfoByName"),
};
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:873-880
// (makeblastdb's streams collect no objects: `m_Objects` is null, checked with makeblastdb
// 2.17.0 on `@0` and `@1`.)
pub(crate) const GET_REGISTERED_OBJECT: ThrowSite = istr_site(
    877,
    "CSerialException::eFormatError",
    "ncbi::CObjectIStream::GetRegisteredObject",
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
/// `ReadInt2`: an `Int4` that does not fit.
pub(crate) const READ_INT2: ThrowSite = istr_site(
    1702,
    "CSerialException::eOverflow",
    "ncbi::CObjectIStream::ReadInt2",
);

/// The source file, line and function of the `END_OBJECT_FRAME` of NCBI's ASN.1 text reader
/// (`CObjectIStreamAsn`) that leaves a frame of `kind` (`skip`: of a skip function): the place of
/// the line that `NCBI_RETHROW_SAME` adds to a `CException` other than a `CSerialException`
/// passing through it.
// NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
// The lines are those of `END_OBJECT_FRAME()` in the functions that push the frames:
// `CObjectIStream::Read` (objistr.cpp:959-969), `CObjectIStreamAsn::ReadClassSequential`
// (objistrasn.cpp:1380-1397, the member frame ends in `ReadClassSequentialContentsEnd` at
// 1393), `CObjectIStream::ReadChoiceSimple` (objistr.cpp:1482-1501),
// `CObjectIStreamAsn::ReadContainer` (objistrasn.cpp:1245-1274), and the skip functions
// `SkipClassSequential` (objistrasn.cpp:1417-1433), `SkipChoiceSimple` (objistr.cpp:1503-1522)
// and `SkipContainer` (objistrasn.cpp:1276-1292).
pub(crate) fn text_frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
    const ASN: &str = "src/serial/objistrasn.cpp";
    const ISTR: &str = "src/serial/objistr.cpp";
    match (kind, skip) {
        (FrameKind::Named, _) => (ISTR, 968, "ncbi::CObjectIStream::Read"),
        (FrameKind::Class, false) => (ASN, 1396, "ncbi::CObjectIStreamAsn::ReadClassSequential"),
        (FrameKind::ClassMember, false) => {
            (ASN, 1393, "ncbi::CObjectIStreamAsn::ReadClassSequential")
        }
        (FrameKind::Choice, false) => (ISTR, 1500, "ncbi::CObjectIStream::ReadChoiceSimple"),
        (FrameKind::ChoiceVariant, false) => (ISTR, 1498, "ncbi::CObjectIStream::ReadChoiceSimple"),
        (FrameKind::Array, false) => (ASN, 1273, "ncbi::CObjectIStreamAsn::ReadContainer"),
        (FrameKind::ArrayElement, false) => (ASN, 1270, "ncbi::CObjectIStreamAsn::ReadContainer"),
        (FrameKind::Class, true) => (ASN, 1432, "ncbi::CObjectIStreamAsn::SkipClassSequential"),
        (FrameKind::ClassMember, true) => {
            (ASN, 1429, "ncbi::CObjectIStreamAsn::SkipClassSequential")
        }
        (FrameKind::Choice, true) => (ISTR, 1521, "ncbi::CObjectIStream::SkipChoiceSimple"),
        (FrameKind::ChoiceVariant, true) => (ISTR, 1519, "ncbi::CObjectIStream::SkipChoiceSimple"),
        (FrameKind::Array, true) => (ASN, 1291, "ncbi::CObjectIStreamAsn::SkipContainer"),
        (FrameKind::ArrayElement, true) => (ASN, 1288, "ncbi::CObjectIStreamAsn::SkipContainer"),
    }
}

/// The `END_OBJECT_FRAME` of the named frame that `CObjectIStream::ReadPointer`, which the text
/// stream does not override, puts around a value of the type that `:Type` names (the frame
/// `GetFrameName` calls by that type): the place of the line that `NCBI_RETHROW_SAME` adds to a
/// `CException` other than a `CSerialException` passing through it.
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1173-1198
pub(crate) const READ_POINTER_FRAME_END: (&str, u32, &str) = (
    "src/serial/objistr.cpp",
    1195,
    "ncbi::CObjectIStream::ReadPointer",
);

/// [`READ_POINTER_FRAME_END`] of `CObjectIStream::SkipPointer`, the skip of a `CRef<>`.
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1256-1270
pub(crate) const SKIP_POINTER_FRAME_END: (&str, u32, &str) = (
    "src/serial/objistr.cpp",
    1267,
    "ncbi::CObjectIStream::SkipPointer",
);
