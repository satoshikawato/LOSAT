//! The frames of NCBI's object stack (`CObjectStackFrame`) as its readers report them: the one
//! owner of a frame's texts (`GetFrameInfo`, which `NCBI_RETHROW_SAME` adds to a `CException`
//! that leaves the frame, and `GetFrameName`, which `CSerialException::AddFrameInfo` adds to the
//! exception's frame stack) and of the places of the `END_OBJECT_FRAME` of NCBI's ASN.1 text
//! reader. The PSSM reader (`report::serial::asn_types::read`) and the reader of makeblastdb's
//! mask files (`writedb::mask_info`) keep their stacks of these frames.

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
// ```c++
// string CObjectStackFrame::GetFrameName(void) const
// {
//     string info;
//     switch ( GetFrameType() ) {
//     case eFrameClassMember:
//     case eFrameChoiceVariant:
//         {
//             if ( m_MemberId ) {
//                 const CMemberId& id = *m_MemberId;
// //                info = '.';
//                 if ( !id.GetName().empty() ) {
//                     info += id.GetName();
//                 }
//                 else {
//                     info += '[';
//                     info += NStr::IntToString(id.GetTag());
//                     info += ']';
//                 }
//             }
//         }
//         break;
//     case eFrameArrayElement:
//         info = "[]";
//         break;
//     case eFrameArray:
//         info = "[]";
//         break;
//     case eFrameNamed:
//         if (!GetNotag()) {
//             info = GetTypeInfo()->GetName();
//         }
//         break;
//     default:
//         {
//             break;
//         }
//     }
//     return info;
// }
// ```
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

/// The source file, line and function of the `END_OBJECT_FRAME` of NCBI's ASN.1 text reader
/// (`CObjectIStreamAsn`) that leaves a frame of `kind` (`skip`: of a skip function): the place of
/// the line that `NCBI_RETHROW_SAME` adds to a `CException` other than a `CSerialException`
/// passing through it.
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
