//! NCBI's object stack (`CObjectStack`) as far as its messages use it: the frames of the value
//! that a stream reads or writes, whose names make `GetStackTraceASN`. One `ObjectStack` is held
//! by each stream that LOSAT reads or writes as NCBI's object streams do: the shared ASN.1 text
//! reader (`asn_text_reader::Reader`: PSSM, Seq-id, rpsbproc archive and mask text), the binary
//! mask reader (`writedb::mask_info::ber`), the writer of definition lines (`writedb::header`)
//! and the ASN.1 text and binary encoders of every serialized document (`asn_text`, `asn_ber`:
//! PSSM, Seq-annot, archive, search strategy, mask files); the trace goes into every note of a
//! replaced byte (`visible::BadChars`).
// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:80-119
// ```c++
// string CObjectStack::GetStackTraceASN(void) const
// {
//     if (!GetStackDepth()) {
//         return "stack is empty";
//     }
//     _ASSERT(FetchFrameFromBottom(0).m_FrameType == TFrame::eFrameNamed);
//     _ASSERT(FetchFrameFromBottom(0).m_TypeInfo);
//     string stack = FetchFrameFromBottom(0).HasTypeInfo() ?
//         FetchFrameFromBottom(0).m_TypeInfo->GetName() : "?";
//     for ( size_t i = 1; i < GetStackDepth(); ++i ) {
//         const TFrame& frame = FetchFrameFromBottom(i);
//         switch ( frame.m_FrameType ) {
//         case TFrame::eFrameClassMember:
//         case TFrame::eFrameChoiceVariant:
//             {
//                 if ( !frame.m_MemberId ) {
//                     _ASSERT(frame.m_FrameType == TFrame::eFrameClassMember);
//                 }
//                 else {
//                     const CMemberId& id = *frame.m_MemberId;
//                     stack += '.';
//                     if ( !id.GetName().empty() ) {
//                         stack += id.GetName();
//                     }
//                     else {
//                         stack += '[';
//                         stack += NStr::IntToString(id.GetTag());
//                         stack += ']';
//                     }
//                 }
//             }
//             break;
//         case TFrame::eFrameArrayElement:
//             stack += ".E";
//             break;
//         default:
//             break;
//         }
//     }
//     return stack;
// }
// ```
// The frames that add nothing (`default`: the class, choice, array frames themselves) are not
// kept; a class-member frame whose member is not known yet (`!frame.m_MemberId`) is the empty
// name.

/// The names of `GetStackTraceASN`: the root type, then the name of each member and variant
/// being read or written (empty for a member frame whose member is not known yet) and `E` for
/// each element.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct ObjectStack {
    frames: Vec<&'static str>,
}

impl ObjectStack {
    /// The stack of a stream that starts the object of type `type_name` (the bottom frame).
    pub(crate) fn with_root(type_name: &'static str) -> Self {
        ObjectStack {
            frames: vec![type_name],
        }
    }

    /// The object of type `type_name` starts: the stack holds its frame only.
    pub(crate) fn set_root(&mut self, type_name: &'static str) {
        self.frames.clear();
        self.frames.push(type_name);
    }

    /// A member, variant or element (`E`) frame; `""` for a member frame whose member is not
    /// known yet (`rename` names it).
    pub(crate) fn push(&mut self, name: &'static str) {
        self.frames.push(name);
    }

    pub(crate) fn pop(&mut self) {
        self.frames.pop();
    }

    /// The number of frames.
    pub(crate) fn depth(&self) -> usize {
        self.frames.len()
    }

    /// The frames above `depth` dropped.
    pub(crate) fn truncate(&mut self, depth: usize) {
        self.frames.truncate(depth);
    }

    /// The frame on top names `name`: a writer's container frame becomes the element frame `E`
    /// at its first element (`BeginContainerElement`; no string is written between two elements,
    /// so one frame serves them all).
    pub(crate) fn set_top(&mut self, name: &'static str) {
        if let Some(top) = self.frames.last_mut() {
            *top = name;
        }
    }

    /// `SetTopMemberId`: the frame at `depth` names the member `name`.
    pub(crate) fn rename(&mut self, depth: usize, name: &'static str) {
        self.frames[depth] = name;
    }

    /// `GetStackTraceASN`.
    pub(crate) fn trace(&self) -> String {
        let mut trace = String::new();
        for name in self.frames.iter().filter(|name| !name.is_empty()) {
            if !trace.is_empty() {
                trace.push('.');
            }
            trace.push_str(name);
        }
        trace
    }
}
