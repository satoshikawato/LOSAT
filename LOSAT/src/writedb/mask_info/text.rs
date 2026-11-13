//! ASN.1 text (`MSerial_AsnText`): the mask reader's [`Input`] over the shared lexical layer of
//! NCBI's ASN.1 text object stream (`report::serial::asn_text_reader`, `CObjectIStreamAsn` over
//! `CIStreamBuffer`). What is particular to the mask reader stays here: the members and variants
//! found by name or by tag in its schema, `m_BlockStart`, and the decision of the frames on the
//! end of the text (`Fail::Eof`).
//!
//! The text is read from a stream in blocks (`StreamSource`): one reader for all objects of a
//! file, as NCBI's object stream gives back the bytes it did not use when an object ends
//! (`CIStreamBuffer::Close`, strbuffer.cpp:125-139), so that the next object starts at the byte
//! after the last one read. Each object gets a new `CObjectIStream` in NCBI (`s_ReadObject`), so
//! its line count starts at 1 again (checked with makeblastdb 2.17.0: a string on the first line
//! of a second object is "at line 1").

use std::io::Read;

use super::reader::{format_err, Fail, Input, StringNote, R};
use super::schema::{ChoiceDef, ClassDef, EnumDef};
use crate::report::serial::asn_text_reader::{
    AsnSource, AsnTextError, AsnTextErrorKind, MemberId, Reader, StreamSource,
};
use crate::report::serial::object_stack::ObjectStack;

/// The text of a mask file.
pub(super) type TextSource = StreamSource<Box<dyn Read>>;

pub(super) struct TextIn<S: AsnSource> {
    r: Reader<S>,
    /// `m_BlockStart`
    block_start: bool,
}

impl From<AsnTextError> for Fail {
    fn from(error: AsnTextError) -> Self {
        match error.kind {
            AsnTextErrorKind::Eof => Fail::Eof,
            AsnTextErrorKind::Format => Fail::Format(error.message),
            AsnTextErrorKind::Unsupported => Fail::Unsupported(error.message),
        }
    }
}

impl<S: AsnSource> TextIn<S> {
    pub(super) fn new(src: S) -> Self {
        TextIn {
            r: Reader::new(src, "mask data", "mask data reader"),
            block_start: false,
        }
    }

    /// The state of a new `CObjectIStreamAsn` at the next byte: line 1, no block started.
    pub(super) fn new_object(&mut self) {
        self.r.line = 1;
        self.block_start = false;
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:389-424
    // ```c++
    // TMemberIndex CObjectIStreamAsn::GetMemberIndex
    //     (const CClassTypeInfo* classType,
    //      const CTempString& id,
    //      const TMemberIndex pos)
    // {
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
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:432-439,462-473
    // ```c++
    // TMemberIndex CItemsInfo::Find(const CTempString& name, TMemberIndex pos) const
    // {
    //     for ( CIterator i(*this, pos); i.Valid(); ++i ) {
    //         if ( name == GetItemInfo(i)->GetId().GetName() )
    //             return *i;
    //     }
    //     return kInvalidMember;
    // }
    // ...
    //         TMemberIndex index = tag + zero_index;
    //         if ( index < pos || index > LastIndex() )
    //             return kInvalidMember;
    //         return index;
    // ```
    // (The members and variants of these types have the automatic tags 0, 1, ... in module
    // order. Their names start with a lower case letter, so `GetAltItemIndex`, which looks up the
    // name with an upper case first letter, finds nothing.)
    /// `ReadMemberId(SkipWhiteSpace())` and the member or variant it names among `names`, searched
    /// (`count` of them, by position) from `pos` on; one that is not there is an error (`UnexpectedMember`).
    fn find_member(
        &mut self,
        count: usize,
        name: impl Fn(usize) -> &'static str,
        pos: usize,
    ) -> R<usize> {
        let (id, found) = match self.r.member_or_tag()? {
            MemberId::Tag(digits) => {
                let id = self.r.bytes(digits);
                let tag = std::str::from_utf8(id)
                    .ok()
                    .and_then(|s| s.parse::<u32>().ok())
                    .ok_or_else(|| Fail::Format("bad member number".into()))?
                    as usize;
                (id, (pos..count).contains(&tag).then_some(tag))
            }
            MemberId::Name(id) => {
                let id = self.r.bytes(id);
                let found = (pos..count).find(|&i| name(i).as_bytes() == id);
                (id, found)
            }
        };
        found.ok_or_else(|| {
            Fail::Format(format!(
                "\"{}\": unexpected member",
                String::from_utf8_lossy(id)
            ))
        })
    }
}

#[cfg(test)]
impl TextIn<TextSource> {
    /// The bytes the stream keeps.
    pub(super) fn buffered(&self) -> usize {
        self.r.src.buffered()
    }
}

impl<S: AsnSource> Input for TextIn<S> {
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        let name = self.r.file_header()?;
        Ok(Some(self.r.bytes(name).to_vec()))
    }

    fn begin_class(&mut self) -> R<()> {
        self.block_start = self.r.start_block()?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1337-1362
    // ```c++
    // TMemberIndex
    // CObjectIStreamAsn::BeginClassMember(const CClassTypeInfo* classType,
    //                                     TMemberIndex pos)
    // {
    //     if ( !NextElement() )
    //         return kInvalidMember;
    //
    //     CTempString id = ReadMemberId(SkipWhiteSpace());
    //     TMemberIndex index = GetMemberIndex(classType, id, pos);
    //     if ( index == kInvalidMember ) {
    //         if (m_TypeAlias && classType->GetMembers().GetItemInfo(pos)->GetId().HasNotag()) {
    //             ...
    //         } else if (CanSkipUnknownMembers()) {
    //             ...
    //         } else {
    //             UnexpectedMember(id, classType->GetMembers());
    //         }
    //     }
    //     return index;
    // }
    // ```
    // (`CanSkipUnknownMembers()` is false by default, `m_TypeAlias` is null.)
    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        if !self.r.next_element(&mut self.block_start)? {
            return Ok(None);
        }
        Ok(Some(self.find_member(
            def.members.len(),
            |i| def.members[i].name,
            pos,
        )?))
    }

    fn end_member(&mut self) -> R<()> {
        Ok(())
    }

    fn end_class(&mut self) -> R<()> {
        Ok(self.r.end_block()?)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1478
    // ```c++
    // void CObjectIStreamAsn::BeginChoice(const CChoiceTypeInfo* choiceType)
    // {
    //     if (choiceType->GetVariantInfo(kFirstMemberIndex)->GetId().IsAttlist()) {
    //         ...
    //     }
    //     m_BlockStart = true;
    // }
    // ...
    // TMemberIndex CObjectIStreamAsn::BeginChoiceVariant(const CChoiceTypeInfo* choiceType)
    // {
    //     bool skipname = !m_BlockStart;
    //     if ( !NextElement() )
    //         return kInvalidMember;
    //     CTempString id = ReadMemberId(SkipWhiteSpace());
    //     if (skipname) {
    //         id = ReadMemberId(SkipWhiteSpace());
    //     }
    //     if ( id.empty() )
    //         ThrowError(fFormatError, "choice variant id expected");
    //
    //     TMemberIndex index = GetChoiceIndex(choiceType, id);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownVariants()) {
    //             SetFailFlags(fUnknownValue);
    //         } else {
    //             UnexpectedMember(id, choiceType->GetVariants());
    //         }
    //     }
    //
    //     return index;
    // }
    // ...
    // void CObjectIStreamAsn::EndChoice(void)
    // {
    //     ...
    //     m_BlockStart = false;
    // }
    // ```
    // (BeginChoice has just set `m_BlockStart`, so `skipname` is false; an empty id is not a
    // name of the module, so it fails as "choice variant id expected" here.)
    fn begin_choice(&mut self) -> R<()> {
        self.block_start = true;
        Ok(())
    }

    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize> {
        if !self.r.next_element(&mut self.block_start)? {
            return format_err("choice variant id expected");
        }
        self.find_member(def.variants.len(), |i| def.variants[i].name, 0)
    }

    fn end_variant(&mut self) -> R<()> {
        Ok(())
    }

    fn end_choice(&mut self) -> R<()> {
        self.block_start = false;
        Ok(())
    }

    fn begin_container(&mut self, _set: bool) -> R<()> {
        self.begin_class()
    }

    fn begin_element(&mut self) -> R<bool> {
        Ok(self.r.next_element(&mut self.block_start)?)
    }

    fn end_container(&mut self) -> R<()> {
        self.end_class()
    }

    fn read_i32(&mut self) -> R<i32> {
        Ok(self.r.int4()?)
    }

    fn read_i64(&mut self) -> R<i64> {
        Ok(self.r.int8()?)
    }

    fn read_u32(&mut self) -> R<u32> {
        Ok(self.r.uint4()?)
    }

    fn read_bool(&mut self) -> R<bool> {
        Ok(self.r.boolean()?)
    }

    fn read_null(&mut self) -> R<()> {
        Ok(self.r.null()?)
    }

    fn read_string(&mut self, notes: &mut Vec<StringNote>) -> R<Vec<u8>> {
        let s = self.r.string();
        notes.extend(self.r.take_bad_chars().into_iter().map(StringNote::Bad));
        Ok(s?)
    }

    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        Ok(if def.integer {
            self.r.named_integer(def.values)?
        } else {
            self.r.enumerated(def.values)?.1
        })
    }

    fn stack(&mut self) -> &mut ObjectStack {
        &mut self.r.stack
    }
}
