//! The typed readers' calls on the object stack of `read`: a caller that builds its own model of
//! a type as it reads it (makeblastdb's mask files, `writedb::mask_info`; the Seq-id,
//! `super::seq_id`) names the type of the table and reads each member, variant or element in a
//! function of its own, while the reader keeps NCBI's frames, members, pointer forms and errors as
//! it keeps them for the generic `Value` (`TypedReader::read_value`, which a typed reader calls
//! for a part it does not build: `value`).
//!
//! A typed reader recurses as the type nests; the types it is given nest a bounded number of
//! levels but for a `Seq-loc` of `mix` and `equiv`, whose depth the stream's limit bounds
//! (`ObjectStream::limit_depth`). It reads (never skips) types that have neither a delayed member
//! nor a hook of NCBI's classes (`super::hooks`).

use super::read::{alias_target, set_flag_error, Pointee, TypedReader, NO_TYPE};
use super::{resolved, type_info, Kind, Member, Value, DELAYED};
use crate::report::serial::asn_input::{Input, R};
use crate::report::serial::object_frame::FrameKind;

impl<I: Input> TypedReader<I> {
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1380-1397 (`ReadClassSequential`;
    // the binary stream's, objistrasnb.cpp:1210-1233, is the same with its own `BeginClass`)
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:84-109
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:871-876
    /// The SEQUENCE `ty` read member by member in module order (`ReadClassSequential`): `read`
    /// reads the value of the member at `index` of the type `member_ty`. A member with a set flag
    /// reads it in `ReadWithSetFlagMember` (`read::set_flag_error`); a member held by `CRef<>`
    /// that is empty when NCBI reads it takes the pointer forms (`member_pointer`): `read` is not
    /// called for a `NULL`, which leaves the member unset.
    pub(crate) fn class(
        &mut self,
        ty: u16,
        mut read: impl FnMut(&mut Self, usize, u16) -> R<()>,
    ) -> R<()> {
        let ty = alias_target(ty);
        let Kind::Sequence(span) = type_info(ty).kind else {
            unreachable!("{} is a SEQUENCE", type_info(ty).name())
        };
        self.push(FrameKind::Class, ty, false)?;
        let result = self.class_contents(span.members(), &mut read);
        self.pop(result)
    }

    /// `BeginClass`, the members in their frame, the mandatory members after the last checked,
    /// `EndClass`.
    #[inline(always)]
    fn class_contents(
        &mut self,
        members: &'static [Member],
        read: &mut impl FnMut(&mut Self, usize, u16) -> R<()>,
    ) -> R<()> {
        self.inp.begin_class()?;
        self.push(FrameKind::ClassMember, NO_TYPE, false)?;
        let read = self.class_members(members, read);
        let pos = self.pop(read)?;
        for member in &members[pos..] {
            self.expected_member(member)?;
        }
        self.inp.end_class()
    }

    /// The members from the first on (`member_begin`), each read by `read`; the position after
    /// the last one read.
    #[inline(always)]
    fn class_members(
        &mut self,
        members: &'static [Member],
        read: &mut impl FnMut(&mut Self, usize, u16) -> R<()>,
    ) -> R<usize> {
        let mut pos = 0;
        while let Some(index) = self.member_begin(members, pos)? {
            let member = &members[index];
            debug_assert!(
                member.flags & DELAYED == 0,
                "a typed reader reads no delayed member"
            );
            match self.member_pointer(member, false)? {
                Pointee::Null => {}
                Pointee::This => {
                    if let Err(fail) = read(self, index, member.ty) {
                        return Err(set_flag_error(fail, member));
                    }
                }
                Pointee::Other => {
                    let result = read(self, index, member.ty);
                    self.end_pointer(result)?;
                }
            }
            pos = index + 1;
            self.inp.end_member()?;
        }
        Ok(pos)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1500
    /// The CHOICE `ty` (`ReadChoiceSimple`): `read` reads the value of the variant at `index`, of
    /// the type `variant_ty`.
    pub(crate) fn choice<T>(
        &mut self,
        ty: u16,
        read: impl FnOnce(&mut Self, usize, u16) -> R<T>,
    ) -> R<T> {
        let ty = alias_target(ty);
        let Kind::Choice(span) = type_info(ty).kind else {
            unreachable!("{} is a CHOICE", type_info(ty).name())
        };
        self.push(FrameKind::Choice, ty, false)?;
        let value = self.choice_contents(span.members(), read);
        self.pop(value)
    }

    /// `BeginChoice`, the variant in its frame, `EndChoice`.
    #[inline(always)]
    fn choice_contents<T>(
        &mut self,
        variants: &'static [Member],
        read: impl FnOnce(&mut Self, usize, u16) -> R<T>,
    ) -> R<T> {
        self.inp.begin_choice()?;
        self.push(FrameKind::ChoiceVariant, NO_TYPE, false)?;
        let value = self.variant(variants, read);
        let value = self.pop(value)?;
        self.inp.end_choice()?;
        Ok(value)
    }

    /// `BeginChoiceVariant`, named in the variant frame, its value, `EndChoiceVariant`.
    #[inline(always)]
    fn variant<T>(
        &mut self,
        variants: &'static [Member],
        read: impl FnOnce(&mut Self, usize, u16) -> R<T>,
    ) -> R<T> {
        let index = self.inp.begin_variant(variants)?;
        let variant = &variants[index];
        self.set_member(variant);
        let value = read(self, index, variant.ty)?;
        self.inp.end_variant()?;
        Ok(value)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1318-1346
    // NCBI reference (598d8ae6): c++/include/serial/impl/stltypes.hpp:246-263
    /// The SEQUENCE OF or SET OF `ty` (`ReadContainer`; a named one is NCBI's implicit class,
    /// read as its list without a class frame): `read` reads one element of the type it is given;
    /// an element of a list of `CRef<>`s (`AddElementIn` pushes an empty one) takes the pointer
    /// forms, and one stated as `NULL` is handed to `read` as `None`, with nothing to read.
    pub(crate) fn list(
        &mut self,
        ty: u16,
        mut read: impl FnMut(&mut Self, Option<u16>) -> R<()>,
    ) -> R<()> {
        let ty = alias_target(ty);
        let Kind::List { element, set, refs } = type_info(ty).kind else {
            unreachable!("{} is a SEQUENCE OF or SET OF", type_info(ty).name())
        };
        self.push(FrameKind::Array, NO_TYPE, false)?;
        let result = self.list_contents(element, set, refs, &mut read);
        self.pop(result)
    }

    /// `BeginContainer`, the elements in their frame, `EndContainer`.
    #[inline(always)]
    fn list_contents(
        &mut self,
        element: u16,
        set: bool,
        refs: bool,
        read: &mut impl FnMut(&mut Self, Option<u16>) -> R<()>,
    ) -> R<()> {
        self.inp.begin_container(set)?;
        self.push(FrameKind::ArrayElement, NO_TYPE, false)?;
        let result = self.elements(element, refs, read);
        self.pop(result)?;
        self.inp.end_container()
    }

    /// Each element (`BeginContainerElement`) read by `read`.
    #[inline(always)]
    fn elements(
        &mut self,
        element: u16,
        refs: bool,
        read: &mut impl FnMut(&mut Self, Option<u16>) -> R<()>,
    ) -> R<()> {
        while self.inp.begin_element()? {
            let pointee = if refs {
                self.begin_pointer(element, false)?
            } else {
                Pointee::This
            };
            match pointee {
                Pointee::Null => read(self, None)?,
                Pointee::This => read(self, Some(element))?,
                Pointee::Other => {
                    let result = read(self, Some(element));
                    self.end_pointer(result)?;
                }
            }
        }
        Ok(())
    }

    /// An INTEGER of the type `ty`, read as the C++ storage of its member (`ReadStd`).
    #[inline(always)]
    pub(crate) fn integer(&mut self, ty: u16) -> R<i64> {
        let Kind::Integer(storage) = resolved(ty).kind else {
            unreachable!("{} is an INTEGER", type_info(ty).name())
        };
        self.read_integer(storage)
    }

    /// An ENUMERATED (or INTEGER with named values) of the type `ty`
    /// (`CEnumeratedTypeInfo::ReadEnum`).
    #[inline(always)]
    pub(crate) fn enumerated(&mut self, ty: u16) -> R<i32> {
        let Kind::Enumerated { integer, values } = resolved(ty).kind else {
            unreachable!("{} is an ENUMERATED", type_info(ty).name())
        };
        self.enum_number(values.values(), integer, false)
    }

    /// A string of the type `ty`, with NCBI's posts (`ReadString`).
    #[inline(always)]
    pub(crate) fn text(&mut self, ty: u16) -> R<Vec<u8>> {
        let Kind::String(kind) = resolved(ty).kind else {
            unreachable!("{} is a string", type_info(ty).name())
        };
        self.string(kind, false)
    }

    /// A BOOLEAN (`ReadBool`).
    #[inline(always)]
    pub(crate) fn boolean(&mut self) -> R<bool> {
        self.inp.read_bool()
    }

    /// A NULL (`ReadNull`).
    #[inline(always)]
    pub(crate) fn null(&mut self) -> R<()> {
        self.inp.read_null()
    }

    /// A value of the type `ty` that the typed reader does not build: the generic reader's
    /// (`read_value`).
    pub(crate) fn value(&mut self, ty: u16) -> R<Value> {
        self.read_value(ty, false)
    }
}
