//! The typed readers' calls on the object stack of `read`: a caller that builds its own model of
//! a type as it reads it (makeblastdb's mask files, `writedb::mask_info`; the Seq-id,
//! `super::seq_id`) names the type of the table and reads each member, variant or element in a
//! function of its own, while the reader keeps NCBI's frames, members, pointer forms and errors as
//! it keeps them for the generic `Value` (`TypedReader::read_value`, which a typed reader calls
//! for a part it does not build: `value`).
//!
//! A typed reader recurses as the type nests; the types it is given nest a bounded number of
//! levels, but for a CHOICE whose variants are lists of it (a `Seq-loc` of `mix` and `equiv`),
//! which `nest` reads to any depth with its levels on the heap, as NCBI reads it with a stack large
//! enough (PD-LOSAT-NCBI-DEFECTS exception 28). It reads (never skips) types that have neither a
//! delayed member nor a hook of NCBI's classes (`super::hooks`). Each call names its type, so the
//! XML stream's frames and tags of a named type (`TypedReader::begin_names`), of a list and of
//! its elements (`Input::ARRAY_FRAMES`, `Input::element_frame`, `Input::end_element`) are those
//! of the generic reader.

use super::read::{
    alias_target, is_missing_value, set_flag_error, Pointee, PostSink, TypedReader, NO_TYPE,
};
use super::{resolved, type_info, Kind, Member, Value, DELAYED};
use crate::report::serial::asn_input::{Input, R};
use crate::report::serial::object_frame::FrameKind;

/// What [`TypedReader::nest`] hands the reader of a nest, in the order of the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NestPart {
    /// The value of the variant at `index` (of the type `ty`), one that is not a list of the
    /// CHOICE: the reader reads it, in its variant's frame.
    Value(usize, u16),
    /// An element of a list of the CHOICE: a value of the CHOICE begins.
    Enter,
    /// The value that the last `Enter` began was read.
    Leave,
    /// An element of the list of the variant at `index` stated as `NULL` (an empty `CRef<>`):
    /// nothing to read.
    Null(usize),
}

/// A level of the nest that [`TypedReader::nest`] reads (`TypedReader::nest`, on the heap): the
/// list variant of the CHOICE whose elements are being read, whether the element being read was
/// stated as `:Type` (in `ReadPointer`'s named frame, which its end leaves), and the named frames
/// of the CHOICE and of the list's type (`begin_names`: the XML stream's, `Input::NAMED_TYPES`;
/// none in the others).
#[derive(Clone, Copy, Debug)]
pub(super) struct NestLevel {
    variant: u16,
    other: bool,
    names: u8,
    list_names: u8,
}

/// Where [`TypedReader::nest`] stands: at the start of a value of the CHOICE, at the next element
/// of the innermost list, or after the end of a value of the CHOICE.
#[derive(Clone, Copy)]
enum NestStep {
    Begin,
    Element,
    End,
}

/// `$read` (an expression of the reader `$r`) in the named frames of the type `$ty`
/// (`begin_names`: the XML stream's `ReadNamedType`, `Input::NAMED_TYPES`); the text and binary
/// streams read it alone, as the code was before the XML stream (no closure: a closure's call
/// left the reads of those streams outlined, +0.8 % instructions in makeblastdb `-mask_data`).
macro_rules! in_names {
    ($r:ident, $ty:expr, $read:expr) => {
        if I::NAMED_TYPES {
            let count = $r.begin_names($ty, false)?;
            let result = $read;
            $r.end_names(count, result)
        } else {
            $read
        }
    };
}

impl<I: Input, P: PostSink> TypedReader<I, P> {
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
        in_names!(self, ty, self.class_frame(ty, &mut read))
    }

    /// The class frame of `class` and its contents.
    #[inline(always)]
    fn class_frame(
        &mut self,
        ty: u16,
        read: &mut impl FnMut(&mut Self, usize, u16) -> R<()>,
    ) -> R<()> {
        let ty = alias_target(ty);
        let Kind::Sequence(span) = type_info(ty).kind else {
            unreachable!("{} is a SEQUENCE", type_info(ty).name())
        };
        self.push(FrameKind::Class, ty, false);
        let result = self.class_contents(span.members(), read);
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
        self.push(FrameKind::ClassMember, NO_TYPE, false);
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
                    let result = read(self, index, member.ty);
                    if I::MISSING_VALUES {
                        Self::missing_value_or_error(result, member)?;
                    } else if let Err(fail) = result {
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

    /// The result of the member `member` read in `ReadWithSetFlagMember` of a stream that throws
    /// `eMissingValue` (`Input::MISSING_VALUES`): a missing value of an OPTIONAL or DEFAULT member
    /// is taken, the member not set (`read::TypedReader::missing_value`); any other error is
    /// rethrown (`set_flag_error`). The other streams test no missing value (`class_members`: a
    /// branch on the const, so their code is the one it was before the XML stream).
    #[inline(always)]
    fn missing_value_or_error(result: R<()>, member: &Member) -> R<()> {
        match result {
            Err(fail) if member.is_ref() || !member.optional() || !is_missing_value(&fail) => {
                Err(set_flag_error(fail, member))
            }
            _ => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1500
    /// The CHOICE `ty` (`ReadChoiceSimple`): `read` reads the value of the variant at `index`, of
    /// the type `variant_ty`.
    pub(crate) fn choice<T>(
        &mut self,
        ty: u16,
        read: impl FnOnce(&mut Self, usize, u16) -> R<T>,
    ) -> R<T> {
        in_names!(self, ty, self.choice_frame(ty, read))
    }

    /// The choice frame of `choice` and its contents.
    #[inline(always)]
    fn choice_frame<T>(
        &mut self,
        ty: u16,
        read: impl FnOnce(&mut Self, usize, u16) -> R<T>,
    ) -> R<T> {
        let ty = alias_target(ty);
        let Kind::Choice(span) = type_info(ty).kind else {
            unreachable!("{} is a CHOICE", type_info(ty).name())
        };
        let variants = span.members();
        let index = self.choice_begin(ty, variants)?;
        let value = read(self, index, variants[index].ty);
        self.choice_end(value)
    }

    /// The start of the CHOICE `ty` (of the variants `variants`) to its variant's value: the
    /// choice frame, `BeginChoice`, the variant frame, `BeginChoiceVariant` named in it; the
    /// variant's index. An error leaves the frames pushed here before it returns. `choice` and
    /// `nest` begin each CHOICE here.
    #[inline(always)]
    fn choice_begin(&mut self, ty: u16, variants: &'static [Member]) -> R<usize> {
        self.push(FrameKind::Choice, ty, false);
        if let Err(fail) = self.inp.begin_choice() {
            return self.pop(Err(fail));
        }
        self.push(FrameKind::ChoiceVariant, NO_TYPE, false);
        match self.inp.begin_variant(variants) {
            Ok(index) => {
                self.set_member(&variants[index]);
                Ok(index)
            }
            Err(fail) => {
                let left = self.pop(Err(fail));
                self.pop(left)
            }
        }
    }

    /// The end of a CHOICE that `choice_begin` began, after its variant's value `value`:
    /// `EndChoiceVariant`, the variant frame left, `EndChoice`, the choice frame left (an error
    /// leaves them instead).
    #[inline(always)]
    fn choice_end<T>(&mut self, value: R<T>) -> R<T> {
        let value = value.and_then(|value| self.inp.end_variant().map(|()| value));
        let value = self.pop(value);
        let value = value.and_then(|value| self.inp.end_choice().map(|()| value));
        self.pop(value)
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
        in_names!(self, ty, self.list_frame(ty, &mut read))
    }

    /// The array frame of `list` (`Input::ARRAY_FRAMES`) and its contents.
    #[inline(always)]
    fn list_frame(
        &mut self,
        ty: u16,
        read: &mut impl FnMut(&mut Self, Option<u16>) -> R<()>,
    ) -> R<()> {
        let ty = alias_target(ty);
        let Kind::List { element, set, refs } = type_info(ty).kind else {
            unreachable!("{} is a SEQUENCE OF or SET OF", type_info(ty).name())
        };
        self.list_begin(element, set)?;
        let result = self.elements(element, refs, read);
        self.list_end(element, result)
    }

    /// The start of a SEQUENCE OF or (`set`) SET OF of `element`s to its elements: the array frame
    /// (`Input::ARRAY_FRAMES`), `BeginContainer`, the element frame (`Input::element_frame`). An
    /// error leaves the frame pushed here before it returns. `list` and `nest` begin each list
    /// here.
    #[inline(always)]
    fn list_begin(&mut self, element: u16, set: bool) -> R<()> {
        if I::ARRAY_FRAMES {
            self.push(FrameKind::Array, NO_TYPE, false);
        }
        if let Err(fail) = self.inp.begin_container(set) {
            if I::ARRAY_FRAMES {
                return self.pop(Err(fail));
            }
            return Err(fail);
        }
        if I::element_frame(element) {
            self.push(FrameKind::ArrayElement, NO_TYPE, false);
        }
        Ok(())
    }

    /// The end of a list of `element`s that `list_begin` began, after its elements (`elements`):
    /// the element frame left, `EndContainer`, the array frame left (an error leaves them
    /// instead).
    #[inline(always)]
    fn list_end(&mut self, element: u16, elements: R<()>) -> R<()> {
        let mut result = elements;
        if I::element_frame(element) {
            result = self.pop(result);
        }
        let result = result.and_then(|()| self.inp.end_container());
        if I::ARRAY_FRAMES {
            return self.pop(result);
        }
        result
    }

    /// Each element (`BeginContainerElement`) read by `read`.
    #[inline(always)]
    fn elements(
        &mut self,
        element: u16,
        refs: bool,
        read: &mut impl FnMut(&mut Self, Option<u16>) -> R<()>,
    ) -> R<()> {
        while let Some(pointee) = self.element_begin(element, refs)? {
            match pointee {
                Pointee::Null => read(self, None)?,
                Pointee::This => read(self, Some(element))?,
                Pointee::Other => {
                    let result = read(self, Some(element));
                    self.end_pointer(result)?;
                }
            }
            self.inp.end_element()?;
        }
        Ok(())
    }

    /// `BeginContainerElement` of a list of `element`s, and where the element's value is: an
    /// element of a list of `CRef<>`s (`refs`) takes the pointer forms (`begin_pointer`); `None`
    /// at the end of the list. `elements` and `nest` begin each element here.
    #[inline(always)]
    fn element_begin(&mut self, element: u16, refs: bool) -> R<Option<Pointee>> {
        if !self.inp.begin_element()? {
            return Ok(None);
        }
        if !refs {
            return Ok(Some(Pointee::This));
        }
        self.begin_pointer(element, false).map(Some)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1482-1500
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1245-1274
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:955-1000
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1173-1198
    /// A value of the CHOICE `ty` whose variants at `nests` are lists of the CHOICE (`Seq-loc`'s
    /// `mix` and `equiv`), nested to any depth: NCBI reads it by recursion, `ReadChoiceSimple`,
    /// `ReadContainer` (text and binary) and, for an element stated as `:Type`, `ReadPointer`'s
    /// named frame, each level a few frames of the stack of its thread. Here the levels are kept
    /// on the heap (`self.nest`) and each frame on the object stack, begun and ended by the steps
    /// of `choice` and `list` (`begin_names`, `choice_begin`, `choice_end`, `list_begin`,
    /// `element_begin`, `end_pointer`, `end_element`, `list_end`, `end_names`), so that a nest is
    /// read to any depth with LOSAT's stack as it is.
    /// `read` reads the value of every other variant (`NestPart::Value`, in its variant's frame)
    /// and is told where a value of the CHOICE begins and ends in a list (`Enter`, `Leave`) and of
    /// an element stated as `NULL` (`Null`). An error leaves the frames of every open level from
    /// the innermost out (`unwind`), as it leaves NCBI's recursion: the frames of the levels above
    /// the error each add what they add to it.
    /// A value that is not a list (most of them: an interval, a `packed-int`) is read by the
    /// steps of `choice` alone (`nest_begin`); the levels begin at its first list.
    pub(crate) fn nest(
        &mut self,
        ty: u16,
        nests: impl Fn(usize) -> bool,
        mut read: impl FnMut(&mut Self, NestPart) -> R<()>,
    ) -> R<()> {
        let named = ty;
        let ty = alias_target(ty);
        let Kind::Choice(span) = type_info(ty).kind else {
            unreachable!("{} is a CHOICE", type_info(ty).name())
        };
        let variants = span.members();
        let depth = self.depth();
        match self.nest_begin(named, ty, variants, &nests, &mut read) {
            Ok(false) => return Ok(()),
            Ok(true) => {}
            Err(fail) => return Err(self.unwind(depth, fail)),
        }
        let levels = self.nest.len() - 1;
        let result = self.nest_steps(named, ty, variants, &nests, &mut read, levels);
        self.nest.truncate(levels);
        result.map_err(|fail| self.unwind(depth, fail))
    }

    /// The start of a value of the CHOICE `ty` (named `named`) in a nest: its named frames and
    /// `choice_begin`; a variant that is not a list read to the end of the value (`false`), or, for
    /// a list, its named frames and `list_begin`, as a new level of `self.nest` (`true`). An error
    /// can leave frames pushed here on the stack, which `nest` leaves (`unwind`).
    #[inline(always)]
    fn nest_begin(
        &mut self,
        named: u16,
        ty: u16,
        variants: &'static [Member],
        nests: &impl Fn(usize) -> bool,
        read: &mut impl FnMut(&mut Self, NestPart) -> R<()>,
    ) -> R<bool> {
        let names = self.names_begin(named)?;
        let index = self.choice_begin(ty, variants)?;
        if !nests(index) {
            let value = read(self, NestPart::Value(index, variants[index].ty));
            let value = self.choice_end(value);
            self.names_end(names, value)?;
            return Ok(false);
        }
        let list_names = self.names_begin(variants[index].ty)?;
        let (element, set, _) = Self::nest_list(ty, variants, index);
        self.list_begin(element, set)?;
        self.nest.push(NestLevel {
            variant: index as u16,
            other: false,
            names,
            list_names,
        });
        Ok(true)
    }

    /// The element type, `set` and `refs` of the list of the variant `variant` of the CHOICE `ty`
    /// (of the variants `variants`) in a nest: a list of the CHOICE itself.
    #[inline(always)]
    fn nest_list(ty: u16, variants: &'static [Member], variant: usize) -> (u16, bool, bool) {
        let list = alias_target(variants[variant].ty);
        let Kind::List { element, set, refs } = type_info(list).kind else {
            unreachable!("{} is a SEQUENCE OF or SET OF", type_info(list).name())
        };
        debug_assert_eq!(alias_target(element), ty, "a list of the CHOICE itself");
        (element, set, refs)
    }

    /// The steps of `nest` above its `base` levels, from the first element of the innermost list to
    /// the end of its value or an error, which leaves the frames of the open levels on the stack.
    #[inline(always)]
    fn nest_steps(
        &mut self,
        named: u16,
        ty: u16,
        variants: &'static [Member],
        nests: &impl Fn(usize) -> bool,
        read: &mut impl FnMut(&mut Self, NestPart) -> R<()>,
        base: usize,
    ) -> R<()> {
        let mut step = NestStep::Element;
        loop {
            step = match step {
                NestStep::Begin => {
                    if self.nest_begin(named, ty, variants, nests, read)? {
                        NestStep::Element
                    } else {
                        NestStep::End
                    }
                }
                NestStep::Element => {
                    let level = self.nest.len() - 1;
                    let variant = usize::from(self.nest[level].variant);
                    let (element, _, refs) = Self::nest_list(ty, variants, variant);
                    match self.element_begin(element, refs)? {
                        Some(Pointee::Null) => {
                            read(self, NestPart::Null(variant))?;
                            self.inp.end_element()?;
                            NestStep::Element
                        }
                        Some(pointee) => {
                            self.nest[level].other = pointee == Pointee::Other;
                            read(self, NestPart::Enter)?;
                            NestStep::Begin
                        }
                        None => {
                            let NestLevel {
                                names, list_names, ..
                            } = self.nest[level];
                            let result = self.list_end(element, Ok(()));
                            self.names_end(list_names, result)?;
                            self.nest.pop();
                            let result = self.choice_end(Ok(()));
                            self.names_end(names, result)?;
                            NestStep::End
                        }
                    }
                }
                NestStep::End => {
                    if self.nest.len() == base {
                        return Ok(());
                    }
                    read(self, NestPart::Leave)?;
                    if self.nest[self.nest.len() - 1].other {
                        self.end_pointer(Ok(()))?;
                    }
                    self.inp.end_element()?;
                    NestStep::Element
                }
            };
        }
    }

    /// `begin_names` of the type `ty` in a stream with named frames (`Input::NAMED_TYPES`, as
    /// `in_names!`): how many it pushed; none in the others.
    #[inline(always)]
    fn names_begin(&mut self, ty: u16) -> R<u8> {
        if I::NAMED_TYPES {
            return self.begin_names(ty, false);
        }
        Ok(0)
    }

    /// `end_names` of the `count` frames of `names_begin`, after `result` (as `in_names!`).
    #[inline(always)]
    fn names_end<T>(&mut self, count: u8, result: R<T>) -> R<T> {
        if I::NAMED_TYPES {
            return self.end_names(count, result);
        }
        result
    }

    /// An INTEGER of the type `ty`, read as the C++ storage of its member (`ReadStd`).
    #[inline(always)]
    pub(crate) fn integer(&mut self, ty: u16) -> R<i64> {
        let Kind::Integer(storage) = resolved(ty).kind else {
            unreachable!("{} is an INTEGER", type_info(ty).name())
        };
        in_names!(self, ty, self.read_integer(storage))
    }

    /// An ENUMERATED (or INTEGER with named values) of the type `ty`
    /// (`CEnumeratedTypeInfo::ReadEnum`).
    #[inline(always)]
    pub(crate) fn enumerated(&mut self, ty: u16) -> R<i32> {
        let Kind::Enumerated { integer, values } = resolved(ty).kind else {
            unreachable!("{} is an ENUMERATED", type_info(ty).name())
        };
        in_names!(self, ty, self.enum_number(values.values(), integer, false))
    }

    /// A string of the type `ty`, with NCBI's posts (`ReadString`).
    #[inline(always)]
    pub(crate) fn text(&mut self, ty: u16) -> R<Vec<u8>> {
        let Kind::String(kind) = resolved(ty).kind else {
            unreachable!("{} is a string", type_info(ty).name())
        };
        in_names!(self, ty, self.string(kind, false))
    }

    /// A string of the type `ty`, with NCBI's posts, kept even while the reader discards other
    /// values ([`super::read::ObjectStream::read_typed_member`]).
    #[inline(always)]
    pub(crate) fn text_kept(&mut self, ty: u16) -> R<Vec<u8>> {
        let Kind::String(kind) = resolved(ty).kind else {
            unreachable!("{} is a string", type_info(ty).name())
        };
        in_names!(self, ty, self.string_kept(kind))
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
