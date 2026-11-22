//! The one generic writer of the values of the table (`Value`, read by `read`) over an
//! `Encoder`: what NCBI's object stream writes of the object that its reader read, so that
//! psiblast's `-out_pssm`, the search strategy's request and the object manager's messages write
//! every member of a PSSM file as NCBI does (line RG item D3-2). The typed writers that own a part
//! of the value keep it: a Seq-id is written by the Seq-id model's writer
//! (`report::seqid::SeqId::Db`), whose model the reader's hook stored (`Value::SeqId`).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:687-714
//! A member is written when it is set (`WriteMember` skips an OPTIONAL member that is not set,
//! and a DEFAULT member only when it was not set: a DEFAULT member that was read is written,
//! whatever its value, member.cpp:973-1013), so the members of a SEQUENCE value are those that
//! the reader set, in module order. NCBI writes an ENUMERATED value by its name and an INTEGER
//! with named values by its name when it has one, else by its number (`WriteEnum`,
//! objostrasn.cpp:113-131); the C++ storage of an INTEGER decides `WriteInt4`, `WriteUint4` or
//! `WriteInt8` (the same decimal digits in ASN.1 text). A value nests to any depth (the reader
//! reads any depth), so it is written from a list on the heap, never by a recursion as deep as
//! the nest.

use super::{resolved, type_info, Kind, Member, Storage, StringKind, Value};
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, Traced};

/// What is left to write: a value of a type, or one call of the encoder.
enum Step<'a> {
    Value(u16, &'a Value),
    Member(&'static Member),
    EndMember,
    EndSequence,
    Element,
    EndList,
    EndVariant,
    EndNamed,
}

/// Writes `value`, a value of the type `ty` of the table, with `enc` (as the value of the member
/// or element that the caller opened, or of a document that it began).
pub(crate) fn write_value<E: Encoder>(enc: &mut E, ty: u16, value: &Value) {
    write_until(enc, ty, value, |_| false);
}

/// `write_value` until `stop(enc)` holds before a step (a value opened or one call of the
/// encoder): whether the value was written whole.
fn write_until<E: Encoder>(enc: &mut E, ty: u16, value: &Value, stop: impl Fn(&E) -> bool) -> bool {
    let mut todo = vec![Step::Value(ty, value)];
    while let Some(step) = todo.pop() {
        if stop(enc) {
            return false;
        }
        match step {
            Step::Value(ty, value) => open(enc, ty, value, &mut todo),
            Step::Member(member) => enc.member(member.name(), member.tag),
            Step::EndMember => enc.end_member(),
            Step::EndSequence => enc.end_sequence(),
            Step::Element => enc.element(),
            Step::EndList => enc.end_list(),
            Step::EndVariant => enc.end_variant(),
            Step::EndNamed => enc.end_named(),
        }
    }
    true
}

/// `value` as one document of its type (an object written to a stream after `MSerial_AsnText`:
/// `WriteFileHeader`, the
/// object, `EndOfWrite` and the separator), for the named type `ty`.
#[cfg(test)]
pub(crate) fn write_document<E: Encoder>(enc: &mut E, ty: u16, value: &Value) {
    write_document_until(enc, ty, value, |_| false);
}

/// `write_document` until `stop(enc)` holds before a step (`write_until`); the document's end
/// only when the value was written whole.
fn write_document_until<E: Encoder>(
    enc: &mut E,
    ty: u16,
    value: &Value,
    stop: impl Fn(&E) -> bool,
) {
    enc.begin_document(type_info(ty).name());
    if write_until(enc, ty, value, stop) {
        enc.end_document();
    }
}

/// What an output stream that has failed gets of `value` written as a document of the named type
/// `ty` (an object written after `MSerial_AsnText` in a message of NCBI's diagnostic stream after
/// it failed, `pssm_search::annot_index`): the bytes that NCBI's buffer held when it first wrote to
/// the stream, which throws (`out_buffer::OutBuffer::first_write`): the document without the end of
/// line after it when it fits the buffer of 4096 bytes, else its start, cut where the buffer first
/// flushed. NCBI's writing stops there: the first pass follows its buffer to find how many bytes it
/// held (`Traced`), the second writes those bytes as they were then
/// (`AsnTextEncoder::failed_after`: a string's wrap that would come later inserts no end of line
/// before the cut). Both stop at the cut, so the text in memory is the cut and at most one step
/// more.
pub(crate) fn failed_stream_text(ty: u16, value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    let mut traced = Traced::new(AsnTextEncoder::new(&mut out), 0);
    write_document_until(&mut traced, ty, value, |enc| enc.first_write().is_some());
    let held = traced.first_write().unwrap_or(0);
    drop(traced);
    out.clear();
    // The same encoder as the first pass (`Traced` writes NCBI's `#` for a byte that the stream
    // replaces), so both passes write the same bytes up to the cut.
    let mut traced = Traced::new(AsnTextEncoder::failed_after(&mut out, held), 0);
    write_document_until(&mut traced, ty, value, |enc| enc.held_len() >= held);
    drop(traced);
    out.truncate(held);
    out
}

/// Writes what `value` starts with and leaves the rest to `todo` (in the reverse order of the
/// calls).
fn open<'a, E: Encoder>(enc: &mut E, ty: u16, value: &'a Value, todo: &mut Vec<Step<'a>>) {
    let info = type_info(ty);
    // NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:716-733 (cited at
    // `AsnTextEncoder::delayed`): the delay buffer is written in place of the member's value,
    // without its type.
    if let Value::Delayed(text) = value {
        enc.delayed(text);
        return;
    }
    // A named type that is not a SEQUENCE or a CHOICE (an alias, a SEQUENCE OF or SET OF, an
    // ENUMERATED, an INTEGER, a string) is written inside `begin_named` (`WriteNamedType`,
    // objostr.cpp:742-756); the SEQUENCE and CHOICE calls carry their names.
    if let Kind::Alias(target) = info.kind {
        enc.begin_named(info.name());
        todo.push(Step::EndNamed);
        todo.push(Step::Value(target, value));
        return;
    }
    if let Value::SeqId(id) = value {
        crate::report::seqid::SeqId::Db(id).write(enc);
        return;
    }
    match (&info.kind, value) {
        (Kind::Sequence(span), Value::Sequence(set)) => {
            let members = span.members();
            enc.begin_sequence(info.name());
            todo.push(Step::EndSequence);
            for (position, value) in set.iter().rev() {
                let member = &members[usize::from(*position)];
                todo.push(Step::EndMember);
                todo.push(Step::Value(member.ty, value));
                todo.push(Step::Member(member));
            }
        }
        (Kind::Choice(span), Value::Choice(index, value)) => {
            let variant = &span.members()[usize::from(*index)];
            enc.variant(info.name(), variant.name(), variant.tag);
            todo.push(Step::EndVariant);
            todo.push(Step::Value(variant.ty, value));
        }
        (kind, value) => {
            let named = !info.name().is_empty();
            if named {
                enc.begin_named(info.name());
            }
            if let Kind::List { element, set } = kind {
                enc.begin_list(*set);
                if named {
                    todo.push(Step::EndNamed);
                }
                todo.push(Step::EndList);
                write_elements(enc, *element, value, todo);
                return;
            }
            write_leaf(enc, kind, value);
            if named {
                enc.end_named();
            }
        }
    }
}

/// The elements of a SEQUENCE OF or SET OF value of `element`s: the numbers of a list of INTEGER
/// or REAL at once, the other values left to `todo`.
fn write_elements<'a, E: Encoder>(
    enc: &mut E,
    element: u16,
    value: &'a Value,
    todo: &mut Vec<Step<'a>>,
) {
    match value {
        Value::List(items) => {
            for item in items.iter().rev() {
                todo.push(Step::Value(element, item));
                todo.push(Step::Element);
            }
        }
        Value::Int4s(numbers) => {
            for &number in numbers {
                enc.element();
                write_number(enc, element, Number::Int(i64::from(number)));
            }
        }
        Value::Ints(numbers) => {
            for &number in numbers {
                enc.element();
                write_number(enc, element, Number::Int(number));
            }
        }
        Value::Reals(numbers) => {
            for &number in numbers {
                enc.element();
                write_number(enc, element, Number::Real(number));
            }
        }
        Value::Texts(runs) => {
            // The elements are Seq-ids (a CHOICE) or strings; each is made, written and dropped.
            let seq_ids = matches!(resolved(element).kind, Kind::Choice(_));
            for run in runs {
                run.for_each(&mut |text| {
                    enc.element();
                    if seq_ids {
                        let id = crate::seqid::parse::seq_id_from_text(
                            text,
                            crate::seqid::accguide::parse_flags::ANY_RAW,
                            &mut |_| {},
                        )
                        .expect("a text that the modifier checked when it made the list");
                        write_value(enc, element, &Value::SeqId(Box::new(id)));
                    } else {
                        write_value(enc, element, &Value::Str(text.to_vec()));
                    }
                });
            }
        }
        _ => unreachable!("a SEQUENCE OF or SET OF value is a list"),
    }
}

/// A number of a list that the reader kept as the numbers (`Value::Int4s`, `Ints`, `Reals`).
#[derive(Clone, Copy)]
enum Number {
    Int(i64),
    Real(f64),
}

/// An element of a list of numbers of the type `ty`, inside the named types of its aliases (the
/// reader read it under `resolved(ty)`).
fn write_number<E: Encoder>(enc: &mut E, mut ty: u16, number: Number) {
    let mut named = 0;
    let kind = loop {
        let info = type_info(ty);
        if !info.name().is_empty() {
            enc.begin_named(info.name());
            named += 1;
        }
        match info.kind {
            Kind::Alias(target) => ty = target,
            ref kind => break kind,
        }
    };
    match (kind, number) {
        (Kind::Integer(storage), Number::Int(number)) => write_int(enc, *storage, number),
        (Kind::Real, Number::Real(number)) => enc.real(number),
        _ => unreachable!("a list of numbers of INTEGER or REAL"),
    }
    for _ in 0..named {
        enc.end_named();
    }
}

/// An INTEGER of the C++ storage `storage` (`WriteInt4`, `WriteUint4`, `WriteInt8`; the reader
/// kept it in the storage's range).
fn write_int<E: Encoder>(enc: &mut E, storage: Storage, number: i64) {
    match storage {
        Storage::Int2 | Storage::Int4 => enc.int(number as i32),
        Storage::Uint4 => match i32::try_from(number) {
            Ok(number) => enc.int(number),
            Err(_) => enc.big_int(number),
        },
        Storage::Int8 => enc.big_int(number),
    }
}

/// A value of a type that holds no other values.
fn write_leaf<E: Encoder>(enc: &mut E, kind: &Kind, value: &Value) {
    match (kind, value) {
        (Kind::Integer(storage), Value::Int(number)) => write_int(enc, *storage, *number),
        // NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:113-131
        // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:235-237 (`ValueToName`)
        // `FindNameEx` takes the last name of the value in module order (each insert replaces
        // the one before); the reader takes an ENUMERATED number only when it has a name.
        (Kind::Enumerated { integer, values }, Value::Enum(number)) => {
            let name = values
                .values()
                .iter()
                .rev()
                .find(|value| value.value == *number)
                .map(|value| value.name());
            match (name, integer) {
                (Some(name), false) => enc.enumerated(name, *number),
                (Some(name), true) => enc.named_int(name, *number),
                (None, true) => enc.int(*number),
                (None, false) => {
                    unreachable!("an ENUMERATED value that the reader read has a name")
                }
            }
        }
        (Kind::Boolean, Value::Bool(value)) => enc.boolean(*value),
        (Kind::Null, Value::Null) => enc.null(),
        (Kind::Real, Value::Real(value)) => enc.real(*value),
        (Kind::String(StringKind::Visible), Value::Str(text)) => enc.visible_str(text),
        (Kind::String(StringKind::Store), Value::Str(text)) => enc.string_store(text),
        (Kind::String(StringKind::Utf8), Value::Str(text)) => enc.utf8_str(text),
        (Kind::Octets, Value::Octets(bytes)) => enc.octets(bytes),
        _ => unreachable!("a value of the shape of its type"),
    }
}
