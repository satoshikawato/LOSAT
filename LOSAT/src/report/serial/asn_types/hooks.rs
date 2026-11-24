//! What NCBI's hand-written classes do when the reader has read one of them (their
//! `PostRead`), for the classes of the table whose `PostRead` changes the object or stops the
//! reader (SURVEY_D3 §2.4 of line RG series F: `CSeq_descr`, `CVariation_ref`; the others keep
//! caches only). (The Seq-id is read into LOSAT's Seq-id model by its own reader, `super::seq_id`.)

use super::read::{Frame, Post};
use super::table::named;
use super::{member_index, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_input::{Fail, R, VECTOR_RESERVE_LENGTH_ERROR};

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:1497-1527
// NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_seg.hpp:181-182 (the hook is global:
// `NCBISER_HAVE_GLOBAL_READ_MEMBER_HOOK` for `CDense_seg` and "starts,lens,strands")
// `numseg` and `dim` are `int` (Dense_seg_.hpp:92-93), so a negative one is a huge `size_t`, and the
// product wraps. `reserve` throws a `std::length_error` ("vector::reserve") above `max_size()`:
// PTRDIFF_MAX / 4 = 2^61 - 1 for the 4-byte elements of `starts` and `lens` (NCBI 2.17.0 probes
// p5: dim 2^30 x numseg 2^31-1 reads on, dim 2^30+1 stops), PTRDIFF_MAX = 2^63 - 1 for the 1-byte
// `ENa_strand`s of `strands` (Na_strand_.hpp:56, `NCBI_PACKED_ENUM_TYPE` of `unsigned char`; the
// mandatory `starts`, reserved first with the same product, stops a larger one before it), which
// no frame catches (they catch
// `CException`s): `CATCH_ALL`'s `std::exception` clause reports it, exit status 255. A smaller
// request that cannot be allocated is a `bad_alloc`, which the hook ignores, as LOSAT does not
// allocate. (`OBJECTS_DENSE_SEG_RESERVE` is NCBI's default, true: LOSAT reads no NCBI parameter.)
// The reserve hook of Seq-graph (`numval` is a TSeqPos: never above `max_size()`) needs nothing
// here; SeqTable-multi-data's is a variant hook (`pre_read_variant`).
/// The hook of NCBI's class `ty` before the value of its member at position `member` is read,
/// with the members `read` so far.
pub(super) fn pre_read_member(ty: u16, member: usize, read: &[(u16, Value)]) -> R<()> {
    if ty != named::DENSE_SEG {
        return Ok(());
    }
    let get = |name: &str| {
        let index = member_index(named::DENSE_SEG, name);
        read.iter()
            .find(|(position, _)| usize::from(*position) == index)
            .map(|(_, value)| int(value))
    };
    // `GetDim()`: the DEFAULT 2 when the member was not read.
    let dim = get("dim").unwrap_or(2) as i32 as i64 as u64;
    let numseg = get("numseg").unwrap_or(0) as i32 as i64 as u64;
    // `max_size()` of a vector of 4-byte and of 1-byte elements.
    const MAX_SIZE_4: u64 = (i64::MAX as u64) / 4;
    const MAX_SIZE_1: u64 = i64::MAX as u64;
    let (size, max_size) = match super::members(ty)[member].name() {
        "starts" => (dim.wrapping_mul(numseg), MAX_SIZE_4),
        "lens" => (numseg, MAX_SIZE_4),
        "strands" => (dim.wrapping_mul(numseg), MAX_SIZE_1),
        _ => return Ok(()),
    };
    if size > max_size {
        return Err(Fail::Std(&VECTOR_RESERVE_LENGTH_ERROR));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:703-765
// NCBI reference (598d8ae6): c++/include/objects/seqtable/SeqTable_multi_data.hpp:231-232 (the hook
// is global: `NCBISER_HAVE_GLOBAL_READ_VARIANT_HOOK` for `CSeqTable_multi_data` and "*")
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1906-1923
// The hook runs when the variant has been named, its frame on top of the stack: the multi-data's
// own choice frame is depth 0, a column's `data` finds its Seq-table at depth 2, an `int-delta`
// nested 1 to 3 times at 3 to 5, 4 times not (no reserve), an `int-scaled`'s data at 4 (NCBI
// BLAST+ 2.17.0, line RG D3-1 `om/` cases t62-t65). `num-rows` is an `int` (Seq_table_.hpp): a
// negative one is a `size_t` above every `max_size()`, and `reserve` throws
// a `std::length_error` ("vector::reserve"), which no frame catches (`CATCH_ALL`'s
// `std::exception` clause, exit status 255: t50, t66, t69); `bit`'s byte count (the sum of `size`
// and 7, divided by 8) wraps below that size, and any other size is at most 2^31-1 elements, whose
// `bad_alloc` the hook ignores (LOSAT allocates nothing). (`OBJECTS_SEQ_TABLE_RESERVE` is NCBI's
// default, true: LOSAT reads no NCBI parameter.)
/// The hook of the CHOICE `ty` before the value of its variant at position `variant` is read,
/// with the reader's frames (the variant's on top).
pub(super) fn pre_read_variant(ty: u16, variant: usize, frames: &[Frame]) -> R<()> {
    if ty != named::SEQTABLE_MULTI_DATA {
        return Ok(());
    }
    let reserved = matches!(
        super::members(ty)[variant].name(),
        "int"
            | "int1"
            | "int2"
            | "int8"
            | "real"
            | "string"
            | "bytes"
            | "common-string"
            | "common-bytes"
            | "loc"
            | "id"
            | "interval"
    );
    if !reserved {
        return Ok(());
    }
    let mut depth = 0;
    for frame in frames.iter().rev().skip(1) {
        if !frame.is_class_or_choice() {
            continue;
        }
        if depth >= 2 && frame.ty() == Some(named::SEQ_TABLE) {
            if frame.num_rows().is_some_and(|rows| rows < 0) {
                return Err(Fail::Std(&VECTOR_RESERVE_LENGTH_ERROR));
            }
            return Ok(());
        }
        depth += 1;
        if depth > 5 {
            break;
        }
    }
    Ok(())
}

/// The hook of the type `ty` after its value was read (not after a skip).
pub(super) fn post_read(ty: u16, value: Value, posts: &mut Vec<Post>) -> R<Value> {
    match ty {
        named::SEQ_DESCR => seq_descr_post_read(value),
        named::VARIATION_REF => Ok(variation_ref_post_read(value, posts)),
        _ => Ok(value),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
// (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless NCBI's environment or configuration file sets
// it, which LOSAT does not read; NCBI 2.17.0 gives the error with the variable set too, c23 of
// SURVEY_D3.)
fn seq_descr_post_read(value: Value) -> R<Value> {
    match &value {
        Value::List(items) if items.is_empty() => Err(Fail::serial(NcbiException::new(
            SEQ_DESCR_POST_READ,
            "empty Seq-descr is not allowed",
        ))),
        _ => Ok(value),
    }
}

/// `NCBI_THROW` of `CSeq_descr::PostRead` (its message has no stream position; NCBI's callers of
/// the PSSM reader print their own line for a `CSerialException`, so the site is never printed).
const SEQ_DESCR_POST_READ: ThrowSite = ThrowSite {
    file: "src/objects/seq/Seq_descr.cpp",
    line: 65,
    kind: "CSerialException::eInvalidData",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Variation_ref.cpp:56-151
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/VariantProperties.cpp:50-55
// `SetVariant_prop()` makes a `VariantProperties` (version 5) when the member is not set. The
// messages are posted when the class has been read, before the members that follow it.
fn variation_ref_post_read(mut value: Value, posts: &mut Vec<Post>) -> Value {
    let member = |name: &str| member_index(named::VARIATION_REF, name);
    let drop = |value: &mut Value, posts: &mut Vec<Post>, name: &str, message: &'static str| {
        if value.take_member(member(name)).is_some() {
            posts.push(Post::Error(message));
        }
    };
    let moved = |value: &mut Value,
                 posts: &mut Vec<Post>,
                 name: &str,
                 property: &str,
                 message: &'static str| {
        if let Some(moved) = value.take_member(member(name)) {
            let prop = variant_prop(value);
            let index = member_index(named::VARIANTPROPERTIES, property);
            if prop.member(index).is_some() {
                posts.push(Post::Error(message));
            } else {
                prop.set_member(index, moved);
            }
        }
    };
    drop(
        &mut value,
        posts,
        "population-data",
        "Variation-ref.population-data is deprecated and will be ignored",
    );
    moved(
        &mut value,
        posts,
        "validated",
        "other-validation",
        "Both Variation-ref.validated and Variation-ref.variant-properties.other-validation are set; ignoring Variation-ref.validated",
    );
    drop(
        &mut value,
        posts,
        "clinical-test",
        "Variation-ref.clinical-test is deprecated and will be ignored",
    );
    moved(
        &mut value,
        posts,
        "allele-origin",
        "allele-origin",
        "Both Variation-ref.allele-origin and Variation-ref.variant-properties.allele-origin are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "allele-state",
        "allele-state",
        "Both Variation-ref.allele-state and Variation-ref.variant-properties.allele-state are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "allele-frequency",
        "allele-frequency",
        "Both Variation-ref.allele-frequency and Variation-ref.variant-properties.allele-frequency are set; ignoring Variation-ref.validated",
    );
    moved(
        &mut value,
        posts,
        "is-ancestral-allele",
        "is-ancestral-allele",
        "Both Variation-ref.is-ancestral-allele and Variation-ref.variant-properties.is-ancestral-allele are set; ignoring Variation-ref.validated",
    );
    drop(
        &mut value,
        posts,
        "pub",
        "Variation-ref.pub is deprecated and will be ignored",
    );
    value
}

/// `SetVariant_prop()`: the member, made (version 5) when it is not set.
fn variant_prop(value: &mut Value) -> &mut Value {
    let index = member_index(named::VARIATION_REF, "variant-prop");
    if value.member(index).is_none() {
        let version = member_index(named::VARIANTPROPERTIES, "version");
        value.set_member(
            index,
            Value::Sequence(vec![(version as u16, Value::Int(5))]),
        );
    }
    match value {
        Value::Sequence(members) => {
            &mut members
                .iter_mut()
                .find(|(position, _)| usize::from(*position) == index)
                .expect("the member was set")
                .1
        }
        _ => unreachable!("a Variation-ref is a SEQUENCE"),
    }
}

fn int(value: &Value) -> i64 {
    match value {
        Value::Int(n) => *n,
        other => unreachable!("an INTEGER of the table, not {other:?}"),
    }
}
