//! Tests of the Seq-table access (`table`): the values of a column at a row as NCBI's
//! `CSeqTable_column`, `CSeqTable_multi_data`, `CSeqTable_single_data` and
//! `CSeqTable_sparse_index` give them (line RG D3-1 `om/SPEC.md` §7).

use super::super::view::V;
use super::super::{Indexer, Stop};
use super::{Caches, Column, Want};
use crate::report::serial::asn_types::{read_document, table::named, Value};

fn column_value(text: &str) -> Value {
    let (value, posts) = read_document(
        format!("SeqTable-column ::= {{ header {{ field-id 3 }}, {text} }}").as_bytes(),
        named::SEQTABLE_COLUMN,
        crate::report::serial::asn_types::PSIBLAST_STACK,
    );
    assert!(posts.is_empty());
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

fn indexer() -> Indexer<'static> {
    Indexer {
        program: "psiblast",
        posts: crate::report::diag_log::Posts::new(),
        muted: false,
        caches: Caches::default(),
        handles: Default::default(),
    }
}

/// `TryGetValue(row, v)` of the column `text` at each of `rows`, as text (`-` for none, `E:` and
/// the message for an exception, `R` for a refusal).
fn values(text: &str, want: Want, rows: std::ops::Range<u64>) -> Vec<String> {
    let value = column_value(text);
    let column = Column(V::new(named::SEQTABLE_COLUMN, &value));
    let mut ix = indexer();
    rows.map(|row| match ix.column_try_get(column, row, want) {
        Ok(Some(v)) => v.to_string(),
        Ok(None) => "-".to_string(),
        Err(Stop::Exception(e)) => format!("E:{}", e.get_msg()),
        Err(Stop::Refuse(r)) => format!("R:{r:?}"),
        Err(other) => format!("{other:?}"),
    })
    .collect()
}

#[test]
fn data_then_default() {
    assert_eq!(
        values("data int { 4, -2 }", Want::Int4, 0..3),
        ["4", "-2", "-"]
    );
    assert_eq!(
        values("data int { 4 }, default int 7", Want::Int4, 0..2),
        ["4", "7"]
    );
    // `bit`: the most significant bit is row 0; rows past the bytes have no value.
    assert_eq!(
        values("data bit 'A0'H", Want::Bool, 0..9),
        ["1", "0", "1", "0", "0", "0", "0", "0", "-"]
    );
    // `int1` is signed bytes.
    assert_eq!(values("data int1 'FF7F'H", Want::Int4, 0..2), ["-1", "127"]);
    // A bool is 0 or 1; another value is too big (E14).
    assert_eq!(
        values("data int { 1, 2 }", Want::Bool, 0..2),
        [
            "1",
            "E:CSeqTable_multi_data::TryGetBool(): value is too big for requested type: 2"
        ]
    );
    // A type without integers (E13), a single value of another type (E15, E16).
    assert_eq!(
        values("data string { \"a\" }", Want::Int8, 0..1),
        ["E:CSeqTable_multi_data::TryGetInt8(): value cannot be converted to requested type"]
    );
    assert_eq!(
        values("default real { 1, 10, 0 }", Want::Int4, 0..1),
        ["E:CSeqTable_single_data::GetValue(int&): real value cannot be converted to requested type"]
    );
    assert_eq!(
        values("default int8 4294967296", Want::Int4, 0..1),
        ["E:CSeqTable_single_data::GetValue(int&): value is too big for requested type: 4294967296"]
    );
    assert_eq!(
        values("default int8 4294967296", Want::Int8, 0..1),
        ["4294967296"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_sparse_index.cpp:533-581 (quoted at
// `Indexer::sparse_index_at`)
#[test]
fn sparse_indexes() {
    // `indexes`: libstdc++'s `lower_bound`, whose halving on the unsorted { 5, 3, 9 } finds 9
    // but neither 3 nor 5.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes { 1, 4, 6 }",
            Want::Int4,
            0..7
        ),
        ["-", "10", "-", "-", "11", "-", "12"]
    );
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes { 5, 3, 9 }",
            Want::Int4,
            3..10
        ),
        ["-", "-", "-", "-", "-", "-", "12"]
    );
    // `indexes-delta`: the rows 2, 3 and 7.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes-delta { 2, 1, 4 }",
            Want::Int4,
            0..8
        ),
        ["-", "-", "10", "11", "-", "-", "-", "12"]
    );
    // `bit-set`: the set bits before the row count the data's index; a skipped row takes
    // `sparse-other`.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse bit-set '4180'H, sparse-other int 0",
            Want::Int4,
            0..10
        ),
        ["0", "10", "0", "0", "0", "0", "0", "11", "12", "0"]
    );
    assert_eq!(
        values(
            "data int { 1 }, sparse bit-set-bvector { size 8, data '00'H }",
            Want::Int4,
            0..1
        ),
        ["R:BitVector"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:138-169 (quoted at
// `Indexer::delta_sum`)
#[test]
fn int_delta_sums_read_whole_blocks() {
    assert_eq!(
        values("data int-delta int { 1, 2, 3 }", Want::Int4, 0..4),
        ["1", "3", "6", "-"]
    );
    // Elements without integers (E13, named after the accumulator's type).
    assert_eq!(
        values("data int-delta string { \"1\", \"2\" }", Want::Int4, 0..1),
        ["E:CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type"]
    );
    // Row 0 reads every element of its block of 128: the second element, too big for the `Int4`
    // accumulator of an `int-scaled` without bounds, throws for row 0.
    assert_eq!(
        values(
            "data int-delta int-scaled { mul 1, add 0, data int8 { 1, 4294967296 } }",
            Want::Int4,
            0..1
        ),
        ["E:CSeqTable_multi_data::TryGetInt4(): value is too big for requested type: 4294967296"]
    );
    // The `Int8` sums of `int8` elements, and an `Int4` sum that overflows (refused).
    assert_eq!(
        values("data int-delta int8 { 4294967296, 1 }", Want::Int8, 0..2),
        ["4294967296", "4294967297"]
    );
    assert_eq!(
        values("data int-delta int { 2147483647, 1 }", Want::Int8, 0..1),
        ["R:DeltaOverflow"]
    );
    // An element of a later block is not read for a row of the first.
    assert_eq!(
        values(
            &format!(
                "data int-delta int {{ {} }}",
                (0..129)
                    .map(|i| if i == 128 { "2147483647" } else { "1" })
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            Want::Int8,
            0..1
        ),
        ["1"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:77-84 (quoted at
// `Indexer::x_try_get_int8`)
#[test]
fn int_scaled_values() {
    assert_eq!(
        values(
            "data int-scaled { mul 3, add -1, data int { 1, 2 } }",
            Want::Int4,
            0..3
        ),
        ["2", "5", "-"]
    );
    assert_eq!(
        values(
            "data int-scaled { mul 2, add 0, data int8 { 9223372036854775807 } }",
            Want::Int8,
            0..1
        ),
        ["R:ScaledOverflow"]
    );
}

/// The multi-data inside multi-data are read from a list on the heap: 50,000 nested `int-scaled`
/// (each adding 1) and 2,000 nested `int-delta` values give their rows on a thread of 256 KiB.
#[test]
fn nested_multi_data_are_read_without_recursion() {
    std::thread::Builder::new()
        .stack_size(256 << 10)
        .spawn(|| {
            let read = |text: String| {
                let (value, posts) = read_document(
                    text.as_bytes(),
                    named::SEQTABLE_COLUMN,
                    crate::report::serial::asn_types::NcbiStack::with_room(false, i64::MAX),
                );
                assert!(posts.is_empty());
                value.expect("read")
            };
            let wrap = |n: usize, open: &str, inner: &str, close: &str| {
                format!("{}{inner}{}", open.repeat(n), close.repeat(n))
            };
            let rows = |value: &Value| {
                let column = Column(V::new(named::SEQTABLE_COLUMN, value));
                let mut ix = indexer();
                (0..3)
                    .map(|row| match ix.column_try_get(column, row, Want::Int8) {
                        Ok(v) => format!("{v:?}"),
                        Err(other) => format!("{other:?}"),
                    })
                    .collect::<Vec<_>>()
            };
            let scaled = read(format!(
                "SeqTable-column ::= {{ header {{ field-id 3 }}, data {} }}",
                wrap(
                    50_000,
                    "int-scaled { mul 1, add 1, data ",
                    "int { 4, 9 }",
                    " }"
                )
            ));
            assert_eq!(rows(&scaled), ["Some(50004)", "Some(50009)", "None"]);
            // `int-delta` of `int-delta`s: the sums of the sums.
            let delta = read(format!(
                "SeqTable-column ::= {{ header {{ field-id 3 }}, data {} }}",
                wrap(2_000, "int-delta ", "int { 4, 9 }", "")
            ));
            assert_eq!(rows(&delta), ["Some(4)", "Some(8009)", "None"]);
        })
        .expect("a thread")
        .join()
        .expect("no overflow");
}
