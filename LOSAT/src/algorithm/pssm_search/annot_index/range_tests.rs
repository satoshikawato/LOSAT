//! Tests of the range model (`range`): `COpenRange`, `CHandleRange` and `CHandleRangeMap` as the
//! object manager's index reads them (line RG D3-1 `om/SPEC.md` §5.8).

use super::super::view::V;
use super::{Handle, HandleRange, HandleRangeMap, Mode, Range, MAX, UNKNOWN};
use crate::report::serial::asn_types::{read_document, table::named, Value};
use crate::seqid::model::ObjectId;
use crate::seqid::ParsedSeqId;

const PLUS: i64 = 1;
const MINUS: i64 = 2;

#[test]
fn ranges_open_at_to_plus_one() {
    assert_eq!(
        Range::new(3, 5),
        Range {
            from: 3,
            to_open: 6
        }
    );
    assert!(!Range::new(3, 3).is_empty());
    assert!(Range::new(5, 3).is_empty());
    // `to` 4294967295 wraps to an open end of 0.
    assert!(Range::new(0, MAX).is_empty());
    assert!(!Range::new(0, MAX - 1).is_empty());
    let mut range = Range::EMPTY;
    range.combine(Range::new(5, 3));
    assert_eq!(range, Range::EMPTY);
    range.combine(Range::new(4, 6));
    range.combine(Range::new(1, 2));
    assert_eq!(
        range,
        Range {
            from: 1,
            to_open: 7
        }
    );
}

// NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:126-206, 338-376 (quoted at
// `HandleRange::add` and `HandleRange::overlapping`)
#[test]
fn a_handle_range_becomes_circular() {
    // Minus strand, a range after (in position) the one before: circular, with an empty second
    // part (f31): the start covers the first range, the end is empty.
    let mut hr = HandleRange::default();
    hr.add(Range::new(0, 5), MINUS, false);
    hr.add(Range::new(10, 9), MINUS, false);
    assert!(hr.circular && hr.has_gaps());
    assert_eq!(hr.overlapping(), Range::WHOLE);
    assert_eq!(
        hr.circular_start(),
        Range {
            from: 0,
            to_open: 6
        }
    );
    assert!(hr.circular_end().is_empty());
    // Plus strand: the end takes the origin (f32).
    let mut hr = HandleRange::default();
    hr.add(Range::new(10, 15), PLUS, false);
    hr.add(Range::new(1, 0), PLUS, false);
    assert!(hr.circular);
    assert!(!hr.circular_start().is_empty() && !hr.circular_end().is_empty());
    // `circular RNA` keeps it linear (f34), and so does a change of strand.
    let mut hr = HandleRange::default();
    hr.add(Range::new(0, 5), MINUS, true);
    hr.add(Range::new(10, 9), MINUS, true);
    assert!(!hr.circular);
    assert_eq!(
        hr.overlapping(),
        Range {
            from: 0,
            to_open: 6
        }
    );
    let mut hr = HandleRange::default();
    hr.add(Range::new(10, 15), PLUS, false);
    hr.add(Range::new(1, 4), PLUS, false);
    hr.add(Range::new(20, 30), UNKNOWN, false);
    assert!(!hr.circular);
    assert_eq!(
        hr.overlapping(),
        Range {
            from: 1,
            to_open: 31
        }
    );
}

fn location(text: &str) -> Value {
    let (value, posts) = read_document(format!("Seq-loc ::= {text}").as_bytes(), named::SEQ_LOC);
    assert!(posts.is_empty());
    value.unwrap()
}

fn map(text: &str, mode: Mode) -> HandleRangeMap {
    let value = location(text);
    HandleRangeMap::of_location(V::new(named::SEQ_LOC, &value), mode)
}

// NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:91-205, 249-281 (quoted at
// `HandleRangeMap::add_location` and `HandleRangeMap::add`)
#[test]
fn a_map_marks_the_ranges_between_seq_ids() {
    // f15: the reversed range of local id 1 after local id 2 gets `more_before`: from 0.
    let m = map(
        "mix { int { from 0, to 5, id local id 2 }, int { from 10, to 5, id local id 1 } }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 2);
    assert_eq!(m.map[1].0, Handle::LocalId(1));
    assert!(!m.map[1].1.overlapping().is_empty());
    // f38: `trans-splicing` sets no such marks.
    let m = map(
        "mix { int { from 0, to 5, id local id 2 }, int { from 10, to 5, id local id 1 } }",
        Mode::TransSplicing,
    );
    assert!(m.map[1].1.overlapping().is_empty());
    // f14: the first Seq-id's reversed range before another one gets `more_after`, which moves
    // the open end of an empty total: still empty.
    let m = map(
        "mix { int { from 10, to 5, id local id 1 }, int { from 0, to 5, id local id 2 } }",
        Mode::None,
    );
    assert!(m.map[0].1.overlapping().is_empty());
    // An `equiv`'s later parts do not follow the Seq-id before them.
    let m = map(
        "equiv { int { from 0, to 5, id local id 2 }, int { from 10, to 5, id local id 1 } }",
        Mode::None,
    );
    assert!(m.map[1].1.overlapping().is_empty());
    // A local string is NCBI's case-insensitive key; two other Seq-ids are ambiguous.
    let m = map(
        "mix { int { from 0, to 5, id local str \"Q\" }, int { from 10, to 12, id local str \"q\" } }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 1);
    assert!(!m.ambiguous());
    assert!(map(
        "mix { int { from 0, to 5, id genbank { accession \"A\" } }, int { from 0, to 5, id embl { accession \"B\" } } }",
        Mode::None
    )
    .ambiguous());
    assert_eq!(
        Handle::of(&ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))),
        Handle::LocalStr(b"query_1".to_vec())
    );
    // `packed-pnt` (f18): the point 4294967295 is empty; `empty` is empty; `whole` is whole.
    assert!(map(
        "packed-pnt { id local id 1, points { 4294967295 } }",
        Mode::None
    )
    .map[0]
        .1
        .overlapping()
        .is_empty());
    assert!(map("empty local id 1", Mode::None).map[0]
        .1
        .overlapping()
        .is_empty());
    assert_eq!(
        map("whole local id 1", Mode::None).map[0].1.overlapping(),
        Range::WHOLE
    );
    assert!(map("feat local id 1", Mode::None).map.is_empty());
}
