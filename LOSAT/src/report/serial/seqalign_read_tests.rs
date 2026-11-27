//! The reader of the Seq-aligns of the tabular formats (`seqalign_read`): the ungapped split,
//! the scores, the Dense-seg of a pair and its ranges.

use super::*;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{read_document, Value};

fn value(text: &str) -> Value {
    let (value, posts) = read_document(text.as_bytes(), named::SEQ_ALIGN, Names::TEST);
    assert!(posts.is_empty(), "{posts:?}");
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

const DENSEG: &str = "Seq-align ::= { type partial, dim 2, score { { id str \"score\", value int 61 }, \
    { id str \"e_value\", value real { 15, 10, -20 } }, { id str \"bit_score\", value real { 1105, 10, -1 } }, \
    { id str \"num_ident\", value int 30 } }, segs denseg { dim 2, numseg 3, ids { local str \"q1\", \
    local str \"s1\" }, starts { 0, 15, 5, -1, 7, 10 }, lens { 5, 2, 3 }, strands { plus, minus, plus, minus, \
    plus, minus } } }";

const DENDIAG: &str = "Seq-align ::= { type diags, dim 2, segs dendiag { { dim 2, ids { local str \"q1\", \
    local str \"s1\" }, starts { 0, 4 }, len 3, scores { { id str \"score\", value int 9 } } }, { dim 2, ids { \
    local str \"q1\", local str \"s1\" }, starts { 10, 20 }, len 4, scores { { id str \"score\", value int 11 } } } } }";

#[test]
fn an_ungapped_search_splits_a_seq_align_of_several_diags() {
    let align = value(DENDIAG);
    let set = [V::new(named::SEQ_ALIGN, &align)];
    assert_eq!(read_set(&set, false).len(), 1);
    let split = read_set(&set, true);
    assert_eq!(split.len(), 2);
    assert_eq!(split[1].scores().unwrap().score, 11);
    let one = value(DENSEG);
    assert_eq!(read_set(&[V::new(named::SEQ_ALIGN, &one)], true).len(), 1);
}

#[test]
fn the_scores_and_the_dense_seg_are_read() {
    let align = value(DENSEG);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let scores = read[0].scores().unwrap();
    assert_eq!((scores.score, scores.num_ident), (61, 30));
    assert_eq!(read[0].named_int("num_ident").unwrap(), Some(30));
    let ds = read[0].dense_seg().unwrap();
    assert_eq!(ds.lens, vec![5, 2, 3]);
    assert_eq!(ds.strands, [Strand::Plus, Strand::Minus]);
    assert!(!ds.first_minus);
    // Each row's extent, gaps skipped; the minus strand's from its last segment.
    assert_eq!(read[0].seq_range(&ds, 0), (0, 9));
    assert_eq!(read[0].seq_range(&ds, 1), (10, 19));
}

#[test]
fn a_std_seg_is_refused_by_name() {
    let text = "Seq-align ::= { type partial, dim 2, segs std { { dim 2, ids { local str \"q1\", \
        local str \"s1\" }, loc { int { from 0, to 2, id local str \"q1\" }, int { from 0, to 8, \
        id local str \"s1\" } } } } }";
    let align = value(text);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    assert!(matches!(read[0].dense_seg(), Err(ReadStop::Unsupported(_))));
}

#[test]
fn a_reversed_dense_seg_takes_its_last_segments_strands() {
    let text =
        "Seq-align ::= { type partial, dim 2, segs denseg { dim 2, numseg 2, ids { local str \
        \"q1\", local str \"s1\" }, starts { 0, 100, 60, 160 }, lens { 60, 60 }, strands { minus, \
        plus, plus, minus } } }";
    let align = value(text);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let ds = read[0].dense_seg().unwrap();
    assert!(ds.first_minus);
    assert_eq!(ds.strands, [Strand::Minus, Strand::Plus]);
    assert_eq!(ds.pair_strands(), [Strand::Plus, Strand::Minus]);
}

/// A Dense-seg Seq-align of one segment whose scores hold the `gi:` score named `name`.
fn with_gi_score(name: &str) -> Value {
    value(&format!(
        "Seq-align ::= {{ type partial, dim 2, score {{ {{ id str \"e_value\", value real {{ 1, 10, \
         -5 }} }}, {{ id str \"{name}\", value int 0 }} }}, segs denseg {{ dim 2, numseg 1, ids {{ \
         local str \"q1\", local str \"s1\" }}, starts {{ 0, 0 }}, lens {{ 5 }} }} }}"
    ))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:262-266
// Frozen tests/blast_formatter/ncbi/of6_gi_ok (`gi:+5`, `gi:1gi:2`, `gi:-0012` and both ends of
// the Int8 range read) and of6_gi_overflow (an exception of `NStr::StringToInt8`).
#[test]
fn a_gi_score_name_is_converted_as_ncbi_converts_it() {
    for name in ["gi:+5", "gi:1gi:2", "gi:gi:7", "gi:-9223372036854775808"] {
        let align = with_gi_score(name);
        let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
        assert!(read[0].scores().is_ok(), "{name}");
    }
    let align = with_gi_score("gi:9223372036854775808");
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    match read[0].scores() {
        Err(ReadStop::Stop(Stop::Exception(exception))) => assert_eq!(
            exception.get_msg(),
            b"Cannot convert string '9223372036854775808' to Int8, overflow"
        ),
        other => panic!("{other:?}"),
    }
}
