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
