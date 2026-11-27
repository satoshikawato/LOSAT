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
    let ds = read[0].dense_seg(false).unwrap();
    assert_eq!(ds.lens, vec![5, 2, 3]);
    assert_eq!(ds.strands, [Strand::Plus, Strand::Minus]);
    assert!(!ds.first_minus);
    // Each row's extent, gaps skipped; the minus strand's from its last segment.
    assert_eq!(read[0].seq_range(&ds, 0), (0, 9));
    assert_eq!(read[0].seq_range(&ds, 1), (10, 19));
}

/// A Std-seg Seq-align of the segments `segs` (each `{ dim 2, loc { ... } }`).
fn std_align(segs: &str) -> Value {
    value(&format!(
        "Seq-align ::= {{ type partial, dim 2, segs std {{ {segs} }} }}"
    ))
}

#[test]
fn a_std_seg_becomes_a_dense_seg_with_the_widths_of_its_rows() {
    // BLASTX: the query's 9 letters on the subject's 3 residues, a gap of the subject, then 6
    // letters on 2 residues.
    let align = std_align(
        "{ dim 2, loc { int { from 0, to 8, strand plus, id local str \"q1\" }, int { from 0, \
         to 2, id local str \"s1\" } } }, { dim 2, loc { int { from 9, to 11, strand plus, id \
         local str \"q1\" }, empty local str \"s1\" } }, { dim 2, loc { int { from 12, to 17, \
         strand plus, id local str \"q1\" }, int { from 3, to 4, id local str \"s1\" } } }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let ds = read[0].dense_seg(false).unwrap();
    assert_eq!(ds.widths, [3, 1]);
    // The gap's segment takes its length from its query row, a row of width 3.
    assert_eq!(ds.lens, vec![3, 1, 2]);
    assert_eq!(ds.starts, vec![0, 0, 9, -1, 12, 3]);
    assert!(ds.translated && !ds.first_minus);
    // `GetSeqRange` of the Std-seg's own rows.
    assert_eq!(read[0].std_range(0).unwrap(), (0, 17));
    assert_eq!(read[0].std_range(1).unwrap(), (0, 4));
    assert_eq!(read[0].row_range(1).unwrap(), (0, 4));
}

/// The line of NCBI's exception that `stop` is, and its message.
fn thrown(stop: ReadStop) -> (u32, String) {
    match stop {
        ReadStop::Stop(Stop::Exception(exception)) => (
            exception.site().line,
            String::from_utf8_lossy(exception.get_msg()).into_owned(),
        ),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_std_seg_that_ncbi_cannot_convert_throws_its_exception() {
    // Another Seq-id in a later segment's row.
    let ids = std_align(
        "{ dim 2, loc { int { from 0, to 8, id local str \"q1\" }, int { from 0, to 2, id local \
         str \"s1\" } } }, { dim 2, loc { int { from 9, to 11, id local str \"q2\" }, int { \
         from 3, to 3, id local str \"s1\" } } }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &ids)], false);
    let (line, message) = thrown(read[0].dense_seg(false).unwrap_err());
    assert_eq!(line, 861);
    assert!(message.starts_with("CreateDensegFromStdseg(): Seq-ids: lcl|q1 and lcl|q2 are not"));
    assert_eq!(
        thrown(read[0].std_range(0).unwrap_err()).1,
        "CSeq_align::GetSeqRange(): Row seqids not consistent. Cannot determine range."
    );
    // Rows whose lengths are neither equal nor in the ratio of a codon.
    let lengths = std_align(
        "{ dim 2, loc { int { from 0, to 8, id local str \"q1\" }, int { from 0, to 4, id local \
         str \"s1\" } } }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &lengths)], false);
    assert_eq!(thrown(read[0].dense_seg(false).unwrap_err()).0, 801);
    // Strands that change in a row.
    let strands = std_align(
        "{ dim 2, loc { int { from 0, to 8, strand plus, id local str \"q1\" }, int { from 0, \
         to 2, id local str \"s1\" } } }, { dim 2, loc { int { from 9, to 17, strand minus, id \
         local str \"q1\" }, int { from 3, to 5, id local str \"s1\" } } }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &strands)], false);
    assert_eq!(thrown(read[0].dense_seg(false).unwrap_err()).0, 827);
    // A row of no letters: `GetSeqRange` throws.
    let empty = std_align(
        "{ dim 2, loc { empty local str \"q1\", int { from 0, to 2, id local str \"s1\" } } }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &empty)], false);
    assert_eq!(
        thrown(read[0].std_range(0).unwrap_err()).1,
        "CSeq_align::GetSeqRange(): Row is empty"
    );
}

#[test]
fn a_gi_score_name_converts_as_ncbis_replace_and_int8() {
    assert_eq!(gi_of_name(b"gi:5"), Some(5));
    assert_eq!(gi_of_name(b"gi:+5"), Some(5));
    assert_eq!(gi_of_name(b"gi:-0012"), Some(-12));
    // Every `gi:` goes, also one that a removal makes.
    assert_eq!(gi_of_name(b"gi:1gi:2"), Some(12));
    assert_eq!(gi_of_name(b"gi:gi:7"), Some(7));
    assert_eq!(gi_of_name(b"gi:9223372036854775807"), Some(i64::MAX));
    assert_eq!(gi_of_name(b"gi:-9223372036854775808"), Some(i64::MIN));
    assert_eq!(gi_of_name(b"gi:9223372036854775808"), None);
    assert_eq!(gi_of_name(b"gi:"), None);
    assert_eq!(gi_of_name(b"gi: 5"), None);
    assert_eq!(gi_of_name(b"gi:-"), None);
    assert_eq!(gi_of_name(b"gi:+-5"), None);
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:644-681
// `CDense_seg::Reverse` flips each strand in place without reversing the strands vector, so the
// first segment's strands remain the pair flipped when `first_minus` is set.
#[test]
fn a_reversed_dense_seg_keeps_its_first_segments_strands_for_flipping() {
    let text =
        "Seq-align ::= { type partial, dim 2, segs denseg { dim 2, numseg 2, ids { local str \
        \"q1\", local str \"s1\" }, starts { 0, 100, 60, 160 }, lens { 60, 60 }, strands { minus, \
        plus, plus, minus } } }";
    let align = value(text);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let ds = read[0].dense_seg(false).unwrap();
    assert!(ds.first_minus);
    assert_eq!(ds.strands, [Strand::Minus, Strand::Plus]);
}
