//! The producer of the tabular rows of Seq-aligns (`tabular_seqalign`): a row's coordinates,
//! counts and letters, and the coverages of a subject's group.

use super::*;
use crate::blastinput::app::FormatChoice;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{read_document, Value};
use crate::report::serial::seqalign_read::read_set;

fn value(text: &str) -> Value {
    let (value, posts) = read_document(text.as_bytes(), named::SEQ_ALIGN, Names::TEST);
    assert!(posts.is_empty(), "{posts:?}");
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

fn spec(fields: &str) -> TabularSpec {
    TabularSpec::of_formatter(&FormatChoice {
        number: 6,
        spec: fields.to_string(),
        delimiter: String::new(),
    })
    .expect("a tabular format")
}

/// The row of `align`, as the driver makes it: its scores, its plan, then its letters.
fn hsp<'r>(
    align: &ReadAlign<'_>,
    ctx: &RowContext<'_>,
    sequences: &dyn Sequences,
    buffers: &'r mut RowBuffers,
) -> Result<TabularHsp<'r>, ReadStop> {
    let scores = row_scores(align, ctx.spec)?;
    let plan = row_plan(align, scores, ctx)?;
    row_hsp(&plan, ctx, sequences, buffers)
}

/// A gapped Seq-align of `q1` (letters 0..=9) on the minus strand of `s1` (10..=19).
const MINUS: &str = "Seq-align ::= { type partial, dim 2, score { { id str \"score\", value int 8 }, \
    { id str \"e_value\", value real { 2, 10, -3 } }, { id str \"bit_score\", value real { 165, 10, -1 } }, \
    { id str \"num_ident\", value int 7 } }, segs denseg { dim 2, numseg 3, ids { local str \"q1\", \
    local str \"s1\" }, starts { 0, 15, 5, -1, 7, 10 }, lens { 5, 2, 3 }, strands { plus, minus, plus, minus, \
    plus, minus } } }";

#[test]
fn a_row_of_the_minus_strand_counts_its_columns() {
    let align = value(MINUS);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec("qstart qend sstart send length gaps gapopen evalue bitscore score qseq sseq");
    let ctx = RowContext {
        spec: &spec,
        no_fetch: true,
        query_length: 10,
    };
    // The subject's letters of the row's range only (10..=19), as the driver reads them.
    let scores = row_scores(&read[0], &spec).unwrap();
    let plan = row_plan(&read[0], scores, &ctx).unwrap();
    assert_eq!(plan.subject_range(&ctx), Some((10, 19)));
    let letters = RowLetters {
        query: b"ACGTACGTAC",
        subject: b"AAAAAAAAAA",
        subject_from: 10,
    };
    let mut buffers = RowBuffers::default();
    let row = hsp(&read[0], &ctx, &letters, &mut buffers).unwrap();
    assert_eq!((row.q_start, row.q_end), (1, 10));
    // The subject's coordinates of a minus strand run down.
    assert_eq!((row.s_start, row.s_end), (20, 11));
    assert_eq!((row.length, row.gaps, row.gap_opens), (10, 2, 1));
    assert_eq!(row.score, 8);
    let (query, subject) = row.rows.expect("the rows' letters");
    assert_eq!((query.len(), subject.len()), (10, 10));
    assert_eq!(subject.iter().filter(|&&b| b == b'-').count(), 2);
}

#[test]
fn a_row_without_counted_fields_reads_no_letters() {
    let align = value(MINUS);
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec("qseqid sseqid evalue");
    assert!(!reads_letters(&spec, true));
    let ctx = RowContext {
        spec: &spec,
        no_fetch: true,
        query_length: 10,
    };
    let empty = RowLetters {
        query: b"",
        subject: b"",
        subject_from: 0,
    };
    let mut buffers = RowBuffers::default();
    let row = hsp(&read[0], &ctx, &empty, &mut buffers).unwrap();
    assert!(row.rows.is_none());
    assert_eq!(row.evalue, 0.002);
}

#[test]
fn qcovs_is_set_on_the_first_row_of_each_subject() {
    let first = value(MINUS);
    let second = value(&MINUS.replace(
        "starts { 0, 15, 5, -1, 7, 10 }",
        "starts { 2, 15, 7, -1, 9, 10 }",
    ));
    let set = [
        V::new(named::SEQ_ALIGN, &first),
        V::new(named::SEQ_ALIGN, &second),
    ];
    let read = read_set(&set, false);
    let scores = subject_scores(&read, Some(20), &spec("qcovs"), "qcovs", false).unwrap();
    assert_eq!(scores.len(), 2);
    // The union of the query's ranges, 0..=11 of 20 letters: 60 %.
    assert_eq!(scores[0].0, 60);
    assert_eq!(scores[1].0, 60);
}
