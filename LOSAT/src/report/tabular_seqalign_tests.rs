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
        spec: fields.as_bytes().to_vec(),
        delimiter: Vec::new(),
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

/// The genetic codes of a BLASTN formatter (`m_QueryGenCode` 1, `m_DbGenCode` 0).
static BLASTN_CODES: std::sync::LazyLock<GeneticCodes> =
    std::sync::LazyLock::new(|| GeneticCodes::new(1, 0));

/// What a BLASTN row of a query of 10 letters and a subject of 20 reads besides its Seq-align.
fn blastn_context(spec: &TabularSpec) -> RowContext<'_> {
    RowContext {
        spec,
        no_fetch: true,
        query_length: 10,
        subject_length: 20,
        nucleotide: [true, true],
        codes: &BLASTN_CODES,
    }
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
    let ctx = blastn_context(&spec);
    // The subject's letters of the row's range only (10..=19), as the driver reads them.
    let scores = row_scores(&read[0], &spec).unwrap();
    let plan = row_plan(&read[0], scores, &ctx).unwrap();
    assert_eq!(plan.subject_range(&ctx), Some((10, 19)));
    // The rows' strands as `CAlnVec` reads them (the Dense-seg reversed: the query plus).
    let [query_plus, subject_plus] = plan.rows_plus().expect("a Dense-seg");
    let letters = RowLetters {
        query: RowSequence {
            letters: RowSource::Held(b"ACGTACGTAC", 0),
            size: 10,
            minus: !query_plus,
        },
        subject: RowSequence {
            letters: RowSource::Held(b"AAAAAAAAAA", 10),
            size: 20,
            minus: !subject_plus,
        },
        codes: &BLASTN_CODES,
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
    let ctx = blastn_context(&spec);
    let empty = RowLetters {
        query: RowSequence::NONE,
        subject: RowSequence::NONE,
        codes: &BLASTN_CODES,
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
    let scores = subject_scores(&read, Some(20), &spec("qcovs"), b"qcovs", false).unwrap();
    assert_eq!(scores.len(), 2);
    // The union of the query's ranges, 0..=11 of 20 letters: 60 %.
    assert_eq!(scores[0].0, 60);
    assert_eq!(scores[1].0, 60);
}

/// A Std-seg Seq-align of one segment whose query row is `query` and subject row `subject`.
fn std_seg(query: &str, subject: &str) -> Value {
    value(&format!(
        "Seq-align ::= {{ type partial, dim 2, score {{ {{ id str \"score\", value int 20 }}, \
         {{ id str \"e_value\", value real {{ 1, 10, -5 }} }}, {{ id str \"bit_score\", value \
         real {{ 120, 10, -1 }} }}, {{ id str \"num_ident\", value int 3 }} }}, segs std {{ {{ \
         dim 2, loc {{ {query}, {subject} }} }} }} }}"
    ))
}

/// What a row of a translated search reads besides its Seq-align: a query of 12 letters and a
/// subject of `subject_length`, of the molecules `nucleotide`, with `codes`.
fn translated_context<'s>(
    spec: &'s TabularSpec,
    subject_length: i64,
    nucleotide: [bool; 2],
    codes: &'s GeneticCodes,
) -> RowContext<'s> {
    RowContext {
        spec,
        no_fetch: false,
        query_length: 12,
        subject_length,
        nucleotide,
        codes,
    }
}

const FIELDS: &str = "qstart qend sstart send qframe sframe frames length nident qseq sseq";

#[test]
fn a_blastx_std_seg_reads_its_query_translated_in_its_frame() {
    let align = std_seg(
        "int { from 0, to 8, strand plus, id local str \"q1\" }",
        "int { from 0, to 2, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec(FIELDS);
    let codes = GeneticCodes::new(1, 0);
    let ctx = translated_context(&spec, 4, [true, false], &codes);
    let letters = RowLetters {
        query: RowSequence {
            letters: RowSource::Held(b"ATGGCCTTTAAA", 0),
            size: 12,
            minus: false,
        },
        subject: RowSequence {
            letters: RowSource::Held(b"MAFK", 0),
            size: 4,
            minus: false,
        },
        codes: &codes,
    };
    let mut buffers = RowBuffers::default();
    let row = hsp(&read[0], &ctx, &letters, &mut buffers).unwrap();
    assert_eq!(
        (row.q_start, row.q_end, row.s_start, row.s_end),
        (1, 9, 1, 3)
    );
    // The query's frame of its start; a protein row's frame 0.
    assert_eq!((row.q_frame, row.s_frame), (1, 0));
    assert_eq!(row.length, 3);
    let (query, subject) = row.rows.expect("the rows' letters");
    assert_eq!((query, subject), (&b"MAF"[..], &b"MAF"[..]));
}

#[test]
fn a_minus_strand_query_row_runs_down_in_a_negative_frame() {
    let align = std_seg(
        "int { from 3, to 11, strand minus, id local str \"q1\" }",
        "int { from 0, to 2, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec(FIELDS);
    let codes = GeneticCodes::new(1, 0);
    let ctx = translated_context(&spec, 4, [true, false], &codes);
    // Letters 3..=11 are AAAGGCCAT, whose reverse complement is ATGGCCTTT (MAF).
    let letters = RowLetters {
        query: RowSequence {
            letters: RowSource::Held(b"CCCAAAGGCCAT", 0),
            size: 12,
            minus: true,
        },
        subject: RowSequence {
            letters: RowSource::Held(b"MAFK", 0),
            size: 4,
            minus: false,
        },
        codes: &codes,
    };
    let mut buffers = RowBuffers::default();
    let row = hsp(&read[0], &ctx, &letters, &mut buffers).unwrap();
    assert_eq!((row.q_start, row.q_end), (12, 4));
    // `GetFrame` of start 11 on the minus strand of 12 letters: -((12 - 11 - 1) % 3 + 1).
    assert_eq!(row.q_frame, -1);
    assert_eq!(row.rows.expect("the rows' letters").0, b"MAF");
}

#[test]
fn a_query_row_without_letters_ends_the_row_at_the_query_block() {
    let align = std_seg(
        "empty local str \"q1\"",
        "int { from 0, to 2, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec("qseqid sseqid evalue qstart length qcovhsp");
    let codes = GeneticCodes::new(1, 0);
    let ctx = translated_context(&spec, 4, [true, false], &codes);
    let scores = row_scores(&read[0], &spec).unwrap();
    let plan = row_plan(&read[0], scores, &ctx).unwrap();
    // `GetSeqRange(0)` throws "Row is empty", which the query block catches.
    assert!(plan.leaves);
    let mut buffers = RowBuffers::default();
    let empty = RowLetters {
        query: RowSequence::NONE,
        subject: RowSequence::NONE,
        codes: &codes,
    };
    let row = row_hsp(&plan, &ctx, &empty, &mut buffers).unwrap();
    assert_eq!((row.q_start, row.length, row.qcovhsp), (0, 0, -1));
    assert_eq!(row.score, 20);
}

#[test]
fn a_tblastx_std_seg_is_a_dense_seg_of_codons_on_both_rows() {
    let align = std_seg(
        "int { from 0, to 8, strand plus, id local str \"q1\" }",
        "int { from 3, to 11, strand minus, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let ds = read[0].dense_seg(true).unwrap();
    assert_eq!((ds.lens.clone(), ds.widths), (vec![3], [3, 3]));
    assert_eq!(ds.strands, [Strand::Plus, Strand::Minus]);
    // A segment of letters that are no whole codons.
    let uneven = std_seg(
        "int { from 0, to 9, strand plus, id local str \"q1\" }",
        "int { from 0, to 9, strand plus, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &uneven)], false);
    match read[0].dense_seg(true) {
        Err(ReadStop::Stop(Stop::Exception(exception))) => {
            assert_eq!(exception.site().line, 985);
            assert!(exception
                .get_msg()
                .ends_with(b"Length of segment 0 is not divisible by 3."));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_translated_row_whose_genetic_code_ncbi_lacks_throws_before_its_letters() {
    // A TBLASTN row (query protein, subject translated) formatted with the db code 0.
    let align = std_seg(
        "int { from 0, to 2, id local str \"q1\" }",
        "int { from 0, to 8, strand plus, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec("qseq sseq");
    let codes = GeneticCodes::new(1, 0);
    let ctx = translated_context(&spec, 9, [false, true], &codes);
    let scores = row_scores(&read[0], &spec).unwrap();
    let plan = row_plan(&read[0], scores, &ctx).unwrap();
    assert!(plan.translation_stop(&ctx, 0).unwrap().is_none());
    let exception = plan
        .translation_stop(&ctx, 1)
        .unwrap()
        .expect("GetTransTable's exception");
    assert_eq!(exception.get_msg(), b"Unable to find genetic code number 0");
    // With code 1 the row reads.
    let codes = GeneticCodes::new(1, 1);
    let ctx = translated_context(&spec, 9, [false, true], &codes);
    assert!(plan.translation_stop(&ctx, 1).unwrap().is_none());
}

#[test]
fn a_translated_row_past_its_sequence_leaves_a_segment_without_whole_codons() {
    // A BLASTX query row of 9 letters on a query of 7 (its Bioseq shorter than the row).
    let align = std_seg(
        "int { from 0, to 8, strand plus, id local str \"q1\" }",
        "int { from 0, to 2, id local str \"s1\" }",
    );
    let read = read_set(&[V::new(named::SEQ_ALIGN, &align)], false);
    let spec = spec("qseq");
    let codes = GeneticCodes::new(1, 0);
    let mut ctx = translated_context(&spec, 3, [true, false], &codes);
    ctx.query_length = 7;
    let scores = row_scores(&read[0], &spec).unwrap();
    let plan = row_plan(&read[0], scores, &ctx).unwrap();
    let exception = plan
        .translation_stop(&ctx, 0)
        .unwrap()
        .expect("TranslateNAToAA's exception");
    assert_eq!(exception.site().line, 898);
}
