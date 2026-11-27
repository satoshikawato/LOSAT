//! Tests of the annotation index (`annot_index`): NCBI BLAST+ 2.17.0's standard error and exit
//! status for the annotations of line RG D3-1's `om/` cases (named by their case ids; the texts are
//! those of `om/ncbi/psiblast/<case>.err`).

use super::{index_annots, Abort, AnnotOutcome, Refusal};
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::{read_document, table::named, Value};

/// The index of the annotations `annots` (each the text of a `Seq-annot`) for psiblast.
fn index(annots: &[&str]) -> Result<AnnotOutcome, Refusal> {
    let values: Vec<Value> = annots
        .iter()
        .map(|text| {
            let (value, posts) = read_document(
                format!("Seq-annot ::= {text}").as_bytes(),
                named::SEQ_ANNOT,
                Names::PSSM,
            );
            assert!(posts.is_empty(), "{posts:?}");
            value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
        })
        .collect();
    index_annots("psiblast", &values, &[])
}

impl Abort {
    /// The `CATCH_ALL` report of NCBI's end (`Abort::Native`); no other end is expected.
    fn message(&self) -> Vec<u8> {
        match self {
            Abort::Native { message, .. } => message.clone(),
            Abort::Refused(text) => panic!("refused: {text}"),
        }
    }
}

/// The messages of an index that goes on.
fn posts(annots: &[&str]) -> String {
    let outcome = index(annots).unwrap();
    assert_eq!(outcome.abort, None);
    String::from_utf8(outcome.posts.text()).unwrap()
}

/// `CATCH_ALL`'s report of an exception at `file`:`line` (one line of the chain): NCBI's text in
/// LOSAT's words (`crate::wording::losat_text`).
fn report(file: &str, line: u32, kind: &str, text: &str) -> Vec<u8> {
    crate::wording::losat_text(&format!("Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/{file}\", line {line}: Error: ({kind}) {text}\n\n")).into_owned().into_bytes()
}

const TABLE: &str = "src/objmgr/seq_table_info.cpp";
const Q1: &str = "local str \"Query_1\"";

/// A feature table of one row and the columns `columns`.
fn table(rows: i32, columns: &str) -> String {
    format!("{{ data seq-table {{ feat-type 1, num-rows {rows}, columns {{ {columns} }} }} }}")
}

/// The minus-strand circular feature whose second range is empty: one "Failed to parse location"
/// (case f31), which fails the stream.
fn failing_feature() -> String {
    format!("{{ data ftable {{ {{ data comment NULL, location mix {{ int {{ from 0, to 5, strand minus, id {Q1} }}, int {{ from 10, to 9, strand minus, id {Q1} }} }} }} }} }}")
}

const FAILED: &str = "Error: (1204.6) [psiblast] Failed to parse location of \n";

#[test]
fn annotations_with_non_empty_ranges_are_silent() {
    for annot in [
        // f01, f10, f12, f15 (the reversed range of another Seq-id after one: `more_before`),
        // f16, f24, f32 (a plus-strand circular location), f34 (`circular RNA`).
        format!("{{ data ftable {{ {{ data comment NULL, location int {{ from 1, to 5, id {Q1} }} }} }} }}"),
        format!("{{ data ftable {{ {{ data comment NULL, location whole {Q1} }} }} }}"),
        "{ data ftable { { data comment NULL, location null NULL } } }".to_string(),
        "{ data ftable { { data comment NULL, location mix { int { from 0, to 5, id local id 2 }, int { from 10, to 5, id local id 1 } } } } }".to_string(),
        format!("{{ data ftable {{ {{ data comment NULL, location packed-int {{ {{ from 0, to 2, id {Q1} }}, {{ from 4, to 6, id {Q1} }} }} }} }} }}"),
        format!("{{ data ftable {{ {{ data comment NULL, location mix {{ int {{ from 10, to 15, id {Q1} }}, int {{ from 0, to 5, id {Q1} }} }} }} }} }}"),
        format!("{{ data ftable {{ {{ data comment NULL, location mix {{ int {{ from 10, to 15, id {Q1} }}, int {{ from 1, to 0, id {Q1} }} }} }} }} }}"),
        format!("{{ data ftable {{ {{ data comment NULL, location mix {{ int {{ from 0, to 5, strand minus, id {Q1} }}, int {{ from 10, to 9, strand minus, id {Q1} }} }}, except-text \"circular RNA\" }} }} }}"),
        // m03, i01: no objects; a graph, an alignment, a location set with its region.
        "{ data ftable { } }".to_string(),
        "{ data ids { local id 1, gi 5 } }".to_string(),
        format!("{{ data graph {{ {{ loc int {{ from 0, to 9, id {Q1} }}, numval 1, graph int {{ max 1, min 0, axis 0, values {{ 1 }} }} }} }} }}"),
        format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 1, ids {{ {Q1}, local id 2 }}, starts {{ 0, 0 }}, lens {{ 5 }} }} }} }} }}"),
        format!("{{ desc {{ region int {{ from 0, to 5, id {Q1} }} }}, data locs {{ whole {Q1} }} }}"),
        // l04: two locations are not indexed.
        format!("{{ desc {{ region int {{ from 9, to 5, id {Q1} }} }}, data locs {{ whole {Q1}, whole {Q1} }} }}"),
    ] {
        assert_eq!(posts(&[&annot]), "", "{annot}");
    }
}

// Probe p4 r0261: a text Seq-id without accession and name is another handle at each range, so
// the empty range of one exon of a Spliced-seg is an empty region (NCBI's `Empty region in`).
#[test]
fn text_seq_ids_without_names_are_never_found_again() {
    let spliced = |id: &str| {
        format!("{{ data align {{ {{ type partial, dim 2, segs spliced {{ product-id {id}, product-type protein, exons {{ {{ product-start nucpos 42, product-end nucpos 0, genomic-start 0, genomic-end 0 }}, {{ product-start nucpos 0, product-end nucpos 42, genomic-start 0, genomic-end 0 }} }} }} }} }} }}")
    };
    assert_eq!(
        index(&[&spliced("pir { }")]).unwrap_err(),
        Refusal::EmptyRegion
    );
    assert_eq!(posts(&[&spliced("pir { name \"p\" }")]), "");
    assert_eq!(posts(&[&spliced(Q1)]), "");
}

// Cases t01, t02, t03, t07, t09, t08, t10: the column checks of `x_Initialize` (phase C), in the
// order of the columns; the NCBI_CATCH pairs.
#[test]
fn seq_table_column_messages() {
    assert_eq!(
        posts(&[&table(1, "{ header { field-id 9999 } }")]),
        "Error: (1213.8) [psiblast] SeqTable-column-info.field-id = 9999\n"
    );
    assert_eq!(
        posts(&[&table(1, "{ header { field-id -5 } }")]),
        "Error: (1213.8) [psiblast] SeqTable-column-info.field-id = -5\n"
    );
    assert_eq!(
        posts(&[&table(1, "{ header { field-name \"\" } }")]),
        "Error: (1213.10) [psiblast] SeqTable-column-info.field-name is empty\n"
    );
    assert_eq!(
        posts(&[&table(1, "{ header { } }")]),
        "Error: (1213.9) [psiblast] SeqTable-column-info: neither field-id nor field-name is set\n"
    );
    let location = format!("{{ header {{ field-name \"Seq-table location\" }}, data loc {{ int {{ from 0, to 9, id {Q1} }} }} }}");
    assert_eq!(
        posts(&[&table(1, &format!("{location}, {location}"))]),
        "Error: [psiblast] Duplicate 'Seq-table location' column\n"
    );
    let sorted = "{ header { field-name \"Sorted, max length\" }, data int { 5 } }";
    assert_eq!(
        posts(&[&table(1, &format!("{sorted}, {sorted}"))]),
        "Error: [psiblast] Duplicate 'Sorted, max length' column\n"
    );
    assert_eq!(
        posts(&[&table(1, "{ header { field-name \"Seq-table location\" }, data int { 1 } }")]),
        // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
        crate::wording::losat_text("Error: (CInvalidChoiceSelection::eFail) [psiblast] CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int\n\
         Error: [psiblast] Bad 'Seq-table location' column (CInvalidChoiceSelection::eFail) CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int\n")
    );
    assert_eq!(
        posts(&[&table(1, "{ header { field-name \"Sorted, max length\" }, data string { \"a\" } }")]),
        crate::wording::losat_text("Error: (CSeqTableException::eIncompatibleValueType) [psiblast] CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type\n\
         Error: [psiblast] Bad 'Sorted, max length' column (CSeqTableException::eIncompatibleValueType) CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type\n")
    );
    // A table that is not a feature table checks only its location (t23).
    assert_eq!(
        posts(&["{ data seq-table { feat-type 0, num-rows 1, columns { { header { field-id 9999 } } } } }"]),
        ""
    );
    // Phase C of every annotation precedes the index (f28: the second annotation's 1213.8 comes
    // before the first one's index).
    let outcome = index(&[
        &failing_feature(),
        &table(1, "{ header { field-id 9999 } }"),
    ])
    .unwrap();
    assert_eq!(
        String::from_utf8(outcome.posts.text()).unwrap(),
        format!("Error: (1213.8) [psiblast] SeqTable-column-info.field-id = 9999\n{FAILED}")
    );
}

// The exceptions of the column checks and of the values (`om/SPEC.md` §7.1: t04 E7, t11 E1, t12 E2,
// t30 E11, t40 E14, t73 E16, a19 E17, a20).
#[test]
fn exceptions_stop_the_run_with_their_report() {
    let abort = |annot: &str| index(&[annot]).unwrap().abort.unwrap();
    let partial = "{ header { field-id partial }, data bit '80'H }";
    assert_eq!(
        abort(&table(1, &format!("{partial}, {partial}"))),
        Abort::Native {
            exit: 255,
            message: report(
                TABLE,
                835,
                "CAnnotException::eOtherError",
                "ncbi::objects::CSeqTableInfo::x_Initialize() - Duplicate partial column"
            )
        }
    );
    let loc = format!(
        "{{ header {{ field-id location }}, data loc {{ int {{ from 0, to 4, id {Q1} }} }} }}"
    );
    assert_eq!(
        abort(&table(1, &format!("{loc}, {loc}"))).message(),
        report(
            TABLE,
            355,
            "CAnnotException::eBadLocation",
            "ncbi::objects::CSeqTableLocColumns::SetColumn() - Duplicate loc column"
        )
    );
    assert_eq!(
        abort(&table(
            1,
            &format!("{loc}, {{ header {{ field-id location-id }}, data id {{ {Q1} }} }}")
        ))
        .message(),
        report(
            TABLE,
            466,
            "CAnnotException::eBadLocation",
            "ncbi::objects::CSeqTableLocColumns::ParseDefaults() - Conflicting loc columns"
        )
    );
    assert_eq!(
        abort(&table(1, "{ header { field-id location }, data int2 { -32768 } }")).message(),
        report(
            "src/objects/seq/../seqtable/SeqTable_multi_data_.cpp",
            207,
            "CInvalidChoiceSelection::eFail",
            "CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int2"
        )
    );
    assert_eq!(
        abort(&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id location-from }}, data int8 {{ 4294967296 }} }}"))).message(),
        report(
            "src/objects/seq/../seqtable/SeqTable_multi_data.cpp",
            441,
            "CSeqTableException::eIncompatibleValueType",
            "CSeqTable_multi_data::TryGetInt4(): value is too big for requested type: 4294967296"
        )
    );
    assert_eq!(
        abort(&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-name \"disabled\" }}, default int 2 }}"))).message(),
        report(
            "src/objects/seq/../seqtable/SeqTable_single_data.cpp",
            58,
            "CSeqTableException::eIncompatibleValueType",
            "ncbi::objects::CSeqTable_single_data::ThrowOverflowError() - CSeqTable_single_data::GetValue(bool&): value is too big for requested type: 2"
        )
    );
    assert_eq!(
        abort("{ data align { { type global, dim 2, segs std { { dim 2, loc { feat local id 1, int { from 0, to 5, id local id 2 } } } } } } }").message(),
        report(
            "src/objects/seq/../seqloc/Seq_loc.cpp",
            1552,
            "CSeqLocException::eUnsupported",
            "ncbi::objects::CSeq_loc_CI_Impl::x_ProcessLocation() - CSeq_loc_CI: unsupported location type: feat"
        )
    );
    assert_eq!(
        abort(&format!("{{ data align {{ {{ type global, dim 2, segs std {{ {{ dim -1, loc {{ int {{ from 0, to 5, id {Q1} }} }} }} }} }} }} }}")),
        Abort::Native {
            exit: 255,
            message: b"Error: vector::_M_default_append\n".to_vec()
        }
    );
}

// Cases a02, a24: the warnings of an alignment's sizes.
#[test]
fn alignment_size_warnings() {
    assert_eq!(
        posts(&[&format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 2, ids {{ {Q1}, local id 2 }}, starts {{ 0, 0 }}, lens {{ 1 }} }} }} }} }}")]),
        "Warning: (1207.4) [psiblast] Invalid 'lens' size in denseg\n"
    );
    assert_eq!(
        posts(&[&format!("{{ data align {{ {{ type global, dim 2, segs sparse {{ rows {{ {{ first-id {Q1}, second-id local id 2, numseg 2, first-starts {{ 0 }}, second-starts {{ 0, 1, 2 }}, lens {{ 1 }} }} }} }} }} }} }}")]),
        "Warning: (1207.9) [psiblast] Invalid size of 'first-starts' in sparse-align\nWarning: (1207.10) [psiblast] Invalid size of 'second-starts' in sparse-align\n"
    );
}

// Cases f31, t34, t58, k01, k02, k05: "Failed to parse location of " fails the stream; every
// later message is its header, and `CATCH_ALL`'s report its line feed.
#[test]
fn the_failed_stream_keeps_the_headers() {
    let outcome = index(&[&failing_feature()]).unwrap();
    assert!(outcome.muted);
    assert_eq!(String::from_utf8(outcome.posts.text()).unwrap(), FAILED);
    let rows = format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id location-from }}, data int {{ 5, 5 }} }}, {{ header {{ field-id location-to }}, data int {{ 1, 1 }} }}");
    assert_eq!(
        posts(&[&table(2, &rows)]),
        format!("{FAILED}Error: (1204.6) [psiblast] \n")
    );
    let lens = format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 2, ids {{ {Q1}, local id 2 }}, starts {{ 0, 0 }}, lens {{ 1 }} }} }} }} }}");
    assert_eq!(
        posts(&[&failing_feature(), &lens]),
        format!("{FAILED}Warning: (1207.4) [psiblast] \n")
    );
    let empty_row = table(
        1,
        &format!(
            "{{ header {{ field-id location }}, data loc {{ int {{ from 9, to 4, id {Q1} }} }} }}"
        ),
    );
    assert_eq!(
        posts(&[&failing_feature(), &empty_row, &lens]),
        format!("{FAILED}Error: (1204.8) [psiblast] \nWarning: (1207.4) [psiblast] \n")
    );
    let outcome = index(&[
        &failing_feature(),
        &table(1, "{ header { field-id location }, data int { 7 } }"),
    ])
    .unwrap();
    assert_eq!(
        outcome.abort,
        Some(Abort::Native {
            exit: 255,
            message: b"\n".to_vec()
        })
    );
    // The null Seq-id handle of a row after the stream failed (k03).
    let outcome = index(&[
        &failing_feature(),
        &table(
            2,
            &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}"),
        ),
    ])
    .unwrap();
    assert_eq!(outcome.abort.unwrap().message(), b"\n");
}

// The values of a table's rows: a sparse bit set (t87), an int-delta column (t52), a sorted table
// (t51), a disabled row (t46).
#[test]
fn table_rows_read_their_values() {
    assert_eq!(
        posts(&[&table(2, &format!("{{ header {{ field-id location-id }}, default id {Q1} }}, {{ header {{ field-id location-from }}, data int {{ 9 }}, sparse bit-set '40'H }}, {{ header {{ field-id location-to }}, data int {{ 5, 5 }} }}"))]),
        FAILED
    );
    assert_eq!(
        posts(&[&table(2, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1}, {Q1} }} }}, {{ header {{ field-id location-from }}, data int-delta int {{ 1, 2 }} }}"))]),
        ""
    );
    assert_eq!(
        posts(&[&table(1, &format!("{{ header {{ field-id location-id }}, default id {Q1} }}, {{ header {{ field-id location-from }}, data int {{ 1 }} }}, {{ header {{ field-id location-to }}, data int {{ 5 }} }}, {{ header {{ field-name \"Seq-table location\" }}, default loc int {{ from 0, to 1599, id {Q1} }} }}, {{ header {{ field-name \"Sorted, max length\" }}, default int 100 }}"))]),
        ""
    );
    // A disabled row is not indexed: its reversed range posts nothing.
    assert_eq!(
        posts(&[&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id location-from }}, data int {{ 5 }} }}, {{ header {{ field-id location-to }}, data int {{ 1 }} }}, {{ header {{ field-name \"disabled\" }}, data bit '80'H }}"))]),
        ""
    );
}

/// f33, k04, k07, k08, k09, k13, k16, k17 (`om/SPEC.md` §3 rule 3): after (1204.6) the empty
/// region of a feature, graph, alignment, location set or table location posts its header and
/// what the failed stream gets of the object's ASN.1 text (without its last end of line; of a
/// longer text than NCBI's buffer of 4096 bytes, the bytes that the buffer held when it first
/// flushed), and the run ends with exit status 0 and nothing more: NCBI BLAST+ 2.17.0's bytes.
#[test]
fn an_empty_region_after_the_failed_stream_writes_the_object() {
    let dumped = |annots: &[&str]| {
        let outcome = index(annots).unwrap();
        assert_eq!(
            outcome.abort,
            Some(Abort::Native {
                exit: 0,
                message: Vec::new()
            })
        );
        assert!(outcome.muted);
        let text = String::from_utf8(outcome.posts.text()).unwrap();
        text.strip_prefix(FAILED)
            .expect("(1204.6) first")
            .to_string()
    };
    let feat = |comment: &str| {
        format!("{{ data ftable {{ {{ data comment NULL,{comment} location int {{ from 10, to 5, id local id 1 }} }} }} }}")
    };
    let feat_text = "Seq-feat ::= {\n  data comment NULL,\n  location int {\n    from 10,\n    to 5,\n    id local id 1\n  }\n}";
    assert_eq!(
        dumped(&[&failing_feature(), &feat("")]),
        format!("Error: (1204.1) [psiblast] {feat_text}\n")
    );
    assert_eq!(
        dumped(&[
            &failing_feature(),
            &format!("{{ data graph {{ {{ loc int {{ from 10, to 5, id {Q1} }}, numval 1, graph int {{ max 1, min 0, axis 0, values {{ 1 }} }} }} }} }}")
        ]),
        format!("Error: (1204.2) [psiblast] Seq-graph ::= {{\n  loc int {{\n    from 10,\n    to 5,\n    id {Q1}\n  }},\n  numval 1,\n  graph int {{\n    max 1,\n    min 0,\n    axis 0,\n    values {{\n      1\n    }}\n  }}\n}}\n")
    );
    assert_eq!(
        dumped(&[
            &failing_feature(),
            &format!("{{ data seq-table {{ feat-type 0, num-rows 1, columns {{ {{ header {{ field-name \"Seq-table location\" }}, data loc {{ int {{ from 9, to 4, id {Q1} }} }} }} }} }} }}")
        ]),
        format!("Error: (1204.7) [psiblast] Seq-loc ::= int {{\n  from 9,\n  to 4,\n  id {Q1}\n}}\n")
    );
    assert_eq!(
        dumped(&[
            &failing_feature(),
            &format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 1, ids {{ {Q1}, local id 2 }}, starts {{ 5, 5 }}, lens {{ 0 }} }} }} }} }}")
        ]),
        format!("Error: (1204.3) [psiblast] Seq-align ::= {{\n  type global,\n  dim 2,\n  segs denseg {{\n    dim 2,\n    numseg 1,\n    ids {{\n      {Q1},\n      local id 2\n    }},\n    starts {{\n      5,\n      5\n    }},\n    lens {{\n      0\n    }}\n  }}\n}}\n")
    );
    // k09: the (1204.6) of a table row first.
    assert_eq!(
        dumped(&[
            &table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id location-from }}, data int {{ 5 }} }}, {{ header {{ field-id location-to }}, data int {{ 1 }} }}")),
            &feat("")
        ]),
        format!("Error: (1204.1) [psiblast] {feat_text}\n")
    );
    // k13: the location set's one location, not the region.
    assert_eq!(
        dumped(&[
            &failing_feature(),
            &format!("{{ desc {{ region int {{ from 9, to 1, id {Q1} }} }}, data locs {{ int {{ from 0, to 4, id {Q1} }} }} }}")
        ]),
        format!("Error: (1204.4) [psiblast] Seq-loc ::= int {{\n  from 0,\n  to 4,\n  id {Q1}\n}}\n")
    );
    // k17: a text of 4063 bytes, whole; k16: the first 4096 bytes of a longer one, the buffer
    // full when it first flushed (a whole line of the string and 31 bytes of the next).
    let comment = |letter: &str, count: usize| format!(" comment \"{}\",", letter.repeat(count));
    let header = "Error: (1204.1) [psiblast] ";
    let k17 = dumped(&[&failing_feature(), &feat(&comment("y", 3900))]);
    assert_eq!(k17.len(), header.len() + 4063 + 1);
    assert!(k17.ends_with("  location int {\n    from 10,\n    to 5,\n    id local id 1\n  }\n}\n"));
    let k16 = dumped(&[&failing_feature(), &feat(&comment("x", 6000))]);
    assert_eq!(k16.len(), header.len() + 4096 + 1);
    assert!(k16.ends_with(&format!("\n{}\n{}\n", "x".repeat(78), "x".repeat(31))));
    // Review 1 F1 (w1): the cut inside a string whose wrap at a word's end would come after it
    // (NCBI's write threw first): no end of line inserted before the cut.
    let words = format!(" comment \"zzzz {}\",", "abcdefghi ".repeat(600));
    let w1 = dumped(&[&failing_feature(), &feat(&words)]);
    assert_eq!(w1.len(), header.len() + 4096 + 1);
    assert!(w1.ends_with(&format!("\n{} abcdefgh\n", " abcdefghi".repeat(7))));
}

// What NCBI does not do the same way on every run: f02, t36, a07, t44 and the bounds.
#[test]
fn refusals_name_their_condition() {
    let refused = |annots: &[&str]| index(annots).unwrap_err();
    assert_eq!(
        refused(&["{ data ftable { { data comment NULL, location int { from 10, to 5, id local id 1 } } } }"]),
        Refusal::EmptyRegion
    );
    assert_eq!(
        refused(&[&table(
            2,
            &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}")
        )]),
        Refusal::NullPointer
    );
    assert_eq!(
        refused(&[&format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 0, ids {{ {Q1}, local id 2 }}, starts {{ 0, 0 }}, lens {{ }} }} }} }} }}")]),
        Refusal::DivideByZero
    );
    // BitMagic data that ends after its first byte: NCBI reads past it (`bm/` case m02).
    assert_eq!(
        refused(&[&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id partial }}, data bit-bvector {{ size 1, data '80'H }} }}"))]),
        Refusal::BitVector
    );
    assert_eq!(
        refused(&[&table(
            10_000_001,
            &format!("{{ header {{ field-id location-id }}, default id {Q1} }}")
        )]),
        Refusal::TableRows
    );
    assert_eq!(
        refused(&[&format!("{{ data align {{ {{ type global, dim 2, segs std {{ {{ dim 1000001, loc {{ int {{ from 0, to 5, id {Q1} }} }} }} }} }} }} }}")]),
        Refusal::StdDim
    );
    // s_i72 unmuted: two equal `giim` Seq-ids are two handles, the reversed one an empty region;
    // s_i73: two different ones of a number are one handle. After a feature that mapped giim 5
    // (its handle lives on), (5, "x") and (5, "y") both find that handle: one.
    let giim_row = |db: &str| {
        table(1, &format!("{{ header {{ field-id location }}, data loc {{ mix {{ int {{ from 10, to 9, id giim {{ id 5, db \"x\" }} }}, int {{ from 0, to 5, id giim {{ id 5, db \"{db}\" }} }} }} }} }}"))
    };
    assert_eq!(refused(&[&giim_row("x")]), Refusal::EmptyRegion);
    assert_eq!(posts(&[&giim_row("y")]), "");
    assert_eq!(
        posts(&[
            "{ data ftable { { data comment NULL, location int { from 0, to 5, id giim { id 5 } } } } }",
            &giim_row("y")
        ]),
        ""
    );
    for refusal in [
        Refusal::EmptyRegion,
        Refusal::NullPointer,
        Refusal::DivideByZero,
        Refusal::BitVector,
        Refusal::PackedOutOfBounds,
        Refusal::DeltaOverflow,
        Refusal::ScaledOverflow,
        Refusal::FieldIdOverflow,
        Refusal::IdOrder,
        Refusal::TableRows,
        Refusal::StdDim,
    ] {
        assert!(
            refusal
                .text(crate::algorithm::pssm_search::in_pssm::PSSM_QUERY)
                .starts_with("a PSSM query whose annot has "),
            "{refusal:?}"
        );
    }
}

// r_ambiguous_ids and the `om/` cases `s_f*`, `s_a*`, `s_t01`, `s_g01` of line RG D3-1 (NCBI BLAST+
// 2.17.0): Seq-ids that NCBI's trees take for one are one handle in every kind of annotation; two
// handles leave the reversed range alone, an empty region.
#[test]
fn seq_ids_that_ncbi_takes_for_one_share_their_ranges() {
    let gb = |acc: &str| format!("genbank {{ accession \"{acc}\" }}");
    let rev_then = |y: &str, x: &str| {
        format!("mix {{ int {{ from 10, to 9, id {y} }}, int {{ from 0, to 5, id {x} }} }}")
    };
    let feature =
        |loc: &str| format!("{{ data ftable {{ {{ data comment NULL, location {loc} }} }} }}");
    let std = |locs: &str| {
        format!("{{ type global, dim 2, segs std {{ {{ dim 2, loc {{ {locs} }} }} }} }}")
    };
    assert_eq!(
        posts(&[&feature(&format!(
            "mix {{ int {{ from 0, to 5, id {} }}, int {{ from 0, to 5, id {} }} }}",
            gb("A"),
            gb("a")
        ))]),
        ""
    );
    assert_eq!(
        posts(&[&feature(&rev_then(&gb("ABCDEFGH"), &gb("abcdefgh")))]),
        ""
    );
    assert_eq!(
        posts(&[&format!("{{ data ftable {{ {{ data comment NULL, product int {{ from 0, to 5, id genbank {{ accession \"ab123456\", version 1 }} }}, location int {{ from 0, to 5, id {} }} }} }} }}", gb("AB123456"))]),
        ""
    );
    assert_eq!(
        index(&[&feature(&rev_then(
            "genbank { accession \"AB123456\", version 1 }",
            &gb("AB123456")
        ))])
        .unwrap_err(),
        Refusal::EmptyRegion
    );
    assert_eq!(
        posts(&[&format!(
            "{{ data align {{ {} }} }}",
            std(&format!(
                "int {{ from 10, to 9, id {} }}, int {{ from 0, to 5, id {} }}",
                gb("AB123456"),
                gb("ab123456")
            ))
        )]),
        ""
    );
    let disc = |y: &str, x: &str| {
        format!(
            "{{ data align {{ {{ type global, dim 2, segs disc {{ {}, {} }} }} }} }}",
            std(&format!("int {{ from 10, to 9, id {y} }}")),
            std(&format!("int {{ from 0, to 5, id {x} }}"))
        )
    };
    assert_eq!(
        posts(&[&disc(
            "general { db \"TESTDB\", tag str \"abc\" }",
            "general { db \"testdb\", tag str \"ABC\" }"
        )]),
        ""
    );
    assert_eq!(
        index(&[&disc(
            "general { db \"TESTDB\", tag str \"5\" }",
            "general { db \"TESTDB\", tag id 5 }"
        )])
        .unwrap_err(),
        Refusal::EmptyRegion
    );
    // A Std-seg's location iterator looks giim 5 up once more; its handle goes with it, so the
    // second part's (5, "x") finds the row's handle of (5): one handle, the union not empty.
    assert_eq!(
        posts(&[&disc("giim { id 5 }", "giim { id 5, db \"x\" }")]),
        ""
    );
    let patent = |country: &str, number: &str| {
        format!("patent {{ seqid 1, cit {{ country \"{country}\", id number \"{number}\" }} }}")
    };
    assert_eq!(
        posts(&[&format!("{{ data seq-table {{ feat-type 0, num-rows 1, columns {{ {{ header {{ field-name \"Seq-table location\" }}, data loc {{ {} }} }} }} }} }}", rev_then(&patent("US", "ABC123"), &patent("us", "abc123")))]),
        ""
    );
    assert_eq!(
        posts(&[&format!("{{ data graph {{ {{ loc {}, numval 5, graph int {{ max 1, min 0, axis 0, values {{ 1, 1, 1, 1, 1 }} }} }} }} }}", rev_then("pdb { mol \"1abc\" }", "pdb { mol \"1ABC\" }"))]),
        ""
    );
}

// The `om/` cases `s_o*` of line RG D3-1 (NCBI BLAST+ 2.17.0, psiblast and tblastn the same): a
// feature table row `mix { int { from 9, to 4, id Y }, X minus 0-5, X minus empty }` (Y empty:
// 1204.8; X circular with an empty end: 1204.6), after the stream failed (`_m`) or before (`_u`,
// where Y first is an "Empty region" with an address). NCBI writes them in `CSeq_id_Handle` order:
// by packed number (gi, accession, `general` tag), unpacked last; one number, the addresses'.
#[test]
fn table_row_messages_follow_the_handle_order() {
    for (case, y, x, muted, expected) in [
        ("o01_gi5F_gi7E_m", "gi 7", "gi 5", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.6) [psiblast] \nError: (1204.8) [psiblast] \n")),
        ("o01_gi5F_gi7E_u", "gi 7", "gi 5", false, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \n")),
        ("o02_gi5F_gi3E_m", "gi 3", "gi 5", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o02_gi5F_gi3E_u", "gi 3", "gi 5", false, Err(Refusal::EmptyRegion)),
        ("o03_acc5F_acc7E_m", "genbank { accession \"AB000007\" }", "genbank { accession \"AB000005\" }", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.6) [psiblast] \nError: (1204.8) [psiblast] \n")),
        ("o03_acc5F_acc7E_u", "genbank { accession \"AB000007\" }", "genbank { accession \"AB000005\" }", false, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \n")),
        ("o04_acc9F_gi7E_m", "gi 7", "genbank { accession \"AB000009\" }", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o04_acc9F_gi7E_u", "gi 7", "genbank { accession \"AB000009\" }", false, Err(Refusal::EmptyRegion)),
        ("o05_localF_gi7E_m", "gi 7", "local str \"Query_1\"", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o05_localF_gi7E_u", "gi 7", "local str \"Query_1\"", false, Err(Refusal::EmptyRegion)),
        ("o06_gi7F_localE_m", "local str \"Query_1\"", "gi 7", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.6) [psiblast] \nError: (1204.8) [psiblast] \n")),
        ("o06_gi7F_localE_u", "local str \"Query_1\"", "gi 7", false, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \n")),
        ("o07_genneg5F_gen0E_m", "general { db \"TESTDB\", tag id 0 }", "general { db \"TESTDB\", tag id -5 }", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.6) [psiblast] \nError: (1204.8) [psiblast] \n")),
        ("o07_genneg5F_gen0E_u", "general { db \"TESTDB\", tag id 0 }", "general { db \"TESTDB\", tag id -5 }", false, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \n")),
        ("o08_genstr5F_genstr3E_m", "general { db \"TESTDB\", tag str \"abc003\" }", "general { db \"TESTDB\", tag str \"abc005\" }", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o08_genstr5F_genstr3E_u", "general { db \"TESTDB\", tag str \"abc003\" }", "general { db \"TESTDB\", tag str \"abc005\" }", false, Err(Refusal::EmptyRegion)),
        ("o09_tie_gi5_acc5_m", "genbank { accession \"AB000005\" }", "gi 5", true, Err(Refusal::IdOrder)),
        ("o09_tie_gi5_acc5_u", "genbank { accession \"AB000005\" }", "gi 5", false, Err(Refusal::IdOrder)),
        ("o10_tie_acc5v1_acc5v2_m", "genbank { accession \"AB000005\", version 1 }", "genbank { accession \"AB000005\", version 2 }", true, Err(Refusal::IdOrder)),
        ("o10_tie_acc5v1_acc5v2_u", "genbank { accession \"AB000005\", version 1 }", "genbank { accession \"AB000005\", version 2 }", false, Err(Refusal::IdOrder)),
        ("o11_gi0F_gi5E_m", "gi 5", "gi 0", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o11_gi0F_gi5E_u", "gi 5", "gi 0", false, Err(Refusal::EmptyRegion)),
        ("o12_gineg_gi5_m", "gi 5", "gi -3", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o12_gineg_gi5_u", "gi 5", "gi -3", false, Err(Refusal::EmptyRegion)),
        ("o13_refseq5F_gb3E_m", "genbank { accession \"AB000003\" }", "other { accession \"NP_000005\" }", true, Ok("Error: (1204.6) [psiblast] Failed to parse location of \nError: (1204.8) [psiblast] \nError: (1204.6) [psiblast] \n")),
        ("o13_refseq5F_gb3E_u", "genbank { accession \"AB000003\" }", "other { accession \"NP_000005\" }", false, Err(Refusal::EmptyRegion)),
        ("o14_tie_local_pdb_m", "pdb { mol \"1ABC\" }", "local str \"Query_1\"", true, Err(Refusal::IdOrder)),
        ("o14_tie_local_pdb_u", "pdb { mol \"1ABC\" }", "local str \"Query_1\"", false, Err(Refusal::IdOrder)),
    ] {
        let row = table(1, &format!("{{ header {{ field-id location }}, data loc {{ mix {{ int {{ from 9, to 4, id {y} }}, int {{ from 0, to 5, strand minus, id {x} }}, int {{ from 10, to 9, strand minus, id {x} }} }} }} }}"));
        let mut annots = vec![row];
        if muted {
            annots.insert(0, failing_feature());
        }
        let annots: Vec<&str> = annots.iter().map(String::as_str).collect();
        let outcome = index(&annots).map(|outcome| {
            assert_eq!(outcome.abort, None, "{case}");
            String::from_utf8(outcome.posts.text()).unwrap()
        });
        assert_eq!(outcome, expected.map(str::to_string), "{case}");
    }
}

// The `bm/` probes p01, p03, p02 of line RG D3-1 (NCBI BLAST+ 2.17.0): the size of a bit-bvector
// is its `size` member, which NCBI reads without decoding BitMagic's data (an int-delta over a
// bit-bvector of size 0 has no value at row 0; an extra fuzz column is only asked `IsSet`); the
// first value read decodes it (here BitMagic's exception for an unknown token).
#[test]
fn a_bit_bvector_decodes_only_for_a_value() {
    let bad = "'02018000000100020005040000000C'H";
    let id = format!("{{ header {{ field-id location-id }}, default id {Q1} }}");
    assert_eq!(
        posts(&[&table(1, &format!("{id}, {{ header {{ field-id location-from }}, data int-delta bit-bvector {{ size 0, data {bad} }} }}"))]),
        ""
    );
    assert_eq!(
        posts(&[&table(1, &format!("{id}, {{ header {{ field-id location-from }}, data int {{ 3 }} }}, {{ header {{ field-id location-fuzz-from-lim }}, data bit-bvector {{ size 5, data {bad} }} }}"))]),
        ""
    );
    assert_eq!(
        index(&[&table(1, &format!("{id}, {{ header {{ field-id location-from }}, data int-delta bit-bvector {{ size 1, data {bad} }} }}"))]).unwrap().abort,
        Some(bm_invalid())
    );
}

/// `CATCH_ALL`'s report of BitMagic's `std::logic_error`: NCBI's text in LOSAT's words
/// (`crate::wording::losat_text`).
fn bm_invalid() -> Abort {
    Abort::Native {
        exit: 255,
        message: crate::wording::losat_text("Error: BM::Invalid serialization format\n")
            .into_owned()
            .into_bytes(),
    }
}

// The `bm/` probes m01-m15 of line RG D3-1 (NCBI BLAST+ 2.17.0, psiblast and tblastn): BitMagic's
// exception as the index reads a value ends the run with `CATCH_ALL`'s report (m03: a 64-bit
// vector, m04: an unknown token, m06: a sync mark that does not match, m07: a super-block of
// length 0, m08: token 25, all in the `disabled` column; m11: in its sparse index; m14: in
// `partial`), its line feed alone after the stream failed. In the column checks `NCBI_CATCH`
// posts it and the run goes on (m12: `Sorted, max length`, m13: the sparse index of `Seq-table
// location`; m15: a good `Sorted, max length`). Data that NCBI reads past or crashes on is refused
// (m01, m10: empty data, which NCBI reads at a null pointer; m09: data without its last byte;
// m05: an XOR block, which reads NCBI's null reference vector).
#[test]
fn bitmagic_exceptions_and_refusals() {
    let location = format!(
        "{{ header {{ field-id location }}, default loc int {{ from 9, to 4, id {Q1} }} }}"
    );
    let disabled = |data: &str| {
        table(4, &format!("{{ header {{ field-name \"disabled\" }}, data bit-bvector {{ size 4, data {data} }} }}, {location}"))
    };
    let sparse = |data: &str| {
        table(4, &format!("{{ header {{ field-name \"disabled\" }}, sparse bit-set-bvector {{ size 4, data {data} }}, default bit FALSE, sparse-other bit TRUE }}, {location}"))
    };
    let unknown = "'02018000000100020005040000000C'H";
    for (case, annot) in [
        ("m03", disabled("'3809'H")),
        ("m04", disabled(unknown)),
        ("m06", disabled("'0201800000010002000504000000320509'H")),
        ("m07", disabled("'02018000000100020005040000003800000009'H")),
        ("m08", disabled("'020180000001000200050400000019010000000000000009'H")),
        ("m11", sparse(unknown)),
        ("m14", table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id partial }}, data bit-bvector {{ size 1, data {unknown} }} }}"))),
    ] {
        assert_eq!(index(&[&annot]).unwrap().abort, Some(bm_invalid()), "{case}");
        assert_eq!(
            index(&[&failing_feature(), &annot]).unwrap().abort,
            Some(Abort::Native {
                exit: 255,
                message: b"\n".to_vec()
            }),
            "{case}"
        );
    }
    let id = format!("{{ header {{ field-id location-id }}, default id {Q1} }}");
    // `STD_CATCH`'s posts (NCBI's text in LOSAT's words, `crate::wording::losat_text`).
    assert_eq!(
        posts(&[&table(1, &format!("{id}, {{ header {{ field-name \"Sorted, max length\" }}, data bit-bvector {{ size 4, data {unknown} }} }}"))]),
        crate::wording::losat_text("Error: [psiblast] [Bad 'Sorted, max length' column] Exception: BM::Invalid serialization format\n")
    );
    assert_eq!(
        posts(&[&table(1, &format!("{id}, {{ header {{ field-name \"Seq-table location\" }}, data loc {{ int {{ from 1, to 2, id {Q1} }} }}, sparse bit-set-bvector {{ size 4, data {unknown} }} }}"))]),
        crate::wording::losat_text("Error: [psiblast] [Bad 'Seq-table location' column] Exception: BM::Invalid serialization format\n")
    );
    assert_eq!(
        posts(&[&table(1, &format!("{id}, {{ header {{ field-name \"Sorted, max length\" }}, data bit-bvector {{ size 1, data '020180000001000200050100000013000009'H }} }}"))]),
        ""
    );
    for (case, annot) in [
        ("m01", disabled("''H")),
        ("m02", disabled("'02'H")),
        ("m05", disabled("'0201800000010002000504000000230000000009'H")),
        ("m09", table(8, &format!("{{ header {{ field-name \"disabled\" }}, data bit-bvector {{ size 8, data '02018000000100020005080000004344720F00'H }} }}, {location}"))),
        ("m10", sparse("''H")),
    ] {
        assert_eq!(index(&[&annot]).unwrap_err(), Refusal::BitVector, "{case}");
    }
}

// The live `giim` handles across annotation objects (NCBI BLAST+ 2.17.0: review 3's gq4-gq6, `om/`
// cases s_i77-s_i80): a mapped handle lives on, a Seq-table column's default is made in phase C,
// a simple row's handle stays in the table's `key`, a Std-seg's location iterator holds its own.
#[test]
fn giim_handles_live_as_ncbis_object_manager_keeps_them() {
    let (g5, gx) = ("giim { id 5 }", "giim { id 5, db \"x\" }");
    let f31 = |a: &str, b: &str| {
        format!("{{ data ftable {{ {{ data comment NULL, location mix {{ int {{ from 0, to 5, strand minus, id {a} }}, int {{ from 10, to 9, strand minus, id {b} }} }} }} }} }}")
    };
    let early = |id: &str| {
        format!("{{ data ftable {{ {{ data comment NULL, location int {{ from 1, to 2, id {id} }} }} }} }}")
    };
    assert_eq!(posts(&[&early(gx), &f31(g5, g5)]), FAILED);
    assert_eq!(
        index(&[&early(g5), &f31(gx, g5)]).unwrap_err(),
        Refusal::EmptyRegion
    );
    assert_eq!(
        index(&[&format!("{{ data ftable {{ {{ data comment NULL, product int {{ from 1, to 2, id {gx} }}, location mix {{ int {{ from 0, to 5, strand minus, id {g5} }}, int {{ from 10, to 9, strand minus, id {g5} }} }} }} }} }}")]).unwrap_err(),
        Refusal::EmptyRegion
    );
    assert_eq!(
        posts(&[
            &f31(g5, g5),
            &table(
                0,
                &format!("{{ header {{ field-id location-id }}, default id {gx} }}")
            )
        ]),
        FAILED
    );
    assert_eq!(
        posts(&[&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {gx} }} }}, {{ header {{ field-id location-from }}, data int {{ 9 }} }}, {{ header {{ field-id location-to }}, data int {{ 4 }} }}, {{ header {{ field-id product }}, data loc {{ mix {{ int {{ from 0, to 5, strand minus, id {g5} }}, int {{ from 10, to 9, strand minus, id {g5} }} }} }} }}"))]),
        format!("{FAILED}Error: (1204.6) [psiblast] \n")
    );
    let std = |locs: &str| {
        format!("{{ data align {{ {{ type global, dim 2, segs std {{ {{ dim 2, loc {{ {locs} }} }} }} }} }} }}")
    };
    assert_eq!(
        index(&[&std(&format!("mix {{ int {{ from 10, to 9, id {gx} }}, int {{ from 0, to 5, id {g5} }} }}, int {{ from 0, to 5, id {g5} }}"))]).unwrap_err(),
        Refusal::EmptyRegion
    );
    assert_eq!(
        posts(&[&std(&format!(
            "int {{ from 10, to 9, id {gx} }}, int {{ from 0, to 5, id {g5} }}"
        ))]),
        ""
    );
}

/// The object manager's step walks the values it reads from lists on the heap: 50,000 nested
/// Bioseq-sets, a feature location of 50,000 nested Seq-loc mixes around 50,000 equivs, and an
/// alignment of 50,000 nested `disc` sets are indexed on a thread of 256 KiB, as on the main
/// thread of every target (the Seq-table's data: `table_tests`).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn nests_of_any_depth_are_walked_without_recursion() {
    std::thread::Builder::new()
        .stack_size(256 << 10)
        .spawn(|| {
            let deep = 50_000;
            let read = |root: u16, text: String| {
                let (value, posts) = read_document(
                    text.as_bytes(),
                    root,
                    Names::PSSM,
                );
                assert!(posts.is_empty(), "{posts:?}");
                value.expect("read")
            };
            let wrap = |open: &str, inner: &str, close: &str| {
                format!("{}{inner}{}", open.repeat(deep), close.repeat(deep))
            };
            let bioseq = format!(
                "seq {{ id {{ {Q1} }}, inst {{ repr raw, mol aa, length 1, seq-data ncbistdaa '01'H }} }}"
            );
            let set = read(
                named::BIOSEQ_SET,
                format!(
                    "Bioseq-set ::= {{ seq-set {{ {} }} }}",
                    wrap("set { seq-set { ", &bioseq, " } }")
                ),
            );
            let outcome = super::index_entry(
                "psiblast",
                &crate::report::serial::pssm::SeqEntry::Set(set),
                true,
                &|_| None,
            )
            .expect("indexed");
            assert!(outcome.abort.is_none());
            let int = format!("int {{ from 0, to 1, id {Q1} }}");
            let location = wrap("mix { ", &wrap("equiv { ", &int, " }"), " }");
            let std = format!(
                "{{ type global, dim 2, segs std {{ {{ dim 2, loc {{ {int}, {int} }} }} }} }}"
            );
            let annots = [
                read(
                    named::SEQ_ANNOT,
                    format!("Seq-annot ::= {{ data ftable {{ {{ data region \"x\", location {location} }} }} }}"),
                ),
                read(
                    named::SEQ_ANNOT,
                    format!(
                        "Seq-annot ::= {{ data align {{ {} }} }}",
                        wrap("{ type disc, segs disc { ", &std, " } }")
                    ),
                ),
            ];
            let outcome = index_annots("psiblast", &annots, &[]).expect("indexed");
            assert!(outcome.abort.is_none());
            assert!(outcome.posts.is_empty());
        })
        .expect("a thread")
        .join()
        .expect("no overflow");
}

// Probe D3c-3b probe6, probe7 (NCBI 2.17.0): an element stated as `NULL` is dereferenced where
// the object manager reads it (rc 255 with a report of memory addresses, or a crash: refused) and
// nowhere else (rc 0: `dbxref`, `qual`, a `locs` of two, the `id`s after an `other` with an
// accession).
#[test]
fn elements_stated_as_null_are_dereferenced_where_ncbi_reads_them() {
    let feature = |more: &str| {
        format!("{{ data ftable {{ {{ data comment NULL, location whole {Q1}{more} }} }} }}")
    };
    for annot in [
        feature(", qual { NULL }, dbxref { NULL }"),
        format!("{{ id {{ other {{ accession \"A\" }}, NULL }}, data ftable {{ {{ data comment NULL, location whole {Q1} }} }} }}"),
        format!("{{ desc {{ user {{ type str \"Other\", data {{ NULL }} }} }}, data locs {{ whole {Q1}, NULL }} }}"),
        "{ data ids { NULL } }".to_string(),
        format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 1, ids {{ NULL, {Q1} }}, starts {{ -1, 0 }}, lens {{ 1 }} }} }} }} }}"),
    ] {
        assert_eq!(posts(&[&annot]), "", "{annot}");
    }
    let refused = |annots: &[&str]| index(annots).expect_err("refused");
    for annot in [
        "{ data ftable { NULL } }".to_string(),
        feature(", ids { NULL }"),
        feature(", xref { NULL }"),
        format!("{{ data ftable {{ {{ data comment NULL, location mix {{ whole {Q1}, NULL }} }} }} }}"),
        format!("{{ data ftable {{ {{ data comment NULL, product packed-int {{ NULL }}, location whole {Q1} }} }} }}"),
        format!("{{ id {{ NULL }}, data ftable {{ {{ data comment NULL, location whole {Q1} }} }} }}"),
        format!("{{ desc {{ name \"n\", NULL }}, data ftable {{ {{ data comment NULL, location whole {Q1} }} }} }}"),
        "{ desc { user { type str \"AnnotationTrack\", data { NULL } } }, data ids { } }".to_string(),
        "{ data align { NULL } }".to_string(),
        "{ data graph { NULL } }".to_string(),
        "{ data align { { type disc, segs disc { NULL } } } }".to_string(),
        format!("{{ data align {{ {{ type global, dim 2, segs denseg {{ dim 2, numseg 1, ids {{ {Q1}, NULL }}, starts {{ 0, 0 }}, lens {{ 1 }} }} }} }} }}"),
        table(1, "NULL"),
    ] {
        assert_eq!(refused(&[&annot]), Refusal::NullElement, "{annot}");
    }
    // The annotation itself (`CBioseq_Base_Info::x_SetAnnot`), after the one before it is made.
    let (first, _) = read_document(
        format!("Seq-annot ::= {}", feature("")).as_bytes(),
        named::SEQ_ANNOT,
        Names::PSSM,
    );
    assert_eq!(
        index_annots("psiblast", &[first.unwrap(), Value::Null], &[]),
        Err(Refusal::NullElement)
    );
    // A row's location stated as `NULL` is the row's empty reference: the table row's null pointer.
    assert_eq!(
        refused(&[&table(
            1,
            "{ header { field-id location }, data loc { NULL } }"
        )]),
        Refusal::NullPointer
    );
    // After the stream failed, the exception's report is its line feed (exit status 255).
    let outcome = index(&[
        &failing_feature(),
        "{ data ftable { { data comment NULL, location mix { NULL } } } }",
    ])
    .unwrap();
    assert_eq!(String::from_utf8(outcome.posts.text()).unwrap(), FAILED);
    assert_eq!(
        outcome.abort,
        Some(Abort::Native {
            exit: 255,
            message: b"\n".to_vec()
        })
    );
}
