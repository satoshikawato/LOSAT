//! Tests of the annotation index (`annot_index`): NCBI BLAST+ 2.17.0's standard error and exit
//! status for the annotations of line RG D3-1's `om/` cases (named by their case ids; the texts are
//! those of `om/ncbi/psiblast/<case>.err`).

use super::{index_annots, Abort, AnnotOutcome, Refusal};
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::report::serial::asn_types::{read_document, table::named, Value};

/// The index of the annotations `annots` (each the text of a `Seq-annot`) for psiblast.
fn index(annots: &[&str]) -> Result<AnnotOutcome, Refusal> {
    let values: Vec<Value> = annots
        .iter()
        .map(|text| {
            let (value, posts) =
                read_document(format!("Seq-annot ::= {text}").as_bytes(), named::SEQ_ANNOT);
            assert!(posts.is_empty(), "{posts:?}");
            value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
        })
        .collect();
    index_annots("psiblast", &values)
}

/// The messages of an index that goes on.
fn posts(annots: &[&str]) -> String {
    let outcome = index(annots).unwrap();
    assert_eq!(outcome.abort, None);
    String::from_utf8(outcome.posts).unwrap()
}

/// `CATCH_ALL`'s report of an exception at `file`:`line` (one line of the chain).
fn report(file: &str, line: u32, kind: &str, text: &str) -> String {
    format!("Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/{file}\", line {line}: Error: ({kind}) {text}\n\n")
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
        "Error: (CInvalidChoiceSelection::eFail) [psiblast] CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int\n\
         Error: [psiblast] Bad 'Seq-table location' column (CInvalidChoiceSelection::eFail) CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int\n"
    );
    assert_eq!(
        posts(&[&table(1, "{ header { field-name \"Sorted, max length\" }, data string { \"a\" } }")]),
        "Error: (CSeqTableException::eIncompatibleValueType) [psiblast] CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type\n\
         Error: [psiblast] Bad 'Sorted, max length' column (CSeqTableException::eIncompatibleValueType) CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type\n"
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
        String::from_utf8(outcome.posts).unwrap(),
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
        Abort {
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
        abort(&table(1, &format!("{loc}, {loc}"))).message,
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
        .message,
        report(
            TABLE,
            466,
            "CAnnotException::eBadLocation",
            "ncbi::objects::CSeqTableLocColumns::ParseDefaults() - Conflicting loc columns"
        )
    );
    assert_eq!(
        abort(&table(1, "{ header { field-id location }, data int2 { -32768 } }")).message,
        report(
            "src/objects/seq/../seqtable/SeqTable_multi_data_.cpp",
            207,
            "CInvalidChoiceSelection::eFail",
            "CSeqTable_multi_data::GetLoc(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.int2"
        )
    );
    assert_eq!(
        abort(&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-id location-from }}, data int8 {{ 4294967296 }} }}"))).message,
        report(
            "src/objects/seq/../seqtable/SeqTable_multi_data.cpp",
            441,
            "CSeqTableException::eIncompatibleValueType",
            "CSeqTable_multi_data::TryGetInt4(): value is too big for requested type: 4294967296"
        )
    );
    assert_eq!(
        abort(&table(1, &format!("{{ header {{ field-id location-id }}, data id {{ {Q1} }} }}, {{ header {{ field-name \"disabled\" }}, default int 2 }}"))).message,
        report(
            "src/objects/seq/../seqtable/SeqTable_single_data.cpp",
            58,
            "CSeqTableException::eIncompatibleValueType",
            "ncbi::objects::CSeqTable_single_data::ThrowOverflowError() - CSeqTable_single_data::GetValue(bool&): value is too big for requested type: 2"
        )
    );
    assert_eq!(
        abort("{ data align { { type global, dim 2, segs std { { dim 2, loc { feat local id 1, int { from 0, to 5, id local id 2 } } } } } } }").message,
        report(
            "src/objects/seq/../seqloc/Seq_loc.cpp",
            1552,
            "CSeqLocException::eUnsupported",
            "ncbi::objects::CSeq_loc_CI_Impl::x_ProcessLocation() - CSeq_loc_CI: unsupported location type: feat"
        )
    );
    assert_eq!(
        abort(&format!("{{ data align {{ {{ type global, dim 2, segs std {{ {{ dim -1, loc {{ int {{ from 0, to 5, id {Q1} }} }} }} }} }} }} }}")),
        Abort {
            exit: 255,
            message: "Error: vector::_M_default_append\n".to_string()
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
    assert_eq!(String::from_utf8(outcome.posts).unwrap(), FAILED);
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
        Some(Abort {
            exit: 255,
            message: "\n".to_string()
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
    assert_eq!(outcome.abort.unwrap().message, "\n");
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

// What NCBI does not do the same way on every run: f02, f33, t36, a07, t44 and the bounds.
#[test]
fn refusals_name_their_condition() {
    let refused = |annots: &[&str]| index(annots).unwrap_err();
    assert_eq!(
        refused(&["{ data ftable { { data comment NULL, location int { from 10, to 5, id local id 1 } } } }"]),
        Refusal::EmptyRegion
    );
    assert_eq!(
        refused(&[&failing_feature(), "{ data ftable { { data comment NULL, location int { from 10, to 5, id local id 1 } } } }"]),
        Refusal::MutedDump
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
    // s_i72: two equal `giim` Seq-ids (NCBI: two handles here, one if another `giim` Seq-id of
    // the number were held); s_i73: two different ones (NCBI: one handle here).
    for db in ["x", "y"] {
        assert_eq!(
            refused(&[&table(1, &format!("{{ header {{ field-id location }}, data loc {{ mix {{ int {{ from 10, to 9, id giim {{ id 5, db \"x\" }} }}, int {{ from 0, to 5, id giim {{ id 5, db \"{db}\" }} }} }} }} }}"))]),
            Refusal::GiimIds
        );
    }
    for refusal in [
        Refusal::EmptyRegion,
        Refusal::MutedDump,
        Refusal::NullPointer,
        Refusal::DivideByZero,
        Refusal::BitVector,
        Refusal::PackedOutOfBounds,
        Refusal::DeltaOverflow,
        Refusal::ScaledOverflow,
        Refusal::FieldIdOverflow,
        Refusal::IdOrder,
        Refusal::GiimIds,
        Refusal::TableRows,
        Refusal::StdDim,
    ] {
        assert!(
            refusal.text().starts_with("a PSSM query whose annot has "),
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
    assert_eq!(
        index(&[&disc("giim { id 5 }", "giim { id 5, db \"x\" }")]).unwrap_err(),
        Refusal::GiimIds
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
            String::from_utf8(outcome.posts).unwrap()
        });
        assert_eq!(outcome, expected.map(str::to_string), "{case}");
    }
}
