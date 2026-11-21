//! Tests of the range model (`range`): `COpenRange`, `CHandleRange` and `CHandleRangeMap` as the
//! object manager's index reads them (line RG D3-1 `om/SPEC.md` §5.8).

use super::super::view::V;
use super::{Handle, HandleRange, HandleRangeMap, LiveHandles, Mode, Range, MAX, UNKNOWN};
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

// NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:126-206, 338-376 (cited at
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

/// The map of the location `text`, as the first annotation object of a run with a database.
fn map(text: &str, mode: Mode) -> HandleRangeMap {
    map_in(text, mode, &mut LiveHandles::default())
}

fn map_in(text: &str, mode: Mode, handles: &mut LiveHandles) -> HandleRangeMap {
    let value = location(text);
    HandleRangeMap::of_location(V::new(named::SEQ_LOC, &value), mode, handles)
}

// NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:91-205, 249-281 (cited at
// `HandleRangeMap::add_location` and `HandleRangeMap::add`)
#[test]
fn a_map_marks_the_ranges_between_seq_ids() {
    // f15: the reversed range of local id 1 after local id 2 gets `more_before`: from 0.
    let m = map(
        "mix { int { from 0, to 5, id local id 2 }, int { from 10, to 5, id local id 1 } }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 2);
    assert_eq!(
        m.map[1].0,
        Handle::of(
            &ParsedSeqId::Local(ObjectId::Id(1)),
            &mut LiveHandles::default()
        )
    );
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
    // A local string is NCBI's case-insensitive key, and so is an accession (r_ambiguous_ids).
    let m = map(
        "mix { int { from 0, to 5, id local str \"Q\" }, int { from 10, to 12, id local str \"q\" } }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 1);
    let m = map(
        "mix { int { from 0, to 5, id genbank { accession \"A\" } }, int { from 0, to 5, id genbank { accession \"a\" } } }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 1);
    assert_eq!(
        Handle::of(
            &ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec())),
            &mut LiveHandles::default()
        ),
        Handle::of(
            &ParsedSeqId::Local(ObjectId::Str(b"query_1".to_vec())),
            &mut LiveHandles::default()
        )
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

// The Seq-id pairs of line RG D3-1's `om/` cases `s_i*` (NCBI BLAST+ 2.17.0, psiblast and tblastn:
// a feature table row `mix { int { from 10, to 9, id Y }, int { from 0, to 5, id X } }` after the
// stream failed posts `Error: (1204.8)` exactly when Y and X are two handles).
#[test]
fn seq_ids_are_one_handle_as_ncbi_trees_take_them() {
    for (case, y, x, one) in [
        (
            "i01_gb_packed_case",
            "genbank { accession \"AB123456\" }",
            "genbank { accession \"ab123456\" }",
            true,
        ),
        (
            "i02_gb_plain_case",
            "genbank { accession \"ABCDEFGH\" }",
            "genbank { accession \"abcdefgh\" }",
            true,
        ),
        (
            "i03_gb_ver_absent",
            "genbank { accession \"AB123456\", version 1 }",
            "genbank { accession \"AB123456\" }",
            false,
        ),
        (
            "i04_gb_ver_1_2",
            "genbank { accession \"AB123456\", version 1 }",
            "genbank { accession \"AB123456\", version 2 }",
            false,
        ),
        (
            "i05_gb_ver_case",
            "genbank { accession \"AB123456\", version 2 }",
            "genbank { accession \"ab123456\", version 2 }",
            true,
        ),
        (
            "i06_gb_name_case",
            "genbank { name \"LOCUSX\" }",
            "genbank { name \"locusx\" }",
            true,
        ),
        (
            "i07_gb_acc_vs_accname",
            "genbank { accession \"AB123456\" }",
            "genbank { name \"LOCUSX\", accession \"AB123456\" }",
            false,
        ),
        (
            "i08_gb_accname_case",
            "genbank { name \"LOCUSX\", accession \"AB123456\" }",
            "genbank { name \"locusx\", accession \"ab123456\" }",
            true,
        ),
        (
            "i09_gb_accname_name_diff",
            "genbank { name \"LOCUSX\", accession \"AB123456\" }",
            "genbank { name \"LOCUSY\", accession \"AB123456\" }",
            false,
        ),
        (
            "i10_gb_name_vs_accname",
            "genbank { name \"LOCUSX\" }",
            "genbank { name \"LOCUSX\", accession \"AB123456\" }",
            false,
        ),
        (
            "i11_gb_embl_packed",
            "genbank { accession \"AB123456\" }",
            "embl { accession \"AB123456\" }",
            true,
        ),
        (
            "i12_gb_ddbj_packed_ver",
            "genbank { accession \"AB123456\", version 1 }",
            "ddbj { accession \"ab123456\", version 1 }",
            true,
        ),
        (
            "i13_gb_embl_plain_name",
            "genbank { name \"N\", accession \"AB123456\" }",
            "embl { name \"N\", accession \"AB123456\" }",
            false,
        ),
        (
            "i14_gb_embl_plain_acc",
            "genbank { accession \"ABCDEFGH\" }",
            "embl { accession \"ABCDEFGH\" }",
            false,
        ),
        (
            "i15_gb_release_case",
            "genbank { accession \"AB123456\", release \"rel\" }",
            "genbank { accession \"AB123456\", release \"REL\" }",
            false,
        ),
        (
            "i16_gb_release_same",
            "genbank { accession \"AB123456\", release \"rel\" }",
            "genbank { accession \"ab123456\", release \"rel\" }",
            true,
        ),
        (
            "i17_other_packed_case",
            "other { accession \"NP_000123\" }",
            "other { accession \"np_000123\" }",
            true,
        ),
        (
            "i18_other_ver_absent",
            "other { accession \"NP_000123\", version 1 }",
            "other { accession \"NP_000123\" }",
            false,
        ),
        (
            "i19_tpg_vs_gb",
            "tpg { accession \"AB123456\" }",
            "genbank { accession \"AB123456\" }",
            false,
        ),
        (
            "i20_sp_case",
            "swissprot { accession \"P12345\" }",
            "swissprot { accession \"p12345\" }",
            true,
        ),
        (
            "i21_sp_release",
            "swissprot { accession \"P12345\", release \"reviewed\" }",
            "swissprot { accession \"P12345\" }",
            false,
        ),
        (
            "i22_sp_name_case",
            "swissprot { name \"ABC_HUMAN\" }",
            "swissprot { name \"abc_human\" }",
            true,
        ),
        (
            "i23_gb_zeros_case",
            "genbank { accession \"AB0001234\" }",
            "genbank { accession \"ab0001234\" }",
            true,
        ),
        (
            "i24_gb_zeros_diff",
            "genbank { accession \"AB0001234\" }",
            "genbank { accession \"AB001234\" }",
            false,
        ),
        (
            "i25_prf_name_case",
            "prf { name \"ABC\" }",
            "prf { name \"abc\" }",
            true,
        ),
        ("i26_pir_empty", "pir { }", "pir { }", false),
        (
            "i30_gen_str_case",
            "general { db \"TESTDB\", tag str \"abc\" }",
            "general { db \"testdb\", tag str \"ABC\" }",
            true,
        ),
        (
            "i31_gen_strpacked_case",
            "general { db \"TESTDB\", tag str \"abc123\" }",
            "general { db \"testdb\", tag str \"ABC123\" }",
            true,
        ),
        (
            "i32_gen_str_vs_id",
            "general { db \"TESTDB\", tag str \"5\" }",
            "general { db \"TESTDB\", tag id 5 }",
            false,
        ),
        (
            "i33_gen_id_db_case",
            "general { db \"TESTDB\", tag id 5 }",
            "general { db \"testdb\", tag id 5 }",
            true,
        ),
        (
            "i34_gen_str_zero_digits",
            "general { db \"TESTDB\", tag str \"0123\" }",
            "general { db \"TESTDB\", tag str \"123\" }",
            false,
        ),
        (
            "i40_local_str_case",
            "local str \"abc\"",
            "local str \"ABC\"",
            true,
        ),
        (
            "i41_local_str_vs_id",
            "local str \"5\"",
            "local id 5",
            false,
        ),
        (
            "i42_local_query_case",
            "local str \"query_1\"",
            "local str \"Query_1\"",
            true,
        ),
        (
            "i50_pdb_mol_case",
            "pdb { mol \"1abc\" }",
            "pdb { mol \"1ABC\" }",
            true,
        ),
        (
            "i51_pdb_mol_first",
            "pdb { mol \"abcd\" }",
            "pdb { mol \"Abcd\" }",
            false,
        ),
        (
            "i52_pdb_chain_vs_chainid",
            "pdb { mol \"1ABC\", chain 65 }",
            "pdb { mol \"1ABC\", chain-id \"A\" }",
            true,
        ),
        (
            "i53_pdb_chainid_case",
            "pdb { mol \"1ABC\", chain-id \"A\" }",
            "pdb { mol \"1ABC\", chain-id \"a\" }",
            false,
        ),
        (
            "i54_pdb_chain_default",
            "pdb { mol \"1ABC\" }",
            "pdb { mol \"1ABC\", chain 32 }",
            false,
        ),
        (
            "i55_pdb_rel_std",
            "pdb { mol \"1ABC\", rel std { year 2000 } }",
            "pdb { mol \"1ABC\" }",
            true,
        ),
        (
            "i56_pdb_rel_str",
            "pdb { mol \"1ABC\", rel str \"x\" }",
            "pdb { mol \"1ABC\" }",
            false,
        ),
        (
            "i60_pat_case",
            "patent { seqid 1, cit { country \"US\", id number \"ABC123\" } }",
            "patent { seqid 1, cit { country \"us\", id number \"abc123\" } }",
            true,
        ),
        (
            "i61_pat_number_vs_app",
            "patent { seqid 1, cit { country \"US\", id number \"123\" } }",
            "patent { seqid 1, cit { country \"US\", id app-number \"123\" } }",
            false,
        ),
        (
            "i62_pat_seqid",
            "patent { seqid 1, cit { country \"US\", id number \"123\" } }",
            "patent { seqid 2, cit { country \"US\", id number \"123\" } }",
            false,
        ),
        (
            "i63_pat_doctype",
            "patent { seqid 1, cit { country \"US\", id number \"123\", doc-type \"x\" } }",
            "patent { seqid 1, cit { country \"US\", id number \"123\" } }",
            true,
        ),
        (
            "i70_gi_vs_acc",
            "gi 5",
            "genbank { accession \"AB000005\" }",
            false,
        ),
        ("i71_gibbsq_vs_gibbmt", "gibbsq 5", "gibbmt 5", false),
    ] {
        let m = map(
            &format!("mix {{ int {{ from 10, to 9, id {y} }}, int {{ from 0, to 5, id {x} }} }}"),
            Mode::None,
        );
        assert_eq!(m.map.len() == 1, one, "{case}");
    }
}

// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-207 (cited at
// `HandleRangeMap::ncbi_order`): `m_Packed - 1` as an unsigned number, then the info's address.
#[test]
fn a_map_iterates_its_seq_ids_by_packed_number() {
    let m = map(
        "mix { whole gi 7, whole local str \"Q\", whole gi 5, whole general { db \"D\", tag id 0 }, whole general { db \"D\", tag id -5 }, whole genbank { accession \"AB000003\" }, whole pdb { mol \"1ABC\" }, whole gi 3, whole genbank { accession \"AB000003\", version 1 }, whole general { db \"D\", tag str \"abc004\" }, whole gi 0, whole gi -3 }",
        Mode::None,
    );
    assert_eq!(m.map.len(), 12);
    // 3: AB000003, gi 3, AB000003.1 (one number: the addresses' order); 4: abc004; 5; 7; then the
    // negative numbers: general -5 (packed -6), gi -3, general 0 (packed -1); unpacked last: local,
    // PDB, gi 0.
    let at = |entry: &(Handle, HandleRange)| {
        m.map
            .iter()
            .position(|other| std::ptr::eq(other, entry))
            .unwrap()
    };
    let order: Vec<Vec<usize>> = m
        .ncbi_order()
        .into_iter()
        .map(|group| group.into_iter().map(at).collect())
        .collect();
    assert_eq!(
        order,
        vec![
            vec![5, 7, 8],
            vec![9],
            vec![2],
            vec![0],
            vec![4],
            vec![11],
            vec![3],
            vec![1, 6, 10]
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2964 (cited at `Handle`):
// the `om/` cases `s_i72`, `s_i73` and review 2's `id13`, `g1`-`g5` (NCBI BLAST+ 2.17.0).
#[test]
fn giim_handles_follow_ncbis_inverted_tree() {
    let pair = |y: &str, x: &str| {
        format!("mix {{ int {{ from 10, to 9, id {y} }}, int {{ from 0, to 5, id {x} }} }}")
    };
    // Equal Seq-ids: two handles; different ones of a number: one; the first made decides.
    assert_eq!(
        map(&pair("giim { id 5 }", "giim { id 5 }"), Mode::None)
            .map
            .len(),
        2
    );
    assert_eq!(
        map(
            &pair("giim { id 5 }", "giim { id 5, db \"x\" }"),
            Mode::None
        )
        .map
        .len(),
        1
    );
    let m = map(
        "mix { int { from 0, to 5, id giim { id 5, db \"x\" } }, int { from 0, to 5, id giim { id 5, db \"y\" } }, int { from 0, to 5, id giim { id 5, db \"x\" } }, int { from 0, to 5, id giim { id 6 } } }",
        Mode::None,
    );
    // (5, x) makes the first handle, (5, y) finds it, the second (5, x) makes another.
    assert_eq!(m.map.len(), 3);
    // The query's Seq-id (5, "q") holds the first handle for the whole run: (5, "q") gets a new
    // handle, any other value of 5 the query's.
    let query = [crate::seqid::ParsedSeqId::Giim(
        crate::seqid::model::GiimportId {
            id: 5,
            db: Some(b"q".to_vec()),
            release: None,
        },
    )];
    let mut handles = LiveHandles::new(&query);
    let m = map_in(
        &pair("giim { id 5, db \"q\" }", "giim { id 5 }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 2);
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5, db \"x\" }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 1);
    // An earlier object's mapped handle lives on (review 3 gq4, gq5; om s_i74); an unmapped one
    // goes with its object.
    let earlier = |mapped: bool| {
        let mut handles = LiveHandles::default();
        let m = map_in(
            "int { from 1, to 2, id giim { id 5 } }",
            Mode::None,
            &mut handles,
        );
        if mapped {
            handles.mapped_all([&m]);
        }
        handles.end_object();
        handles
    };
    let mut handles = earlier(true);
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5, db \"x\" }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 2);
    let m = map_in(
        &pair("giim { id 5, db \"x\" }", "giim { id 5, db \"y\" }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 1);
    let mut handles = earlier(false);
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5, db \"x\" }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 1);
    // A Seq-table's `key` keeps a simple row's handle until the next one (om s_i77); a location
    // iterator's handle is held until it goes (om s_i79).
    let mut handles = LiveHandles::default();
    let key = Handle::of(
        &crate::seqid::ParsedSeqId::Giim(crate::seqid::model::GiimportId {
            id: 5,
            db: Some(b"x".to_vec()),
            release: None,
        }),
        &mut handles,
    );
    handles.set_key(Some(&key));
    handles.end_object();
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5 }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 1);
    handles.end_object();
    handles.end_table();
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5 }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 2);
    let mut handles = LiveHandles::default();
    Handle::look_up_in_iterator(
        &crate::seqid::ParsedSeqId::Giim(crate::seqid::model::GiimportId {
            id: 5,
            db: Some(b"x".to_vec()),
            release: None,
        }),
        &mut handles,
    );
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5 }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 1);
    handles.end_iterator();
    handles.end_object();
    let m = map_in(
        &pair("giim { id 5 }", "giim { id 5 }"),
        Mode::None,
        &mut handles,
    );
    assert_eq!(m.map.len(), 2);
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.hpp:359-388 (F-LH: rows `db.flh.ps.*`
// and `db.flh.tbn.*`, NCBI's "duplicate Bioseq id" message of a PSSM query)
#[test]
fn an_accession_is_written_by_the_seq_id_that_made_the_info_of_its_packed_key() {
    use crate::seqid::model::TextseqId;
    let acc = |text: &str, version: Option<i32>| TextseqId {
        accession: Some(text.as_bytes().to_vec()),
        version,
        name: None,
        release: None,
    };
    let embl = |text: &str, version| ParsedSeqId::Embl(acc(text, version));
    let genbank = |text: &str, version| ParsedSeqId::Genbank(acc(text, version));
    let written = |query: &[ParsedSeqId], id: &ParsedSeqId| {
        let mut handles = LiveHandles::new(query);
        let handle = Handle::of_run(id, &mut handles);
        handles.as_string(&handle, id)
    };
    // Another number of the prefix, digit count and version made the info: its alternative.
    let first = [embl("AAA00002", Some(1))];
    assert_eq!(
        written(&first, &genbank("AAA00001", Some(1))),
        "emb|AAA00001.1|"
    );
    // The accession's own case, the info's alternative; another prefix case shares the info.
    let first = [embl("aaa00002", Some(1))];
    assert_eq!(
        written(&first, &genbank("AAA00001", Some(1))),
        "emb|AAA00001.1|"
    );
    // Another version, digit count or prefix, or no version: an info of its own.
    for other in [
        embl("AAA00002", Some(2)),
        embl("AAA000002", Some(1)),
        embl("AAB00002", Some(1)),
        embl("AAA00002", None),
    ] {
        assert_eq!(
            written(&[other], &genbank("AAA00001", Some(1))),
            "gb|AAA00001.1|"
        );
    }
    // More than 6 digits: the leading zeros move into the prefix, so `AAA0000012` and `AAA0000111`
    // share the prefix `AAA0` and the digit count; `AAA1000012` has 7 real digits (prefix `AAA`).
    assert_eq!(
        written(
            &[embl("AAA0000012", Some(1))],
            &genbank("AAA0000111", Some(1))
        ),
        "emb|AAA0000111.1|"
    );
    assert_eq!(
        written(
            &[embl("AAA1000012", Some(1))],
            &genbank("AAA0000111", Some(1))
        ),
        "gb|AAA0000111.1|"
    );
    // Version 0 is a version: it shares the info with version 0 (written without a version), and
    // is another info than no version, on either side.
    assert_eq!(
        written(&[embl("AAA00002", Some(0))], &genbank("AAA00001", Some(0))),
        "emb|AAA00001|"
    );
    assert_eq!(
        written(&[embl("AAA00002", None)], &genbank("AAA00001", Some(0))),
        "gb|AAA00001|"
    );
    assert_eq!(
        written(&[embl("AAA00002", Some(0))], &genbank("AAA00001", None)),
        "gb|AAA00001|"
    );
}
