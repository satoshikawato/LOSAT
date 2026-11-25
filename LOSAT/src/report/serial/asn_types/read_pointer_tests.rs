//! Unit tests of NCBI's pointer forms in the reader (`super::read`, `super::typed`): where NCBI
//! reads a `CRef<>` that is empty (an OPTIONAL member held by `CRef<>`, an element of a list of
//! `CRef<>`s, any such member in a skip), `NULL` is none, `@n` is NCBI's error, `:Type value`
//! naming the declared type is the value; a mandatory member held by `CRef<>` is read in place;
//! a member with a set flag adds `ReadWithSetFlagMember`'s line to an error of its value. The
//! cases are those of probe D3c-3b (psiblast and tblastn 2.17.0 `-in_pssm`, `probe/SUMMARY*.tsv`
//! of the line's folder; makeblastdb's are the `h_pointer_*` and `b_pointer_*` files of
//! `tests/seqdb/mask_parse`).

use super::read::{read_ber_document, read_document, ReadError};
use super::table::named;
use super::{member_index, members, Value};
use crate::report::serial::asn_input::Names;
use crate::seqid::model::{Date, ParsedSeqId, SEQ_ID_ALTERNATIVES};

fn read(root: u16, text: &str) -> Result<Value, ReadError> {
    let (value, posts) = read_document(text.as_bytes(), root, Names::TEST);
    assert!(posts.is_empty(), "{text}: {posts:?}");
    value
}

fn ok(root: u16, text: &str) -> Value {
    read(root, text).unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

fn ber(root: u16, bytes: &[u8]) -> Result<Value, ReadError> {
    let (value, posts) = read_ber_document(bytes, root, Names::TEST);
    assert!(posts.is_empty(), "{bytes:02x?}: {posts:?}");
    value
}

/// The exception of a `CSerialException` that stops the reader, as its frames and lines show it.
fn serial(result: Result<Value, ReadError>) -> String {
    match result {
        Err(ReadError::Serial(e)) => format!("{e:?}"),
        other => panic!("not a CSerialException: {other:?}"),
    }
}

/// `:name value` of ASN.1 binary: the name in the long form of an application tag, then the
/// value in indefinite length (probe D3c-3b `b_params_other`).
fn other_pointer(name: &str, value: &[u8]) -> Vec<u8> {
    let name = name.as_bytes();
    let mut out = vec![0x7f];
    out.extend(name[..name.len() - 1].iter().map(|b| b | 0x80));
    out.push(name[name.len() - 1]);
    out.push(0x80);
    out.extend_from_slice(value);
    out.extend_from_slice(&[0, 0]);
    out
}

/// A PDB-seq-id `{ mol "1ABC" }` in ASN.1 binary with `rel` (context tag 2) holding `rel`.
fn pdb_ber(rel: Option<&[u8]>) -> Vec<u8> {
    let mut out = vec![
        0x30, 0x80, 0xa0, 0x80, 0x1a, 0x04, b'1', b'A', b'B', b'C', 0, 0,
    ];
    if let Some(rel) = rel {
        out.extend_from_slice(&[0xa2, 0x80]);
        out.extend_from_slice(rel);
        out.extend_from_slice(&[0, 0]);
    }
    out.extend_from_slice(&[0, 0]);
    out
}

// Probe D3c-3b t_params_null, t_ident_null, t_inter_null, t_final_null, t_null_comment,
// t_params_other, t_params_other_space, t_ident_other, b_params_null, b_params_other.
#[test]
fn an_optional_ref_member_reads_null_and_its_own_type() {
    let without = ok(named::PDB_SEQ_ID, "PDB-seq-id ::= { mol \"1ABC\" }");
    for text in [
        "PDB-seq-id ::= { mol \"1ABC\", rel NULL }",
        "PDB-seq-id ::= { mol \"1ABC\", rel NULL -- c --\n }",
    ] {
        assert_eq!(ok(named::PDB_SEQ_ID, text), without, "{text}");
    }
    let with = ok(
        named::PDB_SEQ_ID,
        "PDB-seq-id ::= { mol \"1ABC\", rel str \"x\" }",
    );
    for text in [
        "PDB-seq-id ::= { mol \"1ABC\", rel :Date str \"x\" }",
        "PDB-seq-id ::= { mol \"1ABC\", rel : Date str \"x\" }",
    ] {
        assert_eq!(ok(named::PDB_SEQ_ID, text), with, "{text}");
    }
    let rel = member_index(named::PDB_SEQ_ID, "rel");
    assert!(with.member(rel).is_some() && without.member(rel).is_none());
    // ASN.1 binary: `05 00`, and `7f` with the type's name.
    let date = [0xa0, 0x80, 0x1a, 0x01, b'x', 0, 0];
    assert_eq!(ber(named::PDB_SEQ_ID, &pdb_ber(None)).unwrap(), without);
    assert_eq!(
        ber(named::PDB_SEQ_ID, &pdb_ber(Some(&[0x05, 0x00]))).unwrap(),
        without
    );
    assert_eq!(ber(named::PDB_SEQ_ID, &pdb_ber(Some(&date))).unwrap(), with);
    assert_eq!(
        ber(
            named::PDB_SEQ_ID,
            &pdb_ber(Some(&other_pointer("Date", &date)))
        )
        .unwrap(),
        with
    );
}

// Probe D3c-3b t_params_at, t_params_at0, b_params_at, t_params_nullx, t_null_lower,
// t_params_other_wrong, t_params_other_pssm, t_params_other_unknown, t_params_other_empty,
// b_params_other_wrong.
#[test]
fn an_object_index_and_another_class_are_ncbis_errors() {
    for text in [
        "PDB-seq-id ::= { mol \"1ABC\", rel @1 }",
        "PDB-seq-id ::= { mol \"1ABC\", rel @0 }",
    ] {
        let e = serial(read(named::PDB_SEQ_ID, text));
        assert!(
            e.contains("invalid object index: NO_COLLECT defined")
                && e.contains("GetRegisteredObject"),
            "{text}: {e}"
        );
    }
    // In binary NCBI's `ReadObjectPointer` reads the tag byte `5e` as the length of the index (94
    // bytes): its overflow error (`asn_ber_reader`, the `b_pointer_*_at` mask files).
    let e = serial(ber(named::PDB_SEQ_ID, &pdb_ber(Some(&[0x5e, 0x01, 0x01]))));
    assert!(e.contains("overflow error"), "{e}");
    // Not a pointer form: the value of a Date, which neither is.
    for text in [
        "PDB-seq-id ::= { mol \"1ABC\", rel NULLX }",
        "PDB-seq-id ::= { mol \"1ABC\", rel null }",
    ] {
        serial(read(named::PDB_SEQ_ID, text));
    }
    // A value of another class: NCBI's `CSerialException` whether or not its process knows the
    // class (psiblast and tblastn report it with their own message, `Names::own_serial_message`).
    for (root, text, expected) in [
        (
            named::PSSMWITHPARAMETERS,
            "PssmWithParameters ::= { pssm { numRows 1, numColumns 1 }, params :Object-id id 1 }",
            "incompatible member type",
        ),
        (
            named::PSSMWITHPARAMETERS,
            "PssmWithParameters ::= { pssm { numRows 1, numColumns 1 }, params :Pssm { numRows 1, numColumns 1 } }",
            "incompatible member type",
        ),
        (
            named::PDB_SEQ_ID,
            "PDB-seq-id ::= { mol \"1ABC\", rel :Date-std { year 2026 } }",
            "incompatible member type",
        ),
        (
            named::PDB_SEQ_ID,
            "PDB-seq-id ::= { mol \"1ABC\", rel :No-such-type id 1 }",
            "class not found: No-such-type",
        ),
        // The first `CSerialException` inside the value.
        (
            named::PSSMWITHPARAMETERS,
            "PssmWithParameters ::= { pssm { numRows 1, numColumns 1 }, params :Object-id foo 1 }",
            "unexpected member",
        ),
    ] {
        let e = serial(read(root, text));
        assert!(e.contains(expected), "{text}: {e}");
    }
    let date_std = other_pointer(
        "Date-std",
        &[0x30, 0x80, 0xa0, 0x80, 0x02, 0x01, 0x01, 0, 0, 0, 0],
    );
    let e = serial(ber(named::PDB_SEQ_ID, &pdb_ber(Some(&date_std))));
    assert!(e.contains("incompatible member type"), "{e}");
    let other = other_pointer("Object-id", &[0xa0, 0x80, 0x02, 0x01, 0x05, 0, 0]);
    let e = serial(ber(named::PDB_SEQ_ID, &pdb_ber(Some(&other))));
    assert!(e.contains("class not found: Object-id"), "{e}");
    // A read that posts (a class's `PostRead`) posts only where NCBI knows the class, and so does
    // an error of another kind inside the value (the `CUtilException` of a number in text, which
    // the PSSM reader does not catch): refused.
    let refused = |result: Result<Value, ReadError>| match result {
        Err(ReadError::Unsupported(what)) => what,
        other => panic!("not refused: {other:?}"),
    };
    assert_eq!(
        refused(read(
            named::PSSMWITHPARAMETERS,
            "PssmWithParameters ::= { pssm { numRows 1, numColumns 1 }, params :Object-id id x }"
        )),
        "':Object-id', a value of the type Object-id where a PssmParameters belongs (line 1): \
         whether NCBI's reader reads it depends on the types its process has met before"
    );
    assert_eq!(
        refused(read(
            named::PSSMWITHPARAMETERS,
            "PssmWithParameters ::= { pssm { numRows 1, numColumns 1 }, params :Variation-ref { population-data { }, data unknown NULL } }"
        )),
        "':Variation-ref', a value of the type Variation-ref where a PssmParameters belongs (line 1): \
         whether NCBI's reader reads it depends on the types its process has met before"
    );
    // `:[]`: the name's `std::string` of the largest count of bytes.
    assert_eq!(
        read(
            named::PDB_SEQ_ID,
            "PDB-seq-id ::= { mol \"1ABC\", rel :[] str \"x\" }"
        ),
        Err(ReadError::Std("basic_string::_M_create"))
    );
}

// Probe D3c-3b probe9 (`db.ptr.ps.params_other_outside`: `:Blast4-cutoff`) and PR #1633 R2 probe
// (NCBI 2.17.0 psiblast and tblastn `-in_pssm` on the five `NCBI-BlastDL` types): classes that the
// current read's root does not reach are not in NCBI's class registry during `-in_pssm` reading, so
// NCBI throws "class not found" before reading their value (even when the value holds a bad number
// or a byte outside visible ASCII).
#[test]
fn other_class_lookup_is_scoped_to_the_read_root() {
    for (name, body) in [
        ("Blast4-cutoff", "raw-score x"),
        (
            "Blast-db-mask-info",
            "{ algo-id x, algo-program dust, algo-options \"\", masks { masks { }, more FALSE } }",
        ),
        (
            "Blast-mask-list",
            "{ masks { int { from x, to 1, id local id 1 } }, more FALSE }",
        ),
        (
            "Blast-def-line-set",
            "{ { seqid { local id 1 }, taxid x } }",
        ),
        ("Blast-def-line", "{ seqid { local id 1 }, taxid x }"),
        ("Blast-filter-program", "x"),
    ] {
        let text = format!(
            "PssmWithParameters ::= {{ pssm {{ numRows 1, numColumns 1 }}, params :{name} {body} }}"
        );
        let e = serial(read(named::PSSMWITHPARAMETERS, &text));
        assert!(
            e.contains(&format!("class not found: {name}")),
            "{text}: {e}"
        );
    }
    // ASN.1 binary under `PssmWithParameters`: `Object-id` is reachable (`incompatible member
    // type`), while `Blast4-cutoff` is not (`class not found`).
    let pssm_ber = |params: &[u8]| {
        let mut out = vec![
            0x30, 0x80, // PssmWithParameters
            0xa0, 0x80, // pssm
            0x30, 0x80, // Pssm
            0xa2, 0x80, 0x02, 0x01, 0x01, 0, 0, // numRows 1
            0xa3, 0x80, 0x02, 0x01, 0x01, 0, 0, // numColumns 1
            0, 0, 0, 0, // end Pssm, end pssm
            0xa1, 0x80, // params
        ];
        out.extend_from_slice(params);
        out.extend_from_slice(&[0, 0, 0, 0]);
        out
    };
    let obj_id = other_pointer("Object-id", &[0xa0, 0x80, 0x02, 0x01, 0x05, 0, 0]);
    let e = serial(ber(named::PSSMWITHPARAMETERS, &pssm_ber(&obj_id)));
    assert!(e.contains("incompatible member type"), "{e}");
    let blast4_cutoff = other_pointer("Blast4-cutoff", &[0xa1, 0x80, 0x02, 0x01, 0x05, 0, 0]);
    let e = serial(ber(named::PSSMWITHPARAMETERS, &pssm_ber(&blast4_cutoff)));
    assert!(e.contains("class not found: Blast4-cutoff"), "{e}");
}

// Probe D3c-3b t_inst_null, t_query_null: a mandatory member held by `CRef<>` is allocated by the
// generated class and read in place (its value, which `NULL` is not); an OPTIONAL one is unset.
#[test]
fn a_mandatory_ref_member_is_read_in_place() {
    let e = serial(read(named::DBTAG, "Dbtag ::= { db \"x\", tag NULL }"));
    assert!(e.contains("tag"), "{e}");
    let tag = member_index(named::DBTAG, "tag");
    assert!(ok(named::DBTAG, "Dbtag ::= { db \"x\", tag id 1 }")
        .member(tag)
        .is_some());
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:712-756 (cited at
// `read::set_flag_error`)
#[test]
fn a_member_with_a_set_flag_adds_its_line_and_a_ref_member_does_not() {
    let e = serial(read(named::GIIMPORT_ID, "Giimport-id ::= { id 1, db 5 }"));
    assert!(
        e.contains("error while reading db") && e.contains("ReadWithSetFlagMember"),
        "{e}"
    );
    // Through the typed Seq-id reader too.
    let e = serial(read(named::SEQ_ID, "Seq-id ::= giim { id 1, db 5 }"));
    assert!(e.contains("error while reading db"), "{e}");
    // `rel` is held by `CRef<>` (no set flag) and a variant has none.
    let e = serial(read(
        named::PDB_SEQ_ID,
        "PDB-seq-id ::= { mol \"1ABC\", rel str 5 }",
    ));
    assert!(!e.contains("error while reading"), "{e}");
}

// Probe D3c-3b t_descr_elem_null, t_descr_elem_null2, descr_only_null, descr_empty,
// t_descr_elem_other, t_id_elem_null, t_id_elem_other, t_id_elem_at, b_id_elem_other; probe6 and
// probe7 (NCBI reads every such element; the program that uses the list dereferences it or not).
#[test]
fn an_element_of_a_list_of_refs_takes_the_pointer_forms() {
    let elements = |mut value: Value| match &mut value {
        Value::List(items) => std::mem::take(items),
        other => panic!("not a list: {other:?}"),
    };
    for (text, nulls) in [
        ("Seq-descr ::= { NULL }", vec![true]),
        ("Seq-descr ::= { title \"a\", NULL }", vec![false, true]),
        (
            "Seq-descr ::= { NULL, NULL, title \"a\" }",
            vec![true, true, false],
        ),
    ] {
        let items = elements(ok(named::SEQ_DESCR, text));
        let held: Vec<bool> = items.iter().map(|item| *item == Value::Null).collect();
        assert_eq!(held, nulls, "{text}");
    }
    // The empty list is NCBI's error (`CSeq_descr::PostRead`); one of a `NULL` is not empty.
    serial(read(named::SEQ_DESCR, "Seq-descr ::= { }"));
    // An error after the `NULL` is what NCBI's reader stops at.
    serial(read(named::SEQ_DESCR, "Seq-descr ::= { NULL, title 5 }"));
    let e = serial(read(named::SEQ_DESCR, "Seq-descr ::= { title \"a\", @1 }"));
    assert!(e.contains("NO_COLLECT"), "{e}");
    assert_eq!(
        ok(named::SEQ_DESCR, "Seq-descr ::= { :Seqdesc title \"a\" }"),
        ok(named::SEQ_DESCR, "Seq-descr ::= { title \"a\" }")
    );
    let bioseq = ok(
        named::BIOSEQ,
        "Bioseq ::= { id { local id 1, NULL }, inst { repr raw, mol aa } }",
    );
    let id = bioseq.member(member_index(named::BIOSEQ, "id")).cloned();
    assert!(
        matches!(&id, Some(Value::List(ids)) if ids.len() == 2 && ids[1] == Value::Null),
        "{id:?}"
    );
    assert_eq!(
        ok(
            named::BIOSEQ,
            "Bioseq ::= { id { :Seq-id local id 1 }, inst { repr raw, mol aa } }"
        ),
        ok(
            named::BIOSEQ,
            "Bioseq ::= { id { local id 1 }, inst { repr raw, mol aa } }"
        )
    );
}

// Probe D3c-3b dl_name_null, dl_name_other, dl_elem_null, dl_elem_at: the delayed member is read
// by the skip functions, where every member held by `CRef<>` and every element of a list of them
// takes the pointer forms (`SkipPointer`), and a `NULL` element is no empty reference of the
// object.
#[test]
fn a_skip_reads_the_pointer_forms_of_every_ref() {
    let names = member_index(named::AUTH_LIST, "names");
    for text in [
        "Auth-list ::= { names std { { name NULL } } }",
        "Auth-list ::= { names std { { name :Person-id str \"x\" } } }",
        "Auth-list ::= { names std { { name str \"x\" }, NULL } }",
        "Auth-list ::= { names NULL }",
    ] {
        let value = ok(named::AUTH_LIST, text);
        assert!(
            matches!(value.member(names), Some(Value::Delayed { .. })),
            "{text}: {value:?}"
        );
    }
    let e = serial(read(
        named::AUTH_LIST,
        "Auth-list ::= { names std { { name str \"x\" }, @1 } }",
    ));
    assert!(e.contains("NO_COLLECT"), "{e}");
}

/// The typed Seq-id reader (`super::seq_id`): the variants in the model's order, the members in
/// the table's, and the pointer forms of `PDB-seq-id.rel`.
#[test]
fn the_seq_id_reader_takes_the_tables_order_and_pointer_forms() {
    let variants: Vec<&str> = members(named::SEQ_ID).iter().map(|m| m.name()).collect();
    assert_eq!(variants, SEQ_ID_ALTERNATIVES);
    for (ty, names) in super::seq_id::MEMBER_ORDER {
        let table: Vec<&str> = members(ty).iter().map(|m| m.name()).collect();
        assert_eq!(table, names, "{}", super::type_info(ty).name());
    }
    let pdb = |text: &str| match ok(named::SEQ_ID, text).into_seq_id() {
        Some(ParsedSeqId::Pdb(pdb)) => pdb,
        other => panic!("{text}: {other:?}"),
    };
    assert_eq!(pdb("Seq-id ::= pdb { mol \"1ABC\", rel NULL }").rel, None);
    assert_eq!(
        pdb("Seq-id ::= pdb { mol \"1ABC\", rel :Date str \"x\" }").rel,
        Some(Date::Str(b"x".to_vec()))
    );
    let e = serial(read(
        named::SEQ_ID,
        "Seq-id ::= pdb { mol \"1ABC\", rel @2 }",
    ));
    assert!(e.contains("NO_COLLECT"), "{e}");
}
