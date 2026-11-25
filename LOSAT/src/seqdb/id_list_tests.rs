use super::*;

fn binary(words: &[u32]) -> Vec<u8> {
    words.iter().flat_map(|word| word.to_be_bytes()).collect()
}

#[test]
fn a_text_list_loses_its_last_number_without_a_separator_after_it() {
    assert_eq!(read_gi_list(b"11\n22").unwrap(), vec![11]);
    assert_eq!(read_gi_list(b"11\n22\n").unwrap(), vec![11, 22]);
    assert_eq!(read_gi_list(b"11 0 \r\n\n33\r").unwrap(), vec![11, 33]);
}

#[test]
fn a_text_gi_wraps_as_an_unsigned_32_bit_number() {
    assert_eq!(read_gi_list(b"4294967296\n").unwrap(), Vec::<u32>::new());
    assert_eq!(read_gi_list(b"4294967297\n").unwrap(), vec![1]);
}

#[test]
fn the_first_byte_decides_text_or_binary() {
    assert_eq!(read_gi_list(b"").unwrap_err(), "Specified file is empty.");
    assert_eq!(
        read_gi_list(b" 11\n").unwrap_err(),
        "Specified file is not a valid GI/TI list."
    );
    assert_eq!(
        read_gi_list(b"#c\n11\n").unwrap_err(),
        "Invalid byte in textGI list [35]."
    );
    assert_eq!(
        read_gi_list(b"11\t22\n").unwrap_err(),
        "Invalid byte in textGI list [9]."
    );
    assert_eq!(
        read_gi_list(&[0xFF, 0xFF, 0xFF, 0xFF, 0, 0, 0]).unwrap_err(),
        "Specified file is not a valid GI/TI list."
    );
}

#[test]
fn a_binary_gi_list_needs_its_marker_and_count() {
    assert_eq!(
        read_gi_list(&binary(&[0xFFFF_FFFF, 2, 7, 5])).unwrap(),
        vec![7, 5]
    );
    // A length that is not a whole number of words is rounded down.
    let mut odd = binary(&[0xFFFF_FFFF, 1, 9]);
    odd.push(1);
    assert_eq!(read_gi_list(&odd).unwrap(), vec![9]);
    for bad in [
        binary(&[0xFFFF_FFFF, 3, 7, 5]),
        binary(&[0xFFFF_FFFD, 1, 7]),
    ] {
        assert_eq!(
            read_gi_list(&bad).unwrap_err(),
            "Specified file is not a valid binary GI file."
        );
    }
}

#[test]
fn a_ti_list_is_signed_and_widened_past_int8() {
    assert_eq!(
        read_ti_list(b"9223372036854775807\n9223372036854775808\n5\n").unwrap(),
        vec![
            ListNumber::Value(i64::MAX),
            ListNumber::Past,
            ListNumber::Value(5)
        ]
    );
    assert_eq!(
        read_ti_list(&binary(&[0xFFFF_FFFD, 2, 0xFFFF_FFFF, 3])).unwrap(),
        vec![ListNumber::Value(-1), ListNumber::Value(3)]
    );
    assert_eq!(
        read_ti_list(&binary(&[0xFFFF_FFFC, 1, 1, 2])).unwrap(),
        vec![ListNumber::Value((1 << 32) | 2)]
    );
    for bad in [
        binary(&[0xFFFF_FFFF, 1, 3]),
        binary(&[0xFFFF_FFFD, 2, 3]),
        binary(&[0xFFFF_FFFC, 1, 1, 2, 3]),
    ] {
        assert_eq!(
            read_ti_list(&bad).unwrap_err(),
            "Specified file is not a valid binary GI or TI file."
        );
    }
    assert_eq!(
        read_ti_list(b"1a\n").unwrap_err(),
        "Invalid byte in textTI list [97]."
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:1586 IsBinaryGiList_True
#[test]
fn a_list_whose_first_byte_is_0xff_is_binary() {
    // The first words of NCBI's data/prot345b.gil: the marker, its count of 471, a GI.
    assert_eq!(
        binary_kind(&binary(&[0xFFFF_FFFF, 471, 0x0000_CFC0])),
        Ok(Some(false))
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:1591-1592 IsBinaryGiList_False
#[test]
fn a_text_list_and_an_alias_file_are_not_binary() {
    // The first lines of NCBI's data/prot345t.gil and data/totals.nal.
    assert_eq!(binary_kind(b"43381475\n44389145\n43713354\n"), Ok(None));
    assert_eq!(
        binary_kind(b"# Test case for STATS_NSEQ and STATS_TOTLEN\nDBLIST "),
        Ok(None)
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:1597-1598 IsBinaryGiList_EmptyFile
#[test]
fn an_empty_file_is_no_list() {
    assert_eq!(binary_kind(b"").unwrap_err(), "Specified file is empty.");
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:1603-1606 IsBinaryGiList_InvalidFile
#[test]
fn the_null_device_is_no_list() {
    let null = if cfg!(windows) { "NUL" } else { "/dev/null" };
    if let Ok(bytes) = std::fs::read(null) {
        assert_eq!(binary_kind(&bytes).unwrap_err(), "Specified file is empty.");
    }
}

fn general(db: &str, tag: ObjectId) -> ParsedSeqId {
    ParsedSeqId::General(crate::seqid::model::Dbtag {
        db: db.as_bytes().to_vec(),
        tag,
    })
}

fn text(tag: &str) -> ObjectId {
    ObjectId::Str(tag.as_bytes().to_vec())
}

const GET_ID: &[u8] = b"CObject_id::GetId(): Invalid choice selection: NCBI-General::Object-id.str";

/// `CSeqDBGiList::FindId` converts a `ti` id's string tag whatever the list holds, and its
/// simplification of a `PIG` or `BL_ORD_ID` id with a string tag throws (a GI list asks it too: no
/// string of the id is in its empty Seq-ids); `CSeqDBNegativeList::FindId` asks only the ids of
/// the kind of the id, so a list without TIs converts no `ti` tag, and one without Seq-id strings
/// simplifies nothing (review 1 of line RG D1c-1, F3 and F7).
#[test]
fn the_lookups_throw_where_ncbis_throw() {
    let ti = find_id(&[], &[], &[], &general("ti", text("abc"))).unwrap_err();
    assert!(ti
        .get_msg()
        .starts_with(b"Cannot convert string 'abc' to Int8"));
    assert!(ti.places_known());
    assert_eq!(
        find_id(&[], &[1020], &[], &general("ti", ObjectId::Id(1020))),
        Ok(true)
    );
    for db in ["PIG", "BL_ORD_ID"] {
        let thrown = find_id(&[7001], &[], &[], &general(db, text("abc"))).unwrap_err();
        assert_eq!(thrown.get_msg(), GET_ID);
        assert!(!thrown.places_known());
        let thrown = negative_find_id(&[b"x".to_vec()], &general(db, text("abc"))).unwrap_err();
        assert_eq!(thrown.get_msg(), GET_ID);
        assert_eq!(
            negative_find_id(&[], &general(db, text("abc"))),
            Ok((false, false))
        );
        assert_eq!(
            find_id(&[], &[], &[], &general(db, ObjectId::Id(5))),
            Ok(false)
        );
    }
    assert_eq!(
        negative_find_id(&[b"abc".to_vec()], &general("ti", text("abc"))),
        Ok((false, false))
    );
}

/// `SiToOid` halves the ids in their order: on ids that a second simplification left out of order
/// (sq_a.txt sorted as written, then simplified) it misses one that a search of sorted ids finds.
#[test]
fn find_si_halves_the_ids_in_their_order() {
    let ids = |words: &[&str]| -> Vec<Vec<u8>> {
        words.iter().map(|word| word.as_bytes().to_vec()).collect()
    };
    let sorted = ids(&["a", "b", "c"]);
    for word in ["a", "b", "c"] {
        assert!(find_si(&sorted, word.as_bytes()));
    }
    assert!(!find_si(&sorted, b"d"));
    let unsorted = ids(&["np_000103.1", "np_000118.1", "p10102.1", "loc7"]);
    assert!(!find_si(&unsorted, b"loc7"));
    assert!(find_si(&unsorted, b"p10102.1"));
    assert!(find_si(&unsorted, b"np_000103.1"));
}

#[test]
fn a_taxid_list_is_a_set_of_int4_numbers() {
    let set = read_tax_id_list(b"9606\n10090\n9606\n2147483648\n").unwrap();
    assert_eq!(
        set.into_iter().collect::<Vec<_>>(),
        vec![
            ListNumber::Value(9606),
            ListNumber::Value(10090),
            ListNumber::Past,
        ]
    );
    assert_eq!(
        read_tax_id_list(&binary(&[0xFFFF_FFFF, 1, 562]))
            .unwrap()
            .into_iter()
            .collect::<Vec<_>>(),
        vec![ListNumber::Value(562)]
    );
    assert_eq!(
        read_tax_id_list(&binary(&[0xFFFF_FFFF, 2, 562])).unwrap_err(),
        "Specified file is not a valid binary Tax Id List file."
    );
    assert_eq!(
        read_tax_id_list(b"9606,\n").unwrap_err(),
        "Invalid byte in textTAXID list [44]."
    );
    assert_eq!(
        read_tax_id_list(b"9606").unwrap(),
        std::collections::BTreeSet::new()
    );
}
