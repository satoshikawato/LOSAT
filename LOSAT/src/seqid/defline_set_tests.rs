use super::*;

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:486-522 (cited at
// `ordered_key`)
// The Seq-ids that a map keyed by Seq-ids holds once: the case of compared letters, an integer
// string tag, a Textseq-id's name beside its accession and version, a pdb chain's letter.
#[test]
fn ordered_keys_are_equal_where_compare_ordered_is() {
    let text = |acc: Option<&str>, version: Option<i32>, name: Option<&str>| TextseqId {
        name: name.map(|name| name.as_bytes().to_vec()),
        accession: acc.map(|acc| acc.as_bytes().to_vec()),
        release: None,
        version,
    };
    let same = |a: &ParsedSeqId, b: &ParsedSeqId| ordered_key(a) == ordered_key(b);
    let local = |s: &str| ParsedSeqId::Local(ObjectId::Str(s.as_bytes().to_vec()));
    assert!(same(&local("p5"), &local("P5")));
    assert!(same(&local("5"), &ParsedSeqId::Local(ObjectId::Id(5))));
    assert!(!same(&local("05"), &ParsedSeqId::Local(ObjectId::Id(5))));
    assert!(!same(&local("p5"), &local("p6")));
    let other = |acc, version, name| ParsedSeqId::Other(text(acc, version, name));
    assert!(same(
        &other(Some("NP_1"), Some(1), Some("a")),
        &other(Some("np_1"), Some(1), Some("b"))
    ));
    assert!(!same(
        &other(Some("NP_1"), Some(1), None),
        &other(Some("NP_1"), Some(2), None)
    ));
    assert!(!same(
        &other(Some("NP_1"), None, Some("a")),
        &other(Some("NP_1"), None, Some("b"))
    ));
    assert!(!same(
        &other(Some("NP_1"), Some(1), None),
        &other(Some("NP_1"), None, None)
    ));
    assert!(!same(
        &other(Some("X"), Some(1), None),
        &ParsedSeqId::Genbank(text(Some("X"), Some(1), None))
    ));
    let general = |db: &str, tag: ObjectId| {
        ParsedSeqId::General(Dbtag {
            db: db.as_bytes().to_vec(),
            tag,
        })
    };
    assert!(same(
        &general("BL_ORD_ID", ObjectId::Id(7)),
        &general("bl_ord_id", ObjectId::Str(b"7".to_vec()))
    ));
    assert!(!same(
        &general("BL_ORD_ID", ObjectId::Id(7)),
        &general("BL_ORD_ID", ObjectId::Id(8))
    ));
    let pdb = |chain: Option<i32>, chain_id: Option<&str>| {
        ParsedSeqId::Pdb(PdbSeqId {
            mol: b"1ABC".to_vec(),
            chain,
            rel: None,
            chain_id: chain_id.map(|id| id.as_bytes().to_vec()),
        })
    };
    assert!(same(&pdb(Some(65), None), &pdb(None, Some("A"))));
    assert!(!same(&pdb(None, Some("a")), &pdb(None, Some("A"))));
    assert!(!same(&ParsedSeqId::Gi(5), &ParsedSeqId::Gibbsq(5)));
}

// `CObject_id::GetIdType` (Object_id.cpp:93-161): the strings that are integers, and those
// that stay strings (an error of `StringToInt8`, a `+`, leading zeros, `-0`). The conversion
// ends at a 0 byte, so `5\0x` is the integer 5 (the string's size is still 3).
#[test]
fn object_id_strings_are_typed_as_ncbis_get_id_type() {
    let number = |text: &[u8]| object_id_number(&ObjectId::Str(text.to_vec()));
    let integers: [(&[u8], i64); 7] = [
        (b"5", 5),
        (b"0", 0),
        (b"-12", -12),
        (b"9223372036854775807", i64::MAX),
        (b"-9223372036854775808", i64::MIN),
        (b"5\x00x", 5),
        (b"-7\x00", -7),
    ];
    for (text, value) in integers {
        assert_eq!(number(text), Some(value), "{text:?}");
    }
    let strings: [&[u8]; 14] = [
        b"",
        b"-",
        b"+5",
        b"05",
        b"-05",
        b"-0",
        b"00",
        b"0\x00",
        b" 5",
        b"5 ",
        b"5x",
        b"9223372036854775808",
        b"abc",
        b"\xff5",
    ];
    for text in strings {
        assert_eq!(number(text), None, "{text:?}");
    }
    assert_eq!(object_id_number(&ObjectId::Id(-3)), Some(-3));
}
