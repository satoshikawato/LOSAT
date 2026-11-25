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

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:486-522
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:178-225
// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:145-161
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:60-67
// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:70-91
#[test]
fn compare_ordered_matches_ncbi_seq_id_order() {
    use crate::seqid::model::IdPat;
    use std::cmp::Ordering;

    let text = |acc: Option<&str>, version: Option<i32>, name: Option<&str>| TextseqId {
        name: name.map(|name| name.as_bytes().to_vec()),
        accession: acc.map(|acc| acc.as_bytes().to_vec()),
        release: None,
        version,
    };
    let local_str = |s: &str| ParsedSeqId::Local(ObjectId::Str(s.as_bytes().to_vec()));
    let general = |db: &str, tag: ObjectId| {
        ParsedSeqId::General(Dbtag {
            db: db.as_bytes().to_vec(),
            tag,
        })
    };
    let pdb = |mol: &str, chain: Option<i32>, chain_id: Option<&str>| {
        ParsedSeqId::Pdb(PdbSeqId {
            mol: mol.as_bytes().to_vec(),
            chain,
            rel: None,
            chain_id: chain_id.map(|id| id.as_bytes().to_vec()),
        })
    };
    let patent = |seqid: i32, country: &str, num: &str| {
        ParsedSeqId::Patent(PatentSeqId {
            seqid,
            cit: IdPat {
                country: country.as_bytes().to_vec(),
                id: PatentNumber::Number(num.as_bytes().to_vec()),
                doc_type: None,
            },
        })
    };

    // Alternative tag order: local < gibbsq < genbank < patent < other < general < gi < pdb.
    // Notably `bbs|123` (e_Gibbsq) precedes `gi|1` (e_Gi) even though `gi|1` is packed.
    assert_eq!(
        compare_ordered(&ParsedSeqId::Gibbsq(123), &ParsedSeqId::Gi(1)),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&local_str("z"), &ParsedSeqId::Gibbsq(1)),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(
            &ParsedSeqId::Genbank(text(Some("ZZZ"), Some(1), None)),
            &ParsedSeqId::Other(text(Some("AAA"), Some(1), None))
        ),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&general("BL_ORD_ID", ObjectId::Id(0)), &ParsedSeqId::Gi(1)),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&ParsedSeqId::Gi(1), &pdb("1ABC", Some(65), None)),
        Some(Ordering::Less)
    );

    // Local/General Object-id: integer tags precede string tags; string tags are case-insensitive.
    assert_eq!(
        compare_ordered(&ParsedSeqId::Local(ObjectId::Id(99)), &local_str("a")),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&local_str("5"), &ParsedSeqId::Local(ObjectId::Id(6))),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&local_str("5"), &ParsedSeqId::Local(ObjectId::Id(5))),
        Some(Ordering::Equal)
    );
    assert_eq!(
        compare_ordered(&local_str("alpha"), &local_str("BETA")),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(
            &general("aa", ObjectId::Id(9)),
            &general("AB", ObjectId::Id(1))
        ),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(
            &general("db", ObjectId::Id(10)),
            &general("DB", ObjectId::Str(b"05".to_vec()))
        ),
        Some(Ordering::Less)
    );

    // Textseq-id: accession (None < Some, nocase), version (None < Some), then name unless both
    // accession and version are set.
    let gb = |acc, ver, name| ParsedSeqId::Genbank(text(acc, ver, name));
    assert_eq!(
        compare_ordered(&gb(None, None, Some("z")), &gb(Some("A"), None, None)),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(
            &gb(Some("aaa"), Some(1), Some("z")),
            &gb(Some("AAA"), Some(1), Some("a"))
        ),
        Some(Ordering::Equal)
    );
    assert_eq!(
        compare_ordered(
            &gb(Some("AAA"), None, Some("z")),
            &gb(Some("AAA"), Some(1), None)
        ),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(
            &gb(Some("AAA"), None, Some("a")),
            &gb(Some("AAA"), None, Some("B"))
        ),
        Some(Ordering::Less)
    );

    // PDB: case-insensitive mol, case-sensitive effective chain id (`'A' < 'a'`).
    assert_eq!(
        compare_ordered(&pdb("1abc", Some(65), None), &pdb("1ABC", None, Some("A"))),
        Some(Ordering::Equal)
    );
    assert_eq!(
        compare_ordered(&pdb("1ABC", None, Some("A")), &pdb("1ABC", None, Some("a"))),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&pdb("1ABC", None, None), &pdb("1ABC", Some(65), None)),
        Some(Ordering::Less)
    );

    // Patent: seqid order, then equal if `CId_pat::Match`, else pointer-dependent (`None`).
    assert_eq!(
        compare_ordered(&patent(1, "US", "100"), &patent(2, "EP", "100")),
        Some(Ordering::Less)
    );
    assert_eq!(
        compare_ordered(&patent(1, "us", "100"), &patent(1, "US", "100")),
        Some(Ordering::Equal)
    );
    assert_eq!(
        compare_ordered(&patent(1, "US", "100"), &patent(1, "EP", "100")),
        None
    );
}
