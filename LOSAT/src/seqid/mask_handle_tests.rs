//! The identity table of `mask_handle`: every probe of NCBI BLAST+ 2.17.0 (`tests/seqdb/mask_ids/
//! README.md`) as a mask Seq-id, the database Seq-id that NCBI's FASTA parse gives, and whether
//! NCBI matched them.

use super::*;
use crate::seqid::model::{Dbtag, GiimportId, IdPat, ObjectId, PatentSeqId, PdbSeqId};
use ParsedSeqId as P;

fn b(s: &str) -> Vec<u8> {
    s.as_bytes().to_vec()
}
/// Textseq-id: accession, version, name, release ("" = unset).
fn tx(acc: &str, ver: Option<i32>, name: &str, rel: &str) -> TextseqId {
    let opt = |s: &str| (!s.is_empty()).then(|| b(s));
    TextseqId {
        accession: opt(acc),
        version: ver,
        name: opt(name),
        release: opt(rel),
    }
}
fn lstr(s: &str) -> ParsedSeqId {
    P::Local(ObjectId::Str(b(s)))
}
fn gnl(db: &str, tag: ObjectId) -> ParsedSeqId {
    P::General(Dbtag { db: b(db), tag })
}
fn pat(seqid: i32, country: &str, number: PatentNumber, doc_type: Option<&str>) -> ParsedSeqId {
    P::Patent(PatentSeqId {
        seqid,
        cit: IdPat {
            country: b(country),
            id: number,
            doc_type: doc_type.map(b),
        },
    })
}
fn pdb(mol: &str, chain: Option<i32>, rel: Option<Date>, chain_id: Option<&str>) -> ParsedSeqId {
    P::Pdb(PdbSeqId {
        mol: b(mol),
        chain,
        rel,
        chain_id: chain_id.map(b),
    })
}
fn year(year: i32, month: Option<i32>, season: Option<&str>) -> Option<Date> {
    Some(Date::Std(DateStd {
        year,
        month,
        day: None,
        season: season.map(b),
        hour: None,
        minute: None,
        second: None,
    }))
}
fn giim(id: i32, db: Option<&str>) -> ParsedSeqId {
    P::Giim(GiimportId {
        id,
        db: db.map(b),
        release: None,
    })
}

#[derive(Debug, PartialEq)]
enum E {
    Same,
    Apart,
    /// Found in the Giim tree (`range_set`), not by a key.
    Tree,
}

/// (probe case, mask Seq-id, database Seq-id as NCBI's FASTA parse gives it, NCBI 2.17.0).
fn table() -> Vec<(&'static str, ParsedSeqId, ParsedSeqId, E)> {
    let gb_ay100001 = P::Genbank(tx("AY100001", Some(1), "", ""));
    let sp = P::Swissprot(tx("P12345", Some(1), "NAME_HUMAN", "reviewed"));
    let tr = P::Swissprot(tx("Q9XYZ1", Some(1), "Q9XYZ1_HUMAN", "unreviewed"));
    let pdb_1abc_a = pdb("1ABC", Some(65), None, Some("A"));
    vec![
        ("v01", P::Gi(123), P::Gi(123), E::Same),
        ("v27", P::Gi(124), P::Gi(123), E::Apart),
        (
            "v02",
            P::Other(tx("NC_000002", Some(2), "", "")),
            P::Other(tx("NC_000002", Some(2), "", "")),
            E::Same,
        ),
        (
            "v03",
            P::Other(tx("nc_000002", Some(2), "", "")),
            P::Other(tx("NC_000002", Some(2), "", "")),
            E::Same,
        ),
        (
            "v04",
            P::Other(tx("NC_000002", None, "", "")),
            P::Other(tx("NC_000002", Some(2), "", "")),
            E::Apart,
        ),
        (
            "v12",
            P::Other(tx("NC_000003", Some(1), "", "")),
            P::Other(tx("NC_000003", Some(1), "foo", "")),
            E::Apart,
        ),
        (
            "v13",
            P::Other(tx("NC_000003", Some(1), "foo", "")),
            P::Other(tx("NC_000003", Some(1), "foo", "")),
            E::Same,
        ),
        (
            "n01",
            P::Embl(tx("AY100001", Some(1), "", "")),
            gb_ay100001.clone(),
            E::Same,
        ),
        (
            "n02",
            P::Ddbj(tx("ay100001", Some(1), "", "")),
            gb_ay100001.clone(),
            E::Same,
        ),
        (
            "n03",
            P::Genbank(tx("AY100001", Some(1), "", "x")),
            gb_ay100001.clone(),
            E::Apart,
        ),
        (
            "n04",
            P::Embl(tx("AY000000", Some(1), "", "")),
            P::Genbank(tx("AY000000", Some(1), "", "")),
            E::Apart,
        ),
        (
            "n05",
            P::Genbank(tx("ay000000", Some(1), "", "")),
            P::Genbank(tx("AY000000", Some(1), "", "")),
            E::Same,
        ),
        (
            "n06",
            P::Embl(tx("ABCDE12345", Some(1), "", "")),
            P::Genbank(tx("ABCDE12345", Some(1), "", "")),
            E::Apart,
        ),
        (
            "n07",
            P::Genbank(tx("abcde12345", Some(1), "", "")),
            P::Genbank(tx("ABCDE12345", Some(1), "", "")),
            E::Same,
        ),
        (
            "n08",
            P::Embl(tx("AY100002", Some(1), "LOC2", "")),
            P::Genbank(tx("AY100002", Some(1), "LOC2", "")),
            E::Apart,
        ),
        (
            "n09",
            P::Genbank(tx("ay100002", Some(1), "loc2", "")),
            P::Genbank(tx("AY100002", Some(1), "LOC2", "")),
            E::Same,
        ),
        (
            "n10",
            P::Genbank(tx("AY100002", Some(1), "", "")),
            P::Genbank(tx("AY100002", Some(1), "LOC2", "")),
            E::Apart,
        ),
        (
            "n11",
            P::Genbank(tx("AJ100001", Some(1), "", "")),
            P::Embl(tx("AJ100001", Some(1), "", "")),
            E::Same,
        ),
        (
            "n12",
            P::Embl(tx("AY100003", None, "", "")),
            P::Genbank(tx("AY100003", None, "", "")),
            E::Same,
        ),
        (
            "n13",
            P::Embl(tx("AY100003", Some(1), "", "")),
            P::Genbank(tx("AY100003", None, "", "")),
            E::Apart,
        ),
        (
            "n14",
            P::Genbank(tx("BK100001", Some(1), "", "")),
            P::Tpg(tx("BK100001", Some(1), "", "")),
            E::Apart,
        ),
        (
            "n15",
            P::Tpg(tx("bk100001", Some(1), "", "")),
            P::Tpg(tx("BK100001", Some(1), "", "")),
            E::Same,
        ),
        (
            "n16",
            P::Tpe(tx("BK100001", Some(1), "", "")),
            P::Tpg(tx("BK100001", Some(1), "", "")),
            E::Apart,
        ),
        (
            "n17",
            P::Other(tx("nc_000000", Some(1), "", "")),
            P::Other(tx("NC_000000", Some(1), "", "")),
            E::Same,
        ),
        (
            "n18",
            P::Embl(tx("ay0123456", Some(1), "", "")),
            P::Genbank(tx("AY0123456", Some(1), "", "")),
            E::Same,
        ),
        (
            "n19",
            P::Embl(tx("AY1234567890", Some(1), "", "")),
            P::Genbank(tx("AY1234567890", Some(1), "", "")),
            E::Apart,
        ),
        (
            "n20",
            P::Genbank(tx("ay1234567890", Some(1), "", "")),
            P::Genbank(tx("AY1234567890", Some(1), "", "")),
            E::Same,
        ),
        (
            "n21",
            P::Genbank(tx("AY100001", Some(2), "", "")),
            gb_ay100001.clone(),
            E::Apart,
        ),
        (
            "n22",
            P::Genbank(tx("", None, "LOC2", "")),
            P::Genbank(tx("AY100002", Some(1), "LOC2", "")),
            E::Apart,
        ),
        (
            "v10",
            P::Local(ObjectId::Id(5)),
            P::Local(ObjectId::Id(5)),
            E::Same,
        ),
        ("v11", lstr("5"), P::Local(ObjectId::Id(5)), E::Apart),
        ("n23", P::Local(ObjectId::Id(5)), lstr("05"), E::Apart),
        ("n24", lstr("05"), lstr("05"), E::Same),
        ("v06", lstr("ABC"), lstr("abc"), E::Same),
        ("n26", lstr("caa12345.1"), lstr("CAA12345.1"), E::Same),
        (
            "n27",
            P::Embl(tx("CAA12345", Some(1), "", "")),
            lstr("CAA12345.1"),
            E::Apart,
        ),
        (
            "v08",
            gnl("MYDB", ObjectId::Str(b("tag1"))),
            gnl("mydb", ObjectId::Str(b("tag1"))),
            E::Same,
        ),
        (
            "v09",
            gnl("mydb", ObjectId::Str(b("TAG1"))),
            gnl("mydb", ObjectId::Str(b("tag1"))),
            E::Same,
        ),
        (
            "v22",
            gnl("db2", ObjectId::Id(123)),
            gnl("db2", ObjectId::Id(123)),
            E::Same,
        ),
        (
            "v23",
            gnl("db2", ObjectId::Str(b("123"))),
            gnl("db2", ObjectId::Id(123)),
            E::Apart,
        ),
        (
            "n28",
            gnl("DB3", ObjectId::Str(b("abc12345"))),
            gnl("db3", ObjectId::Str(b("ABC12345"))),
            E::Same,
        ),
        (
            "n29",
            gnl("Db3", ObjectId::Str(b("ab12"))),
            gnl("db3", ObjectId::Str(b("AB12"))),
            E::Same,
        ),
        (
            "n30",
            gnl("db3", ObjectId::Id(123)),
            gnl("db3", ObjectId::Str(b("00123"))),
            E::Apart,
        ),
        (
            "n31",
            gnl("DB3", ObjectId::Str(b("00123"))),
            gnl("db3", ObjectId::Str(b("00123"))),
            E::Same,
        ),
        (
            "n32",
            gnl("db3", ObjectId::Str(b("ABC012345"))),
            gnl("db3", ObjectId::Str(b("ABC12345"))),
            E::Apart,
        ),
        ("n33", P::Gibbsq(123), P::Gibbsq(123), E::Same),
        ("n34", P::Gibbmt(123), P::Gibbsq(123), E::Apart),
        ("n35", P::Gibbmt(77), P::Gibbmt(77), E::Same),
        ("n36", giim(7, None), giim(7, None), E::Tree),
        ("n37", giim(7, Some("x")), giim(7, None), E::Tree),
        ("n39", giim(8, Some("x")), giim(7, None), E::Tree),
        (
            "n46",
            pat(1, "us", PatentNumber::Number(b("re33188")), None),
            pat(1, "US", PatentNumber::Number(b("RE33188")), None),
            E::Same,
        ),
        (
            "n47",
            pat(1, "US", PatentNumber::Number(b("RE33188")), Some("A")),
            pat(1, "US", PatentNumber::Number(b("RE33188")), None),
            E::Same,
        ),
        (
            "n48",
            pat(1, "US", PatentNumber::AppNumber(b("RE33188")), None),
            pat(1, "US", PatentNumber::Number(b("RE33188")), None),
            E::Apart,
        ),
        (
            "n49",
            pat(2, "US", PatentNumber::Number(b("RE33188")), None),
            pat(1, "US", PatentNumber::Number(b("RE33188")), None),
            E::Apart,
        ),
        (
            "n50b",
            pat(7, "ep", PatentNumber::Number(b("0238993")), None),
            pat(7, "EP", PatentNumber::Number(b("0238993")), None),
            E::Same,
        ),
        (
            "p01",
            P::Swissprot(tx("P12345", Some(1), "NAME_HUMAN", "")),
            sp.clone(),
            E::Apart,
        ),
        (
            "p02",
            P::Swissprot(tx("P12345", Some(1), "NAME_HUMAN", "reviewed")),
            sp.clone(),
            E::Same,
        ),
        (
            "p03",
            P::Swissprot(tx("p12345", Some(1), "name_human", "reviewed")),
            sp.clone(),
            E::Same,
        ),
        (
            "p04",
            P::Swissprot(tx("P12345", Some(1), "NAME_HUMAN", "Reviewed")),
            sp.clone(),
            E::Apart,
        ),
        (
            "p05",
            P::Swissprot(tx("P12345", Some(1), "NAME_HUMAN", "unreviewed")),
            sp.clone(),
            E::Apart,
        ),
        (
            "p06",
            P::Swissprot(tx("P12345", Some(1), "", "")),
            sp.clone(),
            E::Apart,
        ),
        (
            "p07",
            P::Swissprot(tx("Q9XYZ1", Some(1), "Q9XYZ1_HUMAN", "unreviewed")),
            tr.clone(),
            E::Same,
        ),
        (
            "p08",
            P::Swissprot(tx("Q9XYZ1", Some(1), "Q9XYZ1_HUMAN", "reviewed")),
            tr.clone(),
            E::Apart,
        ),
        (
            "p09",
            P::Pir(tx("", None, "s12345", "")),
            P::Pir(tx("", None, "S12345", "")),
            E::Same,
        ),
        (
            "p10",
            P::Pir(tx("S12345", None, "", "")),
            P::Pir(tx("", None, "S12345", "")),
            E::Apart,
        ),
        (
            "p11",
            P::Prf(tx("", None, "1234567a", "")),
            P::Prf(tx("", None, "1234567A", "")),
            E::Same,
        ),
        (
            "p35",
            P::Embl(tx("aab10001", Some(1), "", "")),
            P::Genbank(tx("AAB10001", Some(1), "", "")),
            E::Same,
        ),
        (
            "v19",
            pdb("1ABC", None, None, Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "v20",
            pdb("1abc", None, None, Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p12",
            pdb("1abc", Some(65), None, None),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p13",
            pdb("1ABC", Some(66), None, Some("A")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p14",
            pdb("1abc", Some(66), None, Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p15",
            pdb("1ABC", Some(321), None, None),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p16",
            pdb("1ABC", Some(-191), None, None),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p17",
            pdb("1ABC", None, year(2000, None, None), Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p18",
            pdb("1ABC", None, Some(Date::Str(b("2000"))), Some("A")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p19",
            pdb("1ABC", None, year(2000, None, Some("spring")), Some("A")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p20",
            pdb("1ABC", None, year(2000, Some(13), None), Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p21",
            pdb("1ABC", None, year(2000, Some(16), None), Some("A")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p22",
            pdb("1ABC", None, year(0, None, None), Some("A")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p23",
            pdb("1ABC", None, None, Some("a")),
            pdb_1abc_a.clone(),
            E::Apart,
        ),
        (
            "p24",
            pdb("1ABC", Some(66), year(1999, None, None), Some("A")),
            pdb_1abc_a.clone(),
            E::Same,
        ),
        (
            "p25",
            pdb("abcd", None, None, Some("B")),
            pdb("ABCD", Some(66), None, Some("B")),
            E::Apart,
        ),
        (
            "p26",
            pdb("Abcd", None, None, Some("B")),
            pdb("ABCD", Some(66), None, Some("B")),
            E::Same,
        ),
        (
            "p27",
            pdb("1xyzw", None, None, Some("C")),
            pdb("1XYZW", Some(67), None, Some("C")),
            E::Apart,
        ),
        (
            "p28",
            pdb("1xyzW", None, None, Some("C")),
            pdb("1XYZW", Some(67), None, Some("C")),
            E::Same,
        ),
        (
            "p29",
            pdb("2ABC", None, None, Some("AB")),
            pdb("2ABC", None, None, Some("AB")),
            E::Same,
        ),
        (
            "p30",
            pdb("2ABC", Some(65), None, Some("AB")),
            pdb("2ABC", None, None, Some("AB")),
            E::Apart,
        ),
        (
            "p31",
            pdb("2abc", None, None, Some("AB")),
            pdb("2ABC", None, None, Some("AB")),
            E::Same,
        ),
        (
            "p32",
            pdb("3ABC", None, None, None),
            pdb("3ABC", None, None, None),
            E::Same,
        ),
        (
            "p33",
            pdb("3ABC", Some(32), None, None),
            pdb("3ABC", None, None, None),
            E::Apart,
        ),
        (
            "p34",
            pdb("3abc", None, None, None),
            pdb("3ABC", None, None, None),
            E::Same,
        ),
        // source only (no FASTA form): general integer tags around 0, GI 0
        (
            "src-gnl0",
            gnl("db", ObjectId::Id(0)),
            gnl("DB", ObjectId::Id(0)),
            E::Same,
        ),
        (
            "src-gnl0-1",
            gnl("db", ObjectId::Id(0)),
            gnl("db", ObjectId::Id(-1)),
            E::Apart,
        ),
        ("src-gi0", P::Gi(0), P::Gi(0), E::Same),
    ]
}

fn relation(mask: &ParsedSeqId, db: &ParsedSeqId) -> E {
    match (mask_handle_key(mask), mask_handle_key(db)) {
        (Err(_), _) | (_, Err(_)) => E::Tree,
        (Ok(Some(a)), Ok(Some(b))) if a == b => E::Same,
        _ => E::Apart,
    }
}

#[test]
fn ncbi_handle_table() {
    for (case, mask, db, expected) in table() {
        assert_eq!(relation(&mask, &db), expected, "{case}: {mask:?} vs {db:?}");
    }
}

#[test]
fn neither_accession_nor_name_is_a_new_handle() {
    let id = P::Genbank(tx("", Some(1), "", "x"));
    assert_eq!(mask_handle_key(&id), Ok(None));
}

#[test]
fn handle_order_packed() {
    // n56-n58: gi 123 (m_Packed 123) after ref|NC_000001.1| (m_Packed 1);
    // n59-n62: gi 5 = AY000005.1 (5 = 5) and local/general strings (0 = 0) tie: address order.
    let p = |id: &ParsedSeqId| mask_handle_packed(id);
    assert_eq!(p(&P::Gi(123)), 123);
    assert_eq!(p(&P::Other(tx("NC_000001", Some(1), "", ""))), 1);
    assert_eq!(p(&P::Genbank(tx("AY000005", Some(1), "", ""))), 5);
    assert_eq!(p(&P::Genbank(tx("AY0123456", Some(1), "", ""))), 123456);
    assert_eq!(p(&lstr("abc")), 0);
    assert_eq!(p(&gnl("db9", ObjectId::Str(b("xy")))), 0);
    assert_eq!(p(&gnl("db3", ObjectId::Str(b("ABC12345")))), 12345);
    assert_eq!(p(&gnl("db3", ObjectId::Str(b("x000")))), -1);
    assert_eq!(p(&gnl("db", ObjectId::Id(0))), -1);
    assert_eq!(p(&gnl("db", ObjectId::Id(7))), 7);
}
#[test]
fn unpacked_handles_come_last() {
    // `m_Packed - 1` as an unsigned number: ref|NC_000001.1| (1) before gi|123, both before the
    // unpacked ones (0); a negative packed value (a general tag at or below 0 packs to id - 1, a
    // zero digit value to -1) wraps to the top of the range, just below the unpacked ones.
    assert!(mask_handle_order(1) < mask_handle_order(123));
    assert!(mask_handle_order(123) < mask_handle_order(-1));
    assert!(mask_handle_order(-1) < mask_handle_order(0));
}

/// `AsString()` of the handle of `id` whose info `first` made, as NCBI's `duplicate Bioseq id`
/// exception writes it (NCBI BLAST+ 2.17.0 psiblast `-in_pssm`, line RG D3-1 fix round 5
/// `r5/out_ids5/ncbi/`: the case of each pair, `one_*`).
fn text(first: &ParsedSeqId, id: &ParsedSeqId) -> String {
    assert_eq!(mask_handle_key(first), mask_handle_key(id), "one handle");
    let mut out = Vec::new();
    crate::seqid::label::fasta(&handle_seq_id(first, id), &mut out);
    String::from_utf8(out).unwrap()
}

#[test]
fn the_handles_seq_id_is_its_infos_with_the_case_of_its_own() {
    // lcl_case, gnl_db_case, gnl_tag_case, gnl_id_dbcase, gnl_packed_case: the letters' case.
    assert_eq!(text(&lstr("Query_1"), &lstr("QUERY_1")), "lcl|QUERY_1");
    let tag = |s: &str| ObjectId::Str(b(s));
    assert_eq!(text(&gnl("DB", tag("x")), &gnl("db", tag("x"))), "gnl|db|x");
    assert_eq!(
        text(&gnl("DB", tag("Xy")), &gnl("DB", tag("xY"))),
        "gnl|DB|xY"
    );
    assert_eq!(
        text(&gnl("Db", ObjectId::Id(7)), &gnl("dB", ObjectId::Id(7))),
        "gnl|dB|7"
    );
    assert_eq!(
        text(&gnl("DB", tag("ABC123456")), &gnl("db", tag("abc123456"))),
        "gnl|db|abc123456"
    );
    // lcl_long_case_tail: past the 64th letter the info's case stays; gnl_long_case_tail counts
    // the database's letters first.
    let long = "Abcdefghij".repeat(7);
    let tail = format!("{}{}", &long[..66], long[66..].to_ascii_uppercase());
    let all = long.to_ascii_uppercase();
    assert_eq!(text(&lstr(&long), &lstr(&tail)), format!("lcl|{long}"));
    assert_eq!(
        text(&lstr(&long), &lstr(&all)),
        format!("lcl|{}{}", &all[..64], &long[64..])
    );
    assert_eq!(
        text(&gnl("DB", tag(&long)), &gnl("DB", tag(&tail))),
        format!("gnl|DB|{long}")
    );
    // gb_acc_case, gb_ddbj_case, emb_gb, tpg_case, ref_case, sp_acc_only_case: a packed
    // accession takes its own case and the info's alternative.
    let acc = |a: &str| tx(a, Some(1), "", "");
    assert_eq!(
        text(&P::Genbank(acc("AAA00001")), &P::Ddbj(acc("aAa00001"))),
        "gb|aAa00001.1|"
    );
    assert_eq!(
        text(&P::Embl(acc("AAA00001")), &P::Genbank(acc("AAA00001"))),
        "emb|AAA00001.1|"
    );
    assert_eq!(
        text(&P::Other(acc("NP_000549")), &P::Other(acc("np_000549"))),
        "ref|np_000549.1|"
    );
    // gb_plain_case: a plain accession takes its own case; gb_name_case, sp_name_case,
    // prf_name_case: a name keeps the info's.
    assert_eq!(
        text(&P::Genbank(acc("ABCDEFG1")), &P::Genbank(acc("abcdefg1"))),
        "gb|abcdefg1.1|"
    );
    assert_eq!(
        text(
            &P::Genbank(tx("", None, "ABC", "")),
            &P::Genbank(tx("", None, "abc", ""))
        ),
        "gb||ABC"
    );
    assert_eq!(
        text(
            &P::Swissprot(tx("P69905", Some(2), "HBA_HUMAN", "")),
            &P::Swissprot(tx("P69905", Some(2), "hba_human", ""))
        ),
        "sp|P69905.2|HBA_HUMAN"
    );
    // pat_country_case, pat_number_case, pat_doctype: the info's patent.
    let number = |n: &str| PatentNumber::Number(b(n));
    assert_eq!(
        text(
            &pat(1, "US", number("123"), None),
            &pat(1, "us", number("123"), None)
        ),
        "pat|US|123|1"
    );
    assert_eq!(
        text(
            &pat(1, "US", number("RE33188"), None),
            &pat(1, "US", number("re33188"), None)
        ),
        "pat|US|RE33188|1"
    );
    // A PDB id is its own (review 4c `ps.one_pdb_molcase`: `pdb|1abc|A | pdb|1ABC|A`).
    assert_eq!(
        text(
            &pdb("1ABC", Some(65), None, None),
            &pdb("1abc", Some(65), None, None)
        ),
        "pdb|1abc|A"
    );
}

// The packed tree's info (seq_id_tree.hpp:359-388): one for every accession of a prefix (without
// case), digit count and version, across GenBank, EMBL and DDBJ; any other Seq-id's info is its
// handle's.
#[test]
fn a_packed_accession_shares_the_info_of_its_prefix_digits_and_version() {
    let info = |id: &ParsedSeqId| mask_info_key(id).unwrap().unwrap();
    let gb1 = P::Genbank(tx("AB000001", Some(1), "", ""));
    assert_eq!(info(&gb1), info(&P::Ddbj(tx("ab000002", Some(1), "", ""))));
    assert_eq!(info(&gb1), info(&P::Embl(tx("AB999999", Some(1), "", ""))));
    assert_ne!(
        info(&gb1),
        info(&P::Genbank(tx("AB000001", Some(2), "", "")))
    );
    assert_ne!(info(&gb1), info(&P::Genbank(tx("AB000001", None, "", ""))));
    assert_ne!(
        info(&gb1),
        info(&P::Genbank(tx("AB1000001", Some(1), "", "")))
    );
    assert_ne!(
        info(&gb1),
        info(&P::Genbank(tx("AC000001", Some(1), "", "")))
    );
    assert_ne!(info(&gb1), info(&P::Tpg(tx("AB000001", Some(1), "", ""))));
    for id in [
        P::Swissprot(tx("P12345", None, "ABC_HUMAN", "")),
        lstr("rec1"),
        P::Gi(5),
        gnl("db", ObjectId::Id(7)),
    ] {
        assert_eq!(mask_info_key(&id), mask_handle_key(&id), "{id:?}");
    }
}

#[test]
fn which_handles_write_their_own_seq_id() {
    let acc = |a: &str| tx(a, Some(1), "", "");
    let own = [
        lstr("Query_1"),
        lstr(&"a".repeat(64)),
        gnl("DB", ObjectId::Str(b("x"))),
        P::Gi(5),
        P::Gibbsq(7),
        giim(5, Some("x")),
        pdb("1abc", Some(65), None, None),
        P::Other(acc("NP_000549")),
        P::Tpg(acc("DAA00001")),
        P::Genbank(acc("ABCDEFG1")),
        P::Genbank(tx("", Some(1), "", "")),
        pat(1, "12", PatentNumber::Number(b("123")), None),
    ];
    for id in &own {
        assert!(handle_fasta_is_own(id), "{id:?}");
    }
    let shown = [
        lstr(&"a".repeat(65)),
        gnl("DB", ObjectId::Str(b(&"a".repeat(63)))),
        P::Genbank(acc("AAA00001")),
        P::Embl(acc("AAA00001")),
        P::Ddbj(acc("AAA00001")),
        P::Genbank(tx("", None, "ABC", "")),
        P::Swissprot(tx("P69905", Some(2), "HBA_HUMAN", "")),
        pat(1, "US", PatentNumber::Number(b("123")), None),
    ];
    for id in &shown {
        assert!(!handle_fasta_is_own(id), "{id:?}");
    }
}

#[test]
fn the_giim_tree_finds_the_first_handle_not_equal() {
    // seq_id_tree.cpp:2926-2941: two equal `giim` ids are two handles, another `db` finds the first.
    let id = |db: &str| GiimportId {
        id: 5,
        db: Some(b(db)),
        release: None,
    };
    let mut tree = GiimTree::default();
    let (a, new_a) = tree.find_or_create(&id("x"));
    let (b2, new_b) = tree.find_or_create(&id("x"));
    let (c, new_c) = tree.find_or_create(&id("y"));
    assert!(new_a && new_b && !new_c);
    assert_ne!(a, b2);
    assert_eq!(c, a);
    assert_eq!(tree.made_from(5, c), &id("x"));
    tree.drop_handle(5, a);
    assert_eq!(tree.find(&id("y")), Some(b2));
}
