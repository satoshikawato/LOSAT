use super::*;

fn no_warning() -> impl FnMut(SeqIdWarning) {
    |warning| panic!("unexpected warning {warning:?}")
}

fn text(s: &str) -> Option<Vec<u8>> {
    Some(s.as_bytes().to_vec())
}

fn textseq(
    acc: Option<&str>,
    name: Option<&str>,
    release: Option<&str>,
    version: Option<i32>,
) -> TextseqId {
    TextseqId {
        name: name.and_then(text),
        accession: acc.and_then(text),
        release: release.and_then(text),
        version,
    }
}

/// The first parse of makeblastdb (`fParse_PartialOK | fParse_AnyLocal | fParse_RawText`).
const FIRST: u32 = PARTIAL_OK | ANY_LOCAL | RAW_TEXT;
/// `ParseFastaIds` of the second parse (`fParse_RawText | fParse_AnyLocal`).
const SECOND_TAGGED: u32 = RAW_TEXT | ANY_LOCAL;
/// The bare id of the second parse (`fParse_RawText | fParse_ValidLocal`).
const SECOND_BARE: u32 = RAW_TEXT | VALID_LOCAL;

fn ids(s: &str) -> Vec<Option<ParsedSeqId>> {
    parse_ids(s.as_bytes(), FIRST, &mut no_warning()).unwrap()
}

// The ids of survey B fixture p01_types as makeblastdb 2.17.0 wrote them (blastdbcmd %i
// and the BER dump of db.phr).
#[test]
fn fasta_ids_of_every_type() {
    use ParsedSeqId::*;
    assert_eq!(
        ids("sp|P69905|HBA_HUMAN"),
        [Some(Swissprot(textseq(
            Some("P69905"),
            Some("HBA_HUMAN"),
            Some("reviewed"),
            None
        )))]
    );
    assert_eq!(
        ids("tr|A0A024R161|A0A024R161_HUMAN"),
        [Some(Swissprot(textseq(
            Some("A0A024R161"),
            Some("A0A024R161_HUMAN"),
            Some("unreviewed"),
            None
        )))]
    );
    assert_eq!(
        ids("gi|123456|ref|NP_000550.2|"),
        [
            Some(Gi(123456)),
            Some(Other(textseq(Some("NP_000550"), None, None, Some(2))))
        ]
    );
    assert_eq!(
        ids("pir||S12345"),
        [Some(Pir(textseq(None, Some("S12345"), None, None)))]
    );
    assert_eq!(
        ids("pdb|1ABC|A"),
        [Some(Pdb(PdbSeqId {
            mol: b"1ABC".to_vec(),
            chain: Some(65),
            rel: None,
            chain_id: text("A")
        }))]
    );
    assert_eq!(
        ids("pdb|2XYZ|AA"),
        [Some(Pdb(PdbSeqId {
            mol: b"2XYZ".to_vec(),
            chain: None,
            rel: None,
            chain_id: text("AA")
        }))]
    );
    assert_eq!(
        ids("pdb|3DEF|"),
        [Some(Pdb(PdbSeqId {
            mol: b"3DEF".to_vec(),
            chain: None,
            rel: None,
            chain_id: None
        }))]
    );
    assert_eq!(ids("lcl|123"), [Some(Local(ObjectId::Id(123)))]);
    assert_eq!(
        ids("lcl|myseq1"),
        [Some(Local(ObjectId::Str(b"myseq1".to_vec())))]
    );
    assert_eq!(
        ids("gnl|mydb|456"),
        [Some(General(Dbtag {
            db: b"mydb".to_vec(),
            tag: ObjectId::Id(456)
        }))]
    );
    let patent = |number: &str, seqid| {
        Some(Patent(PatentSeqId {
            seqid,
            cit: IdPat {
                country: b"US".to_vec(),
                id: PatentNumber::Number(number.as_bytes().to_vec()),
                doc_type: None,
            },
        }))
    };
    assert_eq!(ids("pat|US|1234567|1"), [patent("1234567", 1)]);
    // pgp: x_Init sets the release to "pat", so the number is not an application number.
    assert_eq!(ids("pgp|US|20010012345|2"), [patent("20010012345", 2)]);
    assert_eq!(
        ids("xyz|ABC"),
        [Some(General(Dbtag {
            db: b"xyz".to_vec(),
            tag: ObjectId::Str(b"ABC".to_vec())
        }))]
    );
    // Bare ids, typed by the accession guide.
    assert_eq!(
        ids("XP_011520000.1"),
        [Some(Other(textseq(
            Some("XP_011520000"),
            None,
            None,
            Some(1)
        )))]
    );
    assert_eq!(
        ids("P69906"),
        [Some(Swissprot(textseq(Some("P69906"), None, None, None)))]
    );
    assert_eq!(
        ids("UniRef50_A0A5A9P0L4"),
        [Some(Local(ObjectId::Str(b"UniRef50_A0A5A9P0L4".to_vec())))]
    );
    assert_eq!(ids("123456"), [Some(Local(ObjectId::Id(123456)))]);
    assert_eq!(
        ids("1ABC_A"),
        [Some(Pdb(PdbSeqId {
            mol: b"1ABC".to_vec(),
            chain: Some(65),
            rel: None,
            chain_id: text("A")
        }))]
    );
    assert_eq!(
        ids("SRA:SRR123.1"),
        [Some(General(Dbtag {
            db: b"SRA".to_vec(),
            tag: ObjectId::Str(b"SRR123.1".to_vec())
        }))]
    );
}

// Exceptions caught by fParse_PartialOK, as makeblastdb 2.17.0 posted them (probe runs of
// makedb-v5 S3a and fixtures p16_unsupported_type, p15_bad_version).
#[test]
fn partial_failures_post_the_exception_and_go_on() {
    let mut posted = Vec::new();
    let result = parse_ids(b"xyz|ABC|DEF", FIRST, &mut |w| posted.push(w));
    let error = result.unwrap_err();
    assert_eq!(error.message, "Unsupported ID type DEF");
    let texts: Vec<_> = posted
        .iter()
        .map(|w| (w.subcode, w.message.clone()))
        .collect();
    assert_eq!(texts.len(), 2);
    assert!(texts[0].1.ends_with("\", line 2663: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::ParseIDs() - Unsupported ID type xyz\n"));
    assert!(texts[0].1.starts_with("NCBI C++ Exception:\n    T0 \""));
    assert!(texts[1].1.ends_with("Unsupported ID type ABC\n"));
    assert_eq!(texts[0].0, 7);

    posted.clear();
    let result = parse_ids(b"ref|NP_1.x|gi|5|", FIRST, &mut |w| posted.push(w)).unwrap();
    assert_eq!(
        result,
        [Some(ParsedSeqId::Other(textseq(
            Some("gi"),
            Some("5"),
            None,
            None
        )))]
    );
    assert_eq!(posted.len(), 1);
    assert!(posted[0].message.contains("seqloc/Textseq_id.cpp\", line 108: Error: (CSeqIdException::eFormat) ncbi::objects::CTextseq_id::Set() - Version embedded in accession NP_1.x is not a positive integer"));

    posted.clear();
    assert!(parse_ids(b"ref|NP_800001.x|", FIRST, &mut |w| posted.push(w)).is_err());
    assert!(posted.is_empty());

    posted.clear();
    assert_eq!(
        parse_ids(b"gi|abc|gi|8|", FIRST, &mut |w| posted.push(w)).unwrap(),
        []
    );
    assert_eq!(posted.len(), 2);
    assert!(posted[1].message.contains("line 2941: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::Set() - Negative, excessively large, or non-numeric gi ID gi"));

    posted.clear();
    // The second x_Init reads "gi" as a PDB molecule again and throws with no piece left.
    assert!(parse_ids(b"pdb|12|gi|6|", FIRST, &mut |w| posted.push(w)).is_err());
    assert_eq!(posted.len(), 1);
    assert!(posted[0].message.contains("line 2813: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::x_Init() - Malformatted PDB ID 12"));

    posted.clear();
    assert!(parse_ids(b"pat|US|12|x|gi|7|", FIRST, &mut |w| posted.push(w)).is_err());
    assert!(posted[0].message.contains("line 2799: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::x_Init() - Bad sequence number x for US patent 12"));

    posted.clear();
    let result = parse_ids(b"gb||gi|10|", FIRST, &mut |w| posted.push(w)).unwrap();
    assert!(posted[0].message.contains("Textseq_id.cpp\", line 129: Error: (CSeqIdException::eFormat) ncbi::objects::CTextseq_id::Set() - Accession and name missing for Textseq-id (but got version 0, release )"));
    assert_eq!(
        result,
        [Some(ParsedSeqId::Genbank(textseq(
            Some("gi"),
            Some("10"),
            None,
            None
        )))]
    );
}

#[test]
fn bare_ids_of_the_second_parse() {
    // NC_000913.3:1-5 has a non-numeric version for the guide: a valid local id.
    assert_eq!(
        seq_id_from_text(b"NC_000913.3:1-5", SECOND_BARE, &mut no_warning()).unwrap(),
        ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3:1-5".to_vec()))
    );
    // Not a valid local id: Malformatted ID (fParse_ValidLocal without fParse_AnyLocal).
    let long = "a".repeat(51);
    let error = seq_id_from_text(long.as_bytes(), SECOND_BARE, &mut no_warning()).unwrap_err();
    assert_eq!(error.message, format!("Malformatted ID {long}"));
    assert_eq!(
        seq_id_from_text(b"", SECOND_BARE, &mut no_warning())
            .unwrap_err()
            .message,
        "Empty bare accession supplied"
    );
    assert_eq!(
        seq_id_from_text(b"1234567A", SECOND_BARE, &mut no_warning()).unwrap(),
        ParsedSeqId::Prf(textseq(None, Some("1234567A"), None, None))
    );
    assert!(parse_ids(b"ref|NP_1.x|gi|5|", SECOND_TAGGED, &mut no_warning()).is_err());
}

// CSeq_id_find_pred (Seq_id.cpp:77-86): a control byte or a byte from 0x80 up makes a local
// id invalid, 0x7F does not; a bare id of the second parse with one is then "Malformatted ID"
// (makeblastdb keeps the first parse's ids; review A1, pr2/review/A/out2/q010, a1/dup.fa).
#[test]
fn control_and_non_ascii_bytes_make_a_local_id_invalid() {
    assert!(!is_valid_local_id(b"a\x02b"));
    assert!(!is_valid_local_id(b"a\x1fb"));
    assert!(!is_valid_local_id(b"ab\xc3\xa9"));
    assert!(!is_valid_local_id(b"abc\xe9def"));
    assert!(!is_valid_local_id(b"M\xc3\xbcller_a"));
    assert!(!is_valid_local_id(b"\x80"));
    assert!(!is_valid_local_id(b"\xff"));
    assert!(is_valid_local_id(b"a\x7fb"));
    assert!(is_valid_local_id(b"a!b~"));
    for c in 0u8..=255 {
        let illegal = c < 0x20 || c >= 0x80 || b" >[]|\"".contains(&c);
        assert_eq!(has_invalid_local_char(&[b'a', c]), illegal, "byte {c:#x}");
    }
    for bare in [&b"a\x02b"[..], b"ab\xc3\xa9", b"M\xc3\xbcller_b"] {
        let error = seq_id_from_text(bare, SECOND_BARE, &mut no_warning()).unwrap_err();
        assert!(
            error.message.starts_with("Malformatted ID "),
            "{}",
            error.message
        );
    }
    assert_eq!(
        seq_id_from_text(b"a\x7fb", SECOND_BARE, &mut no_warning()).unwrap(),
        ParsedSeqId::Local(ObjectId::Str(b"a\x7fb".to_vec()))
    );
    // fParse_AnyLocal (ParseFastaIds of an id with '|') does not check the bytes.
    assert_eq!(
        parse_ids(b"lcl|a\x02b", SECOND_TAGGED, &mut no_warning()).unwrap(),
        vec![Some(ParsedSeqId::Local(ObjectId::Str(b"a\x02b".to_vec())))]
    );
}

// Seq_id.cpp:2534-2545: one "Ignoring extra parts" warning for each non-empty piece left
// with fParse_PartialOK (review A4); without it the first one throws.
#[test]
fn every_extra_piece_is_warned_about() {
    let mut warnings = Vec::new();
    let id = seq_id_from_text(
        b"ref|NP_000001.1|NAME|y||z|",
        PARTIAL_OK | RAW_TEXT,
        &mut |w| warnings.push(w),
    )
    .unwrap();
    assert!(matches!(id, ParsedSeqId::Other(_)), "{id:?}");
    assert_eq!(warnings.len(), 2, "{warnings:?}");
    assert!(warnings.iter().all(|w| w.subcode == 10
        && w.message
            == "Ignoring extra parts (synonyms?) in FASTA-style ID ref|NP_000001.1|NAME|y||z|"));
    let error =
        seq_id_from_text(b"ref|NP_000001.1|NAME|y||z|", RAW_TEXT, &mut no_warning()).unwrap_err();
    assert_eq!(
        error.message,
        "FASTA-style ID ref|NP_000001.1|NAME|y||z| has too many parts."
    );
}

#[test]
fn numbers_parse_as_ncbi_s() {
    assert_eq!(string_to_non_negative_int(b"+12"), 12);
    assert_eq!(string_to_non_negative_int(b"+"), -1);
    assert_eq!(string_to_non_negative_int(b"2147483647"), i32::MAX);
    assert_eq!(string_to_non_negative_int(b"2147483648"), -1);
    assert_eq!(string_to_non_negative_int(b"-1"), -1);
    assert_eq!(string_to_int8(b"4294967296"), Some(4_294_967_296));
    assert_eq!(string_to_int8(b"-5"), Some(-5));
    assert_eq!(string_to_int8(b"9223372036854775808"), None);
    assert_eq!(string_to_int8(b"-9223372036854775808"), Some(i64::MIN));
    assert_eq!(str_or_id(b"0123"), ObjectId::Str(b"0123".to_vec()));
    assert_eq!(str_or_id(b"123a"), ObjectId::Str(b"123a".to_vec()));
    assert_eq!(str_or_id(b"99"), ObjectId::Id(99));
}
