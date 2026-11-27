use super::*;
use crate::report::serial::asn_types::members;
use crate::seqid::model::{ObjectId, TextseqId};

fn parsed(defline: &[u8], ids: Vec<ParsedSeqId>, title: Option<&[u8]>) -> ParsedIds {
    // The DefLine field encodes the defline as an ASN.1 `PrintableString`; for these lines (no
    // backslash, no control byte except ^A and tab) that is the line with \1 and \t escapes.
    let mut field = Vec::new();
    for &c in defline {
        match c {
            0x01 => field.extend_from_slice(b"\\1"),
            b'\t' => field.extend_from_slice(b"\\t"),
            _ => field.push(c),
        }
    }
    ParsedIds {
        ids,
        hist_old_id: None,
        defline_field: field,
        title_descr: title.map(<[u8]>::to_vec),
    }
}

fn refseq(acc: &str, version: i32) -> ParsedSeqId {
    ParsedSeqId::Other(TextseqId {
        accession: Some(acc.as_bytes().to_vec()),
        version: Some(version),
        ..TextseqId::default()
    })
}

fn titles_and_ids(set: &DeflineSet) -> Vec<(Vec<u8>, Vec<ParsedSeqId>)> {
    set.deflines
        .iter()
        .map(|d| (d.title.clone().unwrap(), d.ids.clone()))
        .collect()
}

// survey B fixture p02_multi (makeblastdb 2.17.0, blastdbcmd %i|%t).
#[test]
fn the_second_parse_splits_at_ctrl_a_and_keeps_the_spaces_of_titles() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    let p = parsed(
        b">ref|NP_100003.1| first title\x01ref|NP_100004.1| second title\x01gi|1005| third title",
        vec![refseq("NP_100003", 1)],
        Some(b"first title"),
    );
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(
        titles_and_ids(&set),
        [
            (b"first title".to_vec(), vec![refseq("NP_100003", 1)]),
            (b"second title".to_vec(), vec![refseq("NP_100004", 1)]),
            (b"third title".to_vec(), vec![ParsedSeqId::Gi(1005)]),
        ]
    );
    let p = parsed(
        b">ref|NP_100011.1|\ttab right after id",
        vec![refseq("NP_100011", 1)],
        None,
    );
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(titles_and_ids(&set)[0].0, b"  tab right after id");
    // A space right after '>': the id is empty, NCBI's parse throws.
    let p = parsed(
        b">  ref|NP_100008.1|  leading spaces",
        vec![refseq("NP_100008", 1)],
        Some(b"leading spaces"),
    );
    assert!(second_parse(&p, true, &mut no_warn).is_none());
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: true,
        length: 1,
        seq_data: None,
        parsed: Some(p),
        delta_fix: None,
    };
    let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
    assert_eq!(
        titles_and_ids(&set),
        [(b"leading spaces".to_vec(), vec![refseq("NP_100008", 1)])]
    );
}

// survey B fixture p04_molmismatch and the probe of a nucleotide range.
#[test]
fn ids_of_the_other_molecule_become_local() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    let p = parsed(b">NC_000913.3 raw nucleotide", vec![], None);
    let set = second_parse(&p, true, &mut no_warn).unwrap();
    assert_eq!(
        set.deflines[0].ids,
        [ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))]
    );
    let mut ids = vec![refseq("NC_000913", 3), ParsedSeqId::Gi(5)];
    localize_other_molecule(&mut ids, true, &mut no_warn);
    assert_eq!(
        ids[0],
        ParsedSeqId::Local(ObjectId::Str(b"NC_000913.3".to_vec()))
    );
    assert_eq!(ids[1], ParsedSeqId::Gi(5));
    // A nucleotide id with a range: the Bioseq has lcl|3 and the replaced ref id in
    // its history, so the local second parse is not used.
    let mut p = parsed(
        b">NC_000913.3:1-5 t",
        vec![ParsedSeqId::Local(ObjectId::Id(3))],
        Some(b"t"),
    );
    p.hist_old_id = Some(refseq("NC_000913", 3));
    let bioseq = Bioseq {
        local_id: String::new(),
        title: Vec::new(),
        is_aa: false,
        length: 1,
        seq_data: None,
        parsed: Some(p),
        delta_fix: None,
    };
    let set = extract_bioseq_deflines(&bioseq, &mut no_warn);
    assert_eq!(
        titles_and_ids(&set),
        [(b"t".to_vec(), vec![ParsedSeqId::Local(ObjectId::Id(3))])]
    );
}

#[test]
fn large_gis_skip_the_sequence_and_log_the_last_gi() {
    let mut no_warn = |_: &[u8]| panic!("no warning");
    assert_eq!(
        large_gi_to_skip(&[ParsedSeqId::Gi(4_294_967_295)], &mut no_warn),
        None
    );
    assert_eq!(
        large_gi_to_skip(
            &[ParsedSeqId::Gi(4_294_967_296), ParsedSeqId::Gi(7)],
            &mut no_warn
        ),
        Some(7)
    );
}

// makeprofiledb's query with a descriptor stated as NULL (`Seqdesc::Empty`, line RG item D3c-3b's
// reader): NCBI's `x_GetBioseqBinaryHeader` dereferences each descriptor until it finds the binary
// lines, so without them the definition lines stop at the empty one wherever it stands
// (makeprofiledb 2.17.0, probe D3c-3c `probeB` B01-B03), and so does the taxid loop before the
// first `org` (`first_org_taxid`).
#[test]
fn a_descriptor_stated_as_null_makes_no_definition_lines() {
    for descr in [
        vec![Seqdesc::Empty],
        vec![Seqdesc::Title(b"t".to_vec()), Seqdesc::Empty],
        vec![Seqdesc::Empty, Seqdesc::Title(b"t".to_vec())],
    ] {
        assert_eq!(
            bioseq_deflines(&[], &descr).result,
            Err(BioseqDeflinesError::NullDescriptor),
            "{descr:?}"
        );
        assert_eq!(
            first_org_taxid(&descr),
            Err(BioseqDeflinesError::NullDescriptor),
            "{descr:?}"
        );
    }
    assert_eq!(first_org_taxid(&[Seqdesc::Title(b"t".to_vec())]), Ok(None));
}

// Differential check: the table-driven BER decoder (`decode_binary_deflines_with_names`, root
// `Blast-def-line-set`) matches `seqid::ber::decode_full_set` on valid headers across every
// Seq-id alternative, optional field, empty list, definite/indefinite length, and replaced byte.
#[test]
fn table_driven_blast_def_line_set_matches_ber_decoder() {
    use crate::seqid::model::{
        Date, DateStd, Dbtag, GiimportId, PatentNumber, PatentSeqId, PdbSeqId,
    };
    use crate::writedb::header::encode_defline_set;

    let sample = DeflineSet {
        deflines: vec![
            BlastDefline {
                title: Some(b"first title with ASCII".to_vec()),
                ids: vec![
                    ParsedSeqId::Local(ObjectId::Id(42)),
                    ParsedSeqId::Local(ObjectId::Str(b"loc_str".to_vec())),
                    ParsedSeqId::Gi(123_456),
                    ParsedSeqId::Gibbsq(7),
                    ParsedSeqId::Gibbmt(8),
                    ParsedSeqId::Giim(GiimportId {
                        id: 9,
                        db: Some(b"db".to_vec()),
                        release: Some(b"rel".to_vec()),
                    }),
                    ParsedSeqId::Genbank(TextseqId {
                        name: Some(b"NM1".to_vec()),
                        accession: Some(b"U12345".to_vec()),
                        release: Some(b"r1".to_vec()),
                        version: Some(2),
                    }),
                    ParsedSeqId::Embl(TextseqId {
                        accession: Some(b"X54321".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Ddbj(TextseqId {
                        accession: Some(b"D00001".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    refseq("NP_000001", 3),
                    ParsedSeqId::Pir(TextseqId {
                        name: Some(b"P1".to_vec()),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Prf(TextseqId {
                        name: Some(b"PR1".to_vec()),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Swissprot(TextseqId {
                        name: Some(b"Q1_HUMAN".to_vec()),
                        accession: Some(b"P12345".to_vec()),
                        release: Some(b"reviewed".to_vec()),
                        version: Some(1),
                    }),
                    ParsedSeqId::Tpg(TextseqId {
                        accession: Some(b"BK000001".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Tpe(TextseqId {
                        accession: Some(b"BN000001".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Tpd(TextseqId {
                        accession: Some(b"BR000001".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::Gpipe(TextseqId {
                        accession: Some(b"GP_1".to_vec()),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::NamedAnnotTrack(TextseqId {
                        accession: Some(b"NA000001".to_vec()),
                        version: Some(1),
                        ..TextseqId::default()
                    }),
                    ParsedSeqId::General(Dbtag {
                        db: b"MYDB".to_vec(),
                        tag: ObjectId::Str(b"k1".to_vec()),
                    }),
                    ParsedSeqId::General(Dbtag {
                        db: b"MYDB".to_vec(),
                        tag: ObjectId::Id(11),
                    }),
                    ParsedSeqId::Patent(PatentSeqId {
                        seqid: 1,
                        cit: crate::seqid::model::IdPat {
                            country: b"US".to_vec(),
                            id: PatentNumber::Number(b"1234567".to_vec()),
                            doc_type: Some(b"A".to_vec()),
                        },
                    }),
                    ParsedSeqId::Patent(PatentSeqId {
                        seqid: 2,
                        cit: crate::seqid::model::IdPat {
                            country: b"EP".to_vec(),
                            id: PatentNumber::AppNumber(b"9999".to_vec()),
                            doc_type: None,
                        },
                    }),
                    ParsedSeqId::Pdb(PdbSeqId {
                        mol: b"1ABC".to_vec(),
                        chain: Some(65),
                        rel: Some(Date::Std(DateStd {
                            year: 2020,
                            month: Some(5),
                            day: Some(10),
                            season: None,
                            hour: None,
                            minute: None,
                            second: None,
                        })),
                        chain_id: Some(b"AA".to_vec()),
                    }),
                ],
                taxid: Some(9606),
                memberships: Some(vec![1, 4]),
                links: Some(vec![10090, 9606]),
                other_info: Some(vec![777]),
                ..BlastDefline::default()
            },
            BlastDefline {
                title: None,
                ids: vec![ParsedSeqId::Local(ObjectId::Str(b"second".to_vec()))],
                taxid: None,
                memberships: Some(Vec::new()),
                links: Some(Vec::new()),
                other_info: None,
                ..BlastDefline::default()
            },
        ],
        utf8: false,
    };

    let mut bad = Vec::new();
    let encoded = encode_defline_set(&sample, &mut bad);
    assert!(bad.is_empty());

    // Definite-length BER header with UTF8String tag on title and replaced byte 0x80.
    let utf8_and_bad: Vec<u8> = vec![
        0x30, 0x16, // SEQUENCE OF len 22
        0x30, 0x14, // Blast-def-line len 20
        0xa0, 0x05, 0x0c, 0x03, b'a', 0x80, b'b', // title: UTF8String "a\x80b"
        0xa1, 0x0b, 0x30, 0x09, // seqid: SEQUENCE OF
        0xa0, 0x07, 0xa1, 0x05, 0x1a, 0x03, b'i', b'd', b'1', // local str "id1"
    ];

    for raw in [encoded.as_slice(), utf8_and_bad.as_slice()] {
        let ber_set = crate::seqid::ber::decode_full_set(raw).expect("valid BER");
        let (table_res, posts) = decode_binary_deflines_with_names(raw, Names::DEFLINE_HEADER);
        let decoded = table_res.expect("table decode ok");
        assert!(!decoded.null_defline && decoded.null_seqids.is_empty());
        let table_set = decoded.set;
        assert_eq!(table_set.utf8, ber_set.utf8);
        assert_eq!(table_set.deflines.len(), ber_set.deflines.len());
        for (a, b) in table_set.deflines.iter().zip(&ber_set.deflines) {
            assert_eq!(a.title, b.title);
            assert_eq!(a.ids, b.ids);
            assert_eq!(a.taxid, b.taxid);
            assert_eq!(a.memberships, b.memberships);
            assert_eq!(a.links, b.links);
            assert_eq!(a.other_info, b.other_info);
        }
        if raw == utf8_and_bad.as_slice() {
            assert!(table_set.utf8);
            assert!(posts.iter().any(|p| matches!(p, Post::Utf8Tag(_))));
            assert!(posts.iter().any(|p| matches!(p, Post::BadChar(_))));
        } else {
            assert!(posts.is_empty());
        }
    }

    // Invalid BER header (`t_badhdr`: `30 80 DE AD 00 00`): exact NCBI serial exception and stack.
    let bad_hdr = [0x30, 0x80, 0xde, 0xad, 0x00, 0x00];
    let (res, _) = decode_binary_deflines_with_names(&bad_hdr, Names::DEFLINE_HEADER);
    let ReadError::Serial(err) = res.err().expect("expected CSerialException") else {
        panic!("expected CSerialException");
    };
    let what = String::from_utf8_lossy(&err.ncbi_what()).into_owned();
    assert!(
        what.contains("UnexpectedTagClassByte")
            && what.contains(
                "byte 2: unexpected tag: private/30 (222), should be constructed/None (32)"
            )
            && what.contains("( at Blast-def-line-set.[].[])"),
        "{what}"
    );
}

#[test]
fn blast_def_line_table_members_match_reader_order() {
    let names: Vec<&str> = members(named::BLAST_DEF_LINE)
        .iter()
        .map(|member| member.name())
        .collect();
    assert_eq!(
        names,
        [
            "title",
            "seqid",
            "taxid",
            "memberships",
            "links",
            "other-info"
        ]
    );
}

#[test]
fn decode_binary_deflines_reads_valid_and_reports_ncbi_errors_and_posts() {
    let tlv = |tag: u8, body: &[u8]| -> Vec<u8> {
        let mut out = vec![tag];
        if tag & 0x20 != 0 {
            out.push(0x80);
            out.extend_from_slice(body);
            out.extend_from_slice(&[0, 0]);
        } else {
            assert!(body.len() < 128);
            out.push(body.len() as u8);
            out.extend_from_slice(body);
        }
        out
    };
    let seq = |parts: &[&[u8]]| -> Vec<u8> {
        let body: Vec<u8> = parts.iter().flat_map(|p| p.iter().copied()).collect();
        tlv(0x30, &body)
    };
    let ctx = |n: u8, body: &[u8]| tlv(0xA0 | n, body);
    let vis = |b: &[u8]| tlv(0x1A, b);
    let int = |v: i32| tlv(0x02, &v.to_be_bytes());
    let local_str = |s: &[u8]| ctx(0, &ctx(1, &vis(s)));

    // Full valid Blast-def-line with every optional field + trailing bytes.
    let mut good = seq(&[&seq(&[
        &ctx(0, &vis(b"hello")),
        &ctx(1, &seq(&[&local_str(b"id1")])),
        &ctx(2, &int(9606)),
        &ctx(3, &seq(&[&int(1)])),
        &ctx(4, &seq(&[&int(2), &int(3)])),
        &ctx(5, &seq(&[&int(4)])),
    ])]);
    good.extend_from_slice(&[0xDE, 0xAD]);
    let (res, posts) = decode_binary_deflines(&good);
    assert!(posts.is_empty());
    let decoded = res.unwrap();
    assert!(!decoded.null_defline && decoded.null_seqids.is_empty());
    assert_eq!(decoded.set.deflines.len(), 1);
    let line = &decoded.set.deflines[0];
    assert_eq!(line.title.as_deref(), Some(b"hello".as_slice()));
    assert_eq!(
        line.ids,
        [ParsedSeqId::Local(ObjectId::Str(b"id1".to_vec()))]
    );
    assert_eq!(line.taxid, Some(9606));
    assert_eq!(line.memberships.as_deref(), Some([1].as_slice()));
    assert_eq!(line.links.as_deref(), Some([2, 3].as_slice()));
    assert_eq!(line.other_info.as_deref(), Some([4].as_slice()));

    // RB04 (`0xDEADBEEF`): exact NCBI serial exception text.
    let (res, posts) = decode_binary_deflines(&[0xDE, 0xAD, 0xBE, 0xEF]);
    assert!(posts.is_empty());
    let Err(ReadError::Serial(exc)) = res else {
        panic!("expected ReadError::Serial for RB04");
    };
    assert_eq!(
        exc.get_msg(),
        b"byte 0: unexpected tag: private/30 (222), should be constructed/None (32)"
    );

    // Truncated at root (0x30 -> Eof) vs inside element (0x30 0x80 -> Serial).
    assert!(matches!(
        decode_binary_deflines(&[0x30]).0,
        Err(ReadError::Eof)
    ));
    let Err(ReadError::Serial(exc)) = decode_binary_deflines(&[0x30, 0x80]).0 else {
        panic!("expected ReadError::Serial for 30 80");
    };
    assert_eq!(exc.get_msg(), b"byte 2:  Frame type= eFrameArrayElement");

    // Bad visible chars followed by an error: posts emitted before the error.
    let bad_then_err = seq(&[
        &seq(&[&ctx(0, &vis(b"a\x01b")), &ctx(1, &seq(&[&local_str(b"x")]))]),
        &seq(&[&ctx(9, &int(1))]),
    ]);
    let (res, posts) = decode_binary_deflines(&bad_then_err);
    assert_eq!(posts.len(), 1);
    assert!(matches!(
        &posts[0],
        Post::BadChar(r) if r.trace == "Blast-def-line-set.E.title" && r.byte == 1
    ));
    assert!(matches!(res, Err(ReadError::Serial(_))));

    // NULL Blast-def-line element vs NULL Seq-id elements at various positions.
    let null_def = seq(&[&[0x05, 0x00]]);
    let decoded = decode_binary_deflines(&null_def).0.unwrap();
    assert!(decoded.null_defline);
    let null_id = seq(&[&seq(&[&ctx(1, &seq(&[&[0x05, 0x00]]))])]);
    let decoded = decode_binary_deflines(&null_id).0.unwrap();
    assert!(!decoded.null_defline);
    assert_eq!(decoded.null_seqids, [(0, 0)]);
    let null_multi = seq(&[
        &seq(&[&ctx(1, &seq(&[&local_str(b"id1"), &[0x05, 0x00]]))]),
        &seq(&[&ctx(
            1,
            &seq(&[&[0x05, 0x00], &local_str(b"id2"), &[0x05, 0x00]]),
        )]),
    ]);
    let decoded = decode_binary_deflines(&null_multi).0.unwrap();
    assert!(!decoded.null_defline);
    assert_eq!(decoded.null_seqids, [(0, 1), (1, 0), (1, 2)]);
    assert_eq!(decoded.set.deflines[0].ids.len(), 1);
    assert_eq!(decoded.set.deflines[1].ids.len(), 1);
}
