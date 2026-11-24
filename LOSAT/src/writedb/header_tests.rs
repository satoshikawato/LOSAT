use super::*;

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

// Records that makeblastdb 2.17.0 wrote (2026-10-12).
#[test]
fn the_record_of_a_short_title() {
    let mut bad = Vec::new();
    assert_eq!(
        hex(&defline_set(b"a one", 0, 0, &mut bad)),
        "30803080a0801a0561206f6e650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
    );
    assert_eq!(
        hex(&defline_set(b"c three", 1, 0, &mut bad)),
        "30803080a0801a07632074687265650000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020101000000000000000000000000a280020100000000000000"
    );
    assert_eq!(
        hex(&defline_set(b"", 0, 0, &mut bad)),
        "30803080a0801a000000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
    );
    assert!(bad.is_empty());
}

#[test]
fn bytes_outside_the_visible_range_become_number_signs() {
    let mut bad = Vec::new();
    let record = defline_set(
        b"a   tab  ctl \xc3\xa9 utf8 \x80hi \x7f del",
        0,
        0,
        &mut bad,
    );
    assert_eq!(
        hex(&record),
        "30803080a0801a1e61202020746162202063746c20232320757466382023686920232064656c0000a1803080aa803080a0801a09424c5f4f52445f49440000a180a080020100000000000000000000000000a280020100000000000000"
    );
    // "Bad char [0xC3] in string at byte 21", [0xA9] 22, [0x80] 29, [0x7F] 33
    let reported: Vec<(u8, usize)> = bad.iter().map(|b| (b.byte, b.stream_offset)).collect();
    assert_eq!(reported, [(0xC3, 21), (0xA9, 22), (0x80, 29), (0x7F, 33)]);
}

/// The runs that `write_visible_string` copies whole give the record and the reports of the loop
/// over each byte that it replaced, for texts of every byte and length.
#[test]
fn visible_runs_are_written_as_byte_by_byte() {
    let mut prefix = vec![0x30, 0x80];
    for text in crate::utils::byte_search::tests::search_texts() {
        let earlier = BadTitleChar {
            byte: 0,
            stream_offset: 0,
        };
        let (mut out, mut bad) = (prefix.clone(), vec![earlier]);
        write_visible_string(&mut out, &text, &mut bad);
        let (mut expected, mut expected_bad) = (prefix.clone(), vec![earlier]);
        expected.push(TAG_VISIBLE_STRING);
        write_length(&mut expected, text.len());
        for &byte in &text {
            if (b' '..=b'~').contains(&byte) {
                expected.push(byte);
            } else {
                expected_bad.push(BadTitleChar {
                    byte,
                    stream_offset: expected.len(),
                });
                expected.push(b'#');
            }
        }
        assert_eq!(out, expected, "{text:?}");
        assert_eq!(bad, expected_bad, "{text:?}");
        prefix.push(b'x');
    }
}

#[test]
fn lengths_and_integers_take_the_shortest_form() {
    let mut out = Vec::new();
    write_length(&mut out, 127);
    write_length(&mut out, 128);
    write_length(&mut out, 255);
    write_length(&mut out, 256);
    write_length(&mut out, 5000);
    write_length(&mut out, 0x1_0000);
    assert_eq!(hex(&out), "7f818081ff82010082138883010000");
    let mut out = Vec::new();
    for value in [0, 127, 128, 32767, 32768, 8_388_607, 8_388_608, i32::MAX] {
        write_int4(&mut out, value);
    }
    assert_eq!(
        hex(&out),
        concat!(
            "020100",
            "02017f",
            "02020080",
            "02027fff",
            "0203008000",
            "02037fffff",
            "020400800000",
            "02047fffffff"
        )
    );
}

#[test]
fn the_reader_decodes_what_the_writer_encodes() {
    let path = std::path::Path::new("test.phr");
    for (title, oid) in [
        (&b"sp|P12345|ABC_HUMAN desc"[..], 0u32),
        (&[b'x'; 5000][..], 128),
        (&b""[..], 70_000),
        (&b"t"[..], i32::MAX as u32),
    ] {
        let record = defline_set(title, oid, 0, &mut Vec::new());
        let deflines = crate::seqdb::header::decode_defline_set(&record, path, oid).unwrap();
        assert_eq!(deflines.len(), 1);
        assert_eq!(deflines[0].title, title);
        assert_eq!(deflines[0].volume_ordinal, oid);
    }
}

mod parsed {
    use super::super::*;
    use crate::seqid::model::{DateStd, Dbtag, GiimportId, IdPat, PatentSeqId, PdbSeqId};
    use crate::seqid::BlastDefline;

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|byte| format!("{byte:02x}")).collect()
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

    fn one(title: &str, ids: Vec<ParsedSeqId>) -> DeflineSet {
        DeflineSet {
            deflines: vec![BlastDefline {
                title: text(title),
                ids,
                taxid: Some(0),
                ..BlastDefline::default()
            }],
            utf8: false,
        }
    }

    fn pdb(mol: &str, chain: Option<i32>, chain_id: Option<&str>) -> ParsedSeqId {
        ParsedSeqId::Pdb(PdbSeqId {
            mol: mol.as_bytes().to_vec(),
            chain,
            rel: None,
            chain_id: chain_id.and_then(text),
        })
    }

    fn encode(set: &DeflineSet) -> String {
        let mut bad = Vec::new();
        let record = encode_defline_set(set, &mut bad);
        assert!(bad.is_empty());
        assert_eq!(&crate::seqid::ber::decode_full_set(&record).unwrap(), set);
        hex(&record)
    }

    // Records of the .phr files makeblastdb 2.17.0 wrote with -parse_seqids (survey B
    // fixtures p01_types, p24_gi_range and p02_multi, v4; the v5 files are equal).
    #[test]
    fn records_of_makeblastdb() {
        assert_eq!(
            encode(&one(
                "Hemoglobin subunit alpha OS=Homo sapiens",
                vec![ParsedSeqId::Swissprot(textseq(Some("P69905"), Some("HBA_HUMAN"), Some("reviewed"), None))]
            )),
            "30803080a0801a2848656d6f676c6f62696e20737562756e697420616c706861204f533d486f6d6f2073617069656e730000a1803080a7803080a0801a094842415f48554d414e0000a1801a065036393930350000a2801a08726576696577656400000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one(
                "gi plus ref",
                vec![ParsedSeqId::Gi(123456), ParsedSeqId::Other(textseq(Some("NP_000550"), None, None, Some(2)))]
            )),
            "30803080a0801a0b676920706c7573207265660000a1803080ab80020301e2400000a9803080a1801a094e505f3030303535300000a38002010200000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("pdb chain A", vec![pdb("1ABC", Some(65), Some("A"))])),
            "30803080a0801a0b70646220636861696e20410000a1803080ae803080a0801a04314142430000a1800201410000a3801a014100000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("pdb two-letter chain", vec![pdb("2XYZ", None, Some("AA"))])),
            "30803080a0801a147064622074776f2d6c657474657220636861696e0000a1803080ae803080a0801a043258595a0000a3801a02414100000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("pdb no chain", vec![pdb("3DEF", None, None)])),
            "30803080a0801a0c706462206e6f20636861696e0000a1803080ae803080a0801a043344454600000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("local numeric", vec![ParsedSeqId::Local(ObjectId::Id(123))])),
            "30803080a0801a0d6c6f63616c206e756d657269630000a1803080a080a08002017b0000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one(
                "general numeric tag",
                vec![ParsedSeqId::General(Dbtag { db: b"mydb".to_vec(), tag: ObjectId::Id(456) })]
            )),
            "30803080a0801a1367656e6572616c206e756d65726963207461670000a1803080aa803080a0801a046d7964620000a180a080020201c8000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one(
                "patent",
                vec![ParsedSeqId::Patent(PatentSeqId {
                    seqid: 1,
                    cit: IdPat {
                        country: b"US".to_vec(),
                        id: PatentNumber::Number(b"1234567".to_vec()),
                        doc_type: None,
                    },
                })]
            )),
            "30803080a0801a06706174656e740000a1803080a8803080a0800201010000a1803080a0801a0255530000a180a0801a073132333435363700000000000000000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("gi 2^31", vec![
                ParsedSeqId::Gi(2_147_483_648),
                ParsedSeqId::Other(textseq(Some("NP_000011"), None, None, Some(1))),
            ])),
            "30803080a0801a07676920325e33310000a1803080ab80020500800000000000a9803080a1801a094e505f3030303031310000a38002010100000000000000000000a280020100000000000000"
        );
        assert_eq!(
            encode(&one("gi 2^32-1", vec![ParsedSeqId::Gi(4_294_967_295)])),
            "30803080a0801a09676920325e33322d310000a1803080ab80020500ffffffff000000000000a280020100000000000000"
        );
        let refseq = |acc: &str| ParsedSeqId::Other(textseq(Some(acc), None, None, Some(1)));
        let three = DeflineSet {
            deflines: [
                ("first title", refseq("NP_100003")),
                ("second title", refseq("NP_100004")),
                ("third title", ParsedSeqId::Gi(1005)),
            ]
            .into_iter()
            .map(|(title, id)| BlastDefline {
                title: text(title),
                ids: vec![id],
                taxid: Some(0),
                ..BlastDefline::default()
            })
            .collect(),
            utf8: false,
        };
        assert_eq!(
            encode(&three),
            "30803080a0801a0b6669727374207469746c650000a1803080a9803080a1801a094e505f3130303030330000a38002010100000000000000000000a280020100000000003080a0801a0c7365636f6e64207469746c650000a1803080a9803080a1801a094e505f3130303030340000a38002010100000000000000000000a280020100000000003080a0801a0b7468697264207469746c650000a1803080ab80020203ed000000000000a280020100000000000000"
        );
    }

    #[test]
    fn every_alternative_reads_back() {
        let ts = textseq(Some("A1"), Some("n"), Some("r"), Some(3));
        let ids = vec![
            ParsedSeqId::Local(ObjectId::Str(b"x".to_vec())),
            ParsedSeqId::Gibbsq(7),
            ParsedSeqId::Gibbmt(-8),
            ParsedSeqId::Giim(GiimportId {
                id: 9,
                db: text("d"),
                release: text("r"),
            }),
            ParsedSeqId::Genbank(ts.clone()),
            ParsedSeqId::Embl(ts.clone()),
            ParsedSeqId::Pir(ts.clone()),
            ParsedSeqId::Swissprot(ts.clone()),
            ParsedSeqId::Patent(PatentSeqId {
                seqid: 2,
                cit: IdPat {
                    country: b"EP".to_vec(),
                    id: PatentNumber::AppNumber(b"9".to_vec()),
                    doc_type: text("A"),
                },
            }),
            ParsedSeqId::Other(ts.clone()),
            ParsedSeqId::General(Dbtag {
                db: b"db".to_vec(),
                tag: ObjectId::Str(b"t".to_vec()),
            }),
            ParsedSeqId::Gi(i64::MAX),
            ParsedSeqId::Gi(-1),
            ParsedSeqId::Ddbj(ts.clone()),
            ParsedSeqId::Prf(ts.clone()),
            ParsedSeqId::Pdb(PdbSeqId {
                mol: b"1ABC".to_vec(),
                chain: Some(32),
                rel: Some(Date::Std(DateStd {
                    year: 2020,
                    month: Some(1),
                    day: None,
                    season: text("s"),
                    hour: None,
                    minute: Some(5),
                    second: None,
                })),
                chain_id: None,
            }),
            ParsedSeqId::Pdb(PdbSeqId {
                mol: b"2ABC".to_vec(),
                chain: None,
                rel: Some(Date::Str(b"x".to_vec())),
                chain_id: None,
            }),
            ParsedSeqId::Tpg(ts.clone()),
            ParsedSeqId::Tpe(ts.clone()),
            ParsedSeqId::Tpd(ts.clone()),
            ParsedSeqId::Gpipe(ts.clone()),
            ParsedSeqId::NamedAnnotTrack(ts.clone()),
        ];
        let set = DeflineSet {
            deflines: vec![
                BlastDefline {
                    title: text(&"t".repeat(300)),
                    ids,
                    taxid: Some(9606),
                    memberships: Some(vec![1, 2]),
                    links: Some(vec![]),
                    other_info: Some(vec![-5]),
                    ..BlastDefline::default()
                },
                BlastDefline {
                    title: None,
                    ids: vec![ParsedSeqId::Gi(1)],
                    ..BlastDefline::default()
                },
            ],
            utf8: false,
        };
        encode(&set);
    }

    // makeblastdb 2.17.0 on survey B fixture p23_nonascii (record 1, `>lcl|\xe9abc
    // latin-1 first byte`; record 3, `>ref|NP_000001.1| title with \xe9 byte`).
    #[test]
    fn bad_characters_are_reported_with_their_member() {
        let mut bad = Vec::new();
        let set = one(
            "latin-1 first byte",
            vec![ParsedSeqId::Local(ObjectId::Str(b"\xe9abc".to_vec()))],
        );
        let record = encode_defline_set(&set, &mut bad);
        assert_eq!(
            bad,
            [BadChar {
                path: "Blast-def-line-set.E.seqid.E.local.str".to_string(),
                byte: 0xe9,
                stream_offset: 38,
                string: b"\xe9abc".to_vec(),
            }]
        );
        assert_eq!(record[38], b'#');
        bad.clear();
        let mut set = one(
            "",
            vec![ParsedSeqId::Other(textseq(
                Some("NP_000001"),
                None,
                None,
                Some(1),
            ))],
        );
        set.deflines[0].title = Some(b"title with \xe9 byte".to_vec());
        encode_defline_set(&set, &mut bad);
        assert_eq!(bad.len(), 1);
        assert_eq!(
            (bad[0].path.as_str(), bad[0].stream_offset),
            ("Blast-def-line-set.E.title", 19)
        );
    }
}

/// The title's diagnostics without the record are those of the record, for titles of every byte
/// and of lengths whose BER length takes one to four bytes.
#[test]
fn bad_title_chars_are_those_of_the_record() {
    let mut texts = crate::utils::byte_search::tests::search_texts();
    for len in [127, 128, 255, 256, 65_535, 65_536, 1 << 24] {
        let mut text = vec![b'a'; len];
        text[0] = 0x80;
        text[len / 2] = b'\t';
        text[len - 1] = 0xff;
        texts.push(text);
    }
    for text in texts {
        let (mut alone, mut written) = (vec![], vec![]);
        bad_title_chars(&text, &mut alone);
        defline_set(&text, 7, 9606, &mut written);
        assert_eq!(alone, written, "{} bytes", text.len());
    }
}

/// The record fits the buffer it is written into: at most `DEFLINE_SET_OVERHEAD` bytes besides
/// the title, with the longest lengths and integers.
#[test]
fn the_record_fits_its_buffer() {
    for len in [0, 1, 127, 128, 255, 256, 65_535, 65_536, 1 << 24] {
        let title = vec![b'a'; len];
        for (oid, taxid) in [
            (0, 0),
            (u32::MAX, i32::MIN),
            (0x7fff_ffff, i32::MAX),
            (300, -129),
        ] {
            let record = defline_set(&title, oid, taxid, &mut Vec::new());
            assert!(
                record.len() - len <= 77,
                "{len} {oid} {taxid}: {}",
                record.len() - len
            );
            assert!(record.len() <= DEFLINE_SET_OVERHEAD + len);
        }
    }
}
