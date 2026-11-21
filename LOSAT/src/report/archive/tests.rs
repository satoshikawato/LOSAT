use super::*;
use crate::blastinput::app::Program;

/// `writes_text` and the encoders decide the same bytes: for every byte and each stream, a
/// one-byte VisibleString that `writes_text` accepts is written by the default encoder
/// without an error and by the replacing one with nothing replaced, and one that it refuses
/// is the default encoder's `NonVisibleChar` and the replacing one's single replaced byte.
#[test]
fn writes_text_is_what_each_encoder_writes_as_it_is() {
    use crate::report::serial::asn_ber::AsnBinaryEncoder;
    use crate::report::serial::asn_text::AsnTextEncoder;
    use crate::report::serial::encoder::{Encoder, SerialError};
    use crate::report::serial::xml::XmlEncoder;
    for byte in 0..=u8::MAX {
        for encoding in [
            ArchiveEncoding::Text,
            ArchiveEncoding::Binary,
            ArchiveEncoding::Xml,
        ] {
            let as_is = encoding.writes_text(&[byte]);
            let (error, replaced) = match encoding {
                ArchiveEncoding::Text => {
                    let mut buf = Vec::new();
                    let mut plain = AsnTextEncoder::new(&mut buf);
                    plain.visible_str(&[byte]);
                    let error = plain.error();
                    let mut buf = Vec::new();
                    let mut replacing = AsnTextEncoder::new(&mut buf);
                    replacing.visible_str_replacing(&[byte]);
                    (error, replacing.take_replaced())
                }
                ArchiveEncoding::Binary => {
                    let mut buf = Vec::new();
                    let mut plain = AsnBinaryEncoder::new(&mut buf);
                    plain.visible_str(&[byte]);
                    let error = plain.error();
                    let mut buf = Vec::new();
                    let mut replacing = AsnBinaryEncoder::new(&mut buf);
                    replacing.visible_str_replacing(&[byte]);
                    (error, replacing.take_replaced())
                }
                _ => {
                    let mut buf = Vec::new();
                    let mut plain = XmlEncoder::archive(&mut buf);
                    plain.visible_str(&[byte]);
                    let error = plain.error();
                    let mut buf = Vec::new();
                    let mut replacing = XmlEncoder::archive(&mut buf);
                    replacing.visible_str_replacing(&[byte]);
                    (error, replacing.take_replaced())
                }
            };
            if as_is {
                assert_eq!(error, None, "{encoding:?} 0x{byte:02X}");
                assert!(replaced.is_empty(), "{encoding:?} 0x{byte:02X}");
            } else {
                assert_eq!(
                    error,
                    Some(SerialError::NonVisibleChar(byte)),
                    "{encoding:?} 0x{byte:02X}"
                );
                assert_eq!(replaced.len(), 1, "{encoding:?} 0x{byte:02X}");
                assert_eq!(replaced[0].byte, byte);
            }
        }
    }
}

#[test]
fn archive_format_is_read_like_print_archive() {
    let of = |value: &str| ArchiveEncoding::of(Some(OsStr::new(value)));
    assert_eq!(ArchiveEncoding::of(None), ArchiveEncoding::Text);
    assert_eq!(of(""), ArchiveEncoding::Text);
    assert_eq!(of("xml"), ArchiveEncoding::Xml);
    assert_eq!(of("XML"), ArchiveEncoding::Xml);
    assert_eq!(of("Xml"), ArchiveEncoding::Xml);
    assert_eq!(of("bin"), ArchiveEncoding::Binary);
    assert_eq!(of("BIN"), ArchiveEncoding::Binary);
    assert_eq!(of("binary"), ArchiveEncoding::Binary);
    assert_eq!(of("bin foo"), ArchiveEncoding::Binary);
    assert_eq!(of(" bin"), ArchiveEncoding::Nothing);
    assert_eq!(of("asn"), ArchiveEncoding::Nothing);
    assert_eq!(of("text"), ArchiveEncoding::Nothing);
    assert_eq!(of("foo"), ArchiveEncoding::Nothing);
    assert_eq!(of("xmlx"), ArchiveEncoding::Nothing);
    assert_eq!(of("bi"), ArchiveEncoding::Nothing);
}

/// The strings the archive copies from the command line are written as they are where the
/// encoding does: text and binary write `#` for a byte outside `0x20..=0x7E`, XML writes
/// every byte but a control one as it is (`writes_text`); such a request is written through
/// `Traced` (`writer::write_archive`).
#[test]
fn a_request_tells_whether_its_strings_are_written_as_they_are() {
    let as_is = |encoding, path: Option<&str>, db: Option<&str>| {
        let mut inputs = OptionInputs::new(Program::Blastn, "megablast");
        inputs.window_masker_db = path.map(str::to_string);
        ArchiveRequest::for_inputs(&inputs, db).written_as_is(encoding)
    };
    for encoding in [
        ArchiveEncoding::Text,
        ArchiveEncoding::Binary,
        ArchiveEncoding::Xml,
    ] {
        assert!(as_is(encoding, None, None));
        assert!(as_is(encoding, Some("wm&<\"'q/mg.obinary"), Some("d b&<")));
        // A control byte: NCBI writes `#` and a post, but TAB, LF and CR in XML.
        for bad in ["a\x01b", "a\x1fb", "a\x0bb", "a\x0eb"] {
            assert!(!as_is(encoding, Some(bad), None), "{encoding:?} {bad:?}");
            assert!(!as_is(encoding, None, Some(bad)), "{encoding:?} {bad:?}");
        }
    }
    for good in ["a\tb", "a\nb", "a\rb"] {
        assert!(as_is(ArchiveEncoding::Xml, Some(good), Some(good)));
        assert!(!as_is(ArchiveEncoding::Text, Some(good), None));
    }
    // Bytes from 0x7F up: replaced in text and binary, written as they are in XML.
    for high in ["wm\u{e9}/mg.obinary", "d\u{7f}b", "\u{1f600}"] {
        assert!(
            as_is(ArchiveEncoding::Xml, Some(high), Some(high)),
            "{high:?}"
        );
        for encoding in [ArchiveEncoding::Text, ArchiveEncoding::Binary] {
            assert!(!as_is(encoding, Some(high), None), "{encoding:?} {high:?}");
            assert!(!as_is(encoding, None, Some(high)), "{encoding:?} {high:?}");
        }
    }
}

/// A record's title and Seq-ids: in text and binary a byte outside `0x20..=0x7E`, in XML a
/// control byte (which no FASTA title or id holds: NCBI's reader ends the title there).
#[test]
fn a_bioseq_tells_whether_its_strings_are_written_as_they_are() {
    use requests::RequestRecord;
    let caf = FastaRecord::new("Query_1", "caf\u{e9}".as_bytes(), b"ACGT");
    let del = FastaRecord::new("Query_1", b"ti\x7ftle", b"ACGT");
    let plain = FastaRecord::new("Query_1", b"title", b"ACGT");
    for encoding in [ArchiveEncoding::Text, ArchiveEncoding::Binary] {
        assert!(!caf.bioseq(false).written_as_is(encoding));
        assert!(!del.bioseq(false).written_as_is(encoding));
        assert!(plain.bioseq(false).written_as_is(encoding));
    }
    assert!(caf.bioseq(false).written_as_is(ArchiveEncoding::Xml));
    assert!(del.bioseq(false).written_as_is(ArchiveEncoding::Xml));
    // A parsed local id with DEL (NCBI's `CSeq_id::CheckLocalID` takes it).
    let ids = [crate::seqid::ParsedSeqId::Local(
        crate::seqid::model::ObjectId::Str(b"a\x7fb".to_vec()),
    )];
    let bioseq = bioseq::BioseqInput {
        ids: Some(&ids),
        ..plain.bioseq(false)
    };
    assert!(!bioseq.written_as_is(ArchiveEncoding::Text));
    assert!(bioseq.written_as_is(ArchiveEncoding::Xml));
}

#[test]
fn the_search_space_option_carries_over_unsearched_batches() {
    let handle = ArchiveHandle::new(
        ArchiveEncoding::Text,
        ArchiveRequest {
            program: "blastp",
            application: "blastp",
            service: "plain",
            algorithm: Vec::new(),
            database: None,
            step_number: None,
        },
        OptionsHandle::default(),
    );
    let searched = |space| {
        Some(writer::SearchedContext0 {
            space,
            query_length: 142,
        })
    };
    let space = |handle: &ArchiveHandle, searched| search_space(&handle.request_sets(searched));
    assert_eq!(space(&handle, None), None);
    assert_eq!(space(&handle, searched(8738)), Some(8738));
    assert_eq!(space(&handle, None), Some(8738));
    assert_eq!(space(&handle, searched(0)), Some(0));
    assert_eq!(space(&handle, None), Some(0));
    // A statistics length of 0, and a negative one (probe run2 of line RG D1a-3: the query's
    // length, 142): the run's options handle is set once, so each case has its own.
    let with_length = |length| {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastp",
                application: "blastp",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
                step_number: None,
            },
            OptionsHandle::default(),
        );
        handle.options().set_database_filter(None, length).unwrap();
        handle
    };
    assert_eq!(
        space(&with_length(RequestLength::Lengths(None)), searched(8738)),
        Some(0)
    );
    assert_eq!(
        space(
            &with_length(RequestLength::Lengths(Some((i64::MIN + 1, 23)))),
            searched(15_876)
        ),
        Some(142)
    );
    assert_eq!(
        space(
            &with_length(RequestLength::Lengths(Some((-1, -5)))),
            Some(writer::SearchedContext0 {
                space: 0,
                query_length: 0,
            })
        ),
        Some(0)
    );
}

/// The `EffectiveSearchSpace` option among a request's sets.
pub(super) fn search_space(sets: &[options::Param]) -> Option<i64> {
    sets.iter()
        .find(|param| param.name == "EffectiveSearchSpace")
        .map(|param| match param.value {
            options::ParamValue::BigInteger(space) => space,
            ref other => panic!("EffectiveSearchSpace {other:?}"),
        })
}

/// The application's parameters come first and each batch's set-up's before its search
/// space, by `x_AttachValue` (rpsblast's documents after the first, NCBI BLAST+ 2.17.0
/// `p9.rps.3v.q40.n1.f6`, `p9.rps.err.nofreq`).
#[test]
fn the_handle_sets_fold_in_ncbis_order() {
    use options::{Param, ParamValue};
    let handle = OptionsHandle::default();
    handle
        .set_search_sets(
            vec![Param {
                name: "CompositionBasedStats",
                value: ParamValue::Integer(0),
            }],
            vec![
                Param {
                    name: "DbLength",
                    value: ParamValue::BigInteger(560),
                },
                Param {
                    name: "MatrixName",
                    value: ParamValue::Str("BLOSUM62".to_string()),
                },
            ],
        )
        .unwrap();
    assert!(handle.set_search_sets(Vec::new(), Vec::new()).is_err());
    let mut carry = HandleCarry::default();
    let names = |carry: &HandleCarry| -> Vec<&'static str> {
        handle
            .request_sets(carry)
            .iter()
            .map(|param| param.name)
            .collect()
    };
    assert_eq!(names(&carry), ["CompositionBasedStats"]);
    handle.pass(&mut carry, None);
    assert_eq!(
        names(&carry),
        ["CompositionBasedStats", "DbLength", "MatrixName"]
    );
    let searched = |space| {
        Some(writer::SearchedContext0 {
            space,
            query_length: 142,
        })
    };
    handle.pass(&mut carry, searched(85_680));
    handle.pass(&mut carry, searched(96_628));
    assert_eq!(
        names(&carry),
        [
            "CompositionBasedStats",
            "DbLength",
            "MatrixName",
            "EffectiveSearchSpace"
        ]
    );
    assert_eq!(search_space(&handle.request_sets(&carry)), Some(96_628));
}

#[test]
fn the_tee_saves_lines_for_the_next_archive() {
    let handle = ArchiveHandle::new(
        ArchiveEncoding::Text,
        ArchiveRequest {
            program: "blastn",
            application: "blastn",
            service: "plain",
            algorithm: Vec::new(),
            database: None,
            step_number: None,
        },
        OptionsHandle::default(),
    );
    let mut seen: Vec<u8> = Vec::new();
    {
        let mut tee = DiagTee::new(&mut seen, handle.clone());
        tee.write_all(b"Warning: [blastn] a  b\n").unwrap();
    }
    assert_eq!(seen, b"Warning: [blastn] a  b\n");
    assert_eq!(
        handle.take_messages(),
        vec![(1, b"Warning: [blastn] a b".to_vec())]
    );
    assert!(handle.take_messages().is_empty());
}

/// NCBI 2.17.0 (`PrintErrorArchive`): each text below was run with a missing `-query` and
/// wrote the error archive (`true`) or nothing (`false`).
#[test]
fn the_error_archive_follows_the_outfmt_text_like_print_error_archive() {
    for (text, writes) in [
        ("11", true),
        (" 11", true),
        ("11 ", true),
        ("\t11\n", true),
        ("\x0b11\x0c", true),
        ("011", true),
        ("+11", true),
        ("00000000000011", true),
        ("11 qseqid", false),
        (" 11 qseqid", false),
        ("11qseqid", false),
        ("11.0", false),
        ("1 1", false),
        ("-11", false),
        ("0x0B", false),
        ("1_1", false),
        ("11,0", false),
        ("4294967307", false),
        ("2147483648", false),
        ("", false),
        (" ", false),
        ("+", false),
        ("6", false),
    ] {
        assert_eq!(writes_error_archive(text), writes, "{text:?}");
    }
}

#[test]
fn a_handle_tied_to_a_spec_writes_no_error_archive() {
    let make = || {
        ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "tblastx",
                application: "tblastx",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
                step_number: None,
            },
            OptionsHandle::default(),
        )
    };
    let lines = vec![(ERROR, b"Error: x".to_vec())];
    // No text (the API): the error archive is written.
    let mut sink = Vec::new();
    write_error_archive(&make(), &lines, &mut sink).unwrap();
    assert!(!sink.is_empty());
    // `11`: written; `11 qseqid`: nothing, and the messages stay unread.
    let mut sink = Vec::new();
    write_error_archive(&make().with_outfmt_text(" 11 "), &lines, &mut sink).unwrap();
    assert!(!sink.is_empty());
    let mut sink = Vec::new();
    let spec = make().with_outfmt_text("11 qseqid");
    assert!(!spec.error_archive_enabled());
    write_error_archive(&spec, &lines, &mut sink).unwrap();
    assert!(sink.is_empty());
}

/// PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits (review C of
/// line RG D1a-4, SF1): an archive
/// whose request holds a big integer of -2^63 (`-dbsize -9223372036854775808`), which NCBI's
/// ASN.1 text and XML archives write as `-0` (probes q6, q10b of review C), is written with its
/// digits, not refused.
#[test]
fn an_archive_of_a_dbsize_of_minus_two_to_the_63_writes_its_digits() {
    struct NoAligns;
    impl AlignSource for NoAligns {
        fn write<E: crate::report::serial::encoder::Encoder>(
            &mut self,
            _out: &mut ArchiveOut<'_, E>,
        ) -> anyhow::Result<()> {
            Ok(())
        }
    }
    let mut inputs = options::OptionInputs::new(Program::Blastp, "blastp");
    inputs.dbsize = Some(i64::MIN);
    for (encoding, written) in [
        (
            ArchiveEncoding::Text,
            "name \"DbLength\",\n          value big-integer -9223372036854775808\n",
        ),
        (
            ArchiveEncoding::Xml,
            "<Blast4-value_big-integer>-9223372036854775808</Blast4-value_big-integer>",
        ),
    ] {
        let handle = ArchiveHandle::new(
            encoding,
            ArchiveRequest::for_inputs(&inputs, Some("in/db/bp")),
            OptionsHandle::default(),
        );
        let batch = BatchArchive {
            request: writer::RequestInputs {
                queries: writer::RequestQueries::Bioseqs(&[]),
                subjects: &[],
                required: None,
                lcase_masks: &[],
            },
            masks: &[],
            ka_blocks: [KaBlock::missing(false), KaBlock::missing(true)],
            search_space: 0,
            length_adjustment: 0,
            searched_space: None,
        };
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        write_archive(&handle, &batch, &mut NoAligns, &mut sink, &mut diagnostics).unwrap();
        let text = String::from_utf8(sink).unwrap();
        assert!(text.contains(written), "{encoding:?}: {text}");
        assert!(!text.contains("-0\n") && !text.contains(">-0<"), "{text}");
    }
}
