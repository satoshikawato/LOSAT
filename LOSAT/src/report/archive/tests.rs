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
        inputs.window_masker_db = path.map(|s| s.as_bytes().to_vec());
        ArchiveRequest::for_inputs(&inputs, db.map(str::as_bytes)).written_as_is(encoding)
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
        handle
            .options()
            .set_database_filter(Vec::new(), length)
            .unwrap();
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
        .find(|param| param.name() == "EffectiveSearchSpace")
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
            vec![Param::new(
                fields::COMPOSITION_BASED_STATS,
                ParamValue::Integer(0),
            )],
            vec![
                Param::new(fields::DB_LENGTH, ParamValue::BigInteger(560)),
                Param::new(fields::MATRIX_NAME, ParamValue::Str(b"BLOSUM62".to_vec())),
            ],
            Vec::new(),
        )
        .unwrap();
    assert!(handle
        .set_search_sets(Vec::new(), Vec::new(), Vec::new())
        .is_err());
    let mut carry = HandleCarry::default();
    let names = |carry: &HandleCarry| -> Vec<&'static str> {
        handle
            .request_sets(carry)
            .iter()
            .map(|param| param.name())
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

// psiblast `-phi_pattern` (PS8a fix round 1): psiblast resets the composition-based statistics
// after each batch, which PHI-BLAST's options never set: the reset comes after the batch's search
// space and stays there (NCBI BLAST+ 2.17.0 `st.multi`: requests 2 and 3 hold
// `EffectiveSearchSpace`, then `CompositionBasedStats 0`; `st.multi.p250`, a batch NCBI does not
// search, only the latter).
#[test]
fn the_after_batch_sets_follow_the_search_space() {
    use options::{Param, ParamValue};
    let handle = OptionsHandle::default();
    handle
        .set_search_sets(
            Vec::new(),
            Vec::new(),
            vec![Param::new(
                fields::COMPOSITION_BASED_STATS,
                ParamValue::Integer(0),
            )],
        )
        .unwrap();
    let names = |carry: &HandleCarry| -> Vec<&'static str> {
        handle
            .request_sets(carry)
            .iter()
            .map(|param| param.name())
            .collect()
    };
    let mut carry = HandleCarry::default();
    assert!(names(&carry).is_empty());
    let searched = |space| {
        Some(writer::SearchedContext0 {
            space,
            query_length: 104,
        })
    };
    handle.pass(&mut carry, searched(338_444));
    assert_eq!(
        names(&carry),
        ["EffectiveSearchSpace", "CompositionBasedStats"]
    );
    handle.pass(&mut carry, searched(338_000));
    assert_eq!(
        names(&carry),
        ["EffectiveSearchSpace", "CompositionBasedStats"]
    );
    assert_eq!(search_space(&handle.request_sets(&carry)), Some(338_000));
    let mut unsearched = HandleCarry::default();
    handle.pass(&mut unsearched, None);
    assert_eq!(names(&unsearched), ["CompositionBasedStats"]);
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
            _enc: &mut E,
            _flush: &mut dyn FnMut(&mut E),
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
            ArchiveRequest::for_inputs(&inputs, Some(b"in/db/bp")),
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
            round: None,
        };
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        write_archive(&handle, &batch, &mut NoAligns, &mut sink, &mut diagnostics).unwrap();
        let text = String::from_utf8(sink).unwrap();
        assert!(text.contains(written), "{encoding:?}: {text}");
        assert!(!text.contains("-0\n") && !text.contains(">-0<"), "{text}");
    }
}

/// The error archive after the report stream's failure in outfmt 11 (`error_lines` of
/// `ostream::ReportStream::finish`'s error: the 802.4 post of the archive's `Close` and
/// "BLAST failed to write output"). NCBI 2.17.0 writes it into the failed stream and posts its
/// write fault at byte 356 (blastn), its length (`T/t6/survey` blastn.f11.q1 full2).
#[test]
fn the_error_archive_of_a_failed_report_has_ncbis_length() {
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
    let error: anyhow::Error = crate::cli::NativeError {
        exit: 6,
        message: "Error: (802.4) [blastn] [Cannot close serializing output stream] \
                  Exception: basic_ios::clear: iostream error\nBLAST failed to write output\n"
            .into(),
    }
    .into();
    let lines = error_lines(&error, false).expect("the output error's lines");
    // The fault is at the end of the buffer's first write: the document without the
    // separator, which is a write of its own (`WriteSeparator`).
    assert_eq!(
        writer::error_archive_fault(ArchiveEncoding::Text, &lines),
        356
    );
    let mut sink = Vec::new();
    write_error_archive(&handle, &lines, &mut sink).unwrap();
    assert_eq!(sink.len(), 357);
    assert!(sink.ends_with(b"}\n"));
}

/// psiblast's archive of a round (`ArchiveRound`, psiblast_app.cpp:259-271): `Web_StepNumber` is
/// the round's number with `-db` (none with `-subject`), and a PSSM round records the
/// composition-based statistics of its set-up in place of the handle's.
#[test]
fn a_psiblast_round_writes_its_step_and_composition_mode() {
    struct NoAligns;
    impl AlignSource for NoAligns {
        fn write<E: crate::report::serial::encoder::Encoder>(
            &mut self,
            _enc: &mut E,
            _flush: &mut dyn FnMut(&mut E),
        ) -> anyhow::Result<()> {
            Ok(())
        }
    }
    let mut inputs = options::OptionInputs::new(options::RequestKind::Psiblast, "psiblast");
    inputs.num_iterations = Some(3);
    let archive = |database: Option<&str>, round: Option<writer::ArchiveRound>| {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest::for_inputs(&inputs, database.map(str::as_bytes)),
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
            round,
        };
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        write_archive(&handle, &batch, &mut NoAligns, &mut sink, &mut diagnostics).unwrap();
        String::from_utf8(sink).unwrap()
    };
    let step =
        |number: i32| format!("name \"Web_StepNumber\",\n          value integer {number}\n");
    let cbs =
        |mode: i32| format!("name \"CompositionBasedStats\",\n          value integer {mode}\n");
    let fasta_round = writer::ArchiveRound {
        number: 1,
        comp_based_stats: None,
    };
    let pssm_round = writer::ArchiveRound {
        number: 2,
        comp_based_stats: Some(1),
    };
    let whole = archive(Some("in/db"), None);
    assert!(whole.contains(&step(3)), "{whole}");
    let first = archive(Some("in/db"), Some(fasta_round));
    assert!(
        first.contains(&step(1)) && !first.contains(&step(3)),
        "{first}"
    );
    let second = archive(Some("in/db"), Some(pssm_round));
    assert!(second.contains(&step(2)), "{second}");
    assert!(
        second.contains(&cbs(1)) && !second.contains(&cbs(2)),
        "{second}"
    );
    let subject = archive(None, Some(pssm_round));
    assert!(!subject.contains("Web_StepNumber"), "{subject}");
}

// Review fixes of line PS item PS15b-1, XA (2026-11-23), and line RG-F3's PR 2: what LOSAT's
// archives do not write yet of a PSSM file's members, which NCBI BLAST+ 2.17.0's psiblast writes
// with `-outfmt 11` (`pssm_unwritten`: the run is refused before any output). Since PR 2 every
// stream writes the anonymous types, NULLs, UTF8Strings and REAL specials of these files and a
// delayed member that it copies or reads again without a message or with a replaced byte's post;
// since item O-2 one whose read fails too.
/// The PSSM file `path` under tests/fixtures (gzipped when it ends with `.gz`) as psiblast's
/// `-in_pssm` reads it: ASN.1 binary for `.asnb`, XML for `.xml`, else ASN.1 text.
fn fixture_pssm(path: &str) -> crate::report::serial::pssm::PssmWithParameters {
    let mut bytes = std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(path),
    )
    .unwrap();
    if path.ends_with(".gz") {
        let mut inflated = Vec::new();
        std::io::Read::read_to_end(&mut flate2::read::GzDecoder::new(&bytes[..]), &mut inflated)
            .unwrap();
        bytes = inflated;
    }
    if path.ends_with(".asnb") {
        crate::blastinput::pssm_text::read_pssm_binary(&bytes)
            .result
            .unwrap()
    } else if path.contains(".xml") {
        crate::blastinput::pssm_text::read_pssm_xml(&bytes)
            .result
            .unwrap()
    } else {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&bytes).unwrap()
    }
}

#[test]
fn the_archive_streams_gaps_of_a_pssm_files_members() {
    // (file, the gap of the XML, ASN.1 binary and ASN.1 text archives)
    for (path, xml, binary, text) in [
        // A title only.
        ("db_search/db.ps.asn.rnd.6.asn", None, None, None),
        // A User-object (`User-field.data`, an anonymous CHOICE; a UTF8String).
        ("pssm/d3a_c10_user.asn.gz", None, None, None),
        // A Pubdesc whose `Auth-list.names` the reader keeps in a delay buffer: read again
        // without a message in XML and ASN.1 binary, copied in ASN.1 text.
        ("pssm/d3a_c11_delay.asn.gz", None, None, None),
        // Seq-annot `locs` (`Seq-annot.data`, an anonymous CHOICE) with a NULL Seq-loc.
        ("pssm/d3a_c17_locs.asn.gz", None, None, None),
        // A frequency ratio of PLUS-INFINITY in a PSSM file in ASN.1 binary.
        ("pssm/ps9_ber_real_list_pinf.asnb", None, None, None),
        // A delayed `Auth-list.names` with a control byte in a name: copied by the stream of its
        // own format, read again with NCBI's post (803.7) by the others, which write it through
        // `Traced`.
        ("pssm/ps9_xml_bin_dl_ctrl.asnb", None, None, None),
        ("pssm/ps9_xml_dl_ctrl.xml", None, None, None),
        // One whose number is past its storage: read again with an error, NCBI's exception
        // after the part of the archive written before the member, through `Traced` (line RG-F3
        // item O-2: the XML read's `CUtilException`, the binary read's `CSerialException`).
        ("pssm/ps9_xml_dl_ovf.xml", None, None, None),
        ("pssm/ps9_ber_dl_ovf.asnb", None, None, None),
    ] {
        let pssm = fixture_pssm(path);
        assert_eq!(pssm_unwritten(ArchiveEncoding::Xml, &pssm), xml, "{path}");
        assert_eq!(
            pssm_unwritten(ArchiveEncoding::Binary, &pssm),
            binary,
            "{path}"
        );
        assert_eq!(pssm_unwritten(ArchiveEncoding::Text, &pssm), text, "{path}");
        assert_eq!(
            pssm_unwritten(ArchiveEncoding::Nothing, &pssm),
            None,
            "{path}"
        );
    }
    // An element stated as NULL (a Seq-descr element, a feature's dbxref element: the reader's
    // empty references, line RG item D3c-3b): NCBI's three streams throw at it in the middle of
    // the archive (`NULL element while writing container`), as LOSAT's do (line RG item D3c-3c,
    // `writer::write_archive`): no gap. (The XML stream writes the feature's anonymous SEQUENCE
    // before its dbxref, line RG-F3 PR 2.)
    for path in [
        "pssm/d3b_descr_elem_null.asn.gz",
        "pssm/d3b_f_dbxref_null.asn.gz",
    ] {
        let pssm = fixture_pssm(path);
        assert!(pssm_holds_null_element(&pssm), "{path}");
        for encoding in [
            ArchiveEncoding::Xml,
            ArchiveEncoding::Binary,
            ArchiveEncoding::Text,
        ] {
            assert_eq!(pssm_unwritten(encoding, &pssm), None, "{path} {encoding:?}");
        }
    }
    assert!(!pssm_holds_null_element(&fixture_pssm(
        "db_search/db.ps.asn.rnd.6.asn"
    )));
}

// Line PS item PS7c: the archive block of a Lambda the results do not keep (`stats::tables::ancillary_kept`).
#[test]
fn a_block_of_a_lambda_below_zero_or_not_a_number_is_the_missing_one() {
    let params = |lambda| KarlinParams {
        lambda,
        k: 0.041,
        h: 0.14,
        alpha: 0.0,
        beta: 0.0,
    };
    for lambda in [f64::NAN, f64::NEG_INFINITY, -0.5] {
        let block = ka_block(&params(lambda), true);
        assert_eq!(
            (block.lambda, block.k, block.h, block.gapped),
            (-1.0, -1.0, -1.0, true)
        );
    }
    let block = ka_block(&params(0.267), false);
    assert_eq!(
        (block.lambda, block.k, block.h, block.gapped),
        (0.267, 0.041, 0.14, false)
    );
    // +inf is kept, as NCBI's test keeps it.
    assert_eq!(ka_block(&params(f64::INFINITY), true).lambda, f64::INFINITY);
}

/// A delayed member of another format whose read posts a replaced byte is written through
/// `Traced` in every archive (`pssm_written_as_is`), which posts it while it writes; a buffer of
/// the archive's own format is copied without a read.
#[test]
fn an_archive_writes_a_delayed_member_that_posts_through_traced() {
    let posting = fixture_pssm("pssm/ps9_xml_dl_ctrl.xml");
    assert!(pssm_written_as_is(ArchiveEncoding::Xml, &posting));
    assert!(!pssm_written_as_is(ArchiveEncoding::Binary, &posting));
    let posting = fixture_pssm("pssm/ps9_xml_bin_dl_ctrl.asnb");
    assert!(pssm_written_as_is(ArchiveEncoding::Binary, &posting));
    assert!(!pssm_written_as_is(ArchiveEncoding::Xml, &posting));
    let quiet = fixture_pssm("pssm/ps9_xml_dl_base.xml");
    assert!(pssm_written_as_is(ArchiveEncoding::Binary, &quiet));
}

// NCBI reference (598d8ae6): c++/src/serial/delaybuf.cpp:63-84
/// NCBI's ASN.1 text write of a PSSM file's object (`-out_pssm`, the search strategy's request)
/// reads each delay buffer of another format into the object once (`CDelayBuffer::Update`):
/// `PssmWithParameters::update_delayed` leaves the object that later writes write, the same
/// text without the read's posts, with nothing left that an archive cannot write; a buffer of
/// text, and one whose read fails, stays.
#[test]
fn a_text_write_reads_a_pssm_files_delay_buffers_of_other_formats_into_it() {
    use crate::report::serial::asn_input::Encoding;
    use crate::report::serial::pssm::asn_text_posting;
    // (file, the posts of its first text write, whether the update keeps the object as read)
    for (path, posts, kept) in [
        ("pssm/d3a_c11_delay.asn.gz", 0, true),
        ("pssm/ps9_xml_dl_base.xml", 0, false),
        ("pssm/ps9_xml_dl_ctrl.xml", 1, false),
        ("pssm/ps9_xml_bin_dl_ctrl.asnb", 1, false),
    ] {
        let pssm = fixture_pssm(path);
        let (text, first) = asn_text_posting(&pssm).unwrap();
        assert_eq!(first.len(), posts, "{path}");
        let mut updated = pssm.clone();
        updated.update_delayed(Encoding::Text);
        assert_eq!(updated == pssm, kept, "{path}");
        let (again, second) = asn_text_posting(&updated).unwrap();
        assert!(again == text && second.is_empty(), "{path}");
        for encoding in [
            ArchiveEncoding::Text,
            ArchiveEncoding::Xml,
            ArchiveEncoding::Binary,
        ] {
            assert_eq!(pssm_unwritten(encoding, &updated), None, "{path}");
        }
    }
    // A write of the buffer's own format copies it, so the buffer stays.
    let pssm = fixture_pssm("pssm/ps9_xml_dl_base.xml");
    let mut updated = pssm.clone();
    updated.update_delayed(Encoding::Xml);
    assert!(updated == pssm);
    // A read that fails (an INTEGER past Int4): the buffer stays, which NCBI's text write
    // stops at with its exception.
    let pssm = fixture_pssm("pssm/ps9_xml_dl_ovf.xml");
    let mut updated = pssm.clone();
    updated.update_delayed(Encoding::Text);
    assert!(updated == pssm);
    assert!(asn_text_posting(&updated).is_err());
}

// Line RG-F3 item O-2: a delayed member of another format whose read fails ends NCBI's writer with
// the read's exception, which leaves each frame the writer has open, after the bytes written before
// the member (NCBI BLAST+ 2.17.0, rows db.xml.ps.dl_ovf.outpssm and db.xml.ps.bin.dl_ovf.outpssm:
// psiblast's `-out_pssm` file and standard error, frozen with the official Linux x64 build). The
// XML read's `CUtilException` is rethrown at each frame (`encoder::Traced::rethrown`), the binary
// read's `CSerialException` gets each frame's name (`add_frame_infos`); no post precedes either.
#[test]
fn a_failed_read_of_a_delay_buffer_ends_the_pssm_file_with_ncbis_exception() {
    use crate::report::serial::encoder::Thrown;
    let frozen = |name: &str| {
        std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/db_search")
                .join(name),
        )
        .unwrap()
    };
    for (path, row) in [
        ("pssm/ps9_xml_dl_ovf.xml", "db.xml.ps.dl_ovf.outpssm"),
        ("pssm/ps9_ber_dl_ovf.asnb", "db.xml.ps.bin.dl_ovf.outpssm"),
    ] {
        let pssm = fixture_pssm(path);
        // The plain encoder stops there (`SerialError::DelayedRead`), so the file is written
        // through `Traced`.
        assert!(matches!(
            crate::report::serial::pssm::asn_text_with_karlin(&pssm, None),
            Err(crate::report::serial::encoder::SerialError::DelayedRead)
        ));
        let written = crate::report::serial::pssm::asn_text_written(&pssm, None).unwrap();
        assert_eq!(written.bytes, frozen(&format!("{row}.asn")), "{row}");
        assert!(written.bytes.ends_with(b"authors {"), "{row}");
        assert!(written.replaced.is_empty(), "{row}");
        let thrown = written.thrown.expect("the writer throws");
        assert!(matches!(*thrown, Thrown::DelayedRead(_)), "{row}");
        assert!(thrown.post(Some("psiblast")).is_empty(), "{row}");
        let stderr = [
            b"Error: ".as_slice(),
            &thrown.exception().ncbi_what(),
            b"\n",
        ]
        .concat();
        assert_eq!(
            String::from_utf8(stderr).unwrap(),
            String::from_utf8(frozen(&format!("{row}.err"))).unwrap(),
            "{row}"
        );
    }
}
