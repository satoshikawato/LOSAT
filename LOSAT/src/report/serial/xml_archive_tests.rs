//! The object mode of the XML encoder (`XmlEncoder::archive`, the `Blast4-archive` of
//! outfmt 11 with `ARCHIVE_FORMAT=xml`): the header and the end of a document, the tags of
//! named and anonymous CHOICE types, the `value=` attribute of ENUMERATED and BOOLEAN
//! values, the `_E` tags of primitive list elements, the named-type wrappers, hex, and the
//! streaming. The expected bytes are those of NCBI's `CObjectOStreamXml` as found in its
//! archives (`<Blast4-ka-block_gapped value="false"/>`, `<Dense-seg_starts_E>`,
//! `<Na-strand value="plus"/>`).

use super::encoder::{Encoder, ReplacedChar, SerialError, StreamPosition, Traced};
use super::xml::XmlEncoder;

const HEADER: &str = "<?xml version=\"1.0\"?>\n<!DOCTYPE Blast4-archive PUBLIC \"-//NCBI//NCBI Blast4/EN\" \"NCBI_Blast4.dtd\">";

fn archive_text(write: impl FnOnce(&mut XmlEncoder<'_, true>)) -> String {
    let mut buf = Vec::new();
    let mut enc = XmlEncoder::archive(&mut buf);
    write(&mut enc);
    assert_eq!(enc.error(), None);
    String::from_utf8(buf).unwrap()
}

#[test]
fn a_document_has_the_header_and_an_end_of_line() {
    let text = archive_text(|enc| {
        enc.begin_document("Blast4-archive");
        enc.begin_sequence("Blast4-archive");
        enc.end_sequence();
        enc.end_document();
    });
    assert_eq!(
        text,
        format!("{HEADER}\n<Blast4-archive></Blast4-archive>\n")
    );
}

#[test]
fn documents_are_concatenated_each_with_its_header() {
    let mut buf = Vec::new();
    for round in 0..2 {
        let mut enc = XmlEncoder::archive(&mut buf);
        enc.begin_document("Blast4-archive");
        enc.begin_sequence("Blast4-archive");
        enc.member("request", 0);
        enc.int(round);
        enc.end_member();
        enc.end_sequence();
        enc.end_document();
    }
    let one = |n: i32| {
        format!(
            "{HEADER}\n<Blast4-archive>\n  <Blast4-archive_request>{n}</Blast4-archive_request>\n</Blast4-archive>\n"
        )
    };
    assert_eq!(
        String::from_utf8(buf).unwrap(),
        format!("{}{}", one(0), one(1))
    );
}

#[test]
fn a_named_choice_writes_its_tag_and_an_anonymous_one_only_the_variant() {
    let text = archive_text(|enc| {
        enc.begin_sequence("Score");
        enc.member("id", 0);
        enc.variant("Object-id", "str", 1);
        enc.visible_str(b"score");
        enc.end_variant();
        enc.end_member();
        enc.member("value", 1);
        // `Score_value` is the open member's tag: an anonymous CHOICE.
        enc.variant("Score_value", "int", 1);
        enc.int(-7);
        enc.end_variant();
        enc.end_member();
        enc.end_sequence();
    });
    assert_eq!(
        text,
        "\n<Score>\
\n  <Score_id>\
\n    <Object-id>\
\n      <Object-id_str>score</Object-id_str>\
\n    </Object-id>\
\n  </Score_id>\
\n  <Score_value>\
\n    <Score_value_int>-7</Score_value_int>\
\n  </Score_value>\
\n</Score>"
    );
}

#[test]
fn enumerated_and_boolean_replace_the_close_of_the_innermost_tag() {
    let text = archive_text(|enc| {
        enc.begin_sequence("Seq-align");
        // On a member.
        enc.member("type", 0);
        enc.enumerated("partial", 3);
        enc.end_member();
        // On a member that holds a named ENUMERATED type.
        enc.member("strand", 1);
        enc.begin_named("Na-strand");
        enc.enumerated("minus", 2);
        enc.end_named();
        enc.end_member();
        // On the variant of a named CHOICE.
        enc.member("value", 2);
        enc.variant("Blast4-value", "boolean", 2);
        enc.boolean(true);
        enc.end_variant();
        enc.end_member();
        // On the variant of an anonymous CHOICE.
        enc.member("data", 3);
        enc.variant("Seq-align_data", "bool", 3);
        enc.boolean(false);
        enc.end_variant();
        enc.end_member();
        enc.member("after", 4);
        enc.int(1);
        enc.end_member();
        enc.end_sequence();
    });
    assert_eq!(
        text,
        "\n<Seq-align>\
\n  <Seq-align_type value=\"partial\"/>\
\n  <Seq-align_strand>\
\n    <Na-strand value=\"minus\"/>\
\n  </Seq-align_strand>\
\n  <Seq-align_value>\
\n    <Blast4-value>\
\n      <Blast4-value_boolean value=\"true\"/>\
\n    </Blast4-value>\
\n  </Seq-align_value>\
\n  <Seq-align_data>\
\n    <Seq-align_data_bool value=\"false\"/>\
\n  </Seq-align_data>\
\n  <Seq-align_after>1</Seq-align_after>\
\n</Seq-align>"
    );
}

#[test]
fn primitive_list_elements_have_an_e_tag_and_the_others_do_not() {
    let text = archive_text(|enc| {
        enc.begin_sequence("Dense-seg");
        enc.member("starts", 3);
        enc.begin_list(false);
        for value in [0, -1] {
            enc.element();
            enc.int(value);
        }
        enc.end_list();
        enc.end_member();
        enc.member("stats", 4);
        enc.begin_list(false);
        for value in [&b"a<b"[..], b"c"] {
            enc.element();
            enc.visible_str(value);
        }
        enc.end_list();
        enc.end_member();
        // Elements of a SEQUENCE and of a named ENUMERATED type write their own tags.
        enc.member("ids", 2);
        enc.begin_list(false);
        enc.element();
        enc.variant("Object-id", "id", 0);
        enc.int(5);
        enc.end_variant();
        enc.end_list();
        enc.end_member();
        enc.member("strands", 5);
        enc.begin_list(false);
        for name in ["plus", "minus"] {
            enc.element();
            enc.begin_named("Na-strand");
            enc.enumerated(name, 0);
            enc.end_named();
        }
        enc.end_list();
        enc.end_member();
        // An anonymous BOOLEAN element is `<Parent_E value="x"/>`.
        enc.member("flags", 6);
        enc.begin_list(false);
        enc.element();
        enc.boolean(true);
        enc.element();
        enc.boolean(false);
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
    });
    assert_eq!(
        text,
        "\n<Dense-seg>\
\n  <Dense-seg_starts>\
\n    <Dense-seg_starts_E>0</Dense-seg_starts_E>\
\n    <Dense-seg_starts_E>-1</Dense-seg_starts_E>\
\n  </Dense-seg_starts>\
\n  <Dense-seg_stats>\
\n    <Dense-seg_stats_E>a&lt;b</Dense-seg_stats_E>\
\n    <Dense-seg_stats_E>c</Dense-seg_stats_E>\
\n  </Dense-seg_stats>\
\n  <Dense-seg_ids>\
\n    <Object-id>\
\n      <Object-id_id>5</Object-id_id>\
\n    </Object-id>\
\n  </Dense-seg_ids>\
\n  <Dense-seg_strands>\
\n    <Na-strand value=\"plus\"/>\
\n    <Na-strand value=\"minus\"/>\
\n  </Dense-seg_strands>\
\n  <Dense-seg_flags>\
\n    <Dense-seg_flags_E value=\"true\"/>\
\n    <Dense-seg_flags_E value=\"false\"/>\
\n  </Dense-seg_flags>\
\n</Dense-seg>"
    );
}

#[test]
fn an_empty_container_is_on_one_line_and_a_named_one_has_its_own_tag() {
    let text = archive_text(|enc| {
        enc.begin_sequence("Blast4-get-search-results-reply");
        enc.member("alignments", 0);
        enc.begin_named("Seq-align-set");
        enc.begin_list(true);
        enc.end_list();
        enc.end_named();
        enc.end_member();
        enc.member("masks", 2);
        enc.begin_list(false);
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
    });
    assert_eq!(
        text,
        "\n<Blast4-get-search-results-reply>\
\n  <Blast4-get-search-results-reply_alignments>\
\n    <Seq-align-set></Seq-align-set>\
\n  </Blast4-get-search-results-reply_alignments>\
\n  <Blast4-get-search-results-reply_masks></Blast4-get-search-results-reply_masks>\
\n</Blast4-get-search-results-reply>"
    );
}

#[test]
fn octets_are_unwrapped_upper_case_hex_and_a_string_store_is_escaped_text() {
    let data: Vec<u8> = (0..=255u8).cycle().take(300).collect();
    let text = archive_text(|enc| {
        enc.begin_sequence("Seq-inst");
        enc.member("seq-data", 6);
        enc.variant("Seq-data", "ncbi2na", 2);
        enc.begin_named("NCBI2na");
        enc.octets(&data);
        enc.end_named();
        enc.end_variant();
        enc.end_member();
        enc.member("eaa", 7);
        enc.variant("Seq-data", "ncbieaa", 7);
        enc.begin_named("NCBIeaa");
        enc.string_store(b"AC&*");
        enc.end_named();
        enc.end_variant();
        enc.end_member();
        enc.end_sequence();
    });
    let hex: String = data.iter().map(|byte| format!("{byte:02X}")).collect();
    assert_eq!(hex.len(), 600);
    assert_eq!(
        text,
        format!(
            "\n<Seq-inst>\
\n  <Seq-inst_seq-data>\
\n    <Seq-data>\
\n      <Seq-data_ncbi2na>\
\n        <NCBI2na>{hex}</NCBI2na>\
\n      </Seq-data_ncbi2na>\
\n    </Seq-data>\
\n  </Seq-inst_seq-data>\
\n  <Seq-inst_eaa>\
\n    <Seq-data>\
\n      <Seq-data_ncbieaa>\
\n        <NCBIeaa>AC&amp;*</NCBIeaa>\
\n      </Seq-data_ncbieaa>\
\n    </Seq-data>\
\n  </Seq-inst_eaa>\
\n</Seq-inst>"
        )
    );
}

#[test]
fn a_string_byte_from_0x7f_up_is_written_as_it_is_and_a_control_byte_is_an_error() {
    // The archive's stream has no encoding: `WriteEncodedChar` takes its first branch
    // (objostrxml.cpp:559-565) and passes the byte on (UTF-8 here, 0x7F, a lone 0xFF).
    let mut buf = Vec::new();
    let mut enc = XmlEncoder::archive(&mut buf);
    enc.begin_sequence("T");
    enc.member("m", 0);
    enc.visible_str(b"caf\xc3\xa9 \x7f \xff <&> a\tb");
    enc.end_member();
    enc.end_sequence();
    assert_eq!(enc.error(), None);
    assert_eq!(
        buf,
        b"\n<T>\n  <T_m>caf\xc3\xa9 \x7f \xff &lt;&amp;&gt; a&#x9;b</T_m>\n</T>".to_vec()
    );
    // 0x01-0x1F other than TAB, LF and CR: NCBI writes `#` and a diagnostic, an error of the
    // encoder by itself.
    let mut buf = Vec::new();
    let mut enc = XmlEncoder::archive(&mut buf);
    enc.visible_str(b"a\x01b");
    assert_eq!(enc.error(), Some(SerialError::NonVisibleChar(1)));
    // Through `Traced`: `#` and a `ReplacedChar` at the stream's line, without the string (NCBI
    // 2.17.0 convert2blastmask, probe RG K/k4/probe case c_xml_tab: `Bad char [0x1] in string
    // at line 8` for `-masking_options` on line 8 of a `maskinfo_xml` document).
    let mut buf = Vec::new();
    let mut enc = Traced::new(XmlEncoder::ncbi_object(&mut buf, "NCBI-T"), 0);
    enc.begin_document("T");
    enc.begin_sequence("T");
    enc.member("m", 0);
    enc.visible_str(b"a\tb\x01c\x7f");
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
    assert_eq!(enc.error(), None);
    let replaced = enc.take_replaced_chars();
    assert_eq!(
        replaced,
        vec![ReplacedChar {
            trace: "T.m".to_string(),
            byte: 1,
            position: StreamPosition::Line(4),
            string: std::sync::Arc::from(&[][..]),
        }]
    );
    assert_eq!(
        replaced[0].message(),
        b"T.m\nBad char [0x1] in string at line 4".to_vec()
    );
    assert!(buf.ends_with(b"<T_m>a&#x9;b#c\x7f</T_m>\n</T>\n"));
    // The US-ASCII stream of outfmt 5 decodes a byte >= 0x80 and throws where the character is
    // not ASCII (`report::blast_xml::tests` has the cases).
    let mut buf = Vec::new();
    let mut enc = XmlEncoder::object(&mut buf);
    enc.visible_str(b"caf\xc3\xa9");
    assert_eq!(
        enc.error(),
        Some(SerialError::Utf8(
            crate::ncbi_str::utf8::Utf8Error::SymbolToChar
        ))
    );
    assert_eq!(buf, b"caf".to_vec());
}

#[test]
fn a_string_list_variant_is_written_as_the_default_calls_write_it() {
    let items = ["a&b", "c"];
    let listed = archive_text(|enc| {
        enc.begin_sequence("Holder");
        enc.member("target", 0);
        enc.variant_string_list("Target", "subjects", 1, items);
        enc.end_member();
        enc.end_sequence();
    });
    let by_hand = archive_text(|enc| {
        enc.begin_sequence("Holder");
        enc.member("target", 0);
        enc.variant("Target", "subjects", 1);
        enc.begin_list(false);
        for item in items {
            enc.element();
            enc.visible_str(item.as_bytes());
        }
        enc.end_list();
        enc.end_variant();
        enc.end_member();
        enc.end_sequence();
    });
    assert_eq!(listed, by_hand);
    assert!(
        listed.contains("<Target_subjects_E>a&amp;b</Target_subjects_E>"),
        "{listed}"
    );
}

/// A document that uses every kind of call of the object mode.
fn sample<E: Encoder>(enc: &mut E) {
    enc.begin_document("Blast4-archive");
    enc.begin_sequence("Blast4-archive");
    enc.member("request", 0);
    enc.begin_sequence("Blast4-request");
    enc.member("ident", 0);
    enc.visible_str(b"2.17.0+");
    enc.end_member();
    enc.member("body", 1);
    enc.variant("Blast4-request-body", "queue-search", 8);
    enc.begin_sequence("Blast4-queue-search-request");
    enc.member("algorithm-options", 5);
    enc.begin_named("Blast4-parameters");
    enc.begin_list(false);
    for flag in [true, false] {
        enc.element();
        enc.begin_sequence("Blast4-parameter");
        enc.member("name", 0);
        enc.visible_str(b"x");
        enc.end_member();
        enc.member("value", 1);
        enc.variant("Blast4-value", "boolean", 2);
        enc.boolean(flag);
        enc.end_variant();
        enc.end_member();
        enc.end_sequence();
    }
    enc.element();
    enc.begin_sequence("Blast4-parameter");
    enc.member("value", 1);
    enc.variant("Blast4-value", "cutoff", 3);
    enc.variant("Blast4-cutoff", "e-value", 0);
    enc.real(10.0);
    enc.end_variant();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_named();
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.member("results", 1);
    enc.begin_sequence("Blast4-get-search-results-reply");
    enc.member("starts", 0);
    enc.begin_list(false);
    for value in 0..5 {
        enc.element();
        enc.int(value);
    }
    enc.end_list();
    enc.end_member();
    enc.member("mask", 1);
    enc.begin_sequence("Blast4-mask");
    enc.member("frame", 1);
    enc.begin_named("Blast4-frame-type");
    enc.enumerated("plus1", 1);
    enc.end_named();
    enc.end_member();
    enc.member("seq", 2);
    enc.variant("Seq-data", "ncbi4na", 3);
    enc.begin_named("NCBI4na");
    enc.octets(&[0x12, 0xf8, 0xa0]);
    enc.end_named();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.member("gapped", 2);
    enc.boolean(true);
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
}

#[test]
fn flushing_after_every_call_gives_the_bytes_of_one_final_flush() {
    let whole = archive_text(|enc| sample(enc));
    assert!(whole.contains("<Blast4-value_boolean value=\"false\"/>"));
    assert!(whole.ends_with("</Blast4-archive>\n"));

    let mut buf = Vec::new();
    let mut flushed = Vec::new();
    {
        let mut enc = Flushing {
            inner: XmlEncoder::archive(&mut buf),
            sink: &mut flushed,
        };
        sample(&mut enc);
        assert_eq!(enc.inner.error(), None);
        enc.inner.flush_to(enc.sink, true).unwrap();
        assert_eq!(enc.inner.buffered(), 0);
    }
    assert_eq!(String::from_utf8(flushed).unwrap(), whole);
    assert!(buf.is_empty());
}

#[test]
fn a_flush_keeps_the_last_byte_that_an_attribute_replaces() {
    let mut buf = Vec::new();
    let mut sink = Vec::new();
    let mut enc = XmlEncoder::archive(&mut buf);
    enc.begin_sequence("Seq-inst");
    enc.member("repr", 0);
    enc.flush_to(&mut sink, false).unwrap();
    assert!(enc.buffered() >= 1);
    enc.enumerated("raw", 2);
    enc.end_member();
    enc.end_sequence();
    enc.flush_to(&mut sink, true).unwrap();
    assert_eq!(
        String::from_utf8(sink).unwrap(),
        "\n<Seq-inst>\n  <Seq-inst_repr value=\"raw\"/>\n</Seq-inst>"
    );
}

#[test]
fn a_document_that_starts_after_other_bytes_flushes_only_its_own() {
    let mut buf = b"earlier".to_vec();
    let mut sink = Vec::new();
    let mut enc = XmlEncoder::archive(&mut buf);
    enc.begin_document("Blast4-archive");
    enc.begin_sequence("Blast4-archive");
    enc.end_sequence();
    enc.end_document();
    enc.flush_to(&mut sink, true).unwrap();
    assert_eq!(buf, b"earlier");
    assert_eq!(
        String::from_utf8(sink).unwrap(),
        format!("{HEADER}\n<Blast4-archive></Blast4-archive>\n")
    );
}

#[test]
fn the_other_pieces_still_refuse_what_the_archive_needs() {
    let refused = |build: &dyn Fn(&mut XmlEncoder<'_>)| {
        let mut buf = Vec::new();
        let mut enc = XmlEncoder::object(&mut buf);
        build(&mut enc);
        enc.error()
    };
    assert_eq!(
        refused(&|enc| enc.enumerated("plus", 1)),
        Some(SerialError::XmlUnsupported("ENUMERATED value"))
    );
    assert_eq!(
        refused(&|enc| enc.boolean(true)),
        Some(SerialError::XmlUnsupported("BOOLEAN"))
    );
    assert_eq!(
        refused(&|enc| enc.variant("Seq-id", "local", 0)),
        Some(SerialError::XmlUnsupported("CHOICE"))
    );
    assert_eq!(
        refused(&|enc| enc.octets(&[1])),
        Some(SerialError::OctetString)
    );
    assert_eq!(
        refused(&|enc| enc.variant_string_list("Target", "subjects", 1, ["a"])),
        Some(SerialError::XmlUnsupported("CHOICE"))
    );
    // `begin_named` and `element` write nothing there.
    let mut buf = Vec::new();
    let mut enc = XmlEncoder::object(&mut buf);
    enc.begin_named("Seq-descr");
    enc.element();
    enc.end_named();
    assert_eq!(enc.buffered(), 0);
    assert!(buf.is_empty());
}

/// An encoder that hands the inner encoder's final bytes to a sink after every call.
struct Flushing<'s, E: Encoder> {
    inner: E,
    sink: &'s mut Vec<u8>,
}

impl<E: Encoder> Flushing<'_, E> {
    fn flush(&mut self) {
        self.inner.flush_to(&mut *self.sink, false).unwrap();
    }
}

impl<E: Encoder> Encoder for Flushing<'_, E> {
    fn begin_document(&mut self, type_name: &str) {
        self.inner.begin_document(type_name);
        self.flush();
    }
    fn end_document(&mut self) {
        self.inner.end_document();
        self.flush();
    }
    fn begin_sequence(&mut self, type_name: &str) {
        self.inner.begin_sequence(type_name);
        self.flush();
    }
    fn end_sequence(&mut self) {
        self.inner.end_sequence();
        self.flush();
    }
    fn begin_named(&mut self, type_name: &str) {
        self.inner.begin_named(type_name);
        self.flush();
    }
    fn end_named(&mut self) {
        self.inner.end_named();
        self.flush();
    }
    fn member(&mut self, name: &str, tag: u8) {
        self.inner.member(name, tag);
        self.flush();
    }
    fn end_member(&mut self) {
        self.inner.end_member();
        self.flush();
    }
    fn begin_list(&mut self, set: bool) {
        self.inner.begin_list(set);
        self.flush();
    }
    fn element(&mut self) {
        self.inner.element();
        self.flush();
    }
    fn end_list(&mut self) {
        self.inner.end_list();
        self.flush();
    }
    fn variant(&mut self, type_name: &str, name: &str, tag: u8) {
        self.inner.variant(type_name, name, tag);
        self.flush();
    }
    fn end_variant(&mut self) {
        self.inner.end_variant();
        self.flush();
    }
    fn int(&mut self, value: i32) {
        self.inner.int(value);
        self.flush();
    }
    fn big_int(&mut self, value: i64) {
        self.inner.big_int(value);
        self.flush();
    }
    fn enumerated(&mut self, name: &str, value: i32) {
        self.inner.enumerated(name, value);
        self.flush();
    }
    fn named_int(&mut self, name: &str, value: i32) {
        self.inner.named_int(name, value);
        self.flush();
    }
    fn boolean(&mut self, value: bool) {
        self.inner.boolean(value);
        self.flush();
    }
    fn visible_str(&mut self, value: &[u8]) {
        self.inner.visible_str(value);
        self.flush();
    }
    fn real(&mut self, value: f64) {
        self.inner.real(value);
        self.flush();
    }
    fn octets(&mut self, value: &[u8]) {
        self.inner.octets(value);
        self.flush();
    }
    fn string_store(&mut self, value: &[u8]) {
        self.inner.string_store(value);
        self.flush();
    }
    fn null(&mut self) {
        self.inner.null();
        self.flush();
    }
    fn utf8_str(&mut self, value: &[u8]) {
        self.inner.utf8_str(value);
        self.flush();
    }
    fn delayed(&mut self, text: &[u8]) {
        self.inner.delayed(text);
        self.flush();
    }
    fn error(&self) -> Option<SerialError> {
        self.inner.error()
    }
}

// PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits (review C of line RG
// D1a-4, SF1): `PutInt8` of
// -2^63, which NCBI's `COStreamBuffer` writes as `-0` (a signed overflow, strbuffer.cpp:1002-1010),
// is its digits in both encoders of the archive and of the reports' XML.
#[test]
fn a_big_integer_of_minus_two_to_the_63_is_its_digits() {
    for value in [i64::MIN, i64::MIN + 1, -1, 0, i64::MAX] {
        let text = archive_text(|enc| enc.big_int(value));
        assert_eq!(text, value.to_string());
        let mut buf = Vec::new();
        let mut enc = super::asn_text::AsnTextEncoder::new(&mut buf);
        enc.big_int(value);
        assert_eq!(enc.error(), None);
        assert_eq!(String::from_utf8(buf).unwrap(), value.to_string());
    }
}

/// A fixture of LOSAT's tests (`tests/fixtures/<path>`).
fn test_fixture(path: &str) -> Vec<u8> {
    std::fs::read(format!(
        "{}/tests/fixtures/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

fn gunzip(bytes: Vec<u8>) -> Vec<u8> {
    let mut out = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::GzDecoder::new(&bytes[..]), &mut out).unwrap();
    out
}

/// The `Blast4-archive` frame of NCBI's archive of a PSSM search around its request's PSSM query
/// (`queries pssm`, `pssm::write_value`), as deep as NCBI's (its ident, program and service are
/// NCBI's), closed to the end of the document (`close`) or left after the member `queries`.
fn write_pssm_query<E: Encoder>(enc: &mut E, pssm: &super::pssm::PssmWithParameters, close: bool) {
    enc.begin_document("Blast4-archive");
    enc.begin_sequence("Blast4-archive");
    enc.member("request", 0);
    enc.begin_sequence("Blast4-request");
    enc.member("ident", 0);
    enc.visible_str(b"2.17.0+");
    enc.end_member();
    enc.member("body", 1);
    enc.variant("Blast4-request-body", "queue-search", 8);
    enc.begin_sequence("Blast4-queue-search-request");
    enc.member("program", 0);
    enc.visible_str(b"blastp");
    enc.end_member();
    enc.member("service", 1);
    enc.visible_str(b"psi");
    enc.end_member();
    enc.member("queries", 2);
    enc.variant("Blast4-queries", "pssm", 0);
    super::pssm::write_value(enc, pssm);
    enc.end_variant();
    enc.end_member();
    if !close {
        return;
    }
    enc.end_sequence();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
}

/// The request's PSSM query (`Blast4-queries pssm`) of an archive with `ARCHIVE_FORMAT=xml`, its
/// parts written by the writer of the type table's values (`pssm::write_value`, `asn_types::write`),
/// which names an anonymous SEQUENCE or CHOICE by its empty name: the bytes of NCBI BLAST+ 2.17.0
/// up to the end of that member (`psiblast -db tests/db_search/paja -in_pssm <case>
/// -comp_based_stats 0 -outfmt 11 -num_threads 1` with `ARCHIVE_FORMAT=xml`, line RG-F3 PR 2
/// probe X, 2026-11-23; tests/fixtures/outfmt11/psi_xml_pssm_<case>.head.gz). The cases hold the
/// anonymous CHOICE `Seq-annot.data` in the query's annot (`locs`; `ftable` with
/// `Variation-ref.data`; `align` with `Seq-align.segs`, whose sizes NCBI's object manager warns
/// of before the archive), `User-field.data` in its descr (whose `str` has the byte 0x01: NCBI's
/// 803.7 post, `psi_xml_pssm_d3a_c10_user.err`), and in its inst the anonymous SEQUENCE
/// `Int-fuzz.range` and the anonymous CHOICE `Seq-hist.deleted` (tests/fixtures/pssm/
/// d5x_fuzz_hist.asn.gz: d3a_c17_locs with `fuzz range { max 5, min 1 }` and `hist { deleted
/// bool TRUE }`). The `Blast4-archive` frame around it is this test's, as deep as NCBI's (its
/// ident, program and service are NCBI's lines). LOSAT's psiblast refuses outfmt 11 and its
/// tblastn `-in_pssm` refuses an archive of a search (NCBI's tblastn throws), so no run of the
/// command line writes this member yet.
#[test]
fn a_pssm_query_s_anonymous_types_are_nested_as_ncbi_s_xml_archive_writes_them() {
    let (fixture, gunzip) = (test_fixture, gunzip);
    let mut failures = Vec::new();
    for (case, posts) in [
        ("d3a_c17_locs", Vec::new()),
        ("d3a_c13_variation", Vec::new()),
        ("d3om_a06_denseg_all4", Vec::new()),
        ("d5x_fuzz_hist", Vec::new()),
        (
            "d3a_c10_user",
            fixture("outfmt11/psi_xml_pssm_d3a_c10_user.err"),
        ),
    ] {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(fixture(
            &format!("pssm/{case}.asn.gz"),
        )))
        .unwrap();
        let mut buf = Vec::new();
        let mut enc = Traced::xml(XmlEncoder::archive(&mut buf), 0);
        write_pssm_query(&mut enc, &pssm, true);
        if let Some(error) = enc.error() {
            failures.push(format!("{case}: {error:?}"));
        }
        // The archive's writer takes `Traced`'s path (which posts) for the case with a post.
        let as_is = crate::report::archive::pssm_written_as_is(
            crate::report::archive::ArchiveEncoding::Xml,
            &pssm,
        );
        if as_is != posts.is_empty() {
            failures.push(format!("{case}: written as it is: {as_is}"));
        }
        let written: Vec<u8> = enc
            .take_replaced_chars()
            .iter()
            .flat_map(|replaced| replaced.post(Some("psiblast")))
            .collect();
        drop(enc);
        let ncbi = gunzip(fixture(&format!("outfmt11/psi_xml_pssm_{case}.head.gz")));
        // The first line that differs, with its number.
        let differs = ncbi
            .split_inclusive(|&byte| byte == b'\n')
            .zip(buf.split_inclusive(|&byte| byte == b'\n'))
            .enumerate()
            .find(|(_, (ncbi, losat))| ncbi != losat)
            .map(|(line, (ncbi, losat))| {
                format!(
                    "line {}: NCBI {:?}, LOSAT {:?}",
                    line + 1,
                    String::from_utf8_lossy(ncbi),
                    String::from_utf8_lossy(losat)
                )
            });
        if let Some(differs) = differs {
            failures.push(format!("{case}: {differs}"));
        } else if buf.len() < ncbi.len() {
            failures.push(format!(
                "{case}: {} bytes of NCBI's {}",
                buf.len(),
                ncbi.len()
            ));
        }
        if written != posts {
            failures.push(format!(
                "{case}: posts {:?}, NCBI's {:?}",
                String::from_utf8_lossy(&written),
                String::from_utf8_lossy(&posts)
            ));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The request's PSSM query of an archive with `ARCHIVE_FORMAT=bin` (`AsnBinaryEncoder`, its
/// parts written by `pssm::write_value` and `asn_types::write`): the bytes of NCBI BLAST+ 2.17.0
/// from the start of the archive to the end of the member `queries` (`psiblast -db
/// tests/db_search/paja -in_pssm <case> -comp_based_stats 0 -outfmt 11 -num_threads 1` with
/// `ARCHIVE_FORMAT=bin`, line RG-F3 PR 2 probe XB, 2026-11-23; NCBI's whole archive in
/// tests/fixtures/outfmt11/psi_bin_pssm_<case>.bin.gz). `d3a_c10_user`'s `User-field` `str`
/// is a UTF8String with `é` and the byte 0x01, which NCBI writes with a VisibleString's tag and
/// does not post; `d3a_c17_locs` holds a NULL Seq-loc in the anonymous CHOICE `Seq-annot.data`.
#[test]
fn a_pssm_query_is_written_in_binary_as_ncbi_s_archive_writes_it() {
    let mut failures = Vec::new();
    for case in ["d3a_c10_user", "d3a_c17_locs"] {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(test_fixture(
            &format!("pssm/{case}.asn.gz"),
        )))
        .unwrap();
        let mut buf = Vec::new();
        let mut enc = Traced::new(super::asn_ber::AsnBinaryEncoder::new(&mut buf), 0);
        write_pssm_query(&mut enc, &pssm, false);
        if let Some(error) = enc.error() {
            failures.push(format!("{case}: {error:?}"));
        }
        if !enc.take_replaced_chars().is_empty() {
            failures.push(format!("{case}: a byte replaced"));
        }
        drop(enc);
        if !crate::report::archive::pssm_written_as_is(
            crate::report::archive::ArchiveEncoding::Binary,
            &pssm,
        ) {
            failures.push(format!("{case}: not written as it is"));
        }
        let ncbi = gunzip(test_fixture(&format!(
            "outfmt11/psi_bin_pssm_{case}.bin.gz"
        )));
        let around = |bytes: &[u8], at: usize| bytes[at..(at + 16).min(bytes.len())].to_vec();
        match buf
            .iter()
            .zip(&ncbi)
            .position(|(losat, ncbi)| losat != ncbi)
        {
            Some(at) => failures.push(format!(
                "{case}: from byte {at}: NCBI {:02x?}, LOSAT {:02x?}",
                around(&ncbi, at),
                around(&buf, at)
            )),
            // NCBI's next member is `subject` (`a3 80`).
            None if !ncbi[buf.len().min(ncbi.len())..].starts_with(&[0xa3, 0x80]) => failures.push(
                format!("{case}: {} bytes of NCBI's {}", buf.len(), ncbi.len()),
            ),
            None => {}
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
