//! Unit tests of the generic writer (`super::write`): a document read by the generic reader and
//! written back, against the bytes NCBI BLAST+ 2.17.0 wrote for the same values (the probes
//! `c01`-`c19` of SURVEY_D3 of line RG series F, `F/f1b2/d3/ncbi/*.pssm`: psiblast `-out_pssm`
//! writes the PSSM it read).

use super::read::read_document;
use super::table::named;
use super::write::{write_document, write_value};
use super::Value;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, SerialError};

fn read(root: u16, text: &str) -> Value {
    let (value, posts) = read_document(text.as_bytes(), root, Names::TEST);
    assert!(posts.is_empty(), "{posts:?}");
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

/// The ASN.1 text that NCBI writes of the document `text` read as `root`.
fn written(root: u16, text: &str) -> String {
    let value = read(root, text);
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    write_document(&mut enc, root, &value);
    assert_eq!(enc.error(), None);
    drop(enc);
    String::from_utf8(out).expect("UTF-8")
}

/// c10: a REAL is one line of its normalized mantissa and exponent, a SEQUENCE OF INTEGER one
/// element per line, a UTF8String as read, an odd OCTET STRING with its last digit as a high
/// nibble, `+7` and `007` as their numbers.
#[test]
fn a_user_object_is_written_by_its_types() {
    let text = "Seqdesc ::= user { class \"c\", type str \"t\", data {\n\
        { label str \"real\", data real { 140, 10, -3 } },\n\
        { label str \"ints\", data ints { 1, 10, -3 } },\n\
        { label str \"u\", data str \"\u{e9}\" },\n\
        { label id +7, num 007, data os 'abc'H },\n\
        { label str \"b\", data bool TRUE } } }";
    assert_eq!(
        written(named::SEQDESC, text),
        "Seqdesc ::= user {\n\
        \x20 class \"c\",\n\
        \x20 type str \"t\",\n\
        \x20 data {\n\
        \x20   {\n\
        \x20     label str \"real\",\n\
        \x20     data real { 14, 10, -2 }\n\
        \x20   },\n\
        \x20   {\n\
        \x20     label str \"ints\",\n\
        \x20     data ints {\n\
        \x20       1,\n\
        \x20       10,\n\
        \x20       -3\n\
        \x20     }\n\
        \x20   },\n\
        \x20   {\n\
        \x20     label str \"u\",\n\
        \x20     data str \"\u{e9}\"\n\
        \x20   },\n\
        \x20   {\n\
        \x20     label id 7,\n\
        \x20     num 7,\n\
        \x20     data os 'ABC0'H\n\
        \x20   },\n\
        \x20   {\n\
        \x20     label str \"b\",\n\
        \x20     data bool TRUE\n\
        \x20   }\n\
        \x20 }\n\
        }\n"
    );
}

/// c10: a UTF8String keeps a control byte (`eFNP_Allow`) and loses a zero byte
/// (objostrasn.cpp:436-438), where a VisibleString's reader made `#` of both.
#[test]
fn a_utf8_string_keeps_its_bytes_but_zero() {
    let value = read(
        named::USER_FIELD,
        "User-field ::= { label id 1, data str \"\u{e9}\x01\" }",
    );
    let Value::Sequence(members) = &value else {
        panic!("{value:?}")
    };
    let mut with_zero = value.clone();
    if let Value::Sequence(members) = &mut with_zero {
        if let Value::Choice(_, text) = &mut members[1].1 {
            **text = Value::Str(b"a\0b".to_vec());
        }
    }
    assert_eq!(members.len(), 2);
    let write = |value: &Value| {
        let mut out = Vec::new();
        write_document(&mut AsnTextEncoder::new(&mut out), named::USER_FIELD, value);
        String::from_utf8(out).expect("UTF-8")
    };
    assert_eq!(
        write(&value),
        "User-field ::= {\n  label id 1,\n  data str \"\u{e9}\x01\"\n}\n"
    );
    assert_eq!(
        write(&with_zero),
        "User-field ::= {\n  label id 1,\n  data str \"ab\"\n}\n"
    );
}

/// c01: the numbers of INTEGERs with named values are written by their names; one without a name
/// by its number.
#[test]
fn a_named_integer_is_written_by_its_name() {
    assert_eq!(
        written(
            named::SEQDESC,
            "Seqdesc ::= molinfo { biomol 8, tech +002, completeness 099 }"
        ),
        "Seqdesc ::= molinfo {\n  biomol peptide,\n  tech est,\n  completeness 99\n}\n"
    );
}

/// c02, c03, c17: ENUMERATED numbers by their names, a DEFAULT member that was read written even
/// with its default value, `-0` as 0, a TSeqPos `+29` as 29, an Int8 gi past Int4.
#[test]
fn enumerated_values_and_defaults_read_are_written() {
    assert_eq!(
        written(
            named::SEQ_INST,
            "Seq-inst ::= { repr raw, mol aa, length 3, fuzz pct -0, topology 2, strand 1 }"
        ),
        "Seq-inst ::= {\n  repr raw,\n  mol aa,\n  length 3,\n  fuzz pct 0,\n  topology circular,\n  strand ss\n}\n"
    );
    assert_eq!(
        written(
            named::SEQ_INST,
            "Seq-inst ::= { repr raw, mol aa, topology linear }"
        ),
        "Seq-inst ::= {\n  repr raw,\n  mol aa,\n  topology linear\n}\n"
    );
    assert_eq!(
        written(
            named::SEQ_LOC,
            "Seq-loc ::= int { from 0, to +29, strand 2, id gi 5000000000 }"
        ),
        "Seq-loc ::= int {\n  from 0,\n  to 29,\n  strand minus,\n  id gi 5000000000\n}\n"
    );
    assert_eq!(
        written(named::SEQ_LOC, "Seq-loc ::= null NULL"),
        "Seq-loc ::= null NULL\n"
    );
}

/// c14: an empty SET OF is `{`, an end of line and `}` at its indentation.
#[test]
fn an_empty_list_is_a_block() {
    assert_eq!(
        written(
            named::USER_OBJECT,
            "User-object ::= { type str \"t\", data { } }"
        ),
        "User-object ::= {\n  type str \"t\",\n  data {\n  }\n}\n"
    );
}

/// c11: `Auth-list.names` is written as NCBI's reader kept it (comments, spacing and `level 1`),
/// the members around it as the writer writes them.
#[test]
fn a_delayed_member_is_copied_as_read() {
    let text = "Auth-list ::= { names  std {\n\
        { name name { last \"Smith\", -- note --\n first \"J\" } , level 1 }\n\
        } -- after --\n , affil str \"a\" }";
    assert_eq!(
        written(named::AUTH_LIST, text),
        "Auth-list ::= {\n  names std {\n\
        { name name { last \"Smith\", -- note --\n first \"J\" } , level 1 }\n\
        },\n  affil str \"a\"\n}\n"
    );
}

/// A delayed member read from ASN.1 binary is read again and written by its values
/// (`CDelayBuffer::Update`); where that read posts a message (a VisibleString byte replaced), NCBI
/// reports it while it writes, which LOSAT's writer does not write. (The bytes: `datatool` 2.24.0
/// `-e` of `Auth-list ::= { names str { "J" }, affil str "a" }`, read back by `datatool -d` as the
/// text below; `D/d3/c3a/probe/delay` of line RG.)
#[test]
fn a_delayed_member_read_from_binary_is_written_by_its_values() {
    let bytes: &[u8] = &[
        0x30, 0x80, 0xa0, 0x80, 0xa2, 0x80, 0x30, 0x80, 0x1a, 0x01, 0x4a, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0xa1, 0x80, 0xa0, 0x80, 0x1a, 0x01, 0x61, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    let write = |bytes: &[u8]| {
        let (value, posts) = super::read::read_ber_value(bytes, named::AUTH_LIST);
        assert!(posts.is_empty(), "{posts:?}");
        let value = value.expect("read");
        assert!(matches!(
            &value,
            Value::Sequence(members)
                if matches!(members[0].1, Value::Delayed { encoding: super::Encoding::Binary, .. })
        ));
        let mut out = Vec::new();
        let mut enc = AsnTextEncoder::new(&mut out);
        write_document(&mut enc, named::AUTH_LIST, &value);
        let error = enc.error();
        drop(enc);
        (String::from_utf8(out).expect("UTF-8"), error)
    };
    assert_eq!(
        write(bytes),
        (
            "Auth-list ::= {\n  names str {\n    \"J\"\n  },\n  affil str \"a\"\n}\n".to_owned(),
            None
        )
    );
    // A byte that the read replaces: posted (803.7) under the read's own frames at its offset in
    // the buffer, with the string, before the member is written with `#` (`delayed_posts`); the
    // encoder alone reports it as a replaced byte of its own.
    let mut bad = bytes.to_vec();
    bad[10] = 0xe9;
    assert_eq!(write(&bad).1, Some(SerialError::NonVisibleChar(0xe9)));
    let (value, _) = super::read::read_ber_value(&bad, named::AUTH_LIST);
    let mut out = Vec::new();
    let mut enc = crate::report::serial::encoder::Traced::new(AsnTextEncoder::new(&mut out), 0);
    write_document(&mut enc, named::AUTH_LIST, &value.expect("read"));
    assert_eq!(enc.error(), None);
    let posted = enc.take_replaced_chars();
    drop(enc);
    assert_eq!(
        String::from_utf8(out).expect("ASCII"),
        "Auth-list ::= {\n  names str {\n    \"#\"\n  },\n  affil str \"a\"\n}\n"
    );
    assert_eq!(
        posted,
        vec![crate::report::serial::encoder::ReplacedChar {
            trace: "Auth-list.names.str.E".into(),
            byte: 0xe9,
            position: crate::report::serial::encoder::StreamPosition::Byte(7),
            string: b"\xe9".to_vec().into(),
        }]
    );
    // A UTF8String tag for the name's VisibleString: the read that kept the bytes posts 803.10
    // (once per run, `Utf8TagOnce`), so the re-read before the write posts nothing and the name
    // is written as NCBI's `CDelayBuffer::Update` writes it.
    let mut utf8 = bytes.to_vec();
    utf8[8] = 0x0c;
    let (value, posts) = super::read::read_ber_document(&utf8, named::AUTH_LIST, Names::TEST);
    assert!(value.is_ok(), "{value:?}");
    assert!(
        matches!(posts.as_slice(), [super::read::Post::Utf8Tag(_)]),
        "{posts:?}"
    );
    assert_eq!(write(&utf8), write(bytes));
}

/// The writer writes a value from a list on the heap: 100,000 nested Seq-loc mixes (whose text
/// would hold some 10^10 bytes of indentation) go through an encoder that counts its calls, on a
/// thread of 256 KiB.
#[test]
fn a_nest_of_any_depth_is_written_without_recursion() {
    #[derive(Default)]
    struct Count {
        opened: usize,
        closed: usize,
        leaves: usize,
    }
    impl Encoder for Count {
        fn begin_document(&mut self, _: &str) {}
        fn end_document(&mut self) {}
        fn begin_sequence(&mut self, _: &str) {
            self.opened += 1;
        }
        fn end_sequence(&mut self) {
            self.closed += 1;
        }
        fn member(&mut self, _: &str, _: u8) {}
        fn end_member(&mut self) {}
        fn begin_list(&mut self, _: bool) {
            self.opened += 1;
        }
        fn element(&mut self) {}
        fn end_list(&mut self) {
            self.closed += 1;
        }
        fn variant(&mut self, _: &str, _: &str, _: u8) {}
        fn end_variant(&mut self) {}
        fn int(&mut self, _: i32) {
            self.leaves += 1;
        }
        fn big_int(&mut self, _: i64) {}
        fn enumerated(&mut self, _: &str, _: i32) {}
        fn named_int(&mut self, _: &str, _: i32) {}
        fn boolean(&mut self, _: bool) {}
        fn visible_str(&mut self, _: &[u8]) {}
        fn real(&mut self, _: f64) {}
        fn octets(&mut self, _: &[u8]) {}
        fn null(&mut self) {}
        fn utf8_str(&mut self, _: &[u8]) {}
        fn delayed(&mut self, _: &[u8]) {}
        fn unread_delayed(&mut self) {}
        fn error(&self) -> Option<SerialError> {
            None
        }
    }
    std::thread::Builder::new()
        .stack_size(256 << 10)
        .spawn(|| {
            let depth = 100_000;
            let mut text = String::from("Seq-loc ::= ");
            text.push_str(&"mix { ".repeat(depth));
            text.push_str("int { from 0, to 1, id local id 1 }");
            text.push_str(&" }".repeat(depth));
            let value = read(named::SEQ_LOC, &text);
            let mut count = Count::default();
            write_value(&mut count, named::SEQ_LOC, &value);
            // Each mix is a list; the interval a SEQUENCE with `from` and `to`, the local id one
            // INTEGER more.
            assert_eq!(count.opened, depth + 1);
            assert_eq!(count.closed, depth + 1);
            assert_eq!(count.leaves, 3);
        })
        .expect("thread")
        .join()
        .expect("written");
}

/// The binary stream writes NULL and copies a delay buffer of its own format (the bytes after the
/// member's tag); the JSON stream, whose format no buffer has, writes no delayed member.
#[test]
fn the_binary_stream_copies_its_own_delay_buffers_and_json_writes_none() {
    let mut out = Vec::new();
    let mut enc = crate::report::serial::asn_ber::AsnBinaryEncoder::new(&mut out);
    enc.null();
    enc.delayed(&[0xa0, 0x80, 0x00, 0x00]);
    assert_eq!(enc.error(), None);
    drop(enc);
    assert_eq!(out, [0x05, 0x00, 0xa0, 0x80, 0x00, 0x00]);
    let mut out = Vec::new();
    let mut enc = crate::report::serial::json::JsonEncoder::new(&mut out);
    enc.unread_delayed();
    assert_eq!(
        enc.error(),
        Some(SerialError::Unwritten("a member read into a delay buffer"))
    );
}
