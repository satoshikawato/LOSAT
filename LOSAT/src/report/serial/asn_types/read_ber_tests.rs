//! Unit tests of the binary front end (`read_ber`): what NCBI's `CObjectIStreamAsnBinary` reads
//! of the table's types and where it stops, with its messages (objistrasnb.cpp, cited at each
//! step of `read_ber`), on small values of `Dbtag`, `Object-id`, `Seq-interval`, `Score` and
//! `Seq-loc`.

use super::read::{read_object, Post, Read, ReadError};
use super::read_ber::Ber;
use super::table::named;
use super::write::write_document;
use super::{read_document, Value};
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_text_reader::SliceSource;
use crate::report::serial::encoder::{Encoder, ReplacedChar, StreamPosition};

/// What the binary reader reads of `bytes` as the root type `root`.
fn read(bytes: &[u8], root: u16) -> Read<Value> {
    read_object(&mut Ber::new(SliceSource::new(bytes)), root)
}

/// The ASN.1 text of a value of the root type `root`.
fn text_of(value: &Value, root: u16) -> String {
    let mut out = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut out);
    write_document(&mut enc, root, value);
    assert_eq!(enc.error(), None);
    String::from_utf8(out).expect("ASCII")
}

/// The bytes that NCBI's binary writer gives a text's value (`AsnBinaryEncoder`).
fn ber(root: u16, text: &str) -> Vec<u8> {
    let (value, posts) = read_document(text.as_bytes(), root);
    assert!(posts.is_empty());
    let value = value.expect("a value");
    let mut out = Vec::new();
    let mut enc = AsnBinaryEncoder::new(&mut out);
    write_document(&mut enc, root, &value);
    assert_eq!(enc.error(), None);
    out
}

/// The text of what the binary reader reads of `bytes`, or its error; no message.
fn outcome(bytes: &[u8], root: u16) -> Result<String, ReadError> {
    let (value, posts) = read(bytes, root);
    assert!(posts.is_empty(), "{posts:?}");
    value.map(|value| text_of(&value, root))
}

fn serial(message: &str) -> Result<String, ReadError> {
    Err(ReadError::Serial(message.into()))
}

const DBTAG: &str = "Dbtag ::= {\n  db \"abc\",\n  tag id 5\n}\n";

/// `Dbtag ::= { db "abc", tag id 5 }` with the INTEGER's bytes `int` (`02 01 05` in NCBI's
/// writer).
fn dbtag_with_int(int: &[u8]) -> Vec<u8> {
    [
        &b"\x30\x80\xa0\x80\x1a\x03abc\x00\x00\xa1\x80\xa0\x80"[..],
        int,
        b"\x00\x00\x00\x00\x00\x00",
    ]
    .concat()
}

#[test]
fn a_value_written_by_ncbis_writer_reads_back() {
    let bytes = ber(named::DBTAG, DBTAG);
    assert_eq!(bytes, dbtag_with_int(b"\x02\x01\x05"));
    assert_eq!(outcome(&bytes, named::DBTAG), Ok(DBTAG.into()));
    let interval = "Seq-interval ::= {\n  from 4294967295,\n  to 2147483648,\n  strand other,\n  \
                    id local str \"q\"\n}\n";
    let bytes = ber(named::SEQ_INTERVAL, interval);
    assert_eq!(outcome(&bytes, named::SEQ_INTERVAL), Ok(interval.into()));
}

/// `USE_DEF_LEN`: a definite length is taken wherever an indefinite one may be, its end the end
/// of the value (`HaveMoreElements`, `ExpectEndOfContent`).
#[test]
fn definite_lengths_end_their_values() {
    let definite = b"\x30\x10\xa0\x05\x1a\x03abc\xa1\x07\xa0\x80\x02\x01\x05\x00\x00";
    assert_eq!(definite.len(), 2 + 0x10);
    assert_eq!(outcome(definite, named::DBTAG), Ok(DBTAG.into()));
    // The long form counts from its first length byte (NCBI's build reads the position before
    // the length's bytes: impl D's probes `def_len` and `def_len_plus2`): a length of the content
    // alone ends the class before its last member's end.
    let long = [&b"\x30\x81\x11"[..], &definite[2..]].concat();
    assert_eq!(outcome(&long, named::DBTAG), Ok(DBTAG.into()));
    let content_only = [&b"\x30\x81\x10"[..], &definite[2..]].concat();
    assert_eq!(
        outcome(&content_only, named::DBTAG),
        serial("EndOfContentsByte expected")
    );
    // A member longer than its bytes: the next member starts inside it.
    let longer = b"\x30\x80\xa0\x06\x1a\x03abc\xa1\x80\xa0\x80\x02\x01\x05\x00\x00\x00\x00\x00\x00";
    assert_eq!(
        outcome(longer, named::DBTAG),
        serial("EndOfContentsByte expected")
    );
    // A class shorter than its members: no more member at its limit, then not at the limit.
    let shorter = [&b"\x30\x0f"[..], &definite[2..]].concat();
    assert_eq!(
        outcome(&shorter, named::DBTAG),
        serial("EndOfContentsByte expected")
    );
}

#[test]
fn identifiers_are_checked_as_ncbi_checks_them() {
    let dbtag = dbtag_with_int(b"\x02\x01\x05");
    let with = |at: usize, byte: u8| {
        let mut bytes = dbtag.clone();
        bytes[at] = byte;
        bytes
    };
    assert_eq!(
        outcome(&with(0, 0x31), named::DBTAG),
        serial("Unexpected tag: 17, expected: 16")
    );
    assert_eq!(
        outcome(&with(0, 0x10), named::DBTAG),
        serial("unexpected tag: Sequence (16), should be constructed/None (32)")
    );
    assert_eq!(
        outcome(&with(2, 0x80), named::DBTAG),
        serial(
            "unexpected tag: contextspecific/0 (128), should be \
             contextspecific/constructed/0 (160)"
        )
    );
    // A member that is not one of the class's, or one before the member read last.
    let unexpected = |tag: u32| {
        serial(&format!(
            "unexpected member: [{tag}], should be one of: db[0] tag[1] "
        ))
    };
    assert_eq!(outcome(&with(11, 0xa3), named::DBTAG), unexpected(3));
    assert_eq!(outcome(&with(11, 0xa0), named::DBTAG), unexpected(0));
    assert_eq!(
        outcome(&with(4, 0x04), named::DBTAG),
        serial("unexpected tag: OctetString (4), should be: VisibleString (26)")
    );
    assert_eq!(
        outcome(&with(13, 0xa2), named::DBTAG),
        serial("unexpected member: [2], should be one of: id[0] str[1] ")
    );
    // A long tag number: the bytes after `1f`, 7 bits each, at most 2^24 before the last.
    let long_tag = [&b"\x30\x80\xbf\x88\x80\x80\x80\x00"[..]].concat();
    assert_eq!(
        outcome(&long_tag, named::DBTAG),
        serial("tag number is too big: 16777216")
    );
    let long_tag = b"\x30\x80\xbf\x81\x00\x80\x00\x00";
    assert_eq!(
        outcome(long_tag, named::DBTAG),
        serial("unexpected member: [128], should be one of: db[0] tag[1] ")
    );
    assert_eq!(
        outcome(&with(dbtag.len() - 1, 0x01), named::DBTAG),
        serial("EndOfContentsByte expected")
    );
}

#[test]
fn integers_are_read_as_their_storage() {
    let id = |int: &[u8]| outcome(&dbtag_with_int(int), named::DBTAG);
    let id_is = |value: i64| Ok(DBTAG.replace("id 5", &format!("id {value}")));
    assert_eq!(id(b"\x02\x01\x05"), id_is(5));
    // An application INTEGER (a big integer) is an INTEGER.
    assert_eq!(id(b"\x42\x01\x05"), id_is(5));
    assert_eq!(id(b"\x02\x01\xff"), id_is(-1));
    // More bytes than the storage: the sign's extension only.
    assert_eq!(id(b"\x02\x05\x00\x00\x00\x00\x05"), id_is(5));
    assert_eq!(id(b"\x02\x05\xff\xff\xff\xff\xfe"), id_is(-2));
    assert_eq!(
        id(b"\x02\x05\x00\x80\x00\x00\x00"),
        serial("overflow error")
    );
    assert_eq!(
        id(b"\x02\x05\x01\x00\x00\x00\x05"),
        serial("overflow error")
    );
    assert_eq!(id(b"\x02\x00"), serial("zero length of number"));
    assert_eq!(id(b"\x02\x81\x01\x05"), serial("ShortLength expected"));
    assert_eq!(
        id(b"\x0a\x01\x05"),
        serial("unexpected tag: Enumerated (10), should be: Integer (2)")
    );
    // A Uint4 (`TSeqPos`): fewer bytes than 4 are unsigned, 4 must not set the sign bit.
    let from = |int: &[u8]| {
        let bytes = [
            &b"\x30\x80\xa0\x80"[..],
            int,
            b"\x00\x00\xa1\x80\x02\x01\x01\x00\x00\xa3\x80\xa0\x80\xa0\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00",
        ]
        .concat();
        outcome(&bytes, named::SEQ_INTERVAL)
    };
    let from_is = |value: u64| {
        Ok(format!(
            "Seq-interval ::= {{\n  from {value},\n  to 1,\n  id local id 1\n}}\n"
        ))
    };
    assert_eq!(from(b"\x02\x01\xff"), from_is(255));
    assert_eq!(from(b"\x02\x04\x80\x00\x00\x00"), serial("overflow error"));
    assert_eq!(from(b"\x02\x05\x00\x80\x00\x00\x00"), from_is(2147483648));
    assert_eq!(
        from(b"\x02\x05\x01\x00\x00\x00\x00"),
        serial("overflow error")
    );
}

#[test]
fn enumerated_values_are_named_ones() {
    let strand = |enumerated: &[u8]| {
        let bytes = [
            &b"\x30\x80\xa0\x80\x02\x01\x01\x00\x00\xa1\x80\x02\x01\x01\x00\x00\xa2\x80"[..],
            enumerated,
            b"\x00\x00\xa3\x80\xa0\x80\xa0\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00",
        ]
        .concat();
        outcome(&bytes, named::SEQ_INTERVAL)
    };
    let strand_is = |name: &str| {
        Ok(format!(
            "Seq-interval ::= {{\n  from 1,\n  to 1,\n  strand {name},\n  id local id 1\n}}\n"
        ))
    };
    assert_eq!(strand(b"\x0a\x01\x02"), strand_is("minus"));
    assert_eq!(strand(b"\x0a\x02\x00\xff"), strand_is("other"));
    // `ReadStdSigned`: `ff` alone is -1, no value of Na-strand.
    assert_eq!(
        strand(b"\x0a\x01\xff"),
        serial("invalid value of enumerated type")
    );
    assert_eq!(
        strand(b"\x02\x01\x02"),
        serial("unexpected tag: Integer (2), should be: Enumerated (10)")
    );
}

#[test]
fn reals_are_the_decimal_form_and_the_special_values() {
    let score = |real: &[u8]| {
        let bytes = [
            &b"\x30\x80\xa1\x80\xa0\x80"[..],
            real,
            b"\x00\x00\x00\x00\x00\x00",
        ]
        .concat();
        let (value, posts) = read(&bytes, named::SCORE);
        assert!(posts.is_empty());
        value.map(|value| match &value {
            Value::Sequence(members) => match &members[0].1 {
                Value::Choice(0, real) => match **real {
                    Value::Real(real) => real.to_bits(),
                    _ => panic!("a REAL"),
                },
                other => panic!("{other:?}"),
            },
            other => panic!("{other:?}"),
        })
    };
    assert_eq!(score(b"\x09\x00"), Ok(0.0f64.to_bits()));
    assert_eq!(score(b"\x09\x01\x40"), Ok(f64::INFINITY.to_bits()));
    assert_eq!(score(b"\x09\x01\x41"), Ok(f64::NEG_INFINITY.to_bits()));
    assert_eq!(
        score(b"\x09\x01\x42"),
        Ok(crate::ncbi_str::NCBI_NAN.to_bits())
    );
    assert_eq!(score(b"\x09\x01\x43"), Ok((-0.0f64).to_bits()));
    assert_eq!(
        score(b"\x09\x01\x44"),
        Err(ReadError::Serial("Unrecognized REAL data".into()))
    );
    assert_eq!(score(b"\x09\x04\x001.5"), Ok(1.5f64.to_bits()));
    // The C string ends at its first 0 byte.
    assert_eq!(score(b"\x09\x04\x001\x005"), Ok(1.0f64.to_bits()));
    assert_eq!(
        score(b"\x09\x04\x001x5"),
        Err(ReadError::Serial("bad REAL data string".into()))
    );
    assert_eq!(
        score(b"\x09\x03\x40\x01\x02"),
        Err(ReadError::Serial(
            "Unsupported encoding of REAL data: encoding = 64".into()
        ))
    );
    let long = [&b"\x09\x82\x01\x01\x00"[..], &[b'1'; 256]].concat();
    assert_eq!(
        score(&long),
        Err(ReadError::Serial("too long REAL data: length > 256".into()))
    );
}

#[test]
fn null_and_lengths() {
    let null = |bytes: &[u8]| {
        outcome(
            &[&b"\xa0\x80"[..], bytes, b"\x00\x00"].concat(),
            named::SEQ_LOC,
        )
    };
    assert_eq!(null(b"\x05\x00"), Ok("Seq-loc ::= null NULL\n".into()));
    assert_eq!(
        null(b"\x05\x01\x00"),
        serial("unexpected length: 1, should be: 0")
    );
    let db = |string: &[u8]| {
        let bytes = [
            &b"\x30\x80\xa0\x80"[..],
            string,
            b"\x00\x00\xa1\x80\xa0\x80\x02\x01\x05\x00\x00\x00\x00\x00\x00",
        ]
        .concat();
        outcome(&bytes, named::DBTAG)
    };
    assert_eq!(db(b"\x1a\x81\x03abc"), Ok(DBTAG.into()));
    assert_eq!(db(b"\x1a\x80abc"), serial("unexpected indefinite length"));
    assert_eq!(db(b"\x1a\x81\x00abc"), serial("illegal length start"));
    assert_eq!(db(b"\x1a\x89\x01"), serial("length overflow"));
    // A string longer than the data: the end of the data in the member's frame.
    assert!(matches!(
        db(b"\x1a\x84\x7f\xff\xff\xffabc"),
        Err(ReadError::Serial(message)) if message.starts_with("unexpected end of the data:")
    ));
}

/// `CEofException` leaves the reader from the root's frames; deeper, `HandleEOF` makes it a
/// `CSerialException` (P3.104, P7.001-022).
#[test]
fn the_end_of_the_data_is_an_error_below_the_root() {
    let dbtag = dbtag_with_int(b"\x02\x01\x05");
    let cut = |count: usize| read(&dbtag[..dbtag.len() - count], named::DBTAG).0;
    assert_eq!(cut(1), Err(ReadError::Eof));
    for count in 2..dbtag.len() - 1 {
        assert!(
            matches!(cut(count), Err(ReadError::Serial(_))),
            "{count}: {:?}",
            cut(count)
        );
    }
    // `30` alone: the length's byte is missing in the root class.
    assert_eq!(cut(dbtag.len() - 1), Err(ReadError::Eof));
    assert_eq!(read(b"", named::DBTAG).0, Err(ReadError::Eof));
}

#[test]
fn visible_strings_are_fixed_and_utf8_tags_posted_once() {
    // `0c`, a UTF8String's tag, for a VisibleString: posted once in the stream.
    let bytes = [
        dbtag_with_int(b"\x02\x01\x05").replace_first(b"\x1a", b"\x0c"),
        dbtag_with_int(b"\x02\x01\x05").replace_first(b"\x1a", b"\x0c"),
    ]
    .concat();
    let mut ber = Ber::new(SliceSource::new(&bytes));
    let (first, posts) = read_object(&mut ber, named::DBTAG);
    assert_eq!(text_of(&first.expect("a value"), named::DBTAG), DBTAG);
    assert_eq!(posts, vec![Post::Utf8ForVisible("Dbtag.db".into())]);
    // blast_formatter has no diagnostic prefix (impl D's probe `utf8_for_visible`).
    assert_eq!(
        String::from_utf8(posts[0].line(None)).expect("ASCII"),
        "Warning: (803.10) CObjectIStreamAsnBinary: UTF8String data for \
         VisibleString member Dbtag.db, ASN.1 specification may need an update\n"
    );
    let (second, posts) = read_object(&mut ber, named::DBTAG);
    assert!(second.is_ok());
    assert_eq!(posts, vec![]);
    // A byte outside 0x20..0x7E: `#`, posted with the string as fixed before it and the
    // position after the string.
    let bytes = dbtag_with_int(b"\x02\x01\x05").replace_first(b"abc", b"\x01b\x7f");
    let (value, posts) = read(&bytes, named::DBTAG);
    assert_eq!(
        text_of(&value.expect("a value"), named::DBTAG),
        DBTAG.replace("abc", "#b#")
    );
    let replaced = |byte: u8, string: &[u8]| {
        Post::BadChar(ReplacedChar {
            trace: "Dbtag.db".into(),
            byte,
            position: StreamPosition::Byte(9),
            string: string.into(),
        })
    };
    assert_eq!(
        posts,
        vec![replaced(0x01, b"\x01b\x7f"), replaced(0x7f, b"#b\x7f")]
    );
    assert_eq!(
        posts[1].line(None),
        b"Error: (803.7) Dbtag.db\nBad char [0x7F] in string at byte 9\n#b\x7f\n"
    );
}

/// `Vec::replace_first` for the tests.
trait ReplaceFirst {
    fn replace_first(self, from: &[u8], to: &[u8]) -> Vec<u8>;
}

impl ReplaceFirst for Vec<u8> {
    fn replace_first(self, from: &[u8], to: &[u8]) -> Vec<u8> {
        let at = self
            .windows(from.len())
            .position(|window| window == from)
            .expect("the bytes to replace");
        [&self[..at], to, &self[at + from.len()..]].concat()
    }
}
