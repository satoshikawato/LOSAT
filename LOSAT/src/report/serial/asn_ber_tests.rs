//! The ASN.1 binary encoder (outfmt 9) primitives against the bytes that NCBI BLAST+ 2.17.0
//! writes (outfmt-20261019 S/DESIGN_OF3.md 4.3): the values below are taken from its samples
//! (`tests/fixtures/seqalign/*_fmt9.out`) or follow `objostrasnb.cpp` for a value the samples
//! do not hold (an integer of three or four bytes, a string of 128 bytes or more).

use super::asn_ber::AsnBinaryEncoder;
use super::encoder::{Encoder, SerialError};

/// The bytes that `write` makes the encoder append.
fn bytes(write: impl FnOnce(&mut AsnBinaryEncoder<'_>)) -> Vec<u8> {
    let mut buf = Vec::new();
    let mut enc = AsnBinaryEncoder::new(&mut buf);
    write(&mut enc);
    assert_eq!(enc.error(), None);
    buf
}

/// `WriteNumberValue`: the shortest of one to four two's-complement bytes (230 does not fit
/// the signed byte, so `02 02 00 e6`; the score 329 is `02 02 01 49`).
#[test]
fn integers_use_the_shortest_signed_form() {
    for (value, expected) in [
        (0, &[0x02, 0x01, 0x00][..]),
        (1, &[0x02, 0x01, 0x01]),
        (127, &[0x02, 0x01, 0x7f]),
        (128, &[0x02, 0x02, 0x00, 0x80]),
        (230, &[0x02, 0x02, 0x00, 0xe6]),
        (329, &[0x02, 0x02, 0x01, 0x49]),
        (-1, &[0x02, 0x01, 0xff]),
        (-128, &[0x02, 0x01, 0x80]),
        (-129, &[0x02, 0x02, 0xff, 0x7f]),
        (32767, &[0x02, 0x02, 0x7f, 0xff]),
        (32768, &[0x02, 0x03, 0x00, 0x80, 0x00]),
        (-32769, &[0x02, 0x03, 0xff, 0x7f, 0xff]),
        (0x7f_ffff, &[0x02, 0x03, 0x7f, 0xff, 0xff]),
        (0x80_0000, &[0x02, 0x04, 0x00, 0x80, 0x00, 0x00]),
        (-0x80_0000, &[0x02, 0x03, 0x80, 0x00, 0x00]),
        (-0x80_0001, &[0x02, 0x04, 0xff, 0x7f, 0xff, 0xff]),
        (i32::MAX, &[0x02, 0x04, 0x7f, 0xff, 0xff, 0xff]),
        (i32::MIN, &[0x02, 0x04, 0x80, 0x00, 0x00, 0x00]),
    ] {
        assert_eq!(bytes(|enc| enc.int(value)), expected, "{value}");
    }
}

/// `WriteEnum`: tag 0x0a and the same number form (`diags` = 2, `partial` = 3, `minus` = 2).
#[test]
fn enumerated_values_use_tag_0a() {
    assert_eq!(
        bytes(|enc| enc.enumerated("partial", 3)),
        [0x0a, 0x01, 0x03]
    );
    assert_eq!(
        bytes(|enc| enc.enumerated("unknown", 0)),
        [0x0a, 0x01, 0x00]
    );
    assert_eq!(
        bytes(|enc| enc.enumerated("other", 255)),
        [0x0a, 0x02, 0x00, 0xff]
    );
}

/// `WriteBool`: `01 01 01` and `01 01 00`.
#[test]
fn booleans_are_one_byte() {
    assert_eq!(bytes(|enc| enc.boolean(true)), [0x01, 0x01, 0x01]);
    assert_eq!(bytes(|enc| enc.boolean(false)), [0x01, 0x01, 0x00]);
}

/// `WriteString` of a VisibleString: the short length up to 127 bytes, the long form from 128
/// (`0x80 + count`, then the count bytes big-endian).
#[test]
fn visible_strings_have_short_and_long_lengths() {
    assert_eq!(
        bytes(|enc| enc.visible_str(b"score")),
        [&[0x1a, 0x05][..], b"score"].concat()
    );
    assert_eq!(bytes(|enc| enc.visible_str(b"")), [0x1a, 0x00]);
    let string = vec![b'x'; 127];
    assert_eq!(
        bytes(|enc| enc.visible_str(&string)),
        [&[0x1a, 0x7f][..], &string].concat()
    );
    let string = vec![b'y'; 128];
    assert_eq!(
        bytes(|enc| enc.visible_str(&string)),
        [&[0x1a, 0x81, 0x80][..], &string].concat()
    );
    let string = vec![b'z'; 256];
    assert_eq!(
        bytes(|enc| enc.visible_str(&string)),
        [&[0x1a, 0x82, 0x01, 0x00][..], &string].concat()
    );
    let string = vec![b'z'; 70_000];
    assert_eq!(
        bytes(|enc| enc.visible_str(&string)),
        [&[0x1a, 0x83, 0x01, 0x11, 0x70][..], &string].concat()
    );
    // Above 0xffffff NCBI writes all `sizeof(size_t)` = 8 bytes of the length.
    let string = vec![b'w'; 0x100_0000];
    assert_eq!(
        bytes(|enc| enc.visible_str(&string))[..10],
        [0x1a, 0x88, 0, 0, 0, 0, 0x01, 0, 0, 0]
    );
}

/// A byte outside 0x20..0x7E is an error (the callers refuse such titles before the search),
/// as in the ASN.1 text encoder.
#[test]
fn a_non_visible_byte_is_an_error() {
    let mut buf = Vec::new();
    let mut enc = AsnBinaryEncoder::new(&mut buf);
    enc.visible_str(b"a\tb");
    assert_eq!(enc.error(), Some(SerialError::NonVisibleChar(b'\t')));
}

/// `WriteDouble2(data, DBL_DIG)`: `09`, the length (digits + 1), the decimal form byte 0 and
/// `DoubleToStringPosix` (zero is `09 02 00 30`; the sample's e-value is `09 16 00
/// "2.60916226603389e-178"`; a whole number has no decimal point).
#[test]
fn reals_are_a_decimal_form_byte_and_posix_digits() {
    assert_eq!(bytes(|enc| enc.real(0.0)), [0x09, 0x02, 0x00, b'0']);
    assert_eq!(
        bytes(|enc| enc.real(2.60916226603389e-178)),
        [&[0x09, 0x16, 0x00][..], b"2.60916226603389e-178"].concat()
    );
    assert_eq!(
        bytes(|enc| enc.real(329.0)),
        [0x09, 0x04, 0x00, b'3', b'2', b'9']
    );
    assert_eq!(
        bytes(|enc| enc.real(0.5)),
        [&[0x09, 0x04, 0x00][..], b"0.5"].concat()
    );
}

/// A NaN or an infinity is an error, not a byte string (no BLAST report has one).
#[test]
fn a_non_finite_real_is_an_error() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut buf = Vec::new();
        let mut enc = AsnBinaryEncoder::new(&mut buf);
        enc.real(value);
        assert_eq!(enc.error(), Some(SerialError::NonFiniteReal));
    }
}

/// The structure calls: a SEQUENCE is `30 80 ... 00 00`, a SET OF `31 80`, and every member and
/// variant the explicit tag `a0|i 80 ... 00 00`, here the Score of the sample
/// `blastn-megablast_subject_fmt9.out` (`{ id str "score", value int 329 }`) in a Seq-align's
/// score SET OF (member 2).
#[test]
fn constructed_values_use_indefinite_lengths_and_explicit_tags() {
    let score = bytes(|enc| {
        enc.begin_document("Seq-align");
        enc.begin_sequence("Seq-align");
        enc.member("score", 2);
        enc.begin_list(true);
        enc.element();
        enc.begin_sequence("Score");
        enc.member("id", 0);
        enc.variant("str", 1);
        enc.visible_str(b"score");
        enc.end_variant();
        enc.end_member();
        enc.member("value", 1);
        enc.variant("int", 1);
        enc.int(329);
        enc.end_variant();
        enc.end_member();
        enc.end_sequence();
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
        enc.end_document();
    });
    let expected: Vec<u8> = [
        &[0x30, 0x80, 0xa2, 0x80, 0x31, 0x80][..],
        &[0x30, 0x80, 0xa0, 0x80, 0xa1, 0x80, 0x1a, 0x05],
        b"score",
        &[
            0x00, 0x00, 0x00, 0x00, 0xa1, 0x80, 0xa1, 0x80, 0x02, 0x02, 0x01, 0x49,
        ],
        &[
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ],
    ]
    .concat();
    assert_eq!(score, expected);
}

/// Documents follow each other with nothing between them.
#[test]
fn documents_are_concatenated() {
    let mut buf = Vec::new();
    for _ in 0..2 {
        let mut enc = AsnBinaryEncoder::new(&mut buf);
        enc.begin_document("Seq-annot");
        enc.begin_sequence("Seq-annot");
        enc.end_sequence();
        enc.end_document();
    }
    assert_eq!(buf, [0x30, 0x80, 0x00, 0x00, 0x30, 0x80, 0x00, 0x00]);
}
