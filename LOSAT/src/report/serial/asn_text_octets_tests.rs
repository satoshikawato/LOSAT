//! The OCTET STRING of the ASN.1 text stream (`Seq-data` `ncbi2na` / `ncbi4na` of the BLAST
//! archive) at every column where it can start (NCBI reference (598d8ae6):
//! c++/src/serial/objostrasn.cpp:868-895, c++/include/util/strbuffer.inl:382-393): `'`, then
//! before each byte an end of line (no indentation) when the line holds 78 bytes or more, then
//! two upper-case hex digits; before `'H` the same test.

use super::asn_text::AsnTextEncoder;
use super::encoder::Encoder;

/// An independent model of the rule, byte by byte.
fn model(start_column: usize, data: &[u8]) -> String {
    let mut out = String::new();
    let mut column = start_column;
    out.push('\'');
    column += 1;
    for byte in data {
        if column >= 78 {
            out.push('\n');
            column = 0;
        }
        out.push_str(&format!("{byte:02X}"));
        column += 2;
    }
    if column >= 78 {
        out.push('\n');
    }
    out.push_str("'H");
    out
}

/// The encoder's output after a line of `start_column` bytes (an enumerated name is written
/// as it is).
fn octets_after(start_column: usize, data: &[u8]) -> String {
    let mut buf = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut buf);
    let name = "x".repeat(start_column);
    enc.enumerated(&name, 0);
    enc.octets(data);
    assert_eq!(enc.error(), None);
    String::from_utf8(buf).unwrap()[start_column..].to_string()
}

#[test]
fn octets_wrap_at_every_start_column() {
    let data: Vec<u8> = (0..70u32).map(|i| (i * 37 + 11) as u8).collect();
    for start_column in 0..=100 {
        for length in [0, 1, 2, 3, 38, 39, 40, 41, 70] {
            assert_eq!(
                octets_after(start_column, &data[..length]),
                model(start_column, &data[..length]),
                "column {start_column} length {length}"
            );
        }
    }
}

#[test]
fn a_hex_line_without_prefix_holds_39_bytes_then_wraps() {
    let data = [0x01u8; 80];
    let text = octets_after(0, &data);
    let lines: Vec<&str> = text.split('\n').collect();
    // `'` + 38 pairs = 77 bytes, a 39th pair makes 79 (>= 78): the next pair starts a line.
    assert_eq!(lines[0].len(), 1 + 2 * 39);
    assert_eq!(lines[1].len(), 2 * 39);
    assert_eq!(lines[2], "0101'H");
}

/// `flush_to` hands over the lines that are final and keeps the last end of line, so the
/// string wrap behaves as it does with the whole document in the buffer.
#[test]
fn flushing_does_not_change_the_bytes() {
    let words: Vec<String> = (0..60).map(|i| format!("word{i}")).collect();
    let text = words.join(" ");
    let mut whole = Vec::new();
    {
        let mut enc = AsnTextEncoder::new(&mut whole);
        enc.begin_document("D");
        for _ in 0..5 {
            enc.visible_str(text.as_bytes());
            enc.octets(&[
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ]);
        }
    }
    let mut pieces = Vec::new();
    let mut buf = Vec::new();
    {
        let mut enc = AsnTextEncoder::new(&mut buf);
        enc.begin_document("D");
        for _ in 0..5 {
            enc.visible_str(text.as_bytes());
            enc.flush_to(&mut pieces, false).unwrap();
            enc.octets(&[
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ]);
            enc.flush_to(&mut pieces, false).unwrap();
        }
        enc.flush_to(&mut pieces, true).unwrap();
    }
    assert_eq!(pieces, whole);
    assert!(buf.is_empty());
}
