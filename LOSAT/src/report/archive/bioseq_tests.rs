use super::*;

/// The `DefLine` field is the raw line through `NStr::PrintableString` (fasta.cpp:601) in the
/// default `fNewLine_Quote` mode: a line feed, which only a record the library builds can
/// hold, is written as `\n`.
#[test]
fn the_defline_field_quotes_a_line_feed() {
    let mut buf = Vec::new();
    let mut enc = crate::report::serial::asn_text::AsnTextEncoder::new(&mut buf);
    let bioseq = BioseqInput {
        id: "Query_1",
        ids: None,
        seq_loc: None,
        defline: std::borrow::Cow::Borrowed(b">a\nb"),
        title: b"a",
        letters: b"MKV",
        protein: true,
    };
    write_bioseq(&mut enc, &bioseq);
    assert_eq!(enc.error(), None);
    let text = String::from_utf8(buf).unwrap();
    assert!(text.contains(r#"str ">a\nb""#), "{text}");
}

#[test]
fn nucleotides_pack_into_two_or_four_bits() {
    assert_eq!(pack_nucleotides(b"ACGT"), (true, vec![0b0001_1011]));
    assert_eq!(pack_nucleotides(b"acgta"), (true, vec![0b0001_1011, 0]));
    // N makes the whole sequence ncbi4na: A=1 C=2 N=15 T=8.
    assert_eq!(pack_nucleotides(b"ACNT"), (false, vec![0x12, 0xf8]));
    assert_eq!(pack_nucleotides(b"RYK"), (false, vec![0x5a, 0xc0]));
}
