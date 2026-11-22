use super::{OutBuffer, Piece, Stream};

/// What NCBI had written at each post of `bytes` (one `piece`) after `before` bytes of the
/// same kind of call.
fn written_at(
    stream: Stream,
    before: &[(Piece, &[u8])],
    piece: Piece,
    bytes: &[u8],
    posts: &[usize],
) -> Vec<u64> {
    let mut buffer = OutBuffer::new(stream);
    for (piece, bytes) in before {
        buffer.put(*piece, bytes, &[], &mut Vec::new());
    }
    let mut written = Vec::new();
    buffer.put(piece, bytes, posts, &mut written);
    written
}

#[test]
fn a_full_buffer_keeps_its_back_limit() {
    // 4096 bytes of `PutChar` fill the buffer; the next one flushes it first: ASN.1 text
    // keeps 80 bytes, XML 1, binary none. A post before that byte comes before the flush.
    let chars = vec![b'a'; 4098];
    assert_eq!(
        written_at(Stream::AsnText, &[], Piece::Str, &chars, &[4096, 4097]),
        [0, 4016]
    );
    for (stream, written) in [(Stream::Xml, 4095), (Stream::AsnBinary, 4096)] {
        let mut buffer = OutBuffer::new(stream);
        buffer.put(Piece::Octets, &chars[..4097], &[], &mut Vec::new());
        assert_eq!(buffer.written(), written);
    }
    // The next flush comes when the buffer is full again.
    let chars = vec![b'a'; 8200];
    assert_eq!(
        written_at(Stream::AsnText, &[], Piece::Str, &chars, &[8112, 8113]),
        [4016, 8032]
    );
}

#[test]
fn a_call_that_does_not_fit_flushes_before_it() {
    // 4090 bytes, then a name of 10 bytes: the buffer is flushed before the name (4010 of
    // 4090 bytes written, 80 kept), not at 4096.
    let mut buffer = OutBuffer::new(Stream::AsnText);
    buffer.put(Piece::Str, &[b'a'; 4090], &[], &mut Vec::new());
    buffer.put(Piece::Number, b"1234567890", &[], &mut Vec::new());
    assert_eq!(buffer.written(), 4010);
    assert_eq!(buffer.used, 90);
}

#[test]
fn a_long_string_fills_the_buffer_up() {
    // ASN.1 binary: a run of 1100 bytes of a string after 3004 bytes goes through `Write`,
    // which fills the buffer (4096 written) without growing it; a run of 1000 bytes that
    // does not fit is one `Skip`, which flushes first.
    let mut string = vec![0x1a, 0x82, 0x05, 0xdc];
    string.extend(std::iter::repeat_n(b'L', 1500));
    string[4 + 1100] = b'#';
    let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[0u8; 3000])];
    assert_eq!(
        written_at(Stream::AsnBinary, before, Piece::Str, &string, &[4 + 1100]),
        [4096]
    );
    string[4 + 1100] = b'L';
    string[4 + 1000] = b'#';
    let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[0u8; 3500])];
    assert_eq!(
        written_at(Stream::AsnBinary, before, Piece::Str, &string, &[4 + 1000]),
        [3504]
    );
}

/// The buffer after a number of `before` digits (one `Skip` in ASN.1 text and XML) and the
/// calls of `piece`.
fn after(stream: Stream, before: usize, piece: Piece, bytes: &[u8]) -> OutBuffer {
    let mut buffer = OutBuffer::new(stream);
    buffer.put(Piece::Number, &vec![b'1'; before], &[], &mut Vec::new());
    buffer.put(piece, bytes, &[], &mut Vec::new());
    buffer
}

#[test]
fn asn_text_header_enum_and_real_are_split_as_ncbi_calls_them() {
    // `WriteFileHeader`: `WriteId` of the name (14 bytes fit at 4080) and `PutString` of `" ::= "`,
    // which flushes 4094 - 80 bytes; one call of 19 bytes would flush at 4080.
    let buffer = after(Stream::AsnText, 4080, Piece::Header, b"Blast4-archive ::= ");
    assert_eq!((buffer.written(), buffer.used), (4014, 85));
    // `WriteEnum`: `PutChar` of the first letter (the 4096th byte), then `PutString` of the
    // rest, which flushes.
    let buffer = after(Stream::AsnText, 4095, Piece::Enum, b"plus");
    assert_eq!((buffer.written(), buffer.used), (4016, 83));
    // `WriteDouble2`: `{ ` fits at 4090, the 9 digits flush 4092 - 80 bytes; then `, 10, `,
    // the exponent and ` }`.
    let buffer = after(Stream::AsnText, 4090, Piece::Real, b"{ 123456789, 10, -9 }");
    assert_eq!((buffer.written(), buffer.used), (4012, 99));
    // A negative mantissa puts its sign alone.
    let buffer = after(Stream::AsnText, 4093, Piece::Real, b"{ -12, 10, 3 }");
    assert_eq!((buffer.written(), buffer.used), (4016, 80 + 2 + 6 + 1 + 2));
    // Zero is one `PutString`, which flushes before it.
    let buffer = after(Stream::AsnText, 4090, Piece::Real, b"{ 0, 10, 0 }");
    assert_eq!((buffer.written(), buffer.used), (4010, 92));
}

#[test]
fn asn_text_octets_wrap_with_one_byte_calls() {
    // `BeginBytes` (`'`), two `PutChar`s per byte with the end of line of `WrapAt` at 78 with
    // `false` (`PutEol` with `false`: `Reserve` of 1, no indentation) between them, then
    // `PutString` of `"'H"`: the 97th of the one-byte calls flushes 4096 - 80 bytes.
    let mut octets = b"'".to_vec();
    octets.extend(b"AB".repeat(39));
    octets.push(b'\n');
    octets.extend(b"CD".repeat(10));
    octets.extend(b"'H");
    let buffer = after(Stream::AsnText, 4000, Piece::Octets, &octets);
    assert_eq!((buffer.written(), buffer.used), (4016, 80 + 4 + 2));
}

#[test]
fn a_text_post_after_a_wrap_of_its_string_counts_the_wrap() {
    // The string's wrap (`PutEolAtWordEnd`) inserted one end of line, `Reserve(1)`, before
    // the replaced byte: the post after 112 one-byte calls follows the flush at the 107th.
    let mut string = b"\"".to_vec();
    string.extend([b'a'; 90]);
    string.push(b'\n');
    string.extend([b'b'; 20]);
    string.extend(b"#\"");
    let before: &[(Piece, &[u8])] = &[(Piece::Markup, &[b'x'; 3990])];
    assert_eq!(
        written_at(Stream::AsnText, before, Piece::Str, &string, &[112]),
        [4016]
    );
    // Without the flush: the same string at the start of the buffer.
    assert_eq!(
        written_at(Stream::AsnText, &[], Piece::Str, &string, &[112]),
        [0]
    );
}

#[test]
fn an_xml_character_reference_is_a_put_string_and_one_byte_calls() {
    // `&#x9;`: `PutString` of `"&#x"` fits at 4092, `PutChar` of the digit 9 ends the buffer and
    // `PutChar` of a semicolon flushes 4095 bytes (XML keeps 1); one call of 5 bytes would flush at
    // 4092.
    let buffer = after(Stream::Xml, 4092, Piece::Str, b"&#x9;");
    assert_eq!((buffer.written(), buffer.used), (4095, 2));
    // An entity of a markup character is one `PutString`.
    let buffer = after(Stream::Xml, 4093, Piece::Str, b"&amp;");
    assert_eq!((buffer.written(), buffer.used), (4092, 6));
}

#[test]
fn a_call_larger_than_the_buffer_doubles_it() {
    // `DoReserve`: a `Skip` of 5000 bytes flushes and doubles the buffer to 8192; the next
    // call that does not fit flushes again and keeps that size.
    let mut buffer = after(Stream::AsnText, 100, Piece::Markup, &[b'n'; 5000]);
    assert_eq!(
        (buffer.written(), buffer.used, buffer.capacity),
        (20, 5080, 8192)
    );
    buffer.put(Piece::Markup, &[b'n'; 4000], &[], &mut Vec::new());
    assert_eq!(
        (buffer.written(), buffer.used, buffer.capacity),
        (5020, 4080, 8192)
    );
    // Two doublings for a call of more than twice the buffer.
    let buffer = after(Stream::AsnText, 0, Piece::Markup, &[b'n'; 9000]);
    assert_eq!(
        (buffer.written(), buffer.used, buffer.capacity),
        (0, 9000, 16384)
    );
}

#[test]
fn an_xml_value_attribute_takes_the_tag_end_back() {
    // `<T_b>` ends the buffer; `value` replaces its `>`: one byte of the buffer is given
    // back before the attribute.
    let mut buffer = OutBuffer::new(Stream::Xml);
    buffer.put(Piece::Markup, b"\n<T_b>", &[], &mut Vec::new());
    buffer.put(Piece::Bool, b" value=\"true\"/>", &[], &mut Vec::new());
    assert_eq!(buffer.used, 6 - 1 + 15);
    // An element whose tag opens in the same call.
    let mut buffer = OutBuffer::new(Stream::Xml);
    buffer.put(
        Piece::Enum,
        b"\n  <T_E value=\"plus\"/>",
        &[],
        &mut Vec::new(),
    );
    assert_eq!(buffer.used, 22);
}
