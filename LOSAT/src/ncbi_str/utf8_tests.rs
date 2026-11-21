use super::*;

#[test]
fn the_lenient_decoder_takes_overlong_forms_and_throws_by_ncbi_rule() {
    // Overlong forms of ASCII characters are their characters (no overlong check).
    assert_eq!(decode_to_ascii(b"\xc1\x81"), Ok((b'A', 2)));
    assert_eq!(decode_to_ascii(b"\xc0\x80z"), Ok((0, 2)));
    assert_eq!(decode_to_ascii(b"\xe0\x81\x81b"), Ok((b'A', 3)));
    assert_eq!(decode_to_ascii(b"\xe0\x80\xbc"), Ok((b'<', 3)));
    assert_eq!(decode_to_ascii(b"\xe0\x80\x81"), Ok((0x01, 3)));
    assert_eq!(decode_to_ascii(b"\xf0\x80\x81\x81"), Ok((b'A', 4)));
    // A complete sequence above 0x7F: SymbolToChar (also surrogates and > U+10FFFF).
    for bytes in [
        &b"\xc3\xa9"[..],
        b"\xe2\x82\xac",
        b"\xf0\x9f\x98\x80",
        b"\xed\xa0\x80",
        b"\xf4\x90\x80\x80",
        b"\xf5\x80\x80\x80",
    ] {
        assert_eq!(
            decode_to_ascii(bytes),
            Err(Utf8Error::SymbolToChar),
            "{bytes:x?}"
        );
    }
    // No lead byte: DecodeFirst.
    for bytes in [&b"\x80"[..], b"\xa9s", b"\xbf", b"\xf8\x80", b"\xff"] {
        assert_eq!(
            decode_to_ascii(bytes),
            Err(Utf8Error::DecodeFirst),
            "{bytes:x?}"
        );
    }
    // A missing continuation byte, the end of the string (NUL) included: DecodeNext.
    for bytes in [
        &b"\xc3"[..],
        b"\xe9 musculus",
        b"\xf0\x9f\x98",
        b"\xc3\x00\xa9",
    ] {
        assert_eq!(
            decode_to_ascii(bytes),
            Err(Utf8Error::DecodeNext),
            "{bytes:x?}"
        );
    }
}

#[test]
fn the_exceptions_report_ncbis_lines() {
    // NCBI's reports in LOSAT's words (`crate::wording::losat_text`).
    let ncbi = |text: &str| crate::wording::losat_text(text).into_owned();
    let error = Utf8Error::SymbolToChar.exception();
    assert_eq!(error.site().line, 6901);
    assert_eq!(String::from_utf8(error.what()).unwrap(), ncbi("NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line 6901: Error: (CStringException::eConvert) ncbi::CUtf8::SymbolToChar() - Failed to convert symbol to requested encoding (m_Pos = 0)\n"));
    let error = Utf8Error::DecodeFirst.exception();
    assert_eq!(error.site().line, 7190);
    assert_eq!(String::from_utf8(error.what()).unwrap(), ncbi("NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line 7190: Error: (CStringException::eBadArgs) ncbi::CUtf8::DecodeFirst() - Source string is not in UTF8 format (m_Pos = 0)\n"));
    let error = Utf8Error::DecodeNext.exception();
    assert_eq!(error.site().line, 7202);
    assert_eq!(String::from_utf8(error.what()).unwrap(), ncbi("NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line 7202: Error: (CStringException::eBadArgs) ncbi::CUtf8::DecodeNext() - Source string is not in UTF8 format (m_Pos = 0)\n"));
}
