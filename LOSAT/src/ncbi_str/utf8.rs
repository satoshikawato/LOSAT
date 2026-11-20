//! `CUtf8::Decode` and `CUtf8::SymbolToChar` to US-ASCII: what NCBI's XML object stream does with
//! a byte >= 0x80 of a string when it writes the US-ASCII documents of outfmt 5, 14 and 16
//! (`CObjectOStreamXml::WriteEncodedChar`, `report::serial::xml`), with the three exceptions it
//! throws there ([`Utf8Error`]).
//!
//! NCBI's decoder is lenient: the lead byte alone says how many continuation bytes follow
//! (C0-DF one, E0-EF two, F0-F7 three; C0, C1 and F5-F7 included), each continuation byte only
//! has to be 10xxxxxx, and nothing checks for overlong forms, surrogates or code points past
//! U+10FFFF. So `E0 81 81` is `A`, `C0 80` is 0 and `ED A0 80` is U+D800.

use crate::ncbi_exception::{NcbiException, ThrowSite};

/// The exception that NCBI's decode of a string to US-ASCII throws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Utf8Error {
    /// `CUtf8::DecodeFirst`: a byte where a sequence starts that is no lead byte (80-BF, F8-FF).
    DecodeFirst,
    /// `CUtf8::DecodeNext`: a byte that is no continuation byte where one is due (the string's
    /// terminating NUL included).
    DecodeNext,
    /// `CUtf8::SymbolToChar`: a complete sequence of a code point above 0x7F, which US-ASCII
    /// cannot hold.
    SymbolToChar,
}

impl Utf8Error {
    /// Where NCBI throws the exception (`NCBI_THROW2(CStringException, code, message, 0)`).
    pub(crate) fn site(self) -> ThrowSite {
        let (line, kind, function) = match self {
            Utf8Error::DecodeFirst => (
                7190,
                "CStringException::eBadArgs",
                "ncbi::CUtf8::DecodeFirst",
            ),
            Utf8Error::DecodeNext => (
                7202,
                "CStringException::eBadArgs",
                "ncbi::CUtf8::DecodeNext",
            ),
            Utf8Error::SymbolToChar => (
                6901,
                "CStringException::eConvert",
                "ncbi::CUtf8::SymbolToChar",
            ),
        };
        ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line,
            kind,
            function: Some(function),
        }
    }

    /// The exception, at position 0 (`m_Pos = 0`).
    pub(crate) fn exception(self) -> NcbiException {
        let message = match self {
            Utf8Error::DecodeFirst | Utf8Error::DecodeNext => "Source string is not in UTF8 format",
            Utf8Error::SymbolToChar => "Failed to convert symbol to requested encoding",
        };
        NcbiException::with_pos(self.site(), message, 0)
    }
}

/// `CUtf8::Decode(src)` followed by `CUtf8::SymbolToChar(sym, eEncoding_Ascii)`: the character
/// that the sequence at the start of `src` (a byte >= 0x80 first) stands for, and how many bytes
/// it takes. The end of `src` reads as the string's terminating NUL.
///
/// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:5665-5674
/// ```c++
/// template <typename TIterator> inline TUnicodeSymbol
/// CUtf8::x_Decode(TIterator& src)
/// {
///     SIZE_TYPE more=0;
///     TUnicodeSymbol sym = DecodeFirst(*src,more);
///     while (more--) {
///         sym = DecodeNext(sym, *(++src));
///     }
///     return sym;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7174-7206
/// ```c++
/// TUnicodeSymbol CUtf8::DecodeFirst(char ch, SIZE_TYPE& more)
/// {
///     TUnicodeSymbol chRes = 0;
///     more = 0;
///     if ((ch & 0x80) == 0) {
///         chRes = ch;
///     } else if ((ch & 0xE0) == 0xC0) {
///         chRes = (ch & 0x1F);
///         more = 1;
///     } else if ((ch & 0xF0) == 0xE0) {
///         chRes = (ch & 0x0F);
///         more = 2;
///     } else if ((ch & 0xF8) == 0xF0) {
///         chRes = (ch & 0x07);
///         more = 3;
///     } else {
///         NCBI_THROW2(CStringException, eBadArgs,
///             "Source string is not in UTF8 format", 0);
///     }
///     return chRes;
/// }
///
/// TUnicodeSymbol CUtf8::DecodeNext(TUnicodeSymbol chU, char ch)
/// {
///     if ((ch & 0xC0) == 0x80) {
///         return (chU << 6) | (ch & 0x3F);
///     } else {
///         NCBI_THROW2(CStringException, eBadArgs,
///             "Source string is not in UTF8 format", 0);
///     }
///     return 0;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6884-6904
/// ```c++
/// char CUtf8::SymbolToChar(TUnicodeSymbol cp, EEncoding encoding)
/// {
///     ...
///     if ( cp <= 0x7F) {
///         return (char)cp;
///     }
///     ...
///     if (cp > 0xFF || (encoding == eEncoding_Ascii && cp > 0x7F)) {
///         NCBI_THROW2(CStringException, eConvert,
///                     "Failed to convert symbol to requested encoding", 0);
///     }
/// ```
pub(crate) fn decode_to_ascii(src: &[u8]) -> Result<(u8, usize), Utf8Error> {
    let first = src[0];
    let (mut symbol, more) = if first & 0x80 == 0 {
        (u32::from(first), 0)
    } else if first & 0xE0 == 0xC0 {
        (u32::from(first & 0x1F), 1)
    } else if first & 0xF0 == 0xE0 {
        (u32::from(first & 0x0F), 2)
    } else if first & 0xF8 == 0xF0 {
        (u32::from(first & 0x07), 3)
    } else {
        return Err(Utf8Error::DecodeFirst);
    };
    for next in 1..=more {
        let byte = src.get(next).copied().unwrap_or(0);
        if byte & 0xC0 != 0x80 {
            return Err(Utf8Error::DecodeNext);
        }
        symbol = (symbol << 6) | u32::from(byte & 0x3F);
    }
    if symbol > 0x7F {
        return Err(Utf8Error::SymbolToChar);
    }
    Ok((symbol as u8, 1 + more))
}

#[cfg(test)]
mod tests {
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
}
