//! The reader of NCBI's XML object stream (`CObjectIStreamXml` over `CIStreamBuffer`) for types
//! of an ASN.1 module, as NCBI's serial library writes them (not the "standard XML" of a DTD or
//! schema module): the file header (byte order mark, XML declaration and processing
//! instructions, `DOCTYPE`), comments, element tags (`<Type>`, `<Type_member>`, self-closed
//! `<.../>`), character data with entity and character references and CDATA sections in the
//! document's encoding, and the values of VisibleString, INTEGER, BigInt and REAL. The grammar of
//! each type stays with its reader (rpsbproc's `BlastOutput`): this module reads a SEQUENCE's
//! members in the module's order (`class`) and a SEQUENCE OF named elements (`container`).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:68-78
//! ```c++
//! CObjectIStreamXml::CObjectIStreamXml(CNcbiIstream& in, EOwnership deleteIn)
//!     : CObjectIStream(eSerial_Xml),
//!       m_TagState(eTagOutside), m_LeadingWs(0), m_Attlist(false),
//!       m_StdXml(false), m_Doctype_found(false), m_IsNil(false),
//!       m_Encoding( eEncoding_Unknown ),
//!       m_StringEncoding( eEncoding_UTF8 ),
//!       m_SkipNextTag(false)
//! ```
//!
//! The bytes come from an `AsnSource` (the stream of `asn_text_reader`, which keeps only the
//! bytes not yet read): one reader reads documents one after another from one stream, each from
//! the byte after the previous one's last `>` (`CIStreamBuffer::Close` gives the bytes it did
//! not use back to the stream, quoted at `StreamSource`).
//!
//! Where NCBI's reader throws (a `CSerialException`, the `CUtilException` of a bad number, the
//! `CEofException` of the end of the input), LOSAT's gives `XmlError::Format`: the document is
//! not read. NCBI's reader catches the end of the input in a few places (`HasMoreElements`,
//! `ReadTagData`, `ReadWord`) and then always reads a closing tag, which meets the end again: an
//! end of the input inside a document is therefore its failure, as LOSAT reads it.
//!
//! What LOSAT does not read is `XmlError::Unsupported`, refused where the reader meets it, also in
//! a document that NCBI would then reject for another reason: XML attributes (NCBI skips unknown
//! ones, with its own handling of `xmlns` and `nil`), names with a namespace prefix (NCBI strips
//! the prefix and keeps namespace maps), control characters in character data (NCBI replaces
//! them with `#` and posts a warning with its stack of the read). NCBI's BLAST XML writers produce
//! none of them.

use super::asn_text_reader::{string_to_double_posix, AsnSource, RealRefusal};
use crate::ncbi_str::is_space;

/// A failure while a document is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum XmlError {
    /// NCBI's reader throws: the document is not read (the message is LOSAT's, for tests).
    Format(String),
    /// What LOSAT does not read (an explicit refusal, the whole message).
    Unsupported(String),
}

pub(crate) type Result<T> = std::result::Result<T, XmlError>;

fn fail<T>(message: &str) -> Result<T> {
    Err(XmlError::Format(message.to_string()))
}

fn unsupported(what: &str) -> XmlError {
    XmlError::Unsupported(format!(
        "{what} in an XML document is not supported by LOSAT's rpsbproc"
    ))
}

/// `m_TagState`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TagState {
    Outside,
    InsideOpening,
    InsideClosing,
    SelfClosed,
}

/// `EEncoding` of a document (`m_Encoding`): `Unknown` until a declaration names one.
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:198-209
// ```c++
// enum EEncoding {
//     eEncoding_Unknown,
//     eEncoding_UTF8,
//     eEncoding_Ascii,
//     eEncoding_ISO8859_1,    ///< Note:  From the point of view of the C++
//     ...
//     eEncoding_Windows_1252,
//     eEncoding_CESU8
// };
// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Encoding {
    Unknown,
    Utf8,
    Ascii,
    Iso8859_1,
    Windows1252,
}

/// `CUtf8::StringToEncoding` (`None` is `eEncoding_Unknown`).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6815-6853
// ```c++
// EEncoding CUtf8::StringToEncoding(const CTempString& str)
// {
//     if (NStr::CompareNocase(str,"UTF-8")==0) {
//         return eEncoding_UTF8;
//     }
//     if (NStr::CompareNocase(str,"windows-1252")==0) {
//         return eEncoding_Windows_1252;
//     }
//     int i;
//     const char* ascii[] = {
//     "ANSI_X3.4-1968","iso-ir-6","ANSI_X3.4-1986","ISO_646.irv:1991",
//     "ASCII","ISO646-US","US-ASCII","us","IBM367","cp367","csASCII", NULL};
//     for (i=0; ascii[i]; ++i) {
//         if (NStr::CompareNocase(str,ascii[i])==0) {
//             return eEncoding_Ascii;
//         }
//     }
//     const char* iso8859_1[] = {
//     "ISO_8859-1:1987","iso-ir-100","ISO_8859-1","ISO-8859-1",
//     "latin1","l1","IBM819","CP819","csISOLatin1", NULL};
//     for (i=0; iso8859_1[i]; ++i) {
//         if (NStr::CompareNocase(str,iso8859_1[i])==0) {
//             return eEncoding_ISO8859_1;
//         }
//     }
//     const char* cesu[] = {
//     "CESU-8","csCESU8","csCESU-8",NULL};
//     for (i=0; cesu[i]; ++i) {
//         if (NStr::CompareNocase(str,iso8859_1[i])==0) {
//             return eEncoding_CESU8;
//         }
//     }
//     return eEncoding_Unknown;
// }
// ```
// (The CESU-8 loop compares with the Latin-1 names, which the loop before has matched: no name
// gives `eEncoding_CESU8`, which `SkipQDecl` would reject as it rejects `Unknown`.)
fn string_to_encoding(name: &[u8]) -> Option<Encoding> {
    let is = |names: &[&str]| {
        names
            .iter()
            .any(|n| n.as_bytes().eq_ignore_ascii_case(name))
    };
    if is(&["UTF-8"]) {
        return Some(Encoding::Utf8);
    }
    if is(&["windows-1252"]) {
        return Some(Encoding::Windows1252);
    }
    if is(&[
        "ANSI_X3.4-1968",
        "iso-ir-6",
        "ANSI_X3.4-1986",
        "ISO_646.irv:1991",
        "ASCII",
        "ISO646-US",
        "US-ASCII",
        "us",
        "IBM367",
        "cp367",
        "csASCII",
    ]) {
        return Some(Encoding::Ascii);
    }
    if is(&[
        "ISO_8859-1:1987",
        "iso-ir-100",
        "ISO_8859-1",
        "ISO-8859-1",
        "latin1",
        "l1",
        "IBM819",
        "CP819",
        "csISOLatin1",
    ]) {
        return Some(Encoding::Iso8859_1);
    }
    None
}

/// The code points of cp1252's bytes 0x80 to 0x9F.
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6850-6856
// ```c++
// static const TUnicodeSymbol s_cp1252_table[] = {
//     0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021,
//     0x02C6, 0x2030, 0x0160, 0x2039, 0x0152, 0x003F, 0x017D, 0x003F,
//     0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022, 0x2013, 0x2014,
//     0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178
// };
// ```
const CP1252: [u32; 32] = [
    0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160, 0x2039,
    0x0152, 0x003F, 0x017D, 0x003F, 0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022, 0x2013, 0x2014,
    0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178,
];

/// `CUtf8::x_AppendChar`: the UTF-8 bytes of a code point (any 32-bit value: the bits above 21
/// are cut by the byte conversion).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7040-7060
// ```c++
// CStringUTF8& CUtf8::x_AppendChar( CStringUTF8& self, TUnicodeSymbol c)
// {
//     Uint4 ch = c;
//     if (ch < 0x80) {
//         self.append(1, Uint1(ch));
//     }
//     else if (ch < 0x800) {
//         self.append(1, Uint1( (ch >>  6)         | 0xC0));
//         self.append(1, Uint1( (ch        & 0x3F) | 0x80));
//     } else if (ch < 0x10000) {
//         self.append(1, Uint1( (ch >> 12)         | 0xE0));
//         self.append(1, Uint1(((ch >>  6) & 0x3F) | 0x80));
//         self.append(1, Uint1(( ch        & 0x3F) | 0x80));
//     } else {
//         self.append(1, Uint1( (ch >> 18)         | 0xF0));
//         self.append(1, Uint1(((ch >> 12) & 0x3F) | 0x80));
//         self.append(1, Uint1(((ch >>  6) & 0x3F) | 0x80));
//         self.append(1, Uint1( (ch        & 0x3F) | 0x80));
//     }
//     return self;
// }
// ```
// (`CUtf8::AsUTF8(&chU, 1)` of one `TUnicodeSymbol`, ncbistr.hpp:5718-5758, appends it with
// this; the surrogate pair test of `x_TCharToUnicodeSymbol` is for 2-byte characters only.)
pub(crate) fn utf8_of(c: u32) -> Vec<u8> {
    if c < 0x80 {
        vec![c as u8]
    } else if c < 0x800 {
        vec![((c >> 6) | 0xC0) as u8, ((c & 0x3F) | 0x80) as u8]
    } else if c < 0x10000 {
        vec![
            ((c >> 12) | 0xE0) as u8,
            (((c >> 6) & 0x3F) | 0x80) as u8,
            ((c & 0x3F) | 0x80) as u8,
        ]
    } else {
        vec![
            ((c >> 18) | 0xF0) as u8,
            (((c >> 12) & 0x3F) | 0x80) as u8,
            (((c >> 6) & 0x3F) | 0x80) as u8,
            ((c & 0x3F) | 0x80) as u8,
        ]
    }
}

/// `CUtf8::AsUTF8(CTempString(&ch, 1), enc_in)` of one byte of 0x80 and above in a document of
/// a single-byte encoding: the byte itself for ASCII (`x_Append` copies ASCII and UTF-8), its
/// code point in UTF-8 for Latin-1 and cp1252 (`CharToSymbol`).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7079
// ```c++
// CStringUTF8& CUtf8::x_Append( CStringUTF8& self, const CTempString& src,
//     EEncoding encoding, EValidate validate)
// {
//     ...
//     if (encoding == eEncoding_UTF8 || encoding == eEncoding_Ascii) {
//         self.append(src);
//         return self;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6859-6881
// ```c++
// TUnicodeSymbol CUtf8::CharToSymbol(char c, EEncoding encoding)
// {
//     Uint1 ch = c;
//     switch (encoding)
//     {
//     ...
//     case eEncoding_Ascii:
//     case eEncoding_ISO8859_1:
//         break;
//     case eEncoding_Windows_1252:
//         if (ch > 0x7F && ch < 0xA0) {
//             return s_cp1252_table[ ch - 0x80 ];
//         }
//         break;
//     ...
//     }
//     return (TUnicodeSymbol)ch;
// }
// ```
fn utf8_of_byte(c: u8, encoding: Encoding) -> Vec<u8> {
    match encoding {
        Encoding::Unknown | Encoding::Utf8 | Encoding::Ascii => vec![c],
        Encoding::Iso8859_1 => utf8_of(u32::from(c)),
        Encoding::Windows1252 => {
            if (0x80..0xA0).contains(&c) {
                utf8_of(CP1252[usize::from(c - 0x80)])
            } else {
                utf8_of(u32::from(c))
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:156-221
// ```c++
// static inline
// bool IsBaseChar(char c)
// {
//     return
//         (c >= 'A' && c <='Z') ||
//         (c >= 'a' && c <= 'z') ||
//         (c >= '\xC0' && c <= '\xD6') ||
//         (c >= '\xD8' && c <= '\xF6') ||
//         (c >= '\xF8' && c <= '\xFF');
// }
// ...
// static inline
// bool IsFirstNameChar(char c)
// {
//     return IsLetter(c) || c == '_' || c == ':';
// }
// ...
// static inline
// bool IsExtender(char c)
// {
//     return c == '\xB7';
// }
//
// static inline
// bool IsNameChar(char c)
// {
//     return IsFirstNameChar(c) ||
//         IsDigit(c) || c == '.' || c == '-' ||
//         IsCombiningChar(c) || IsExtender(c);
// }
//
// static inline
// bool IsWhiteSpace(char c)
// {
//     return c == ' ' || c == '\t' || c == '\n' || c == '\r';
// }
//
// static inline
// bool IsEndOfTagChar(char c)
// {
//     return c == '>' || c == '/';
// }
// ```
// (`char` is signed in NCBI's build; the ranges of bytes from 0xC0 read as the unsigned ranges.)
fn is_base_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || (0xC0..=0xD6).contains(&c) || (0xD8..=0xF6).contains(&c) || c >= 0xF8
}

fn is_first_name_char(c: u8) -> bool {
    is_base_char(c) || c == b'_' || c == b':'
}

fn is_name_char(c: u8) -> bool {
    is_first_name_char(c) || c.is_ascii_digit() || c == b'.' || c == b'-' || c == 0xB7
}

fn is_white_space(c: u8) -> bool {
    c == b' ' || c == b'\t' || c == b'\n' || c == b'\r'
}

fn is_end_of_tag_char(c: u8) -> bool {
    c == b'>' || c == b'/'
}

/// `BAD_CHAR`: the control characters that XML 1.0 does not allow (NCBI's `x_VerifyChar`
/// replaces them).
// NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:811-818
// ```c++
// inline bool BAD_CHAR(int x) {
//     return (x < 0x20 && x > 0x0 && x != 0x9 && x != 0xA && x != 0xD);
// }
// inline int CObjectIStreamXml::x_VerifyChar(int x) {
//     return BAD_CHAR(x) ?
//         ReplaceVisibleChar((char)x, x_FixCharsMethod(), this, kEmptyStr, x_FixCharsSubst()) : x;
// }
// ```
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1886,302
// ```c++
// char ReplaceVisibleChar(char c, EFixNonPrint fix_method,
//     const CObjectStack* io, const CTempString& str, char subst)
// {
//     ...
//     if ( fix_method != eFNP_Replace ) {
//         string message;
//         if (io != NULL) {
//             message += io->GetStackTrace() + "\n";
//         }
//         message += "Bad char [0x" +
//                    NStr::NumericToString((unsigned char)c,0,16)+
//                    "] in string";
//         ...
//         case eFNP_ReplaceAndWarn:
//             ...
//             ERR_POST_X(7, message);
// ...
// NCBI_PARAM_ENUM_DEF(EFixNonPrint, SERIAL, WRONG_CHARS_READ, eFNP_ReplaceAndWarn);
// ```
fn bad_char(c: i32) -> bool {
    c < 0x20 && c > 0 && c != 0x9 && c != 0xA && c != 0xD
}

/// A member of a SEQUENCE: its name in the module and whether it is OPTIONAL.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Member {
    pub(crate) name: &'static str,
    pub(crate) optional: bool,
}

/// A mandatory member.
pub(crate) const fn member(name: &'static str) -> Member {
    Member {
        name,
        optional: false,
    }
}

/// An OPTIONAL member.
pub(crate) const fn optional(name: &'static str) -> Member {
    Member {
        name,
        optional: true,
    }
}

/// The reader of XML documents from one source.
pub(crate) struct XmlReader<S: AsnSource> {
    src: S,
    /// The position of the next byte to read (`m_CurrentPos`).
    pos: usize,
    state: TagState,
    /// `m_RejectedTag`: a tag name read and given back (`UndoClassMember`).
    rejected: Option<Vec<u8>>,
    /// `m_LastTag`.
    last_tag: Vec<u8>,
    /// `m_Encoding`.
    encoding: Encoding,
    /// `m_Utf8Buf` and `m_Utf8Pos`: the UTF-8 bytes of a character being returned one by one.
    utf8: Vec<u8>,
    utf8_at: usize,
}

impl<S: AsnSource> XmlReader<S> {
    pub(crate) fn new(src: S) -> Self {
        XmlReader {
            src,
            pos: 0,
            state: TagState::Outside,
            rejected: None,
            last_tag: Vec::new(),
            encoding: Encoding::Unknown,
            utf8: Vec::new(),
            utf8_at: 0,
        }
    }

    /// `PeekChar(offset)`: the byte, or the failure of the end of the input (`CEofException`).
    #[inline]
    fn peek(&mut self, offset: usize) -> Result<u8> {
        match self.src.byte_at(self.pos + offset) {
            Some(c) => Ok(c),
            None => fail("end of the input"),
        }
    }

    /// `PeekCharNoEOF(offset)`: the byte, 0 past the end of the input.
    #[inline]
    fn peek_no_eof(&mut self, offset: usize) -> u8 {
        self.src.byte_at(self.pos + offset).unwrap_or(0)
    }

    #[inline]
    fn skip(&mut self, count: usize) {
        self.pos += count;
    }

    /// `GetChar()`.
    #[inline]
    fn get_char(&mut self) -> Result<u8> {
        let c = self.peek(0)?;
        self.pos += 1;
        Ok(c)
    }

    /// `CIStreamBuffer::SkipSpaces`: skips `' '` and gives the next byte (not read).
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:213-252
    // ```c++
    // char CIStreamBuffer::SkipSpaces(void)
    //     THROWS1((CIOException))
    // {
    //     ...
    //         do {
    //             // cache current char
    //             char c = *pos;
    //             if ( c != ' ' ) { // it's not space (' ')
    //                 // point m_CurrentPos to first non space char
    //                 m_CurrentPos = pos;
    //                 // return char value
    //                 return c;
    //             }
    //             // skip space char
    //         } while ( ++pos < end );
    // ```
    fn skip_spaces(&mut self) -> Result<u8> {
        loop {
            let c = self.peek(0)?;
            if c != b' ' {
                return Ok(c);
            }
            self.pos += 1;
        }
    }

    /// `CIStreamBuffer::SkipEndOfLine(lastChar)`: a `\r\n` or `\n\r` pair is one end of line.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-537
    // ```c++
    // void CIStreamBuffer::SkipEndOfLine(char lastChar)
    //     THROWS1((CIOException))
    // {
    //     _ASSERT(lastChar == '\n' || lastChar == '\r');
    //     _ASSERT(m_CurrentPos > m_Buffer && m_CurrentPos[-1] == lastChar);
    //     m_Line++;
    //     char nextChar = PeekCharNoEOF();
    //     // lastChar either '\r' or \n'
    //     // if nextChar is compliment, skip it
    //     if ( (lastChar + nextChar) == ('\r' + '\n') )
    //         SkipChar();
    // }
    // ```
    fn skip_end_of_line(&mut self, last: u8) {
        let next = self.peek_no_eof(0);
        if i32::from(last) + i32::from(next as i8) == i32::from(b'\r' + b'\n') {
            self.pos += 1;
        }
    }

    /// `FindChar(c)`: moves to the next `c`.
    fn find_char(&mut self, c: u8) -> Result<()> {
        while self.peek(0)? != c {
            self.pos += 1;
        }
        Ok(())
    }

    /// `PeekFindChar(c, limit)`: the offset of `c` among the next `limit` bytes (fewer at the end
    /// of the input), `limit` when it is not there.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:292-306
    // ```c++
    // size_t CIStreamBuffer::PeekFindChar(char c, size_t limit)
    //     THROWS1((CIOException))
    // {
    //     _ASSERT(limit > 0);
    //     PeekCharNoEOF(limit - 1);
    //     // cache pointers
    //     const char* pos = m_CurrentPos;
    //     size_t bufferSize = m_DataEndPos - pos;
    //     if ( bufferSize != 0 ) {
    //         const void* found = memchr(pos, c, min(limit, bufferSize));
    //         if ( found )
    //             return static_cast<const char*>(found) - pos;
    //     }
    //     return limit;
    // }
    // ```
    fn peek_find_char(&mut self, c: u8, limit: usize) -> usize {
        for i in 0..limit {
            match self.src.byte_at(self.pos + i) {
                Some(b) if b == c => return i,
                Some(_) => {}
                None => break,
            }
        }
        limit
    }

    /// `SkipWS`: white space inside a tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:223-241
    // ```c++
    // char CObjectIStreamXml::SkipWS(void)
    // {
    // //    _ASSERT(InsideTag());
    //     for ( ;; ) {
    //         char c = m_Input.SkipSpaces();
    //         switch ( c ) {
    //         case '\t':
    //             m_Input.SkipChar();
    //             continue;
    //         case '\r':
    //         case '\n':
    //             m_Input.SkipChar();
    //             m_Input.SkipEndOfLine(c);
    //             continue;
    //         default:
    //             return c;
    //         }
    //     }
    // }
    // ```
    fn skip_ws(&mut self) -> Result<u8> {
        loop {
            let c = self.skip_spaces()?;
            match c {
                b'\t' => self.skip(1),
                b'\r' | b'\n' => {
                    self.skip(1);
                    self.skip_end_of_line(c);
                }
                _ => return Ok(c),
            }
        }
    }

    /// `SkipWSAndComments`: white space, processing instructions and comments between tags.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:243-306
    // ```c++
    // char CObjectIStreamXml::SkipWSAndComments(void)
    // {
    //     _ASSERT(OutsideTag());
    //     for ( ;; ) {
    //         Int8 before = m_Input.GetStreamPosAsInt8();
    //         char c = m_Input.SkipSpaces();
    //         m_LeadingWs += m_Input.GetStreamPosAsInt8() - before;
    //         switch ( c ) {
    //         case '\t':
    //             m_Input.SkipChar();
    //             continue;
    //         case '\r':
    //         case '\n':
    //             m_Input.SkipChar();
    //             m_Input.SkipEndOfLine(c);
    //             continue;
    //         case '<':
    // // https://www.w3.org/TR/xml11/#sec-pi
    //             if ( m_Input.PeekChar(1) == '?') {
    //                 m_Input.SkipChar();
    //                 Found_lt();
    //                 SkipQDecl();
    //                 return SkipWSAndComments();
    //             }
    // // https://www.w3.org/TR/xml11/#sec-comments
    //             if ( m_Input.PeekChar(1) == '!' &&
    //                  m_Input.PeekChar(2) == '-' &&
    //                  m_Input.PeekChar(3) == '-' ) {
    //                 // start of comment
    //                 m_Input.SkipChars(4);
    //                 if (m_Input.PeekChar(0) == '-' &&
    //                     m_Input.PeekChar(1) == '-') {
    //                         ThrowError(fFormatError,
    //                                    "double-hyphen '--' is not allowed in XML comments");
    //                 }
    //                 for ( ;; ) {
    //                     m_Input.FindChar('-');
    //                     if ( m_Input.PeekChar(1) == '-' ) {
    //                         // --
    //                         if ( m_Input.PeekChar(2) == '>' ) {
    //                             // -->
    //                             m_Input.SkipChars(3);
    //                             break;
    //                         }
    //                         else {
    //                             // --[^>]
    //                             ThrowError(fFormatError,
    //                                        "double-hyphen '--' is not allowed in XML comments");
    //                         }
    //                     }
    //                     else {
    //                         // -[^-]
    //                         m_Input.SkipChars(2);
    //                     }
    //
    //                 }
    //                 continue; // skip the next WS or comment
    //             }
    //             return '<';
    //         default:
    //             return c;
    //         }
    //     }
    // }
    // ```
    // (`m_LeadingWs`, the spaces counted here, is given to the next `ReadTagData`; every tag
    // read resets it (`Found_lt`, `Found_gt`), and the character data of a member is read just
    // after its opening tag, so it is always 0 there: LOSAT does not keep it.)
    fn skip_ws_and_comments(&mut self) -> Result<u8> {
        self.src.release_before(self.pos);
        loop {
            let c = self.skip_spaces()?;
            match c {
                b'\t' => self.skip(1),
                b'\r' | b'\n' => {
                    self.skip(1);
                    self.skip_end_of_line(c);
                }
                b'<' => {
                    if self.peek(1)? == b'?' {
                        self.skip(1);
                        self.state = TagState::InsideOpening;
                        self.skip_qdecl()?;
                        continue;
                    }
                    if self.peek(1)? == b'!' && self.peek(2)? == b'-' && self.peek(3)? == b'-' {
                        self.skip(4);
                        if self.peek(0)? == b'-' && self.peek(1)? == b'-' {
                            return fail("double-hyphen '--' is not allowed in XML comments");
                        }
                        loop {
                            self.find_char(b'-')?;
                            if self.peek(1)? == b'-' {
                                if self.peek(2)? == b'>' {
                                    self.skip(3);
                                    break;
                                }
                                return fail("double-hyphen '--' is not allowed in XML comments");
                            }
                            self.skip(2);
                        }
                        continue;
                    }
                    return Ok(b'<');
                }
                _ => return Ok(c),
            }
        }
    }

    /// `EndTag`: the end (`>`) of the tag being read.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:308-340
    // ```c++
    // void CObjectIStreamXml::EndTag(void)
    // {
    //     char c = SkipWS();
    //     if (m_Attlist) {
    //         ...
    //     }
    //     if ( c != '>' ) {
    //         c = ReadUndefinedAttributes();
    //         if ( c != '>' ) {
    //             ThrowError(fFormatError, "'>' expected");
    //         }
    //     }
    //     m_Input.SkipChar();
    //     Found_gt();
    // }
    // ```
    fn end_tag(&mut self) -> Result<()> {
        let c = self.skip_ws()?;
        if c != b'>' && self.read_undefined_attributes()? != b'>' {
            return fail("'>' expected");
        }
        self.skip(1);
        self.state = TagState::Outside;
        Ok(())
    }

    /// `ReadUndefinedAttributes`: the end of a tag, or an attribute (refused).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:914-931
    // ```c++
    // char CObjectIStreamXml::ReadUndefinedAttributes(void)
    // {
    //     char c;
    //     m_Attlist = true;
    //     for (;;) {
    //         c = SkipWS();
    //         if (IsEndOfTagChar(c)) {
    //             m_Attlist = false;
    //             break;
    //         }
    //         CTempString tagName = ReadName(c);
    //         if (!tagName.empty()) {
    //             string value;
    //             ReadAttributeValue(value, true);
    //         }
    //     }
    //     return c;
    // }
    // ```
    // (`ReadName` throws for a byte that cannot start a name.)
    fn read_undefined_attributes(&mut self) -> Result<u8> {
        let c = self.skip_ws()?;
        if is_end_of_tag_char(c) {
            return Ok(c);
        }
        if is_first_name_char(c) {
            return Err(unsupported("an XML attribute"));
        }
        fail("Name begins with an invalid character")
    }

    /// `EndOpeningTagSelfClosed`: ends an opening tag being read; whether it was self-closed.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:342-378
    // ```c++
    // bool CObjectIStreamXml::EndOpeningTagSelfClosed(void)
    // {
    //     if (!StackIsEmpty() && TopFrame().GetNotag()) {
    //         if (SelfClosedTag()) {
    //             return true;
    //         }
    //     }
    //     if( InsideOpeningTag() ) {
    //         char c = SkipWS();
    //         if (m_Attlist) {
    //             return false;
    //         }
    //         if ( c == '/' && m_Input.PeekChar(1) == '>' ) {
    //             // end of self closed tag
    //             m_Input.SkipChars(2);
    //             Found_slash_gt();
    //             return true;
    //         }
    //
    //         if ( c != '>' ) {
    //             c = ReadUndefinedAttributes();
    //             if ( c == '/' && m_Input.PeekChar(1) == '>' ) {
    //                 // end of self closed tag
    //                 m_Input.SkipChars(2);
    //                 Found_slash_gt();
    //                 return true;
    //             }
    //             if ( c != '>' )
    //                 ThrowError(fFormatError, "end of tag expected");
    //         }
    //
    //         // end of open tag
    //         m_Input.SkipChar(); // '>'
    //         Found_gt();
    //     }
    //     return SelfClosedTag();
    // }
    // ```
    // (No frame of these types is read without its tag (`GetNotag`).)
    fn end_opening_tag_self_closed(&mut self) -> Result<bool> {
        if self.state == TagState::InsideOpening {
            let mut c = self.skip_ws()?;
            if c == b'/' && self.peek(1)? == b'>' {
                self.skip(2);
                self.state = TagState::SelfClosed;
                return Ok(true);
            }
            if c != b'>' {
                c = self.read_undefined_attributes()?;
                if c == b'/' && self.peek(1)? == b'>' {
                    self.skip(2);
                    self.state = TagState::SelfClosed;
                    return Ok(true);
                }
                if c != b'>' {
                    return fail("end of tag expected");
                }
            }
            self.skip(1);
            self.state = TagState::Outside;
        }
        Ok(self.state == TagState::SelfClosed)
    }

    /// `BeginData`: ends an opening tag being read.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrxml.inl:111-131
    // ```c++
    // inline
    // void CObjectIStreamXml::EndOpeningTag(void)
    // {
    //     _ASSERT(InsideOpeningTag());
    //     EndTag();
    // }
    // ...
    // inline
    // void CObjectIStreamXml::BeginData(void)
    // {
    //     if ( InsideOpeningTag() )
    //         EndOpeningTag();
    //     _ASSERT(OutsideTag());
    // }
    // ```
    fn begin_data(&mut self) -> Result<()> {
        if self.state == TagState::InsideOpening {
            self.end_tag()?;
        }
        Ok(())
    }

    /// `BeginOpeningTag`: reads `<` and gives the byte after it (not read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:405-418
    // ```c++
    // char CObjectIStreamXml::BeginOpeningTag(void)
    // {
    //     BeginData();
    //     // find beginning '<'
    //     char c = SkipWSAndComments();
    //     if ( c != '<' )
    //         ThrowError(fFormatError, "'<' expected");
    //     c = m_Input.PeekChar(1);
    //     if ( c == '/' )
    //         ThrowError(fFormatError, "unexpected '</'");
    //     m_Input.SkipChar();
    //     Found_lt();
    //     return c;
    // }
    // ```
    fn begin_opening_tag(&mut self) -> Result<u8> {
        self.begin_data()?;
        if self.skip_ws_and_comments()? != b'<' {
            return fail("'<' expected");
        }
        let c = self.peek(1)?;
        if c == b'/' {
            return fail("unexpected '</'");
        }
        self.skip(1);
        self.state = TagState::InsideOpening;
        Ok(c)
    }

    /// `BeginClosingTag`: reads `</` and gives the byte after it (not read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:420-430
    // ```c++
    // char CObjectIStreamXml::BeginClosingTag(void)
    // {
    //     BeginData();
    //     // find beginning '<'
    //     char c = SkipWSAndComments();
    //     if ( c != '<' || m_Input.PeekChar(1) != '/' )
    //         ThrowError(fFormatError, "'</' expected");
    //     m_Input.SkipChars(2);
    //     Found_lt_slash();
    //     return m_Input.PeekChar();
    // }
    // ```
    fn begin_closing_tag(&mut self) -> Result<u8> {
        self.begin_data()?;
        if self.skip_ws_and_comments()? != b'<' || self.peek(1)? != b'/' {
            return fail("'</' expected");
        }
        self.skip(2);
        self.state = TagState::InsideClosing;
        self.peek(0)
    }

    /// `ReadName(c)`: the name that starts with `c` (not read yet), into `last_tag`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:432-469
    // ```c++
    // CTempString CObjectIStreamXml::ReadName(char c)
    // {
    //     _ASSERT(InsideTag());
    //     if ( !IsFirstNameChar(c) )
    //         ThrowError(fFormatError,
    //             "Name begins with an invalid character: #"
    //             +NStr::UIntToString((unsigned int)c));
    //
    //     // find end of tag name
    //     size_t i = 1, iColon = 0;
    //     while ( IsNameChar(c = m_Input.PeekChar(i)) ) {
    //         if (!m_Doctype_found && c == ':') {
    //             iColon = i+1;
    //         }
    //         ++i;
    //     }
    //
    //     // save beginning of tag name
    //     const char* ptr = m_Input.GetCurrentPos();
    //     m_LastTag = string(ptr+iColon, i-iColon);
    //     string ns_prefix;
    //     if (iColon > 1) {
    //         ns_prefix = string(ptr, iColon-1);
    //     }
    //     m_Input.SkipChars(i);
    //
    //     // check end of tag name
    //     if (c == '\n' || c == '\r') {
    //         m_Input.SkipChar();
    //         m_Input.SkipEndOfLine(c);
    //     }
    //     if (iColon > 1) {
    //         ...
    // ```
    // (`m_Doctype_found` is false for the types of an ASN.1 module (ReadFileHeader:
    // `GetDataSpec() == EDataSpec::eDTD`). A `:` after the first byte makes a namespace prefix,
    // refused; the `xmlns` and `nil` names of `m_Attlist` are attributes, refused before.)
    fn read_name(&mut self, c: u8) -> Result<()> {
        if !is_first_name_char(c) {
            return fail("Name begins with an invalid character");
        }
        let mut i = 1;
        let next = loop {
            let next = self.peek(i)?;
            if !is_name_char(next) {
                break next;
            }
            if next == b':' {
                return Err(unsupported("a name with a namespace prefix"));
            }
            i += 1;
        };
        self.last_tag.clear();
        let name = self.src.bytes(self.pos..self.pos + i);
        self.last_tag.extend_from_slice(name);
        self.skip(i);
        if next == b'\n' || next == b'\r' {
            self.skip(1);
            self.skip_end_of_line(next);
        }
        Ok(())
    }

    /// `UndoClassMember`: gives back the name of the opening tag being read.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2379-2389
    // ```c++
    // void CObjectIStreamXml::UndoClassMember(void)
    // {
    //     if (InsideOpeningTag()) {
    //         m_RejectedTag = m_LastTag;
    //         m_TagState = eTagOutside;
    //         m_LeadingWs = 0;
    //     ...
    // ```
    fn undo_class_member(&mut self) {
        if self.state == TagState::InsideOpening {
            self.rejected = Some(self.last_tag.clone());
            self.state = TagState::Outside;
        }
    }

    /// `SkipQDecl`: the rest of a processing instruction after `<`; an `encoding` pseudo-attribute
    /// sets the document's encoding.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:545-580
    // ```c++
    // void CObjectIStreamXml::SkipQDecl(void)
    // {
    //     _ASSERT(InsideOpeningTag());
    //     m_Input.SkipChar();
    //
    //     CTempString tagName;
    //     tagName = ReadName( SkipWS());
    // //    _ASSERT(tagName == "xml");
    //     for (;;) {
    //         char ch = SkipWS();
    //         if (ch == '?') {
    //             break;
    //         }
    //         tagName = ReadName(ch);
    //         string value;
    //         ReadAttributeValue(value);
    //         if (tagName == "encoding") {
    //             m_Encoding = CUtf8::StringToEncoding(value);
    //             if (m_Encoding == eEncoding_Unknown || m_Encoding > eEncoding_Windows_1252) {
    //                 ThrowError(fFormatError, "unsupported encoding: " + value);
    //             }
    //             break;
    //         }
    //     }
    //     for ( ;; ) {
    //         m_Input.FindChar('?');
    //         if ( m_Input.PeekChar(1) == '>' ) {
    //             // ?>
    //             m_Input.SkipChars(2);
    //             Found_gt();
    //             return;
    //         }
    //         else
    //             m_Input.SkipChar();
    //     }
    // }
    // ```
    fn skip_qdecl(&mut self) -> Result<()> {
        self.skip(1);
        let c = self.skip_ws()?;
        self.read_name(c)?;
        loop {
            let ch = self.skip_ws()?;
            if ch == b'?' {
                break;
            }
            self.read_name(ch)?;
            let is_encoding = self.last_tag == b"encoding";
            let value = self.read_attribute_value()?;
            if is_encoding {
                match string_to_encoding(&value) {
                    Some(encoding) => self.encoding = encoding,
                    None => return fail("unsupported encoding"),
                }
                break;
            }
        }
        loop {
            self.find_char(b'?')?;
            if self.peek(1)? == b'>' {
                self.skip(2);
                self.state = TagState::Outside;
                return Ok(());
            }
            self.skip(1);
        }
    }

    /// `ReadAttributeValue(value)` of a processing instruction (`m_Attlist` false: the closing
    /// quote is read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:891-912
    // ```c++
    // void CObjectIStreamXml::ReadAttributeValue(string& value, bool skipClosing)
    // {
    //     if ( SkipWS() != '=' )
    //         ThrowError(fFormatError, "'=' expected");
    //     m_Input.SkipChar(); // '='
    //     char startChar = SkipWS();
    //     if ( startChar != '\'' && startChar != '\"' )
    //         ThrowError(fFormatError, "attribute value must start with ' or \"");
    //     m_Input.SkipChar();
    //     bool encoded = false;
    //     for ( ;; ) {
    //         int c = ReadEncodedChar(startChar,eStringTypeUTF8,encoded);
    //         if ( c < 0 )
    //             break;
    //         if (c != 0) {
    //             value += char(c);
    //         }
    //     }
    //     if (!m_Attlist || skipClosing) {
    //         m_Input.SkipChar();
    //     }
    // }
    // ```
    fn read_attribute_value(&mut self) -> Result<Vec<u8>> {
        if self.skip_ws()? != b'=' {
            return fail("'=' expected");
        }
        self.skip(1);
        let start = self.skip_ws()?;
        if start != b'\'' && start != b'"' {
            return fail("attribute value must start with ' or \"");
        }
        self.skip(1);
        let mut value = Vec::new();
        loop {
            let c = self.read_encoded_char(start)?;
            if c < 0 {
                break;
            }
            if c != 0 {
                value.push(c as u8);
            }
        }
        self.skip(1);
        Ok(value)
    }

    /// `SkipAttributeValue(c)`: a quoted string of the `DOCTYPE`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:537-543
    // ```c++
    // void CObjectIStreamXml::SkipAttributeValue(char c)
    // {
    //     _ASSERT(InsideOpeningTag());
    //     m_Input.SkipChar();
    //     m_Input.FindChar(c);
    //     m_Input.SkipChar();
    // }
    // ```
    fn skip_attribute_value(&mut self, c: u8) -> Result<()> {
        self.skip(1);
        self.find_char(c)?;
        self.skip(1);
        Ok(())
    }

    /// `ReadEscapedChar(endingChar, &encoded)`: the next character (a byte, or the value of an
    /// entity or character reference, which is `encoded`), or -1 at `ending` (not read). A
    /// character reference of 2^31 and more is negative, as NCBI's `int` of its `unsigned`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:722-800
    // ```c++
    // int CObjectIStreamXml::ReadEscapedChar(char endingChar, bool* encoded)
    // {
    //     char c = m_Input.PeekChar();
    //     if (encoded) {
    //         *encoded = false;
    //     }
    //     if ( c == '&' ) {
    //         if (encoded) {
    //             *encoded = true;
    //         }
    //         m_Input.SkipChar();
    //         const size_t limit = 32;
    //         size_t offset = m_Input.PeekFindChar(';', limit);
    //         if ( offset >= limit )
    //             ThrowError(fFormatError, "entity reference is too long");
    //         const char* p = m_Input.GetCurrentPos(); // save entity string pointer
    //         m_Input.SkipChars(offset + 1); // skip it
    //         if ( offset == 0 )
    //             ThrowError(fFormatError, "invalid entity reference");
    //         if ( *p == '#' ) {
    //             const char* end = p + offset;
    //             ++p;
    //             // char ref
    //             if ( p == end )
    //                 ThrowError(fFormatError, "invalid char reference");
    //             unsigned v = 0;
    //             if ( *p == 'x' ) {
    //                 // hex
    //                 if ( ++p == end )
    //                     ThrowError(fFormatError, "invalid char reference");
    //                 do {
    //                     c = *p++;
    //                     if ( c >= '0' && c <= '9' )
    //                         v = v * 16 + (c - '0');
    //                     else if ( c >= 'A' && c <='F' )
    //                         v = v * 16 + (c - 'A' + 0xA);
    //                     else if ( c >= 'a' && c <='f' )
    //                         v = v * 16 + (c - 'a' + 0xA);
    //                     else
    //                         ThrowError(fFormatError,
    //                             "invalid symbol in char reference");
    //                 } while ( p < end );
    //             }
    //             else {
    //                 // dec
    //                 if ( p == end )
    //                     ThrowError(fFormatError, "invalid char reference");
    //                 do {
    //                     c = *p++;
    //                     if ( c >= '0' && c <= '9' )
    //                         v = v * 10 + (c - '0');
    //                     else
    //                         ThrowError(fFormatError,
    //                             "invalid symbol in char reference");
    //                 } while ( p < end );
    //             }
    //             return v;
    //         }
    //         else {
    //             CTempString e(p, offset);
    //             if ( e == "lt" )
    //                 return '<';
    //             if ( e == "gt" )
    //                 return '>';
    //             if ( e == "amp" )
    //                 return '&';
    //             if ( e == "apos" )
    //                 return '\'';
    //             if ( e == "quot" )
    //                 return '"';
    //             ThrowError(fFormatError, "unknown entity name: " + string(e));
    //         }
    //     }
    //     else if ( c == endingChar ) {
    //         return -1;
    //     }
    //     m_Input.SkipChar();
    //     return c & 0xFF;
    // }
    // ```
    fn read_escaped_char(&mut self, ending: u8) -> Result<(i32, bool)> {
        let c = self.peek(0)?;
        if c == b'&' {
            self.skip(1);
            const LIMIT: usize = 32;
            let offset = self.peek_find_char(b';', LIMIT);
            if offset >= LIMIT {
                return fail("entity reference is too long");
            }
            let start = self.pos;
            self.skip(offset + 1);
            if offset == 0 {
                return fail("invalid entity reference");
            }
            let entity = self.src.bytes(start..start + offset);
            if entity[0] == b'#' {
                let digits = &entity[1..];
                if digits.is_empty() {
                    return fail("invalid char reference");
                }
                let mut v: u32 = 0;
                if digits[0] == b'x' {
                    let digits = &digits[1..];
                    if digits.is_empty() {
                        return fail("invalid char reference");
                    }
                    for &d in digits {
                        let value = match d {
                            b'0'..=b'9' => d - b'0',
                            b'A'..=b'F' => d - b'A' + 0xA,
                            b'a'..=b'f' => d - b'a' + 0xA,
                            _ => return fail("invalid symbol in char reference"),
                        };
                        v = v.wrapping_mul(16).wrapping_add(u32::from(value));
                    }
                } else {
                    for &d in digits {
                        if !d.is_ascii_digit() {
                            return fail("invalid symbol in char reference");
                        }
                        v = v.wrapping_mul(10).wrapping_add(u32::from(d - b'0'));
                    }
                }
                return Ok((v as i32, true));
            }
            let value = match entity {
                b"lt" => b'<',
                b"gt" => b'>',
                b"amp" => b'&',
                b"apos" => b'\'',
                b"quot" => b'"',
                _ => return fail("unknown entity name"),
            };
            return Ok((i32::from(value), true));
        }
        if c == ending {
            return Ok((-1, false));
        }
        self.skip(1);
        Ok((i32::from(c), false))
    }

    /// `ReadEncodedChar(endingChar, type, encoded)` of a VisibleString (or of a UTF-8 attribute
    /// value: the output encoding is UTF-8 either way, `m_StringEncoding`): the bytes of the
    /// next character in UTF-8, one per call, or a negative value at `ending`; a control
    /// character is refused (`x_VerifyChar`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:819-868
    // ```c++
    // inline
    // int CObjectIStreamXml::ReadEncodedChar(char endingChar, EStringType type, bool& encoded)
    // {
    //     return x_VerifyChar(x_ReadEncodedChar(endingChar,type,encoded));
    // }
    //
    // int CObjectIStreamXml::x_ReadEncodedChar(char endingChar, EStringType type, bool& encoded)
    // {
    //     EEncoding enc_out( type == eStringTypeUTF8 ? eEncoding_UTF8 : m_StringEncoding);
    //     EEncoding enc_in(m_Encoding == eEncoding_Unknown ? eEncoding_UTF8 : m_Encoding);
    //
    //     if (enc_out == eEncoding_UTF8 &&
    //         !m_Utf8Buf.empty() && m_Utf8Pos != m_Utf8Buf.end()) {
    //         if (++m_Utf8Pos != m_Utf8Buf.end()) {
    //             return *m_Utf8Pos & 0xFF;
    //         } else {
    //             m_Utf8Buf.clear();
    //         }
    //     }
    //     int c = ReadEscapedChar(endingChar, &encoded);
    //     if (c < 0) {
    //         return c;
    //     }
    //     if (enc_out != eEncoding_Unknown) {
    //         if (encoded) {
    //             TUnicodeSymbol chU = c;
    //             if (enc_out == eEncoding_UTF8) {
    //                 m_Utf8Buf = CUtf8::AsUTF8( &chU, 1);
    //                 m_Utf8Pos = m_Utf8Buf.begin();
    //                 return *m_Utf8Pos & 0xFF;
    //             } else {
    //                 return CUtf8::SymbolToChar( chU, enc_out);
    //             }
    //         }
    //         if (enc_in != enc_out) {
    //             if (enc_out != eEncoding_UTF8) {
    //                 ...
    //             }
    //             if ((c & 0x80) == 0) {
    //                 return c;
    //             }
    //             char ch = (char)c;
    //             m_Utf8Buf = CUtf8::AsUTF8( CTempString(&ch,1), enc_in);
    //             m_Utf8Pos = m_Utf8Buf.begin();
    //             return *m_Utf8Pos & 0xFF;
    //         }
    //     }
    //     return c;
    // }
    // ```
    fn read_encoded_char(&mut self, ending: u8) -> Result<i32> {
        if !self.utf8.is_empty() && self.utf8_at != self.utf8.len() {
            self.utf8_at += 1;
            if self.utf8_at != self.utf8.len() {
                return Ok(i32::from(self.utf8[self.utf8_at]));
            }
            self.utf8.clear();
        }
        let (c, encoded) = self.read_escaped_char(ending)?;
        if c < 0 {
            return Ok(c);
        }
        let c = if encoded {
            self.utf8 = utf8_of(c as u32);
            self.utf8_at = 0;
            i32::from(self.utf8[0])
        } else if !matches!(self.encoding, Encoding::Unknown | Encoding::Utf8) && c & 0x80 != 0 {
            self.utf8 = utf8_of_byte(c as u8, self.encoding);
            self.utf8_at = 0;
            i32::from(self.utf8[0])
        } else {
            c
        };
        if bad_char(c) {
            return Err(XmlError::Unsupported(format!(
                "a control character (0x{c:02x}) in the character data of an XML document is not supported by LOSAT's rpsbproc (NCBI replaces it and posts a warning)"
            )));
        }
        Ok(c)
    }

    /// `ReadCDSection(str)`: a CDATA section where character data meets `<`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1271-1295
    // ```c++
    // bool CObjectIStreamXml::ReadCDSection(string& str)
    // // http://www.w3.org/TR/2000/REC-xml-20001006#dt-cdsection
    // // must begin with <![CDATA[
    // // must end with   ]]>
    // {
    //     if (m_Input.PeekChar() != '<' || m_Input.PeekChar(1) != '!') {
    //         return false;
    //     }
    //     m_Input.SkipChars(2);
    //     const char* open = "[CDATA[";
    //     for ( ; *open; ++open) {
    //         if (m_Input.PeekChar() != *open) {
    //             ThrowError(fFormatError, "CDATA section expected");
    //         }
    //         m_Input.SkipChar();
    //     }
    //     while ( m_Input.PeekChar(0) != ']' ||
    //             m_Input.PeekChar(1) != ']' ||
    //             m_Input.PeekChar(2) != '>') {
    //         str += m_Input.PeekChar();
    //         m_Input.SkipChar();
    //     }
    //     m_Input.SkipChars(3);
    //     return true;
    // }
    // ```
    fn read_cd_section(&mut self, out: &mut Vec<u8>) -> Result<bool> {
        if self.peek(0)? != b'<' || self.peek(1)? != b'!' {
            return Ok(false);
        }
        self.skip(2);
        for &c in b"[CDATA[" {
            if self.peek(0)? != c {
                return fail("CDATA section expected");
            }
            self.skip(1);
        }
        while self.peek(0)? != b']' || self.peek(1)? != b']' || self.peek(2)? != b'>' {
            let c = self.peek(0)?;
            out.push(c);
            self.skip(1);
        }
        self.skip(3);
        Ok(true)
    }

    /// `ReadTagData(str, eStringTypeVisible)`: the character data up to the next tag, with
    /// NCBI's end-of-line handling (a first `\r` is dropped; from it to the next `\n`, every
    /// `\r` is a `\n`) and NUL characters dropped.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1297-1354
    // ```c++
    // void CObjectIStreamXml::ReadTagData(string& str, EStringType type)
    // {
    //     if (m_LeadingWs != 0) {
    //         str.append(m_LeadingWs, ' ');
    //         m_LeadingWs = 0;
    //     }
    //     BeginData();
    //     bool encoded = false;
    //     bool CR = false;
    //     try {
    //         for ( ;; ) {
    //             int c = ReadEncodedChar(m_Attlist ? '\"' : '<', type, encoded);
    //             if ( c < 0 ) {
    //                 if (m_Attlist || !ReadCDSection(str)) {
    //                     break;
    //                 }
    //                 CR = false;
    //                 continue;
    //             }
    //             if (c == 0) {
    //                 continue;
    //             }
    //             if (CR) {
    //                 if (c == '\n') {
    //                     CR = false;
    //                 } else if (c == '\r') {
    //                     c = '\n';
    //                 }
    //             } else if (c == '\r') {
    //                 CR = true;
    //                 continue;
    //             }
    //             if (m_Attlist && !encoded && IsWhiteSpace((char)c)) {
    //                 c = ' ';
    //             }
    //             str += (char)c;
    //             ...
    //         }
    //     } catch (CEofException&) {
    //     }
    //     str.reserve(str.size());
    // }
    // ```
    fn read_tag_data(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.begin_data()?;
        let mut cr = false;
        loop {
            let c = self.read_encoded_char(b'<')?;
            if c < 0 {
                if !self.read_cd_section(out)? {
                    break;
                }
                cr = false;
                continue;
            }
            if c == 0 {
                continue;
            }
            let mut c = c as u8;
            if cr {
                if c == b'\n' {
                    cr = false;
                } else if c == b'\r' {
                    c = b'\n';
                }
            } else if c == b'\r' {
                cr = true;
                continue;
            }
            out.push(c);
        }
        Ok(())
    }

    /// `ReadWord(str)`: the character data up to white space or the next tag, after white space.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1356-1374
    // ```c++
    // void CObjectIStreamXml::ReadWord(string& str, EStringType type)
    // {
    //     BeginData();
    //     bool encoded = false;
    //     SkipWS();
    //     try {
    //         for ( ;; ) {
    //             int c = ReadEncodedChar(m_Attlist ? '\"' : '<', type, encoded);
    //             if ( c < 0  || IsWhiteSpace((char)c)) {
    //                 break;
    //             }
    //             if (c != 0) {
    //                 str += (char)c;
    //             }
    //         }
    //     } catch (CEofException&) {
    //     }
    //     str.reserve(str.size());
    // }
    // ```
    fn read_word(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.begin_data()?;
        self.skip_ws()?;
        loop {
            let c = self.read_encoded_char(b'<')?;
            if c < 0 || is_white_space(c as u8) {
                break;
            }
            if c != 0 {
                out.push(c as u8);
            }
        }
        Ok(())
    }

    /// `CIStreamBuffer::GetInt4` and `GetInt8` up to their sign: whether the number is
    /// negative and its magnitude, at most `max` (`kMax_I4` or `kMax_I8`) or one more when
    /// negative.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:617-670
    // ```c++
    // char CIStreamBuffer::SkipWs(void)
    // {
    //     char c;
    //     while (isspace((unsigned char)(c = GetChar())))
    //     ;
    //     return c;
    // }
    //
    //
    // Int4 CIStreamBuffer::GetInt4(void)
    //     THROWS1((CIOException,CUtilException))
    // {
    //     bool sign;
    //     char c = SkipWs();
    //     switch ( c ) {
    //     case '-':
    //         sign = true;
    //         c = GetChar();
    //         break;
    //     case '+':
    //         sign = false;
    //         c = GetChar();
    //         break;
    //     default:
    //         sign = false;
    //         break;
    //     }
    //
    //     Uint4 n = c - '0';
    //     if ( n > 9 )
    //         BadNumber();
    //
    //     // overflow limits
    //     const Uint4 kMaxBeforeMul = kMax_I4/10;
    //     const Uint1 kMaxLimitAdd = Uint1(kMax_I4%10 + sign);
    //
    //     for ( ;; ) {
    //         Uint1 d = (Uint1)(PeekCharNoEOF() - '0');
    //         if ( d > 9 )
    //             break;
    //         SkipChar();
    //
    //         // check multiplication overflow
    //         if ( n > kMaxBeforeMul || (n == kMaxBeforeMul && d > kMaxLimitAdd) )
    //             NumberOverflow();
    //
    //         n = n * 10 + d;
    //     }
    //     if ( sign )
    //         return -Int4(n);
    //     else
    //         return n;
    // }
    // ```
    // (`GetInt8`, strbuffer.cpp:703-745, is the same with `Uint8` and `kMax_I8`.)
    fn get_number(&mut self, max: u64) -> Result<(bool, u64)> {
        let mut c = loop {
            let c = self.get_char()?;
            if !is_space(c) {
                break c;
            }
        };
        let negative = match c {
            b'-' => {
                c = self.get_char()?;
                true
            }
            b'+' => {
                c = self.get_char()?;
                false
            }
            _ => false,
        };
        if !c.is_ascii_digit() {
            return fail("bad number");
        }
        let mut n = u64::from(c - b'0');
        let max_before_mul = max / 10;
        let max_limit_add = (max % 10) as u8 + u8::from(negative);
        loop {
            let d = self.peek_no_eof(0).wrapping_sub(b'0');
            if d > 9 {
                break;
            }
            self.skip(1);
            if n > max_before_mul || (n == max_before_mul && d > max_limit_add) {
                return fail("number overflow");
            }
            n = n * 10 + u64::from(d);
        }
        Ok((negative, n))
    }

    /// `ReadInt4`: an INTEGER.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:987-994
    // ```c++
    // Int4 CObjectIStreamXml::ReadInt4(void)
    // {
    //     if (ExpectSpecialCase()!=0 && UseSpecialCaseRead()) {
    //         return x_UseMemberDefault<Int4>();
    //     }
    //     BeginData();
    //     return m_Input.GetInt4();
    // }
    // ```
    // (No member of these types has a default: `ExpectSpecialCase` is 0. For `-2147483648`,
    // `-Int4(n)` negates `INT_MIN`, a signed overflow; LOSAT gives the number's value, `INT_MIN`,
    // which both NCBI builds give too (PS16d oracle `p_intmin`).)
    pub(crate) fn int4(&mut self) -> Result<i32> {
        self.begin_data()?;
        let (negative, n) = self.get_number(i32::MAX as u64)?;
        let n = n as i64;
        Ok((if negative { -n } else { n }) as i32)
    }

    /// `ReadInt8`: a BigInt.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1005-1012
    // ```c++
    // Int8 CObjectIStreamXml::ReadInt8(void)
    // {
    //     if (ExpectSpecialCase()!=0 && UseSpecialCaseRead()) {
    //         return x_UseMemberDefault<Int8>();
    //     }
    //     BeginData();
    //     return m_Input.GetInt8();
    // }
    // ```
    // (`-9223372036854775808` is read as its value, as `int4` reads `-2147483648`.)
    pub(crate) fn int8(&mut self) -> Result<i64> {
        self.begin_data()?;
        let (negative, n) = self.get_number(i64::MAX as u64)?;
        let n = i128::from(n);
        Ok((if negative { -n } else { n }) as i64)
    }

    /// `ReadDouble`: a REAL.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1027-1042
    // ```c++
    // double CObjectIStreamXml::ReadDouble(void)
    // {
    //     if (ExpectSpecialCase()!=0 && UseSpecialCaseRead()) {
    //         return x_UseMemberDefault<double>();
    //     }
    //     string s;
    //     ReadWord(s);
    //     char* endptr;
    //     double result = NStr::StringToDoublePosix(s.c_str(), &endptr, NStr::fDecimalPosixFinite);
    //     while (IsWhiteSpace(*endptr)) {
    //         ++endptr;
    //     }
    //     if ( *endptr != 0 )
    //         ThrowError(fFormatError, "invalid float number");
    //     return result;
    // }
    // ```
    pub(crate) fn real(&mut self) -> Result<f64> {
        let mut word = Vec::new();
        self.read_word(&mut word)?;
        let (value, mut end) = string_to_double_posix(&word, true).map_err(|refusal| {
            unsupported(match refusal {
                RealRefusal::Digits => "a REAL of more than 18 significant digits",
                RealRefusal::ExponentRange => {
                    "a REAL whose exponent is outside the range of NCBI's int arithmetic"
                }
            })
        })?;
        while end < word.len() && is_white_space(word[end]) {
            end += 1;
        }
        if end != word.len() {
            return fail("invalid float number");
        }
        Ok(value)
    }

    /// `ReadString(str, eStringTypeVisible)`: a VisibleString (empty for a self-closed tag).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1233-1258
    // ```c++
    // void CObjectIStreamXml::ReadString(string& str, EStringType type)
    // {
    //     str.erase();
    //     if (ExpectSpecialCase()!=0 && UseSpecialCaseRead()) {
    //         ...
    //     }
    //     if (EndOpeningTagSelfClosed()) {
    //         return;
    //     }
    //     ...
    //     ReadTagData(str, type);
    // }
    // ```
    pub(crate) fn string(&mut self, out: &mut Vec<u8>) -> Result<()> {
        out.clear();
        if self.end_opening_tag_self_closed()? {
            return Ok(());
        }
        self.read_tag_data(out)
    }

    /// `ReadFileHeader`: the byte order mark, the XML declaration, processing instructions,
    /// comments and the `DOCTYPE` before the first element, whose name it gives back
    /// (`UndoClassMember`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:582-670
    // ```c++
    // string CObjectIStreamXml::ReadFileHeader(void)
    // {
    // // check for UTF8 Byte Order Mark (EF BB BF)
    // // http://unicode.org/faq/utf_bom.html#BOM
    //     {
    //         char c = m_Input.PeekChar();
    //         if ((unsigned char)c == 0xEF) {
    //             if ((unsigned char)m_Input.PeekChar(1) == 0xBB &&
    //                 (unsigned char)m_Input.PeekChar(2) == 0xBF) {
    //                 m_Input.SkipChars(3);
    //                 m_Encoding = eEncoding_UTF8;
    //             }
    //         }
    //     }
    //
    //     m_Doctype_found = false;
    //     ...
    //     for ( ;; ) {
    //         switch ( BeginOpeningTag() ) {
    //         case '?':
    //             SkipQDecl();
    //             break;
    //         case '!':
    //             {
    //                 m_Input.SkipChar();
    //                 CTempString tagName = ReadName(m_Input.PeekChar());
    //                 if ( tagName == "DOCTYPE" ) {
    // //                    m_Doctype_found = true;
    //                     ReadName(SkipWS());
    //                     // skip the rest of !DOCTYPE
    //                     for ( ;; ) {
    //                         char c = SkipWS();
    //                         if ( c == '>' ) {
    //                             m_Input.SkipChar();
    //                             Found_gt();
    //                             break;
    //                         }
    //                         else if ( c == '"' || c == '\'' ) {
    //                             SkipAttributeValue(c);
    //                         }
    //                         else {
    //                             ReadName(c);
    //                         }
    //                     }
    //                 }
    //                 else {
    //                     // unknown tag
    //                     ThrowError(fFormatError,
    //                         "unknown tag in file header: "+string(tagName));
    //                 }
    //             }
    //             break;
    //         default:
    //             {
    //                 string typeName = ReadName(m_Input.PeekChar());
    //                 if (!m_Doctype_found && !StackIsEmpty()) {
    //                     // verify typename
    //                     ...
    //                         if ( !typeName.empty() && !tname.empty() && typeName != tname ) {
    //                             string tmp = m_CurrNsPrefix + ":" + typeName;
    //                             ...
    //                 }
    //                 UndoClassMember();
    //                 return typeName;
    //             }
    //         }
    //     }
    //     return NcbiEmptyString;
    // }
    // ```
    // (The namespace prefix of the type name is refused by `read_name`; without one,
    // `m_CurrNsPrefix` is empty and `":" + typeName` never equals a type's name.)
    fn read_file_header(&mut self) -> Result<Vec<u8>> {
        if self.peek(0)? == 0xEF && self.peek(1)? == 0xBB && self.peek(2)? == 0xBF {
            self.skip(3);
            self.encoding = Encoding::Utf8;
        }
        loop {
            match self.begin_opening_tag()? {
                b'?' => self.skip_qdecl()?,
                b'!' => {
                    self.skip(1);
                    let c = self.peek(0)?;
                    self.read_name(c)?;
                    if self.last_tag != b"DOCTYPE" {
                        return fail("unknown tag in file header");
                    }
                    let c = self.skip_ws()?;
                    self.read_name(c)?;
                    loop {
                        let c = self.skip_ws()?;
                        if c == b'>' {
                            self.skip(1);
                            self.state = TagState::Outside;
                            break;
                        } else if c == b'"' || c == b'\'' {
                            self.skip_attribute_value(c)?;
                        } else {
                            self.read_name(c)?;
                        }
                    }
                }
                _ => {
                    let c = self.peek(0)?;
                    self.read_name(c)?;
                    let name = self.last_tag.clone();
                    self.undo_class_member();
                    return Ok(name);
                }
            }
        }
    }

    /// The start of a document whose root is `root` (`SkipFileHeader`): a new object stream's
    /// state, the file header, and the check of the root's name.
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:883-900
    // ```c++
    // void CObjectIStream::SkipFileHeader(TTypeInfo typeInfo)
    // {
    //     ResetState();
    //     ...
    //     string name = ReadFileHeader();
    //     const string& tname = typeInfo->GetName();
    //     if ( !name.empty() && !tname.empty() && name != tname ) {
    //         ThrowError(fFormatError,
    //                    "incompatible type "+name+"<>"+typeInfo->GetName());
    //     }
    // ```
    // (rpsbproc opens one object stream per document, `ObjStreamIn`: its encoding starts
    // unknown, as `m_Encoding` of a new `CObjectIStreamXml`.)
    pub(crate) fn begin_document(&mut self, root: &str) -> Result<()> {
        self.state = TagState::Outside;
        self.rejected = None;
        self.last_tag.clear();
        self.encoding = Encoding::Unknown;
        self.utf8.clear();
        self.utf8_at = 0;
        let name = self.read_file_header()?;
        if name != root.as_bytes() {
            return fail("incompatible type");
        }
        Ok(())
    }

    /// `OpenTag(e)`: the opening tag `e` (its name given back by `UndoClassMember`, or read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:525-535,1563-1573
    // ```c++
    // CTempString CObjectIStreamXml::RejectedName(void)
    // {
    //     _ASSERT(!m_RejectedTag.empty());
    //     m_LastTag = m_RejectedTag;
    //     m_RejectedTag.erase();
    //     m_TagState = eTagInsideOpening;
    //     ...
    //     return m_LastTag;
    // }
    // ...
    // void CObjectIStreamXml::OpenTag(const string& e)
    // {
    //     CTempString tagName;
    //     if (m_RejectedTag.empty()) {
    //         tagName = ReadName(BeginOpeningTag());
    //     } else {
    //         tagName = RejectedName();
    //     }
    //     if ( tagName != e )
    //         ThrowError(fFormatError, "tag '"+e+"' expected: "+string(tagName));
    // }
    // ```
    fn open_tag(&mut self, e: &str) -> Result<()> {
        match self.rejected.take() {
            Some(name) => {
                self.last_tag = name;
                self.state = TagState::InsideOpening;
            }
            None => {
                let c = self.begin_opening_tag()?;
                self.read_name(c)?;
            }
        }
        if self.last_tag != e.as_bytes() {
            return fail("unexpected tag");
        }
        Ok(())
    }

    /// The closing tag of a class (`CloseTag(e)`) or of a member (`CloseStackTag(0)`, whose name
    /// is the class's name, `_` and the member's): nothing after a self-closed opening tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1575-1586,1604-1623
    // ```c++
    // void CObjectIStreamXml::CloseTag(const string& e)
    // {
    //     if ( SelfClosedTag() ) {
    //         EndSelfClosedTag();
    //     }
    //     else {
    //         CTempString tagName = ReadName(BeginClosingTag());
    //         if ( tagName != e )
    //             ThrowError(fFormatError, "tag '"+e+"' expected: "+string(tagName));
    //         EndClosingTag();
    //     }
    // }
    // ...
    // void CObjectIStreamXml::CloseStackTag(size_t level)
    // {
    //     if ( SelfClosedTag() ) {
    //         EndSelfClosedTag();
    //     }
    //     else {
    //         if (m_Attlist) {
    //             m_TagState = eTagInsideClosing;
    //         } else {
    //             CTempString tagName = ReadName(BeginClosingTag());
    //             if (!x_IsStdXml()) {
    //                 CTempString rest = SkipStackTagName(tagName, level);
    //                 if ( !rest.empty() )
    //                     ThrowError(fFormatError,
    //                         "unexpected tag: "+string(tagName)+string(rest));
    //             }
    //         }
    //         EndClosingTag();
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1508-1561
    // ```c++
    // CTempString CObjectIStreamXml::SkipTagName(CTempString tag,
    //                                             const char* str, size_t length)
    // {
    //     if ( tag.size() < length ||
    //          memcmp(tag.data(), str, length) != 0 )
    //         ThrowError(fFormatError, "invalid tag name: "+string(tag));
    //     return CTempString(tag.data() + length, tag.size() - length);
    // }
    //
    // CTempString CObjectIStreamXml::SkipStackTagName(CTempString tag,
    //                                                  size_t level)
    // {
    //     const TFrame& frame = FetchFrameFromTop(level);
    //     switch ( frame.GetFrameType() ) {
    //     ...
    //     case TFrame::eFrameClass:
    //     ...
    //         {
    //             const string& name = frame.GetTypeInfo()->GetName();
    //             if ( !name.empty() )
    //                 return SkipTagName(tag, name);
    //     ...
    //     case TFrame::eFrameClassMember:
    //     case TFrame::eFrameChoiceVariant:
    //         {
    //             tag = SkipStackTagName(tag, level + 1, '_');
    //             return SkipTagName(tag, frame.GetMemberId().GetName());
    //         }
    // ...
    // CTempString CObjectIStreamXml::SkipStackTagName(CTempString tag,
    //                                                  size_t level, char c)
    // {
    //     tag = SkipStackTagName(tag, level);
    //     if ( tag.empty() || tag[0] != c )
    //         ThrowError(fFormatError, "invalid tag name: "+string(tag));
    //     return CTempString(tag.data() + 1, tag.size() - 1);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrxml.inl:104-124
    // ```c++
    // inline
    // void CObjectIStreamXml::EndSelfClosedTag(void)
    // {
    //     _ASSERT(SelfClosedTag());
    //     m_TagState = eTagOutside;
    //     m_LeadingWs = 0;
    // }
    // ...
    // inline
    // void CObjectIStreamXml::EndClosingTag(void)
    // {
    //     _ASSERT(InsideClosingTag());
    //     EndTag();
    // }
    // ```
    fn close_tag(&mut self, class: &str, member: Option<&str>) -> Result<()> {
        if self.state == TagState::SelfClosed {
            self.state = TagState::Outside;
            return Ok(());
        }
        let c = self.begin_closing_tag()?;
        self.read_name(c)?;
        if !is_tag(&self.last_tag, class, member) {
            return fail("unexpected tag");
        }
        self.end_tag()
    }

    /// `ThisTagIsSelfClosed`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1672-1678
    // ```c++
    // bool CObjectIStreamXml::ThisTagIsSelfClosed(void)
    // {
    //     if (InsideOpeningTag()) {
    //         return EndOpeningTagSelfClosed();
    //     }
    //     return SelfClosedTag();
    // }
    // ```
    fn this_tag_is_self_closed(&mut self) -> Result<bool> {
        if self.state == TagState::InsideOpening {
            return self.end_opening_tag_self_closed();
        }
        Ok(self.state == TagState::SelfClosed)
    }

    /// `NextTagIsClosing` and `NextIsTag`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1658-1670
    // ```c++
    // bool CObjectIStreamXml::NextIsTag(void)
    // {
    //     BeginData();
    //     return SkipWSAndComments() == '<' &&
    //         m_Input.PeekChar(1) != '/' &&
    //         m_Input.PeekChar(1) != '!';
    // }
    //
    // bool CObjectIStreamXml::NextTagIsClosing(void)
    // {
    //     BeginData();
    //     return SkipWSAndComments() == '<' && m_Input.PeekChar(1) == '/';
    // }
    // ```
    fn next_tag_is_closing(&mut self) -> Result<bool> {
        self.begin_data()?;
        Ok(self.skip_ws_and_comments()? == b'<' && self.peek(1)? == b'/')
    }

    fn next_is_tag(&mut self) -> Result<bool> {
        self.begin_data()?;
        Ok(self.skip_ws_and_comments()? == b'<' && self.peek(1)? != b'/' && self.peek(1)? != b'!')
    }

    /// `HasAttlist` of an opening tag: whether something other than its end follows.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1650-1656
    // ```c++
    // bool CObjectIStreamXml::HasAttlist(void)
    // {
    //     if (InsideTag()) {
    //         return !IsEndOfTagChar( SkipWS() );
    //     }
    //     return false;
    // }
    // ```
    fn has_attlist(&mut self) -> Result<bool> {
        if matches!(
            self.state,
            TagState::InsideOpening | TagState::InsideClosing
        ) {
            return Ok(!is_end_of_tag_char(self.skip_ws()?));
        }
        Ok(false)
    }

    /// `BeginClassMember(classType, pos)` of a SEQUENCE of an ASN.1 module (not standard XML,
    /// no member without a tag, no attribute list, unknown members not skipped): the index of
    /// the next member, `None` at the end of the members.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2185-2367
    // ```c++
    // TMemberIndex
    // CObjectIStreamXml::BeginClassMember(const CClassTypeInfo* classType,
    //                                     TMemberIndex pos)
    // {
    //     CTempString tagName;
    //     TMemberIndex first = classType->GetMembers().FirstIndex();
    //     if (m_RejectedTag.empty()) {
    //         if (m_Attlist && InsideTag()) {
    //             ...
    //         } else {
    //             if (!m_Attlist) {
    //                 if (pos == first) {
    //                     if (classType->GetMemberInfo(first)->GetId().IsAttlist()) {
    //                         ...
    //                     }
    // // if class spec defines no attributes, but there are some - skip them
    //                     if (HasAttlist()) {
    //                         ReadUndefinedAttributes();
    //                     }
    //                 }
    //             }
    //             ...
    //             if ( SelfClosedTag() || ThisTagIsSelfClosed()) {
    //                 m_Attlist = false;
    //                 TMemberIndex last = classType->GetMembers().LastIndex();
    //                 if (pos == last) {
    //                     if (classType->GetMemberInfo(pos)->GetId().HasNotag() &&
    //                         ...
    //                 }
    //                 return kInvalidMember;
    //             }
    //             if ( NextTagIsClosing() )
    //                 return kInvalidMember;
    //             if (pos <= classType->GetItems().LastIndex()) {
    //                 const CMemberInfo* mem_info = classType->GetMemberInfo(pos);
    //                 if (mem_info->GetId().HasNotag() &&
    //                     ...
    //             } else {
    //                 if (CanSkipUnknownMembers()) {
    //                     ...
    //                 }
    //                 return kInvalidMember;
    //             }
    //             if (!NextIsTag()) {
    //                 return kInvalidMember;
    //             }
    //             tagName = ReadName(BeginOpeningTag());
    //         }
    //     } else {
    //         tagName = RejectedName();
    //     }
    //
    //     TMemberIndex ind = classType->GetMembers().Find(tagName);
    //     if (ind == kInvalidMember) {
    //         ind = classType->GetMembers().FindDeep(tagName, pos);
    //         if (ind != kInvalidMember && ind >= pos) {
    //             TopFrame().SetNotag();
    //             UndoClassMember();
    //             return ind;
    //         }
    //     } else {
    //         ...
    //     }
    //     ...
    //     CTempString id = SkipStackTagName(tagName, 1, '_');
    //     TMemberIndex index = classType->GetMembers().Find(id, pos);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownMembers()) {
    //             ...
    //         } else {
    //            UnexpectedMember(id, classType->GetMembers());
    //         }
    //     }
    //     return index;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:294-313,432-439
    // ```c++
    // TMemberIndex CItemsInfo::FindDeep(const CTempString& name, TMemberIndex pos) const
    // {
    //     TMemberIndex ind = Find(name, pos);
    //     if (ind != kInvalidMember) {
    //         return ind;
    //     }
    //     for (CIterator item(*this, pos); item.Valid(); ++item) {
    //         const CItemInfo* info = GetItemInfo(item);
    //         const CClassTypeInfoBase* classType =
    //             dynamic_cast<const CClassTypeInfoBase*>(
    //                 FindRealTypeInfo(info->GetTypeInfo()));
    //         if (classType) {
    //             if (classType->GetItems().FindDeep(name)
    //                 != kInvalidMember) {
    //                 return *item;
    //             }
    //         }
    //     }
    //     return kInvalidMember;
    // }
    // ...
    // TMemberIndex CItemsInfo::Find(const CTempString& name, TMemberIndex pos) const
    // {
    //     for ( CIterator i(*this, pos); i.Valid(); ++i ) {
    //         if ( name == GetItemInfo(i)->GetId().GetName() )
    //             return *i;
    //     }
    //     return kInvalidMember;
    // }
    // ```
    // (Unknown members are not skipped: `SERIAL_SKIP_UNKNOWN_MEMBERS` is a registry parameter,
    // which LOSAT rejects. A tag that is a member's name without the class prefix, or that
    // `FindDeep` finds among the members of a member's class, ends in an error: the first
    // goes to `SkipStackTagName`, which needs the prefix; the second is read as that member
    // without its tag, whose class then needs its own name as the tag (`OpenTag`), and no member
    // is named like a class. The tag of a member is therefore `<Class_member>` of a member at
    // `pos` or after, or the document fails.)
    fn begin_class_member(
        &mut self,
        class: &str,
        members: &[Member],
        pos: usize,
    ) -> Result<Option<usize>> {
        if pos == 0 && self.has_attlist()? {
            self.read_undefined_attributes()?;
        }
        if self.state == TagState::SelfClosed || self.this_tag_is_self_closed()? {
            return Ok(None);
        }
        if self.next_tag_is_closing()? || pos >= members.len() || !self.next_is_tag()? {
            return Ok(None);
        }
        let c = self.begin_opening_tag()?;
        self.read_name(c)?;
        let Some(id) = self
            .last_tag
            .strip_prefix(class.as_bytes())
            .and_then(|rest| rest.strip_prefix(b"_"))
        else {
            return fail("invalid tag name");
        };
        match members[pos..]
            .iter()
            .position(|member| member.name.as_bytes() == id)
        {
            Some(offset) => Ok(Some(pos + offset)),
            None => fail("unexpected member"),
        }
    }

    /// A SEQUENCE of `class` (`ReadClassSequential`): its opening tag, each member that the
    /// text holds in the module's order, read by `read` with the member's index (its opening
    /// tag read, its value not), and its closing tag; a mandatory member that the text lacks
    /// fails the document.
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1398-1428
    // ```c++
    // void CObjectIStream::ReadClassSequential(const CClassTypeInfo* classType,
    //                                          TObjectPtr classPtr)
    // {
    //     TMemberIndex prevIndex = kInvalidMember;
    //     BEGIN_OBJECT_FRAME3(eFrameClass, classType, classPtr);
    //     BeginClass(classType);
    //
    //     ReadClassSequentialContentsBegin(classType);
    //
    //     TMemberIndex index;
    //     while ( (index = BeginClassMember(classType, *pos)) != kInvalidMember ) {
    //
    //         if ((prevIndex != kInvalidMember) && (prevIndex >= index)) {
    //             const CMemberInfo *mem_info = classType->GetMemberInfo(index);
    //             if (mem_info->GetId().HaveNoPrefix()) {
    //                 UndoClassMember();
    //                 break;
    //             }
    //         }
    //         prevIndex = index;
    //
    //         ReadClassSequentialContentsMember(classPtr);
    //
    //         EndClassMember();
    //     }
    //
    //     ReadClassSequentialContentsEnd(classPtr);
    //
    //     EndClass();
    //     END_OBJECT_FRAME();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:85-111
    // ```c++
    // #define ClassSequentialContentsMember(Func, Args) \
    //     { \
    //         const CMemberInfo* memberInfo = classType->GetMemberInfo(index); \
    //         SetTopMemberId(memberInfo->GetId()); \
    //         for ( TMemberIndex i = *pos; i < index; ++i ) { \
    //             classType->GetMemberInfo(i)->NCBI_NAME2(Func,MissingMember)Args; \
    //         } \
    //         { \
    //             memberInfo->NCBI_NAME2(Func,Member)Args; \
    //         } \
    //         pos.SetIndex(index + 1); \
    //     }
    //
    // #define ClassSequentialContentsEnd(Func, Args) \
    //     } \
    //     END_OBJECT_FRAME(); \
    //     for ( ; pos.Valid(); ++pos ) { \
    //         classType->GetMemberInfo(pos)->NCBI_NAME2(Func,MissingMember)Args; \
    //     } \
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
    // ```c++
    // bool CObjectIStream::ExpectedMember(const CMemberInfo* memberInfo)
    // {
    //     const CItemInfo* info = CItemsInfo::FindNextMandatory(memberInfo);
    //     if (info) {
    //         if (GetVerifyData() == eSerialVerifyData_Yes) {
    //             ThrowError(fFormatError,
    //                     "member "+info->GetId().ToString()+" expected");
    //         } else {
    //             ...
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2069-2103,2369-2377
    // ```c++
    // void CObjectIStreamXml::BeginClass(const CClassTypeInfo* classInfo)
    // {
    //     CheckStdXml(classInfo);
    //     if (m_SkipNextTag) {
    //         ...
    //     }
    //     if (x_IsStdXml()) {
    //         ...
    //     } else {
    //         OpenTagIfNamed(classInfo);
    //     }
    // }
    //
    // void CObjectIStreamXml::EndClass(void)
    // {
    //     if (TopFrame().GetNotag()) {
    //         TopFrame().SetNotag(false);
    //     } else {
    //         CloseTagIfNamed(TopFrame().GetTypeInfo());
    //     }
    //     x_EndTypeNamespace();
    // }
    // ...
    // void CObjectIStreamXml::EndClassMember(void)
    // {
    //     m_SkipNextTag = false;
    //     if (TopFrame().GetNotag()) {
    //         TopFrame().SetNotag(false);
    //     } else {
    //         CloseStackTag(0);
    //     }
    // }
    // ```
    // (The verification of the read is `eSerialVerifyData_Yes` unless `SERIAL_VERIFY_DATA_READ`
    // turns it off, which LOSAT's rpsbproc rejects; `FindNextMandatory` of a member with a tag is
    // the member when it is not OPTIONAL. Members come in the module's order only
    // (`begin_class_member`), so `prevIndex >= index` never holds.)
    pub(crate) fn class(
        &mut self,
        class: &str,
        members: &[Member],
        mut read: impl FnMut(&mut Self, usize) -> Result<()>,
    ) -> Result<()> {
        self.open_tag(class)?;
        let mut pos = 0;
        while let Some(index) = self.begin_class_member(class, members, pos)? {
            if members[pos..index].iter().any(|member| !member.optional) {
                return fail("member expected");
            }
            read(self, index)?;
            pos = index + 1;
            self.close_tag(class, Some(members[index].name))?;
        }
        if members[pos..].iter().any(|member| !member.optional) {
            return fail("member expected");
        }
        self.close_tag(class, None)
    }

    /// A SEQUENCE OF a named type as a member's value (`ReadContainerContents` of a container
    /// without a name): `element` is called for each element while the member has more (it
    /// reads one with `class`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1765-1776,1921-1979
    // ```c++
    // bool CObjectIStreamXml::HasMoreElements(TTypeInfo elementType)
    // {
    //     bool no_more=false;
    //     try {
    //         no_more = ThisTagIsSelfClosed() || NextTagIsClosing();
    //     } catch (CEofException&) {
    //         no_more = true;
    //     }
    //     if (no_more) {
    //         m_LastPrimitive.erase();
    //         return false;
    //     }
    //     if (x_IsStdXml()) {
    //         ...
    //     }
    //     return true;
    // }
    // ...
    // void CObjectIStreamXml::ReadContainerContents(const CContainerTypeInfo* cType,
    //                                               TObjectPtr containerPtr)
    // {
    //     int count = 0;
    //     TTypeInfo elementType = cType->GetElementType();
    //     if ( !WillHaveName(elementType) ) {
    //         ...
    //     }
    //     else {
    //         CContainerTypeInfo::CIterator iter;
    //         bool old_element = cType->InitIterator(iter, containerPtr);
    //         while ( HasMoreElements(elementType) ) {
    //             if ( old_element ) {
    //                 ...
    //             }
    //             else {
    //                 cType->AddElement(containerPtr, *this);
    //             }
    //             ++count;
    //         }
    //         ...
    //     }
    //     ...
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1860-1871
    // ```c++
    // void CObjectIStreamXml::ReadContainer(const CContainerTypeInfo* containerType,
    //                                       TObjectPtr containerPtr)
    // {
    //     if ( m_StdXml || containerType->GetName().empty() ) {
    //         ReadContainerContents(containerType, containerPtr);
    //     }
    // ```
    // (The end of the input in `HasMoreElements` ends the elements; the member's closing tag
    // then meets it again: the failure of the document.)
    pub(crate) fn container(
        &mut self,
        mut element: impl FnMut(&mut Self) -> Result<()>,
    ) -> Result<()> {
        while !(self.this_tag_is_self_closed()? || self.next_tag_is_closing()?) {
            element(self)?;
        }
        Ok(())
    }
}

/// Whether `name` is the tag of `class` or of its `member` (`Class_member`).
fn is_tag(name: &[u8], class: &str, member: Option<&str>) -> bool {
    match member {
        None => name == class.as_bytes(),
        Some(member) => {
            name.len() == class.len() + 1 + member.len()
                && name.starts_with(class.as_bytes())
                && name[class.len()] == b'_'
                && name.ends_with(member.as_bytes())
        }
    }
}
