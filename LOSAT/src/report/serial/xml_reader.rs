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
//!
//! The bytes come from an `AsnSource` (the stream of `asn_text_reader`, which keeps only the
//! bytes not yet read): one reader reads documents one after another from one stream, each from
//! the byte after the previous one's last `>` (`CIStreamBuffer::Close` gives the bytes it did
//! not use back to the stream, cited at `StreamSource`).
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

use super::asn_text_reader::{buffer_number, AsnSource, NumberStop};
use crate::ncbi_str::{is_space, string_to_double_posix};

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

/// The message of the end of the input (`CEofException`, `CIStreamBuffer::FillBuffer`): the
/// object stream (`asn_types::read_xml`) tells it from the other failures by it.
pub(super) const EOF_MESSAGE: &str = "end of the input";

fn unsupported(what: &str) -> XmlError {
    XmlError::Unsupported(format!(
        "{what} in an XML document is not supported by LOSAT's rpsbproc"
    ))
}

/// `m_TagState`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TagState {
    Outside,
    InsideOpening,
    InsideClosing,
    SelfClosed,
}

/// `EEncoding` of a document (`m_Encoding`): `Unknown` until a declaration names one.
// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:198-209
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
const CP1252: [u32; 32] = [
    0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160, 0x2039,
    0x0152, 0x003F, 0x017D, 0x003F, 0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022, 0x2013, 0x2014,
    0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178,
];

/// `CUtf8::x_AppendChar`: the UTF-8 bytes of a code point (any 32-bit value: the bits above 21
/// are cut by the byte conversion).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7040-7060
// (`CUtf8::AsUTF8` of one `TUnicodeSymbol` (`chU`, length 1), ncbistr.hpp:5718-5758, appends it
// with this; the surrogate pair test of `x_TCharToUnicodeSymbol` is for 2-byte characters only.)
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

/// `CUtf8::AsUTF8` of a `CTempString` holding `ch` (length 1) and `enc_in`, for one byte of 0x80
/// and above in a document of a single-byte encoding: the byte itself for ASCII (`x_Append` copies
/// ASCII and UTF-8), its code point in UTF-8 for Latin-1 and cp1252 (`CharToSymbol`).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7079
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6859-6881
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
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1886,302
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
    pub(super) src: S,
    /// The position of the next byte to read (`m_CurrentPos`).
    pub(super) pos: usize,
    pub(super) state: TagState,
    /// `m_RejectedTag`: a tag name read and given back (`UndoClassMember`).
    pub(super) rejected: Option<Vec<u8>>,
    /// `m_LastTag`.
    pub(super) last_tag: Vec<u8>,
    /// `m_Encoding`.
    pub(super) encoding: Encoding,
    /// `m_Utf8Buf` and `m_Utf8Pos`: the UTF-8 bytes of a character being returned one by one.
    utf8: Vec<u8>,
    utf8_at: usize,
    /// The object stream of a type table (`asn_types::read_xml`, `CObjectIStream::Open` with
    /// `eSerial_Xml` for blast_formatter's archives) rather than rpsbproc's reader of a BLAST XML document:
    /// attributes are skipped as NCBI skips the ones it does not know, and a control character
    /// is replaced by `#` with a post (`x_VerifyChar`), where rpsbproc's refuses both.
    object: bool,
    /// `m_Input.GetLine()`: the line counter of NCBI's input buffer, which only
    /// `SkipEndOfLine` advances.
    pub(super) line: usize,
    /// The bytes replaced by `#` (`ReplaceVisibleChar` of `x_VerifyChar`) and the line of each,
    /// for the object stream to post.
    pub(super) replaced: Vec<(u8, usize)>,
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
            object: false,
            line: 1,
            replaced: Vec::new(),
        }
    }

    /// The reader of an object stream (see `object`).
    pub(super) fn object_stream(src: S) -> Self {
        XmlReader {
            object: true,
            ..Self::new(src)
        }
    }

    /// `UngetChar` of `count` bytes just read.
    pub(super) fn skip_back(&mut self, count: usize) {
        self.pos -= count;
    }

    /// `ResetState` of an object stream at its root frame: no tag open, none given back.
    pub(super) fn reset_state(&mut self) {
        self.state = TagState::Outside;
        self.rejected = None;
        self.last_tag.clear();
        self.utf8.clear();
        self.utf8_at = 0;
    }

    /// `PeekChar` at `offset`: the byte, or the failure of the end of the input (`CEofException`).
    #[inline]
    pub(super) fn peek(&mut self, offset: usize) -> Result<u8> {
        match self.src.byte_at(self.pos + offset) {
            Some(c) => Ok(c),
            None => fail(EOF_MESSAGE),
        }
    }

    /// `PeekCharNoEOF` at `offset`: the byte, 0 past the end of the input.
    #[inline]
    pub(super) fn peek_no_eof(&mut self, offset: usize) -> u8 {
        self.src.byte_at(self.pos + offset).unwrap_or(0)
    }

    #[inline]
    pub(super) fn skip(&mut self, count: usize) {
        self.pos += count;
    }

    /// `GetChar()`.
    #[inline]
    pub(super) fn get_char(&mut self) -> Result<u8> {
        let c = self.peek(0)?;
        self.pos += 1;
        Ok(c)
    }

    /// `CIStreamBuffer::SkipSpaces`: skips `' '` and gives the next byte (not read).
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:213-252
    fn skip_spaces(&mut self) -> Result<u8> {
        loop {
            let c = self.peek(0)?;
            if c != b' ' {
                return Ok(c);
            }
            self.pos += 1;
        }
    }

    /// `CIStreamBuffer::SkipEndOfLine` of `lastChar`: a `\r\n` or `\n\r` pair is one end of line.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-537
    pub(super) fn skip_end_of_line(&mut self, last: u8) {
        self.line += 1;
        let next = self.peek_no_eof(0);
        if i32::from(last) + i32::from(next as i8) == i32::from(b'\r' + b'\n') {
            self.pos += 1;
        }
    }

    /// `FindChar` of `c`: moves to the next `c`.
    fn find_char(&mut self, c: u8) -> Result<()> {
        while self.peek(0)? != c {
            self.pos += 1;
        }
        Ok(())
    }

    /// `PeekFindChar` of `c` and `limit`: the offset of `c` among the next `limit` bytes (fewer at
    /// the end of the input), `limit` when it is not there.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:292-306
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
    pub(super) fn skip_ws(&mut self) -> Result<u8> {
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
    // (`m_LeadingWs`, the spaces counted here, is given to the next `ReadTagData`; every tag
    // read resets it (`Found_lt`, `Found_gt`), and the character data of a member is read just
    // after its opening tag, so it is always 0 there: LOSAT does not keep it.)
    pub(super) fn skip_ws_and_comments(&mut self) -> Result<u8> {
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
    pub(super) fn end_tag(&mut self) -> Result<()> {
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
    // (`ReadName` throws for a byte that cannot start a name.)
    // (An object stream skips each attribute: its name, `=` and the quoted value. `xmlns` and `nil`
    // would also set namespace maps and `m_IsNil`, which a type of an ASN.1 module never
    // reads; the `nil` value's own check (`NStr::StringToBool`) is not made.)
    pub(super) fn read_undefined_attributes(&mut self) -> Result<u8> {
        if self.object {
            loop {
                let c = self.skip_ws()?;
                if is_end_of_tag_char(c) {
                    return Ok(c);
                }
                self.read_name(c)?;
                self.read_attribute_value()?;
            }
        }
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
    // (No frame of these types is read without its tag (`GetNotag`).)
    pub(super) fn end_opening_tag_self_closed(&mut self) -> Result<bool> {
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
    pub(super) fn begin_data(&mut self) -> Result<()> {
        if self.state == TagState::InsideOpening {
            self.end_tag()?;
        }
        Ok(())
    }

    /// `BeginOpeningTag`: reads `<` and gives the byte after it (not read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:405-418
    pub(super) fn begin_opening_tag(&mut self) -> Result<u8> {
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
    pub(super) fn begin_closing_tag(&mut self) -> Result<u8> {
        self.begin_data()?;
        if self.skip_ws_and_comments()? != b'<' || self.peek(1)? != b'/' {
            return fail("'</' expected");
        }
        self.skip(2);
        self.state = TagState::InsideClosing;
        self.peek(0)
    }

    /// `ReadName` of `c`: the name that starts with `c` (not read yet), into `last_tag`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:432-469
    // (`m_Doctype_found` is false for the types of an ASN.1 module (ReadFileHeader:
    // `GetDataSpec()` is `EDataSpec::eDTD`). A `:` after the first byte makes a namespace prefix,
    // refused; the `xmlns` and `nil` names of `m_Attlist` are attributes, refused before.)
    pub(super) fn read_name(&mut self, c: u8) -> Result<()> {
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
    pub(super) fn undo_class_member(&mut self) {
        if self.state == TagState::InsideOpening {
            self.rejected = Some(self.last_tag.clone());
            self.state = TagState::Outside;
        }
    }

    /// `SkipQDecl`: the rest of a processing instruction after `<`; an `encoding` pseudo-attribute
    /// sets the document's encoding.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:545-580
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

    /// `ReadAttributeValue` of `value` for a processing instruction (`m_Attlist` false: the closing
    /// quote is read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:891-912
    pub(super) fn read_attribute_value(&mut self) -> Result<Vec<u8>> {
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

    /// `SkipAttributeValue` of `c`: a quoted string of the `DOCTYPE`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:537-543
    fn skip_attribute_value(&mut self, c: u8) -> Result<()> {
        self.skip(1);
        self.find_char(c)?;
        self.skip(1);
        Ok(())
    }

    /// `ReadEscapedChar` of `endingChar` (and `encoded`, by address): the next character (a byte,
    /// or the value of an entity or character reference, which is `encoded`), or -1 at `ending`
    /// (not read). A character reference of 2^31 and more is negative, as NCBI's `int` of its
    /// `unsigned`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:722-800
    pub(super) fn read_escaped_char(&mut self, ending: u8) -> Result<(i32, bool)> {
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

    /// `ReadEncodedChar` of `endingChar`, `type` and `encoded`, for a VisibleString (or of a UTF-8
    /// attribute value: the output encoding is UTF-8 either way, `m_StringEncoding`): the bytes of
    /// the next character in UTF-8, one per call, or a negative value at `ending`; a control
    /// character is refused (`x_VerifyChar`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:818-868
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
            if self.object {
                // `ReplaceVisibleChar` of `c` with `eFNP_ReplaceAndWarn`: the post at the line now,
                // and `#` (`x_FixCharsSubst`).
                self.replaced.push((c as u8, self.line));
                return Ok(i32::from(b'#'));
            }
            return Err(XmlError::Unsupported(format!(
                "a control character (0x{c:02x}) in the character data of an XML document is not supported by LOSAT's rpsbproc (NCBI replaces it and posts a warning)"
            )));
        }
        Ok(c)
    }

    /// `ReadCDSection` of `str`: a CDATA section where character data meets `<`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1271-1295
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

    /// `ReadTagData` of `str` with `eStringTypeVisible`: the character data up to the next tag,
    /// with NCBI's end-of-line handling (a first `\r` is dropped; from it to the next `\n`, every
    /// `\r` is a `\n`) and NUL characters dropped.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1297-1354
    pub(super) fn read_tag_data(&mut self, out: &mut Vec<u8>) -> Result<()> {
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

    /// `ReadWord` of `str`: the character data up to white space or the next tag, after white
    /// space.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1356-1374
    pub(super) fn read_word(&mut self, out: &mut Vec<u8>) -> Result<()> {
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

    /// `CIStreamBuffer::GetInt4`, `GetUint4` (`signed` false) and `GetInt8` (`max`): `SkipWs`, then
    /// the one port of the buffer's numbers (`asn_text_reader::buffer_number`). The end of the
    /// input in `SkipWs` is the lexer's failure; the number's own stop is given to the caller,
    /// whose stream reports it in its form.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:617-624
    pub(super) fn get_number(
        &mut self,
        max: u64,
        signed: bool,
    ) -> Result<std::result::Result<i64, NumberStop>> {
        let c = loop {
            let c = self.get_char()?;
            if !is_space(c) {
                break c;
            }
        };
        let start = self.pos - 1;
        let number = buffer_number(&mut self.src, &mut self.pos, start, c, max, signed);
        if number == Err(NumberStop::Eof) {
            return fail(EOF_MESSAGE);
        }
        Ok(number)
    }

    /// `get_number` in rpsbproc's reader: a number's stop is a failure of the document.
    fn number(&mut self, max: u64) -> Result<i64> {
        match self.get_number(max, true)? {
            Ok(n) => Ok(n),
            Err(NumberStop::Bad) => fail("bad number"),
            Err(_) => fail("number overflow"),
        }
    }

    /// `ReadInt4`: an INTEGER (`get_number`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:987-997
    // (No member of these types has a default: `ExpectSpecialCase` is 0. For `-2147483648`,
    // the negation of `n` cast to `Int4` negates `INT_MIN`, a signed overflow; LOSAT gives the
    // number's value, `INT_MIN`, which both NCBI builds give too (PS16d oracle `p_intmin`).)
    pub(crate) fn int4(&mut self) -> Result<i32> {
        self.begin_data()?;
        Ok(self.number(i32::MAX as u64)? as i32)
    }

    /// `ReadInt8`: a BigInt.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1005-1015
    // (`-9223372036854775808` is read as its value, as `int4` reads `-2147483648`.)
    pub(crate) fn int8(&mut self) -> Result<i64> {
        self.begin_data()?;
        self.number(i64::MAX as u64)
    }

    /// `ReadDouble`: a REAL.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1027-1042
    pub(crate) fn real(&mut self) -> Result<f64> {
        let mut word = Vec::new();
        self.read_word(&mut word)?;
        let read = string_to_double_posix(&word, true);
        let (value, mut end) = (read.value, read.end);
        while end < word.len() && is_white_space(word[end]) {
            end += 1;
        }
        if end != word.len() {
            return fail("invalid float number");
        }
        Ok(value)
    }

    /// `ReadString` of `str` with `eStringTypeVisible`: a VisibleString (empty for a self-closed
    /// tag).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1233-1258
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
    // (The namespace prefix of the type name is refused by `read_name`; without one,
    // `m_CurrNsPrefix` is empty and a colon followed by `typeName` never equals a type's name.)
    pub(super) fn read_file_header(&mut self) -> Result<Vec<u8>> {
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

    /// `OpenTag` of `e`: the opening tag `e` (its name given back by `UndoClassMember`, or read).
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:525-535,1563-1573
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

    /// The closing tag of a class (`CloseTag` of `e`) or of a member (`CloseStackTag` of 0, whose
    /// name is the class's name, `_` and the member's): nothing after a self-closed opening tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1575-1586,1604-1623
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1508-1561
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrxml.inl:104-124
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
    pub(super) fn this_tag_is_self_closed(&mut self) -> Result<bool> {
        if self.state == TagState::InsideOpening {
            return self.end_opening_tag_self_closed();
        }
        Ok(self.state == TagState::SelfClosed)
    }

    /// `NextTagIsClosing` and `NextIsTag`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1658-1670
    pub(super) fn next_tag_is_closing(&mut self) -> Result<bool> {
        self.begin_data()?;
        Ok(self.skip_ws_and_comments()? == b'<' && self.peek(1)? == b'/')
    }

    pub(super) fn next_is_tag(&mut self) -> Result<bool> {
        self.begin_data()?;
        Ok(self.skip_ws_and_comments()? == b'<' && self.peek(1)? != b'/' && self.peek(1)? != b'!')
    }

    /// `HasAttlist` of an opening tag: whether something other than its end follows.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1650-1656
    pub(super) fn has_attlist(&mut self) -> Result<bool> {
        if matches!(
            self.state,
            TagState::InsideOpening | TagState::InsideClosing
        ) {
            return Ok(!is_end_of_tag_char(self.skip_ws()?));
        }
        Ok(false)
    }

    /// `BeginClassMember` of `classType` at `pos` for a SEQUENCE of an ASN.1 module (not standard
    /// XML, no member without a tag, no attribute list, unknown members not skipped): the index of
    /// the next member, `None` at the end of the members.
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2185-2367
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:294-313,432-439
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
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:85-111
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1045
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:2069-2103,2369-2377
    // (The verification of the read is `eSerialVerifyData_Yes` unless `SERIAL_VERIFY_DATA_READ`
    // turns it off, which LOSAT's rpsbproc rejects; `FindNextMandatory` of a member with a tag is
    // the member when it is not OPTIONAL. Members come in the module's order only
    // (`begin_class_member`), so `prevIndex` is never at or above `index`.)
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
    // NCBI reference (598d8ae6): c++/src/serial/objistrxml.cpp:1860-1871
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
