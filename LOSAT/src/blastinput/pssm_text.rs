//! The reader of a PSSM in ASN.1 text (`PssmWithParameters`, the file of psiblast's `-out_pssm`
//! and of `-in_pssm`): NCBI's ASN.1 text object stream (`CObjectIStreamAsn`) for the types of
//! `report::serial::pssm`, whose writer gives the same bytes back for a PSSM that psiblast wrote.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1209-1216
//! ```c++
//!     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
//!         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
//!         _ASSERT(m_Pssm.Empty());
//!         m_Pssm.Reset(new CPssmWithParameters);
//!         try {
//!             switch (CFormatGuess().Format(in)) {
//!             case CFormatGuess::eBinaryASN:
//!                 in >> MSerial_AsnBinary >> *m_Pssm;
//! ```
//!
//! The values are read as NCBI's stream reads them: white space is spaces, tabs and ends of line
//! (`\n`, `\r`, `\r\n` or `\n\r`, each one line), and a comment runs from `--` to the next `--` or
//! the end of the line; the members of a SEQUENCE come in the order of the module; a REAL is
//! `{ mantissa, base, exponent }` or a word (`ReadDouble`, `read_real`); a VisibleString joins its
//! lines (an end of line inside it is dropped) and `""` is one quote; an OCTET STRING is hex digits
//! across lines. The members and variants that LOSAT does not keep are skipped as values
//! (`skip_value`): the descriptors of the query other than its titles, its `annot`, the members
//! of its `Seq-inst` other than `repr`, `mol`, `length` and `seq-data`, and the `constraints` and
//! `annotatedSites` of the parameters; `read_pssm_with_parameters_and_skips` names them, since NCBI
//! writes them back where it writes the PSSM it read (psiblast `-in_pssm` with `-out_pssm`, item
//! PS5b). An empty `descr` is NCBI's error (`CSeq_descr::PostRead`). A string keeps its bytes as
//! read (NCBI replaces a byte outside the visible range by `#` with an error message; a PSSM that
//! LOSAT writes has none, `algorithm::psiblast::pssm`). LOSAT reads the forms below and refuses the others with an error
//! naming its line: a query `Seq-entry` that is a `set`, a `Seq-id` other than `local`, `gi`,
//! `general` and the `Textseq-id` ones, `Seq-data` other than `iupacaa`, `ncbieaa` and
//! `ncbistdaa`, and a REAL mantissa of more than 18 digits.

use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmSeqId, PssmWithParameters, SeqData, SeqInst, TextseqId, SEQ_INST_MOL,
    SEQ_INST_REPR, TEXTSEQ_VARIANTS,
};
use crate::utils::x87::X87;

/// A PSSM that the reader cannot read: the line (from 1) where it stopped and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PssmTextError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for PssmTextError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for PssmTextError {}

type Result<T> = std::result::Result<T, PssmTextError>;

/// Reads the `PssmWithParameters` document at the start of `text` (what follows it is not read,
/// as `stream >> object` reads one object).
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:709-721
// ```c++
// string CObjectIStreamAsn::ReadFileHeader()
// {
//     CTempString id = ReadTypeId(SkipWhiteSpace());
//     string s(id);
//     if ( SkipWhiteSpace() == ':' &&
//          m_Input.PeekCharNoEOF(1) == ':' &&
//          m_Input.PeekCharNoEOF(2) == '=' ) {
//         m_Input.SkipChars(3);
//     }
//     else
//         ThrowError(fFormatError, "'::=' expected");
//     return s;
// }
// ```
pub fn read_pssm_with_parameters(text: &[u8]) -> Result<PssmWithParameters> {
    read_pssm_with_parameters_and_skips(text).map(|(pssm, _)| pssm)
}

/// `read_pssm_with_parameters`, with the names of the members that the reader skipped (in the
/// order read, `"Seqdesc"` for a descriptor of the query other than a title), which NCBI's object
/// keeps and writes back.
pub fn read_pssm_with_parameters_and_skips(
    text: &[u8],
) -> Result<(PssmWithParameters, Vec<&'static str>)> {
    let mut reader = Reader {
        text,
        pos: 0,
        line: 1,
        skipped: Vec::new(),
    };
    let c = reader.skip_white_space()?;
    let id = reader.scan_end_of_id(is_first_id_char(c));
    if id != b"PssmWithParameters" {
        return Err(reader.error(format!(
            "the document is a {}, not a PssmWithParameters",
            String::from_utf8_lossy(id)
        )));
    }
    if reader.skip_white_space()? != b':' || reader.peek_at(1) != b':' || reader.peek_at(2) != b'='
    {
        return Err(reader.error("'::=' expected"));
    }
    reader.pos += 3;
    let pssm = reader.pssm_with_parameters()?;
    Ok((pssm, reader.skipped))
}

/// `FirstIdChar` (objistrasn.cpp:118-121).
fn is_first_id_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:123-127
// ```c++
// bool CObjectIStreamAsn::IdChar(char c)
// {
//     return isalnum((unsigned char) c) || c == '_' || c == '.';
// }
// ```
fn is_id_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// The C locale's `isspace`.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
}

struct Reader<'a> {
    text: &'a [u8],
    pos: usize,
    /// `CIStreamBuffer::m_Line`.
    line: usize,
    /// The members skipped as values (`skip_value`), by name.
    skipped: Vec<&'static str>,
}

/// The members of a SEQUENCE as they are read: each name in turn, in the order of the module.
struct Members {
    block_start: bool,
    /// The position in the module of the last member read.
    last: Option<usize>,
}

impl<'a> Reader<'a> {
    fn error(&self, message: impl Into<String>) -> PssmTextError {
        PssmTextError {
            line: self.line,
            message: message.into(),
        }
    }

    /// The byte at `pos + offset`, 0 past the end (`PeekCharNoEOF`).
    fn peek_at(&self, offset: usize) -> u8 {
        self.text.get(self.pos + offset).copied().unwrap_or(0)
    }

    /// `SkipEndOfLine(c)` after the end-of-line byte `c` was taken.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-537
    // ```c++
    //     m_Line++;
    //     char nextChar = PeekCharNoEOF();
    //     // lastChar either '\r' or \n'
    //     // if nextChar is compliment, skip it
    //     if ( (lastChar + nextChar) == ('\r' + '\n') )
    //         SkipChar();
    // ```
    fn skip_end_of_line(&mut self, c: u8) {
        self.line += 1;
        if c as u32 + self.peek_at(0) as u32 == b'\r' as u32 + b'\n' as u32 {
            self.pos += 1;
        }
    }

    /// `SkipWhiteSpace`: the next byte that is not white space or in a comment, not taken. The
    /// end of the text is an error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243
    // ```c++
    //         for ( ;; ) {
    //             char c = m_Input.SkipSpaces();
    //             switch ( c ) {
    //             case '\t':
    //                 m_Input.SkipChar();
    //                 continue;
    //             case '\r':
    //             case '\n':
    //                 m_Input.SkipChar();
    //                 SkipEndOfLine(c);
    //                 continue;
    //             case '-':
    //                 // check for comments
    //                 if ( m_Input.PeekChar(1) != '-' ) {
    //                     return '-';
    //                 }
    //                 m_Input.SkipChars(2);
    //                 // skip comments
    //                 SkipComments();
    //                 continue;
    //             default:
    //                 return c;
    //             }
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:245-274
    // ```c++
    //         for ( ;; ) {
    //             char c = GetChar();
    //             switch ( c ) {
    //             case '\r':
    //             case '\n':
    //                 SkipEndOfLine(c);
    //                 return;
    //             case '-':
    //                 c = GetChar();
    //                 switch ( c ) {
    //                 case '\r':
    //                 case '\n':
    //                     SkipEndOfLine(c);
    //                     return;
    //                 case '-':
    //                     return;
    //                 }
    //                 continue;
    //             default:
    //                 continue;
    //             }
    //         }
    // ```
    fn skip_white_space(&mut self) -> Result<u8> {
        loop {
            let Some(&c) = self.text.get(self.pos) else {
                return Err(self.error("unexpected end of the PSSM"));
            };
            match c {
                b' ' | b'\t' => self.pos += 1,
                b'\r' | b'\n' => {
                    self.pos += 1;
                    self.skip_end_of_line(c);
                }
                b'-' if self.peek_at(1) == b'-' => {
                    self.pos += 2;
                    loop {
                        let Some(&c) = self.text.get(self.pos) else {
                            return Err(self.error("unexpected end of the PSSM"));
                        };
                        self.pos += 1;
                        match c {
                            b'\r' | b'\n' => {
                                self.skip_end_of_line(c);
                                break;
                            }
                            b'-' => {
                                let Some(&c) = self.text.get(self.pos) else {
                                    return Err(self.error("unexpected end of the PSSM"));
                                };
                                self.pos += 1;
                                match c {
                                    b'\r' | b'\n' => {
                                        self.skip_end_of_line(c);
                                        break;
                                    }
                                    b'-' => break,
                                    _ => {}
                                }
                            }
                            _ => {}
                        }
                    }
                }
                _ => return Ok(c),
            }
        }
    }

    /// `Expect(c, true)`.
    fn expect(&mut self, c: u8) -> Result<()> {
        if self.skip_white_space()? != c {
            return Err(self.error(format!("'{}' expected", c as char)));
        }
        self.pos += 1;
        Ok(())
    }

    /// `ScanEndOfId(is_id)` at the current byte: the byte, then the identifier bytes and the `-`
    /// before one; empty when `is_id` is false.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:276-290
    // ```c++
    // CTempString CObjectIStreamAsn::ScanEndOfId(bool isId)
    // {
    //     if ( isId ) {
    //         for ( size_t i = 1; ; ++i ) {
    //             char c = m_Input.PeekCharNoEOF(i);
    //             if ( !IdChar(c) &&
    //                  (c != '-' || !IdChar(m_Input.PeekChar(i + 1))) ) {
    //                 const char* ptr = m_Input.GetCurrentPos();
    //                 m_Input.SkipChars(i);
    //                 return CTempString(ptr, i);
    //             }
    //         }
    //     }
    //     return CTempString();
    // }
    // ```
    fn scan_end_of_id(&mut self, is_id: bool) -> &'a [u8] {
        if !is_id {
            return &[];
        }
        let mut i = 1;
        loop {
            let c = self.peek_at(i);
            if !is_id_char(c) && (c != b'-' || !is_id_char(self.peek_at(i + 1))) {
                let start = self.pos;
                self.pos += i;
                return &self.text[start..self.pos];
            }
            i += 1;
        }
    }

    /// `ReadMemberId(SkipWhiteSpace())` (objistrasn.cpp:343-364): a member or variant name.
    fn member_id(&mut self) -> Result<&'a [u8]> {
        let c = self.skip_white_space()?;
        Ok(self.scan_end_of_id(c.is_ascii_lowercase()))
    }

    /// `StartBlock` (objistrasn.cpp:1199-1203).
    fn start_block(&mut self) -> Result<bool> {
        self.expect(b'{')?;
        Ok(true)
    }

    /// `NextElement`: whether another element follows (its `,` taken), or the block ends.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1205-1223
    // ```c++
    //     char c = SkipWhiteSpace();
    //     if ( m_BlockStart ) {
    //         // first element
    //         m_BlockStart = false;
    //         return c != '}';
    //     }
    //     else {
    //         // next element
    //         if ( c == ',' ) {
    //             m_Input.SkipChar();
    //             return true;
    //         }
    //         else if ( c != '}' )
    //             ThrowError(fFormatError, "',' or '}' expected");
    //         return false;
    //     }
    // ```
    fn next_element(&mut self, block_start: &mut bool) -> Result<bool> {
        let c = self.skip_white_space()?;
        if *block_start {
            *block_start = false;
            return Ok(c != b'}');
        }
        if c == b',' {
            self.pos += 1;
            return Ok(true);
        }
        if c != b'}' {
            return Err(self.error("',' or '}' expected"));
        }
        Ok(false)
    }

    /// `EndBlock`: the `}` that `next_element` stopped at.
    fn end_block(&mut self) -> Result<()> {
        self.expect(b'}')
    }

    fn begin_members(&mut self) -> Result<Members> {
        Ok(Members {
            block_start: self.start_block()?,
            last: None,
        })
    }

    /// The next member of a SEQUENCE whose members are `names` (in module order), as its
    /// position; `None` at the end of the SEQUENCE (its `}` taken). A member out of order or not
    /// in the module is an error.
    fn next_member(&mut self, members: &mut Members, names: &[&str]) -> Result<Option<usize>> {
        if !self.next_element(&mut members.block_start)? {
            self.end_block()?;
            return Ok(None);
        }
        let id = self.member_id()?;
        let Some(index) = names.iter().position(|name| name.as_bytes() == id) else {
            return Err(self.error(format!("unexpected member {}", String::from_utf8_lossy(id))));
        };
        if members.last.is_some_and(|last| index <= last) {
            return Err(self.error(format!(
                "member {} out of order",
                String::from_utf8_lossy(id)
            )));
        }
        members.last = Some(index);
        Ok(Some(index))
    }

    /// The elements of a SEQUENCE OF or SET OF, each read by `element`.
    fn list<T>(&mut self, mut element: impl FnMut(&mut Self) -> Result<T>) -> Result<Vec<T>> {
        let mut block_start = self.start_block()?;
        let mut out = Vec::new();
        while self.next_element(&mut block_start)? {
            out.push(element(self)?);
        }
        self.end_block()?;
        Ok(out)
    }

    /// `GetInt8` after `SkipWhiteSpace` (`ReadInt4` / `ReadInt8`): an optional sign and digits,
    /// in `min..=max`.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-668
    // ```c++
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
    //     ...
    //     Uint4 n = c - '0';
    //     if ( n > 9 )
    //         BadNumber();
    //     ...
    //         // check multiplication overflow
    //         if ( n > kMaxBeforeMul || (n == kMaxBeforeMul && d > kMaxLimitAdd) )
    //             NumberOverflow();
    // ```
    fn integer(&mut self, min: i128, max: i128) -> Result<i128> {
        let mut c = self.skip_white_space()?;
        let negative = c == b'-';
        if c == b'-' || c == b'+' {
            self.pos += 1;
            c = self.peek_at(0);
        }
        if !c.is_ascii_digit() {
            return Err(self.error("bad number"));
        }
        let mut n: i128 = 0;
        while self.peek_at(0).is_ascii_digit() {
            n = n * 10 + i128::from(self.peek_at(0) - b'0');
            if n > max.max(-min) {
                return Err(self.error("number overflow"));
            }
            self.pos += 1;
        }
        let n = if negative { -n } else { n };
        if n < min || n > max {
            return Err(self.error("number overflow"));
        }
        Ok(n)
    }

    fn int4(&mut self) -> Result<i32> {
        Ok(self.integer(i32::MIN.into(), i32::MAX.into())? as i32)
    }

    fn uint4(&mut self) -> Result<u32> {
        Ok(self.integer(0, u32::MAX.into())? as u32)
    }

    fn int8(&mut self) -> Result<i64> {
        Ok(self.integer(i64::MIN.into(), i64::MAX.into())? as i64)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:739-764
    // ```c++
    //     switch ( SkipWhiteSpace() ) {
    //     case 'T':
    //         if ( m_Input.PeekCharNoEOF(1) == 'R' &&
    //              m_Input.PeekCharNoEOF(2) == 'U' &&
    //              m_Input.PeekCharNoEOF(3) == 'E' &&
    //              !IdChar(m_Input.PeekCharNoEOF(4)) ) {
    //             m_Input.SkipChars(4);
    //             return true;
    //         }
    //         break;
    //     case 'F':
    //         ...
    //     ThrowError(fFormatError, "TRUE or FALSE expected");
    // ```
    fn boolean(&mut self) -> Result<bool> {
        let c = self.skip_white_space()?;
        let rest = &self.text[self.pos..];
        if c == b'T' && rest.starts_with(b"TRUE") && !is_id_char(self.peek_at(4)) {
            self.pos += 4;
            return Ok(true);
        }
        if c == b'F' && rest.starts_with(b"FALSE") && !is_id_char(self.peek_at(5)) {
            self.pos += 5;
            return Ok(false);
        }
        Err(self.error("TRUE or FALSE expected"))
    }

    /// `ReadEnum` of an ENUMERATED type with the values `values`: a name or a number.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    // ```c++
    //     CTempString id = ReadLCaseId(SkipWhiteSpace());
    //     if ( !id.empty() ) {
    //         // enum element by name
    //         return values.FindValue(id);
    //     }
    //     // enum element by value
    //     TEnumValueType value = m_Input.GetInt4();
    //     if ( !values.IsInteger() ) // check value
    //         values.FindName(value, false);
    // ```
    fn enumerated(&mut self, values: &[(&'static str, i32)]) -> Result<(&'static str, i32)> {
        let c = self.skip_white_space()?;
        let id = self.scan_end_of_id(c.is_ascii_lowercase());
        let found = if id.is_empty() {
            let value = self.int4()?;
            values.iter().find(|&&(_, number)| number == value)
        } else {
            values.iter().find(|&&(name, _)| name.as_bytes() == id)
        };
        found
            .copied()
            .ok_or_else(|| self.error("invalid enum value"))
    }

    /// `ReadStringValue`: a quoted string, its lines joined, `""` read as one quote.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-986
    // ```c++
    //     Expect('\"', true);
    //     ...
    //         for (;;) {
    //             char c = m_Input.PeekChar(i);
    //             switch ( c ) {
    //             case '\r':
    //             case '\n':
    //                 // flush string
    //                 AppendLongStringData(s, i, fix_method, startLine);
    //                 m_Input.SkipChar(); // '\r' or '\n'
    //                 i = 0;
    //                 // skip end of line
    //                 SkipEndOfLine(c);
    //                 break;
    //             case '\"':
    //                 s.reserve(s.size() + i);
    //                 AppendStringData(s, i, fix_method, startLine);
    //                 m_Input.SkipChar(); // quote
    //                 if ( m_Input.PeekCharNoEOF() != '\"' ) {
    //                     // end of string
    //                     return;
    //                 }
    //                 else {
    //                     // double quote -> one quote
    //                     i = 1;
    //                 }
    //                 break;
    // ```
    fn string(&mut self) -> Result<Vec<u8>> {
        self.expect(b'"')?;
        let start_line = self.line;
        let mut out = Vec::new();
        loop {
            let Some(&c) = self.text.get(self.pos) else {
                return Err(PssmTextError {
                    line: self.line,
                    message: format!("unclosed string starts at line {start_line}"),
                });
            };
            self.pos += 1;
            match c {
                b'\r' | b'\n' => self.skip_end_of_line(c),
                b'"' => {
                    if self.peek_at(0) != b'"' {
                        return Ok(out);
                    }
                    out.push(b'"');
                    self.pos += 1;
                }
                _ => out.push(c),
            }
        }
    }

    /// An OCTET STRING: `'`, hex digits across lines, `'H`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1477-1538
    // ```c++
    // void CObjectIStreamAsn::BeginBytes(ByteBlock& )
    // {
    //     Expect('\'', true);
    // }
    //
    // int CObjectIStreamAsn::GetHexChar(void)
    // {
    //     for ( ;; ) {
    //         char c = GetChar();
    //         if ( c >= '0' && c <= '9' ) {
    //             return c - '0';
    //         }
    //         else if ( c >= 'A' && c <= 'F' ) {
    //             return c - 'A' + 10;
    //         }
    //         else if ( c >= 'a' && c <= 'f' ) {
    //             return c - 'a' + 10;
    //         }
    //         switch ( c ) {
    //         case '\'':
    //             return -1;
    //         case '\r':
    //         case '\n':
    //             SkipEndOfLine(c);
    //             break;
    //         default:
    //             m_Input.UngetChar(c);
    //             ThrowError(fFormatError, "bad char in octet string: #"
    //                 + NStr::IntToString(c));
    //         }
    //     }
    // }
    // ...
    //         int c1 = GetHexChar();
    //         if ( c1 < 0 ) {
    //             block.EndOfBlock();
    //             return count;
    //         }
    //         int c2 = GetHexChar();
    //         if ( c2 < 0 ) {
    //             *dst++ = char(c1 << 4);
    //             count++;
    //             block.EndOfBlock();
    //             return count;
    //         }
    //         else {
    //             *dst++ = char((c1 << 4) | c2);
    //     ...
    // void CObjectIStreamAsn::EndBytes(const ByteBlock& )
    // {
    //     Expect('H');
    // }
    // ```
    fn octets(&mut self) -> Result<Vec<u8>> {
        self.expect(b'\'')?;
        let mut out = Vec::new();
        let mut high: Option<u8> = None;
        loop {
            let Some(&c) = self.text.get(self.pos) else {
                return Err(self.error("unexpected end of the PSSM"));
            };
            self.pos += 1;
            let digit = match c {
                b'0'..=b'9' => c - b'0',
                b'A'..=b'F' => c - b'A' + 10,
                b'a'..=b'f' => c - b'a' + 10,
                b'\'' => break,
                b'\r' | b'\n' => {
                    self.skip_end_of_line(c);
                    continue;
                }
                _ => {
                    self.pos -= 1;
                    return Err(self.error(format!("bad char in octet string: #{}", c as i8)));
                }
            };
            match high.take() {
                Some(h) => out.push((h << 4) | digit),
                None => high = Some(digit),
            }
        }
        if let Some(h) = high {
            out.push(h << 4);
        }
        if self.peek_at(0) != b'H' {
            return Err(self.error("'H' expected"));
        }
        self.pos += 1;
        Ok(out)
    }

    /// `ReadNumber`: a sign or a digit, then digits.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:315-329
    // ```c++
    //     char c = SkipWhiteSpace();
    //     if ( c != '-' && c != '+' && !isdigit((unsigned char) c) )
    //         ThrowError(fFormatError, "invalid number");
    //     for ( size_t i = 1; ; ++i ) {
    //         c = m_Input.PeekChar(i);
    //         if ( !isdigit((unsigned char) c) ) {
    //             const char* ptr = m_Input.GetCurrentPos();
    //             m_Input.SkipChars(i);
    //             return CTempString(ptr, i);
    //         }
    //     }
    // ```
    fn number(&mut self) -> Result<&'a [u8]> {
        let c = self.skip_white_space()?;
        if c != b'-' && c != b'+' && !c.is_ascii_digit() {
            return Err(self.error("invalid number"));
        }
        let start = self.pos;
        let mut i = 1;
        while self.peek_at(i).is_ascii_digit() {
            i += 1;
        }
        self.pos += i;
        Ok(&self.text[start..self.pos])
    }

    /// `CObjectIStreamAsn::ReadDouble`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857
    // ```c++
    //     double result = 0;
    //     if (PeekChar(true) != '{') {
    //         CTempString tmp( ScanEndOfId(true) );
    //         if (NStr::strncasecmp(tmp.data(), "PLUS-INFINITY", 13) == 0) {
    //             return HUGE_VAL;
    //         } else if (NStr::strncasecmp(tmp.data(),"MINUS-INFINITY", 14) == 0) {
    //             return -HUGE_VAL;
    //         } else if (NStr::strncasecmp(tmp.data(),"NOT-A-NUMBER", 12) == 0) {
    // 	    return HUGE_VAL/HUGE_VAL; /* NCBI_FAKE_WARNING */
    //         }
    //         char* endptr;
    //         return NStr::StringToDoublePosix( string(tmp).c_str(), &endptr, NStr::fDecimalPosixFinite );
    //     }
    //     Expect('{', true);
    //     bool is_negative = PeekChar(true) == '-';
    //     CTempString mantissaStr = ReadNumber();
    //     size_t mantissaLength = mantissaStr.size();
    //     char buffer[128];
    //     if ( mantissaLength >= sizeof(buffer) - 1 )
    //         ThrowError(fOverflow, "buffer overflow");
    //     memcpy(buffer, mantissaStr.data(), mantissaLength);
    //     buffer[mantissaLength] = '\0';
    //     char* endptr;
    //     double mantissa = NStr::StringToDoublePosix(buffer, &endptr);
    //     if ( *endptr != 0 )
    //         ThrowError(fFormatError, "bad double in line "
    //             + NStr::SizetToString(m_Input.GetLine()));
    //     Expect(',', true);
    //     unsigned base = ReadUint4();
    //     Expect(',', true);
    //     int exp = ReadInt4();
    //     Expect('}', true);
    //     if ( base != 2 && base != 10 ) {
    //         ThrowError(fFormatError, "illegal REAL base (must be 2 or 10)");
    //     }
    //     if (mantissa == 0.) {
    //         return mantissa;
    //     }
    //     if (is_negative) {
    //         mantissa = -mantissa;
    //     }
    //
    //     if ( base == 10 ) {
    //         result = mantissa * pow((double)10.0, exp);
    //     } else {
    //         result = ldexp( mantissa, exp);
    //     }
    //
    //     if (result >= 0 && result <= DBL_MIN) {
    //         result = DBL_MIN;
    //     } else if (!finite(result)) {
    //         result = DBL_MAX;
    //     }
    //     return is_negative ? -result : result;
    // ```
    // `strncasecmp` compares the stream's bytes from the word's start (the word is not
    // terminated); the comparison stops at the first difference, so a byte past the text's end
    // compares as the terminating 0.
    fn real(&mut self) -> Result<f64> {
        if self.skip_white_space()? != b'{' {
            let start = self.pos;
            let word = self.scan_end_of_id(true);
            let rest = &self.text[start..];
            let starts_with = |name: &[u8]| {
                rest.len() >= name.len() && rest[..name.len()].eq_ignore_ascii_case(name)
            };
            if starts_with(b"PLUS-INFINITY") {
                return Ok(f64::INFINITY);
            }
            if starts_with(b"MINUS-INFINITY") {
                return Ok(f64::NEG_INFINITY);
            }
            if starts_with(b"NOT-A-NUMBER") {
                return Ok(f64::NAN);
            }
            return string_to_double_posix(word, true)
                .map(|(value, _)| value)
                .map_err(|message| self.error(message));
        }
        self.pos += 1;
        let is_negative = self.skip_white_space()? == b'-';
        let mantissa_text = self.number()?;
        if mantissa_text.len() >= 128 - 1 {
            return Err(self.error("buffer overflow"));
        }
        let (mut mantissa, used) =
            string_to_double_posix(mantissa_text, false).map_err(|message| self.error(message))?;
        if used != mantissa_text.len() {
            return Err(self.error(format!("bad double in line {}", self.line)));
        }
        self.expect(b',')?;
        let base = self.uint4()?;
        self.expect(b',')?;
        let exp = self.int4()?;
        self.expect(b'}')?;
        if base != 2 && base != 10 {
            return Err(self.error("illegal REAL base (must be 2 or 10)"));
        }
        if mantissa == 0.0 {
            return Ok(mantissa);
        }
        if is_negative {
            mantissa = -mantissa;
        }
        let mut result = if base == 10 {
            mantissa * libm_pow10(exp)
        } else {
            ldexp(mantissa, exp)
        };
        if (0.0..=f64::MIN_POSITIVE).contains(&result) {
            result = f64::MIN_POSITIVE;
        } else if !result.is_finite() {
            result = f64::MAX;
        }
        Ok(if is_negative { -result } else { result })
    }

    /// `skip_value` of the member `name`, noted in `skipped`.
    fn skip_member(&mut self, name: &'static str) -> Result<()> {
        self.skipped.push(name);
        self.skip_value()
    }

    /// Skips one value of any type: a block with what it holds, a string, an OCTET STRING or
    /// BIT STRING, or the words and numbers up to the `,` or `}` that ends it (a CHOICE is its
    /// variant's name and value).
    fn skip_value(&mut self) -> Result<()> {
        let mut depth = 0usize;
        let mut any = false;
        loop {
            let c = self.skip_white_space()?;
            match c {
                b'{' => {
                    self.pos += 1;
                    depth += 1;
                }
                b'}' | b',' if depth == 0 => {
                    if !any {
                        return Err(self.error("value expected"));
                    }
                    return Ok(());
                }
                b'}' => {
                    self.pos += 1;
                    depth -= 1;
                }
                b'"' => {
                    self.string()?;
                }
                b'\'' => {
                    self.pos += 1;
                    loop {
                        let Some(&c) = self.text.get(self.pos) else {
                            return Err(self.error("unexpected end of the PSSM"));
                        };
                        self.pos += 1;
                        match c {
                            b'\'' => break,
                            b'\r' | b'\n' => self.skip_end_of_line(c),
                            _ => {}
                        }
                    }
                    if matches!(self.peek_at(0), b'H' | b'B') {
                        self.pos += 1;
                    }
                }
                _ => self.pos += 1,
            }
            any = true;
        }
    }

    fn object_id(&mut self) -> Result<ObjectId> {
        let variant = self.member_id()?;
        match variant {
            b"id" => Ok(ObjectId::Id(self.int4()?)),
            b"str" => Ok(ObjectId::Str(self.string()?)),
            _ => Err(self.error(format!(
                "unexpected Object-id variant {}",
                String::from_utf8_lossy(variant)
            ))),
        }
    }

    fn pssm_with_parameters(&mut self) -> Result<PssmWithParameters> {
        let mut members = self.begin_members()?;
        let mut pssm = None;
        let mut params = None;
        while let Some(index) = self.next_member(&mut members, &["pssm", "params"])? {
            match index {
                0 => pssm = Some(self.pssm()?),
                _ => params = Some(self.params()?),
            }
        }
        Ok(PssmWithParameters {
            pssm: pssm.ok_or_else(|| self.error("member pssm missing"))?,
            params,
        })
    }

    fn pssm(&mut self) -> Result<Pssm> {
        const NAMES: [&str; 9] = [
            "isProtein",
            "identifier",
            "numRows",
            "numColumns",
            "rowLabels",
            "byRow",
            "query",
            "intermediateData",
            "finalData",
        ];
        let mut members = self.begin_members()?;
        let mut pssm = Pssm {
            is_protein: None,
            identifier: None,
            num_rows: 0,
            num_columns: 0,
            row_labels: None,
            by_row: None,
            query: None,
            intermediate_data: None,
            final_data: None,
        };
        let (mut rows, mut columns) = (None, None);
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => pssm.is_protein = Some(self.boolean()?),
                1 => pssm.identifier = Some(self.object_id()?),
                2 => rows = Some(self.int4()?),
                3 => columns = Some(self.int4()?),
                4 => pssm.row_labels = Some(self.list(Self::string)?),
                5 => pssm.by_row = Some(self.boolean()?),
                6 => pssm.query = Some(self.seq_entry()?),
                7 => pssm.intermediate_data = Some(self.intermediate_data()?),
                _ => pssm.final_data = Some(self.final_data()?),
            }
        }
        pssm.num_rows = rows.ok_or_else(|| self.error("member numRows missing"))?;
        pssm.num_columns = columns.ok_or_else(|| self.error("member numColumns missing"))?;
        Ok(pssm)
    }

    /// The query's `Seq-entry`: its `seq` variant, a Bioseq.
    fn seq_entry(&mut self) -> Result<PssmQuery> {
        let variant = self.member_id()?;
        if variant != b"seq" {
            return Err(self.error(format!(
                "a query Seq-entry {} is not supported by LOSAT's PSSM reader",
                String::from_utf8_lossy(variant)
            )));
        }
        let mut members = self.begin_members()?;
        let (mut ids, mut titles, mut inst) = (None, None, None);
        while let Some(index) = self.next_member(&mut members, &["id", "descr", "inst", "annot"])? {
            match index {
                0 => ids = Some(self.list(Self::seq_id)?),
                1 => {
                    // NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
                    // ```c++
                    // void CSeq_descr::PostRead(void) const
                    // {
                    //     if ( !s_SeqDescrAllowEmpty->Get() && Get().empty() ) {
                    //         NCBI_THROW(CSerialException, eInvalidData,
                    //                    "empty Seq-descr is not allowed");
                    // ```
                    // (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless set in NCBI's environment or
                    // configuration file, which LOSAT does not read.)
                    let descriptors = self.list(Self::seqdesc_title)?.into_iter();
                    if descriptors.len() == 0 {
                        return Err(self.error("empty Seq-descr is not allowed"));
                    }
                    let found: Vec<Vec<u8>> = descriptors.flatten().collect();
                    titles = (!found.is_empty()).then_some(found);
                }
                2 => inst = Some(self.seq_inst()?),
                _ => self.skip_member("annot")?,
            }
        }
        Ok(PssmQuery {
            ids: ids.ok_or_else(|| self.error("member id missing"))?,
            titles,
            inst: inst.ok_or_else(|| self.error("member inst missing"))?,
        })
    }

    /// One `Seqdesc`: its text when it is a `title`; another variant is skipped.
    fn seqdesc_title(&mut self) -> Result<Option<Vec<u8>>> {
        let variant = self.member_id()?;
        if variant.is_empty() {
            return Err(self.error("choice variant id expected"));
        }
        if variant == b"title" {
            return Ok(Some(self.string()?));
        }
        self.skip_member("Seqdesc")?;
        Ok(None)
    }

    fn seq_id(&mut self) -> Result<PssmSeqId> {
        let variant = self.member_id()?;
        match variant {
            b"local" => Ok(PssmSeqId::Local(self.object_id()?)),
            b"gi" => Ok(PssmSeqId::Gi(self.int8()?)),
            b"general" => {
                let mut members = self.begin_members()?;
                let (mut db, mut tag) = (None, None);
                while let Some(index) = self.next_member(&mut members, &["db", "tag"])? {
                    match index {
                        0 => db = Some(self.string()?),
                        _ => tag = Some(self.object_id()?),
                    }
                }
                Ok(PssmSeqId::General {
                    db: db.ok_or_else(|| self.error("member db missing"))?,
                    tag: tag.ok_or_else(|| self.error("member tag missing"))?,
                })
            }
            _ => {
                let Some(&(name, tag)) = TEXTSEQ_VARIANTS
                    .iter()
                    .find(|(name, _)| name.as_bytes() == variant)
                else {
                    return Err(self.error(format!(
                        "a query Seq-id {} is not supported by LOSAT's PSSM reader",
                        String::from_utf8_lossy(variant)
                    )));
                };
                let mut members = self.begin_members()?;
                let mut id = TextseqId::default();
                while let Some(index) =
                    self.next_member(&mut members, &["name", "accession", "release", "version"])?
                {
                    match index {
                        0 => id.name = Some(self.string()?),
                        1 => id.accession = Some(self.string()?),
                        2 => id.release = Some(self.string()?),
                        _ => id.version = Some(self.int4()?),
                    }
                }
                Ok(PssmSeqId::Textseq {
                    variant: name,
                    tag,
                    id,
                })
            }
        }
    }

    fn seq_inst(&mut self) -> Result<SeqInst> {
        const NAMES: [&str; 9] = [
            "repr", "mol", "length", "fuzz", "topology", "strand", "seq-data", "ext", "hist",
        ];
        let mut members = self.begin_members()?;
        let (mut repr, mut mol, mut length, mut seq_data) = (None, None, None, None);
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => repr = Some(self.enumerated(&SEQ_INST_REPR)?),
                1 => mol = Some(self.enumerated(&SEQ_INST_MOL)?),
                2 => length = Some(self.int4()?),
                6 => seq_data = Some(self.seq_data()?),
                other => self.skip_member(NAMES[other])?,
            }
        }
        Ok(SeqInst {
            repr: repr.ok_or_else(|| self.error("member repr missing"))?,
            mol: mol.ok_or_else(|| self.error("member mol missing"))?,
            length,
            seq_data,
        })
    }

    fn seq_data(&mut self) -> Result<SeqData> {
        let variant = self.member_id()?;
        match variant {
            b"iupacaa" => Ok(SeqData::Iupacaa(self.string()?)),
            b"ncbieaa" => Ok(SeqData::Ncbieaa(self.string()?)),
            b"ncbistdaa" => Ok(SeqData::Ncbistdaa(self.octets()?)),
            _ => Err(self.error(format!(
                "query Seq-data {} is not supported by LOSAT's PSSM reader",
                String::from_utf8_lossy(variant)
            ))),
        }
    }

    fn intermediate_data(&mut self) -> Result<PssmIntermediateData> {
        const NAMES: [&str; 9] = [
            "resFreqsPerPos",
            "weightedResFreqsPerPos",
            "freqRatios",
            "informationContent",
            "gaplessColumnWeights",
            "sigma",
            "intervalSizes",
            "numMatchingSeqs",
            "numIndeptObsr",
        ];
        let mut members = self.begin_members()?;
        let mut data = PssmIntermediateData::default();
        let mut freq_ratios = None;
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => data.res_freqs_per_pos = Some(self.list(Self::int4)?),
                1 => data.weighted_res_freqs_per_pos = Some(self.list(Self::real)?),
                2 => freq_ratios = Some(self.list(Self::real)?),
                3 => data.information_content = Some(self.list(Self::real)?),
                4 => data.gapless_column_weights = Some(self.list(Self::real)?),
                5 => data.sigma = Some(self.list(Self::real)?),
                6 => data.interval_sizes = Some(self.list(Self::int4)?),
                7 => data.num_matching_seqs = Some(self.list(Self::int4)?),
                _ => data.num_indept_obsr = Some(self.list(Self::real)?),
            }
        }
        data.freq_ratios = freq_ratios.ok_or_else(|| self.error("member freqRatios missing"))?;
        Ok(data)
    }

    fn final_data(&mut self) -> Result<PssmFinalData> {
        const NAMES: [&str; 9] = [
            "scores",
            "lambda",
            "kappa",
            "h",
            "scalingFactor",
            "lambdaUngapped",
            "kappaUngapped",
            "hUngapped",
            "wordScoreThreshold",
        ];
        let mut members = self.begin_members()?;
        let (mut scores, mut lambda, mut kappa, mut h) = (None, None, None, None);
        let mut data = PssmFinalData {
            scores: Vec::new(),
            lambda: 0.0,
            kappa: 0.0,
            h: 0.0,
            scaling_factor: None,
            lambda_ungapped: None,
            kappa_ungapped: None,
            h_ungapped: None,
            word_score_threshold: None,
        };
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => scores = Some(self.list(Self::int4)?),
                1 => lambda = Some(self.real()?),
                2 => kappa = Some(self.real()?),
                3 => h = Some(self.real()?),
                4 => data.scaling_factor = Some(self.int4()?),
                5 => data.lambda_ungapped = Some(self.real()?),
                6 => data.kappa_ungapped = Some(self.real()?),
                7 => data.h_ungapped = Some(self.real()?),
                _ => data.word_score_threshold = Some(self.real()?),
            }
        }
        data.scores = scores.ok_or_else(|| self.error("member scores missing"))?;
        data.lambda = lambda.ok_or_else(|| self.error("member lambda missing"))?;
        data.kappa = kappa.ok_or_else(|| self.error("member kappa missing"))?;
        data.h = h.ok_or_else(|| self.error("member h missing"))?;
        Ok(data)
    }

    fn params(&mut self) -> Result<PssmParameters> {
        const NAMES: [&str; 6] = [
            "pseudocount",
            "rpsdbparams",
            "constraints",
            "bitScoreThresh",
            "bitScoreReportingThresh",
            "annotatedSites",
        ];
        let mut members = self.begin_members()?;
        let mut params = PssmParameters::default();
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => params.pseudocount = Some(self.int4()?),
                1 => params.rpsdbparams = Some(self.rpsdbparams()?),
                3 => params.bit_score_thresh = Some(self.real()?),
                4 => params.bit_score_reporting_thresh = Some(self.real()?),
                other => self.skip_member(NAMES[other])?,
            }
        }
        Ok(params)
    }

    fn rpsdbparams(&mut self) -> Result<FormatRpsDbParameters> {
        let mut members = self.begin_members()?;
        let (mut matrix_name, mut gap_open, mut gap_extend) = (None, None, None);
        while let Some(index) =
            self.next_member(&mut members, &["matrixName", "gapOpen", "gapExtend"])?
        {
            match index {
                0 => matrix_name = Some(self.string()?),
                1 => gap_open = Some(self.int4()?),
                _ => gap_extend = Some(self.int4()?),
            }
        }
        Ok(FormatRpsDbParameters {
            matrix_name: matrix_name.ok_or_else(|| self.error("member matrixName missing"))?,
            gap_open,
            gap_extend,
        })
    }
}

/// glibc's `pow(10.0, exp)` (the C library NCBI 2.17.0 for x86-64 Linux runs with, version 2.28
/// on): the power of ten rounded to the nearest double, except at the exponents 23 and 210,
/// where glibc's result is the double above the nearest one (its error bound is 0.52 units in
/// the last place, and 10^23 and 10^210 lie almost halfway between two doubles). The table keeps
/// the bytes the same on every platform; `libm_pow10_is_glibcs_pow` checks it against the C
/// library on Linux.
pub(crate) fn libm_pow10(exp: i32) -> f64 {
    match exp {
        23 => f64::from_bits(0x44b5_2d02_c7e1_4af7),
        210 => f64::from_bits(0x6b88_557f_3132_6bbc),
        _ => format!("1e{exp}")
            .parse()
            .expect("a power of ten parses as a double"),
    }
}

/// The C library's `ldexp(x, exp)`: `x * 2^exp` rounded once to the nearest double (ties to even).
fn ldexp(x: f64, exp: i32) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let bits = x.to_bits();
    let sign = bits & (1 << 63);
    let mut biased = ((bits >> 52) & 0x7ff) as i64;
    let mut significand = bits & ((1u64 << 52) - 1);
    if biased == 0 {
        // A subnormal x: normalize its significand.
        let shift = significand.leading_zeros() as i64 - 11;
        significand <<= shift;
        biased = 1 - shift;
    } else {
        significand |= 1 << 52;
    }
    // x = significand * 2^(biased - 1075), significand in [2^52, 2^53).
    let new_biased = biased + i64::from(exp);
    if new_biased >= 0x7ff {
        return f64::from_bits(sign | f64::INFINITY.to_bits());
    }
    if new_biased >= 1 {
        return f64::from_bits(
            sign | ((new_biased as u64) << 52) | (significand & ((1 << 52) - 1)),
        );
    }
    // A subnormal or zero result: significand * 2^(new_biased - 1075) in units of 2^-1074.
    let drop = (1 - new_biased) as u32;
    if drop > 54 {
        return f64::from_bits(sign);
    }
    let kept = significand >> drop;
    let rest = significand & ((1u64 << drop) - 1);
    let half = 1u64 << (drop - 1);
    let up = rest > half || (rest == half && (kept & 1) != 0);
    f64::from_bits(sign | (kept + up as u64))
}

/// The refusal of a decimal REAL whose exponent leaves the range of `int` in
/// `StringToDoublePosix`'s arithmetic (`dot_position - digits`, `exponent - expvalue`): NCBI's
/// result is then undefined behaviour, which LOSAT does not reproduce.
fn exponent_out_of_range() -> String {
    "a REAL whose exponent is outside the range of NCBI's int arithmetic is not supported by LOSAT's PSSM reader"
        .to_string()
}

/// `NStr::StringToDoublePosix(text, &endptr, flags)`: the value and the bytes of `text` that it
/// read (`endptr - text`); `finite` is `fDecimalPosixFinite`. The 80-bit `long double` arithmetic
/// of NCBI's x86-64 build is `X87`'s. A number of more than 18 significant digits, which NCBI
/// adds up in `long double` (`third`), is refused, and so is one whose exponent overflows NCBI's
/// `int` (`exponent_out_of_range`).
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:985-1278
// ```c++
//     const char* start = ptr;
//     char c = *ptr++;
//
//     // skip leading whitespace
//     while ( isspace((unsigned char)c) ) {
//         c = *ptr++;
//     }
//
//     int sign = 0;
//     if ( c == '-' ) {
//         sign = -1;
//         c = *ptr++;
//     }
//     else if ( c == '+' ) {
//         sign = +1;
//         c = *ptr++;
//     }
//
//     if (c == 0) {
//         if (endptr) {
//             *endptr = (char*)start;
//         }
//         err_guard.Set(EINVAL);
//         return 0.;
//     }
//
//     // short-cut - single digit
//     if ( !*ptr && c >= '0' && c <= '9' ) {
//         if (endptr) {
//             *endptr = (char*)ptr;
//         }
//         double result = c-'0';
//         // some compilers fail to negate zero
//         return sign < 0 ? (c == '0' ? -0. : -result) : result;
//     }
//
//     bool         dot = false, expn = false, anydigits = false;
//     int          digits = 0, dot_position = 0;
//     unsigned int first=0, second=0, first_mul=1;
//     long double  second_mul = NCBI_CONST_LONGDOUBLE(1.),
//                  third      = NCBI_CONST_LONGDOUBLE(0.);
//
//     // up to exponent
//     for ( ; ; c = *ptr++ ) {
//         if (c >= '0' && c <= '9') {
//             // digits: accumulate
//             c = (char)(c - '0');
//             anydigits = true;
//             ++digits;
//             if (first == 0) {
//                 first = c;
//                 if ( first == 0 ) {
//                     // omit leading zeros
//                     --digits;
//                     if (dot) {
//                         --dot_position;
//                     }
//                 }
//             } else if (digits <= 9) {
//                 // first 9 digits come to 'first'
//                 first = first*10 + c;
//             } else if (digits <= 18) {
//                 // next 9 digits come to 'second'
//                 first_mul *= 10;
//                 second = second*10 + c;
//             } else {
//                 // other digits come to 'third'
//                 second_mul *= NCBI_CONST_LONGDOUBLE(10.);
//                 third = third * NCBI_CONST_LONGDOUBLE(10.) + c;
//             }
//         }
//         else if (c == '.') {
//             // dot
//             // if second dot, stop
//             if (dot) {
//                 --ptr;
//                 break;
//             }
//             dot_position = digits;
//             dot = true;
//         }
//         else if (c == 'e' || c == 'E') {
//             // if exponent, stop
//             if (!anydigits) {
//                 --ptr;
//                 break;
//             }
//             expn = true;
//             break;
//         }
//         else {
//             --ptr;
//             if (!anydigits) {
//                 if ( !dot && (c == 'n' || c == 'N') &&
//                      NStr::strncasecmp(ptr,"nan",3)==0) {
//                     if (endptr) {
//                         *endptr = (char*)(ptr+3);
//                     }
//                     return HUGE_VAL/HUGE_VAL; /* NCBI_FAKE_WARNING */
//                 }
//                 if ( (c == 'i' || c == 'I') ) {
//                     if ( NStr::strncasecmp(ptr,"inf",3)==0) {
//                         ptr += 3;
//                         if ( NStr::strncasecmp(ptr,"inity",5)==0) {
//                             ptr += 5;
//                         }
//                         if (endptr) {
//                             *endptr = (char*)ptr;
//                         }
//                         return sign < 0 ? -HUGE_VAL : HUGE_VAL;
//                     }
//                 }
//             }
//             break;
//         }
//     }
//     // if no digits, stop now - error
//     if (!anydigits) {
//         if (endptr) {
//             *endptr = (char*)start;
//         }
//         err_guard.Set(EINVAL);
//         return 0.;
//     }
//     int exponent = dot ? dot_position - digits : 0;
//
//     // read exponent
//     if (expn && *ptr) {
//         int expvalue = 0;
//         bool expsign = false, expnegate= false;
//         int expdigits= 0;
//         for( ; ; ++ptr) {
//             c = *ptr;
//             // sign: should be no digits at this point
//             if (c == '-' || c == '+') {
//                 // if there was sign or digits, stop
//                 if (expsign || expdigits) {
//                     break;
//                 }
//                 expsign = true;
//                 expnegate = c == '-';
//             }
//             // digits: accumulate
//             else if (c >= '0' && c <= '9') {
//                 ++expdigits;
//                 int newexpvalue = expvalue*10 + (c-'0');
//                 if (newexpvalue > expvalue) {
//                     expvalue = newexpvalue;
//                 }
//             }
//             else {
//                 break;
//             }
//         }
//         // if no digits, rollback
//         if (!expdigits) {
//             // rollback sign
//             if (expsign) {
//                 --ptr;
//             }
//             // rollback exponent
//             if (expn) {
//                 --ptr;
//             }
//         }
//         else {
//             exponent = expnegate ? exponent - expvalue : exponent + expvalue;
//         }
//     }
//     long double ret;
//     if ( first_mul > 1 ) {
//         _ASSERT(first);
//         ret = ((long double)first * first_mul + second)* second_mul + third;
//     }
//     else {
//         ...
//         ret = first;
//     }
//     // calculate exponent
//     if ( first && exponent ) {
//         // multiply by power of 10 only non-zero mantissa
//         if (exponent > 2*DBL_MAX_10_EXP) {
//             ret = (flags & fDecimalPosixFinite) ? DBL_MAX :  HUGE_VAL;
//             err_guard.Set(ERANGE);
//         } else if (exponent < 2*DBL_MIN_10_EXP) {
//             ret = (flags & fDecimalPosixFinite) ? DBL_MIN : 0.;
//             err_guard.Set(ERANGE);
//         } else {
//             if ( exponent > 0 ) {
//                 static const double mul1[16] = {
//                     1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7,
//                     1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15
//                 };
//                 ret *= mul1[exponent&15];
//                 if ( exponent >>= 4 ) {
//                     static const long double mul2[16] = {
//                         NCBI_CONST_LONGDOUBLE(1e0),
//                         NCBI_CONST_LONGDOUBLE(1e16),
//                         ...
//                         NCBI_CONST_LONGDOUBLE(1e240)
//                     };
//                     ret *= mul2[exponent&15];
//                     for ( exponent >>= 4; exponent; --exponent ) {
//                         ret *= NCBI_CONST_LONGDOUBLE(1e256);
//                     }
//                 }
//                 if (!finite(double(ret))) {
//                     if (flags & fDecimalPosixFinite) {
//                         ret = DBL_MAX;
//                     }
//                     err_guard.Set(ERANGE);
//                 }
//             }
//             else {
//                 exponent = -exponent;
//                 static const long double mul1[16] = {
//                     NCBI_CONST_LONGDOUBLE(1e-0),
//                     ...
//                     NCBI_CONST_LONGDOUBLE(1e-15)
//                 };
//                 ret *= mul1[exponent&15];
//                 if ( exponent >>= 4 ) {
//                     static const long double mul2[16] = {
//                         NCBI_CONST_LONGDOUBLE(1e-0),
//                         NCBI_CONST_LONGDOUBLE(1e-16),
//                         ...
//                         NCBI_CONST_LONGDOUBLE(1e-240)
//                     };
//                     ret *= mul2[exponent&15];
//                     for ( exponent >>= 4; exponent; --exponent ) {
//                         ret *= NCBI_CONST_LONGDOUBLE(1e-256);
//                     }
//                 }
//                 if ( ret < DBL_MIN ) {
//                     if (flags & fDecimalPosixFinite) {
//                         ret = DBL_MIN;
//                     }
//                     err_guard.Set(ERANGE);
//                 }
//             }
//         }
//     }
//     if ( sign < 0 ) {
//         ret = -ret;
//     }
//     // done
//     if (endptr) {
//         *endptr = (char*)ptr;
//     }
//     return (double)ret;
// ```
pub(crate) fn string_to_double_posix(
    text: &[u8],
    finite: bool,
) -> std::result::Result<(f64, usize), String> {
    // `*p` of the NUL-terminated string.
    let at = |p: usize| text.get(p).copied().unwrap_or(0);
    let mut ptr = 0usize;
    let mut c = at(ptr);
    ptr += 1;
    while is_c_space(c) {
        c = at(ptr);
        ptr += 1;
    }
    let mut sign = 0;
    if c == b'-' {
        sign = -1;
        c = at(ptr);
        ptr += 1;
    } else if c == b'+' {
        sign = 1;
        c = at(ptr);
        ptr += 1;
    }
    if c == 0 {
        return Ok((0.0, 0));
    }
    if at(ptr) == 0 && c.is_ascii_digit() {
        let result = f64::from(c - b'0');
        let value = if sign < 0 {
            if c == b'0' {
                -0.0
            } else {
                -result
            }
        } else {
            result
        };
        return Ok((value, ptr));
    }
    let (mut dot, mut expn, mut anydigits) = (false, false, false);
    let (mut digits, mut dot_position): (i32, i32) = (0, 0);
    let (mut first, mut second, mut first_mul): (u32, u32, u32) = (0, 0, 1);
    loop {
        if c.is_ascii_digit() {
            let d = u32::from(c - b'0');
            anydigits = true;
            digits += 1;
            if first == 0 {
                first = d;
                if first == 0 {
                    digits -= 1;
                    if dot {
                        dot_position = dot_position
                            .checked_sub(1)
                            .ok_or_else(exponent_out_of_range)?;
                    }
                }
            } else if digits <= 9 {
                first = first * 10 + d;
            } else if digits <= 18 {
                first_mul *= 10;
                second = second * 10 + d;
            } else {
                return Err(
                    "a REAL of more than 18 significant digits is not supported by LOSAT's PSSM reader"
                        .to_string(),
                );
            }
        } else if c == b'.' {
            if dot {
                ptr -= 1;
                break;
            }
            dot_position = digits;
            dot = true;
        } else if c == b'e' || c == b'E' {
            if !anydigits {
                ptr -= 1;
                break;
            }
            expn = true;
            break;
        } else {
            ptr -= 1;
            if !anydigits {
                let rest = &text[ptr.min(text.len())..];
                let starts_with = |name: &[u8]| {
                    rest.len() >= name.len() && rest[..name.len()].eq_ignore_ascii_case(name)
                };
                if !dot && (c == b'n' || c == b'N') && starts_with(b"nan") {
                    return Ok((f64::NAN, ptr + 3));
                }
                if (c == b'i' || c == b'I') && starts_with(b"inf") {
                    let mut end = ptr + 3;
                    let after = &text[end.min(text.len())..];
                    if after.len() >= 5 && after[..5].eq_ignore_ascii_case(b"inity") {
                        end += 5;
                    }
                    let value = if sign < 0 {
                        f64::NEG_INFINITY
                    } else {
                        f64::INFINITY
                    };
                    return Ok((value, end));
                }
            }
            break;
        }
        c = at(ptr);
        ptr += 1;
    }
    if !anydigits {
        return Ok((0.0, 0));
    }
    let mut exponent: i32 = if dot {
        dot_position
            .checked_sub(digits)
            .ok_or_else(exponent_out_of_range)?
    } else {
        0
    };
    if expn && at(ptr) != 0 {
        let (mut expvalue, mut expsign, mut expnegate, mut expdigits) = (0i32, false, false, 0);
        loop {
            let c = at(ptr);
            if c == b'-' || c == b'+' {
                if expsign || expdigits != 0 {
                    break;
                }
                expsign = true;
                expnegate = c == b'-';
            } else if c.is_ascii_digit() {
                expdigits += 1;
                let newexpvalue = expvalue.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
                if newexpvalue > expvalue {
                    expvalue = newexpvalue;
                }
            } else {
                break;
            }
            ptr += 1;
        }
        if expdigits == 0 {
            if expsign {
                ptr -= 1;
            }
            if expn {
                ptr -= 1;
            }
        } else {
            exponent = if expnegate {
                exponent.checked_sub(expvalue)
            } else {
                exponent.checked_add(expvalue)
            }
            .ok_or_else(exponent_out_of_range)?;
        }
    }
    // The integer of up to 18 digits is exact in the 64-bit significand.
    let mut ret = if first_mul > 1 {
        X87::from_u64(u64::from(first) * u64::from(first_mul) + u64::from(second))
    } else {
        X87::from_u64(u64::from(first))
    };
    let ld = |n: i32| {
        crate::report::serial::real::long_double_pow10(n)
            .expect("the long double literals of StringToDoublePosix")
    };
    // DBL_MAX_10_EXP 308, DBL_MIN_10_EXP -307.
    if first != 0 && exponent != 0 {
        if exponent > 2 * 308 {
            ret = if finite {
                X87::from_f64(f64::MAX)
            } else {
                return Ok((
                    if sign < 0 {
                        f64::NEG_INFINITY
                    } else {
                        f64::INFINITY
                    },
                    ptr,
                ));
            };
        } else if exponent < 2 * -307 {
            ret = if finite {
                X87::from_f64(f64::MIN_POSITIVE)
            } else {
                X87::ZERO
            };
        } else if exponent > 0 {
            const MUL1: [f64; 16] = [
                1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14,
                1e15,
            ];
            let mut e = exponent;
            ret = ret.mul(X87::from_f64(MUL1[(e & 15) as usize]));
            e >>= 4;
            if e != 0 {
                ret = ret.mul(ld(16 * (e & 15)));
                e >>= 4;
                while e != 0 {
                    ret = ret.mul(ld(256));
                    e -= 1;
                }
            }
            if !ret.to_f64().is_finite() && finite {
                ret = X87::from_f64(f64::MAX);
            }
        } else {
            let mut e = -exponent;
            ret = ret.mul(ld(-(e & 15)));
            e >>= 4;
            if e != 0 {
                ret = ret.mul(ld(-16 * (e & 15)));
                e >>= 4;
                while e != 0 {
                    ret = ret.mul(ld(-256));
                    e -= 1;
                }
            }
            if ret.lt(X87::from_f64(f64::MIN_POSITIVE)) && finite {
                ret = X87::from_f64(f64::MIN_POSITIVE);
            }
        }
    }
    if sign < 0 {
        ret = ret.neg();
    }
    Ok((ret.to_f64(), ptr))
}

#[cfg(test)]
#[path = "pssm_text_tests.rs"]
mod tests;
