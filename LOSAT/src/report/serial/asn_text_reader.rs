//! The lexical layer of NCBI's ASN.1 text object stream (`CObjectIStreamAsn` over
//! `CIStreamBuffer`), shared by LOSAT's readers of ASN.1 text: the PSSM reader
//! (`blastinput::pssm_text`, `PssmWithParameters` from a byte slice) and the Blast4-archive reader
//! of `losat rpsbproc` (`rpsbproc::archive`, objects one after another from a stream). The
//! grammar of each type stays with its reader; this module reads what every type is made of:
//! white space and comments, identifiers and `::=`, blocks and their elements, the members of a
//! SEQUENCE in the order of the module, INTEGER, BOOLEAN, ENUMERATED, VisibleString, OCTET
//! STRING and REAL values, and the skipping of one value of any type.
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243
//! ```c++
//! char CObjectIStreamAsn::SkipWhiteSpace(void)
//! {
//!     try { // catch CEofException
//!         for ( ;; ) {
//!             char c = m_Input.SkipSpaces();
//! ```
//!
//! The text comes from an `AsnSource`: a slice (`SliceSource`, the PSSM file read whole) or a
//! stream read in blocks (`StreamSource`, an archive of any size: only the bytes not yet read
//! and those of the value being read are kept). The readers are generic over the source, so the
//! slice reader is what it was before the stream existed.

use crate::ncbi_str::string_to_double_posix;
use std::io::Read;
use std::ops::Range;

/// What kind of failure an `AsnTextError` is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsnTextErrorKind {
    /// NCBI's object stream fails too (`CSerialException`): a format error, a missing or
    /// unexpected member, a number out of range, the end of the text inside a value.
    Format,
    /// A form that NCBI reads and LOSAT's reader does not (an explicit LOSAT refusal).
    Unsupported,
}

/// A document that the reader cannot read: the line (from 1) where it stopped and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AsnTextError {
    pub line: usize,
    pub message: String,
    pub(crate) kind: AsnTextErrorKind,
}

impl std::fmt::Display for AsnTextError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for AsnTextError {}

pub(crate) type Result<T> = std::result::Result<T, AsnTextError>;

/// The bytes that an ASN.1 text reader reads, by position from the start of the text.
pub(crate) trait AsnSource {
    /// The byte at `pos`, `None` past the end of the text. A position before one given to
    /// `release_before` is never asked for again.
    fn byte_at(&mut self, pos: usize) -> Option<u8>;
    /// The bytes `range`, all of which `byte_at` has given.
    fn bytes(&self, range: Range<usize>) -> &[u8];
    /// The reader will not ask for a byte before `pos` again (a stream may drop them).
    #[inline]
    fn release_before(&mut self, _pos: usize) {}
}

/// A text held whole in memory.
pub(crate) struct SliceSource<'a> {
    text: &'a [u8],
}

impl<'a> SliceSource<'a> {
    pub(crate) fn new(text: &'a [u8]) -> Self {
        SliceSource { text }
    }

    /// The bytes `range` with the text's lifetime.
    pub(crate) fn text(&self, range: Range<usize>) -> &'a [u8] {
        &self.text[range]
    }
}

impl AsnSource for SliceSource<'_> {
    #[inline]
    fn byte_at(&mut self, pos: usize) -> Option<u8> {
        self.text.get(pos).copied()
    }

    #[inline]
    fn bytes(&self, range: Range<usize>) -> &[u8] {
        &self.text[range]
    }
}

/// A text read from a stream in blocks. It keeps the bytes from the last position given to
/// `release_before` on (dropping the bytes before it once they make up most of the buffer), so
/// its memory is about the largest value read between two releases, not the whole text.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:125-139
// ```c++
// void CIStreamBuffer::Close(void)
// {
//     if ( m_Input ) {
//         size_t unused = m_DataEndPos - m_CurrentPos;
//         if ( unused ) {
//             m_Input->Pushback(m_CurrentPos, unused);
//         }
//         m_Input.reset();
//     }
// ```
// NCBI's object stream reads its input in blocks too and gives the bytes it did not use back
// to the stream when an object ends, so the next object is read from the byte after the last
// one. The reader of a stream of objects keeps one `StreamSource` for all of them.
pub(crate) struct StreamSource<R: Read> {
    input: R,
    buf: Vec<u8>,
    /// The position of `buf[0]`.
    base: usize,
    eof: bool,
    /// The block each read fills before its bytes are appended to `buf` (allocated once, so a
    /// short read costs only the bytes it gives).
    scratch: Box<[u8]>,
}

/// The size of one read from the stream.
const STREAM_BLOCK: usize = 64 * 1024;

impl<R: Read> StreamSource<R> {
    /// A source whose text starts with `first` (bytes already taken from `input`) and goes on
    /// with what `input` gives.
    pub(crate) fn new(input: R, first: &[u8]) -> Self {
        StreamSource {
            input,
            buf: first.to_vec(),
            base: 0,
            eof: false,
            scratch: vec![0; STREAM_BLOCK].into_boxed_slice(),
        }
    }

    /// The bytes kept.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.buf.len()
    }

    /// Reads blocks until `pos` is in the buffer or the stream ends. A read error ends the
    /// text, as a failed read ends an `istream` (its bad bit, which the object stream reads as
    /// the end of its input).
    #[cold]
    fn fill_to(&mut self, pos: usize) -> Option<u8> {
        assert!(pos >= self.base, "a released byte was asked for");
        while pos - self.base >= self.buf.len() {
            if self.eof {
                return None;
            }
            let read = loop {
                match self.input.read(&mut self.scratch) {
                    Ok(read) => break read,
                    Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                    Err(_) => break 0,
                }
            };
            self.buf.extend_from_slice(&self.scratch[..read]);
            if read == 0 {
                self.eof = true;
            }
        }
        Some(self.buf[pos - self.base])
    }
}

impl<R: Read> AsnSource for StreamSource<R> {
    #[inline]
    fn byte_at(&mut self, pos: usize) -> Option<u8> {
        let index = pos.wrapping_sub(self.base);
        if index < self.buf.len() {
            return Some(self.buf[index]);
        }
        self.fill_to(pos)
    }

    #[inline]
    fn bytes(&self, range: Range<usize>) -> &[u8] {
        &self.buf[range.start - self.base..range.end - self.base]
    }

    #[inline]
    fn release_before(&mut self, pos: usize) {
        let dead = pos - self.base;
        // Drop the read bytes once they are most of a buffer of more than a few blocks: each
        // byte is moved at most once on average.
        if dead > 4 * STREAM_BLOCK && dead * 2 > self.buf.len() {
            self.buf.drain(..dead);
            self.base = pos;
        }
    }
}

/// `FirstIdChar` (objistrasn.cpp:118-121).
pub(crate) fn is_first_id_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:123-127
// ```c++
// bool CObjectIStreamAsn::IdChar(char c)
// {
//     return isalnum((unsigned char) c) || c == '_' || c == '.';
// }
// ```
pub(crate) fn is_id_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// The reader of one ASN.1 text document (or of several from one stream).
pub(crate) struct Reader<S: AsnSource> {
    pub(crate) src: S,
    /// The position of the next byte to read.
    pub(crate) pos: usize,
    /// `CIStreamBuffer::m_Line`.
    pub(crate) line: usize,
    /// The members skipped as values (`skip_member`), by name.
    pub(crate) skipped: Vec<&'static str>,
    /// What the text is, for the error at its end ("unexpected end of the PSSM").
    document: &'static str,
    /// The reader, for the refusals ("... is not supported by LOSAT's PSSM reader").
    reader: &'static str,
    /// Whether a string with a byte outside the visible range is refused. NCBI replaces such a
    /// byte by `#` and posts an error with the reader's stack trace (`ReplaceVisibleChar`,
    /// objistr.cpp:1854-1895, `eFNP_ReplaceAndWarn`, the default of `SERIAL_WRONG_CHARS_READ`),
    /// which a reader that sets this does not reproduce.
    visible_only: bool,
    /// Whether the refusals of the lexical layer (a REAL or a string that LOSAT does not read)
    /// are noted in `refusal` and the reading goes on (a reader that must know whether NCBI's
    /// reader would read the whole object before it refuses), rather than returned at once.
    defer_refusals: bool,
    /// The first refusal noted (`defer_refusals`).
    pub(crate) refusal: Option<AsnTextError>,
}

/// The members of a SEQUENCE as they are read: each name in turn, in the order of the module.
pub(crate) struct Members {
    block_start: bool,
    /// The position in the module of the last member read.
    last: Option<usize>,
}

impl<S: AsnSource> Reader<S> {
    /// A reader at the start of `src`; `document` and `reader` name the text and the reader in
    /// error messages.
    pub(crate) fn new(src: S, document: &'static str, reader: &'static str) -> Self {
        Reader {
            src,
            pos: 0,
            line: 1,
            skipped: Vec::new(),
            document,
            reader,
            visible_only: false,
            defer_refusals: false,
            refusal: None,
        }
    }

    /// Note the refusals of the lexical layer in `refusal` and go on (`defer_refusals`).
    pub(crate) fn defer_refusals(mut self) -> Self {
        self.defer_refusals = true;
        self
    }

    /// `error`, returned, or noted when the refusals are deferred (then `Ok(value)`).
    fn refuse<T>(&mut self, error: AsnTextError, value: T) -> Result<T> {
        if !self.defer_refusals {
            return Err(error);
        }
        if self.refusal.is_none() {
            self.refusal = Some(error);
        }
        Ok(value)
    }

    /// Refuse strings with bytes outside the visible range (`visible_only`).
    pub(crate) fn refuse_invisible_bytes(mut self) -> Self {
        self.visible_only = true;
        self
    }

    pub(crate) fn error(&self, message: impl Into<String>) -> AsnTextError {
        AsnTextError {
            line: self.line,
            message: message.into(),
            kind: AsnTextErrorKind::Format,
        }
    }

    /// An explicit LOSAT refusal at the current line.
    pub(crate) fn unsupported(&self, message: impl Into<String>) -> AsnTextError {
        AsnTextError {
            line: self.line,
            message: message.into(),
            kind: AsnTextErrorKind::Unsupported,
        }
    }

    /// "`<what>` is not supported by LOSAT's `<reader>`".
    pub(crate) fn not_supported(&self, what: &str) -> AsnTextError {
        self.unsupported(format!(
            "{what} is not supported by LOSAT's {}",
            self.reader
        ))
    }

    fn end_of_text(&self) -> AsnTextError {
        self.error(format!("unexpected end of the {}", self.document))
    }

    /// The byte at `pos + offset`, 0 past the end (`PeekCharNoEOF`).
    #[inline]
    pub(crate) fn peek_at(&mut self, offset: usize) -> u8 {
        self.src.byte_at(self.pos + offset).unwrap_or(0)
    }

    /// The bytes `range` that the reader has read.
    #[inline]
    pub(crate) fn bytes(&self, range: Range<usize>) -> &[u8] {
        self.src.bytes(range)
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
    pub(crate) fn skip_white_space(&mut self) -> Result<u8> {
        loop {
            let Some(c) = self.src.byte_at(self.pos) else {
                return Err(self.end_of_text());
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
                        let Some(c) = self.src.byte_at(self.pos) else {
                            return Err(self.end_of_text());
                        };
                        self.pos += 1;
                        match c {
                            b'\r' | b'\n' => {
                                self.skip_end_of_line(c);
                                break;
                            }
                            b'-' => {
                                let Some(c) = self.src.byte_at(self.pos) else {
                                    return Err(self.end_of_text());
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
    pub(crate) fn expect(&mut self, c: u8) -> Result<()> {
        if self.skip_white_space()? != c {
            return Err(self.error(format!("'{}' expected", c as char)));
        }
        self.pos += 1;
        Ok(())
    }

    /// `ScanEndOfId(is_id)` at the current byte: the byte, then the identifier bytes and the `-`
    /// before one; empty when `is_id` is false. Returns the positions of the identifier.
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
    pub(crate) fn scan_end_of_id(&mut self, is_id: bool) -> Range<usize> {
        if !is_id {
            return self.pos..self.pos;
        }
        let mut i = 1;
        loop {
            let c = self.peek_at(i);
            if !is_id_char(c) && (c != b'-' || !is_id_char(self.peek_at(i + 1))) {
                let start = self.pos;
                self.pos += i;
                return start..self.pos;
            }
            i += 1;
        }
    }

    /// `ReadMemberId(SkipWhiteSpace())`: a member or variant name, bare (from a lower-case letter)
    /// or between `[` and `]` on one line (the positions of the text between them).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:343-363
    // ```c++
    // CTempString CObjectIStreamAsn::ReadMemberId(char c)
    // {
    //     if ( c == '[' ) {
    //         for ( size_t i = 1; ; ++i ) {
    //             switch ( m_Input.PeekChar(i) ) {
    //             case '\r':
    //             case '\n':
    //                 ThrowError(fFormatError, "end of line: expected ']'");
    //                 break;
    //             case ']':
    //                 {
    //                     const char* ptr = m_Input.GetCurrentPos();
    //                     m_Input.SkipChars(++i);
    //                     return CTempString(ptr + 1, i - 2);
    //                 }
    //             }
    //         }
    //     }
    //     else {
    //         return ScanEndOfId(islower((unsigned char) c) != 0);
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:394-397
    // ```c++
    //     if (!id.empty()  &&  isdigit((unsigned char) id[0])) {
    //         idx = classType->GetMembers().Find
    //             (CMemberId::TTag(NStr::StringToInt(id)), CAsnBinaryDefs::eContextSpecific);
    //     }
    // ```
    // An id in brackets that starts with a digit names the member or variant by its tag
    // (`GetMemberIndex`, and `GetChoiceIndex` at :431-433), which depends on the tags NCBI's
    // generated type information gives each member; LOSAT's readers know the members by name
    // only, so that form is an explicit refusal.
    pub(crate) fn member_id_range(&mut self) -> Result<Range<usize>> {
        let c = self.skip_white_space()?;
        if c != b'[' {
            return Ok(self.scan_end_of_id(c.is_ascii_lowercase()));
        }
        let mut i = 1;
        loop {
            let Some(c) = self.src.byte_at(self.pos + i) else {
                return Err(self.end_of_text());
            };
            match c {
                b'\r' | b'\n' => return Err(self.error("end of line: expected ']'")),
                b']' => {
                    let id = self.pos + 1..self.pos + i;
                    self.pos += i + 1;
                    if self
                        .bytes(id.clone())
                        .first()
                        .is_some_and(u8::is_ascii_digit)
                    {
                        return Err(self.not_supported("a member or variant named by its tag"));
                    }
                    return Ok(id);
                }
                _ => i += 1,
            }
        }
    }

    /// The variant of a CHOICE whose variants are `names` (in module order), as its position;
    /// a name that is not one of them is an error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1475
    // ```c++
    // void CObjectIStreamAsn::BeginChoice(const CChoiceTypeInfo* choiceType)
    // {
    //     if (choiceType->GetVariantInfo(kFirstMemberIndex)->GetId().IsAttlist()) {
    //         ...
    //     }
    //     m_BlockStart = true;
    // }
    // ...
    // TMemberIndex CObjectIStreamAsn::BeginChoiceVariant(const CChoiceTypeInfo* choiceType)
    // {
    //     bool skipname = !m_BlockStart;
    //     if ( !NextElement() )
    //         return kInvalidMember;
    //     CTempString id = ReadMemberId(SkipWhiteSpace());
    //     ...
    //     if ( id.empty() )
    //         ThrowError(fFormatError, "choice variant id expected");
    //
    //     TMemberIndex index = GetChoiceIndex(choiceType, id);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownVariants()) {
    //             SetFailFlags(fUnknownValue);
    //         } else {
    //             UnexpectedMember(id, choiceType->GetVariants());
    //         }
    //     }
    // ```
    // (`NextElement` at the choice's start fails at a `}`, which reads as no id here; the variant
    // names of the modules LOSAT reads are their ASN.1 names, `GetChoiceIndex` (:426-442).)
    pub(crate) fn choice(&mut self, names: &[&str]) -> Result<usize> {
        let id = self.member_id_range()?;
        if id.is_empty() {
            return Err(self.error("choice variant id expected"));
        }
        let text = self.bytes(id);
        match names.iter().position(|name| name.as_bytes() == text) {
            Some(index) => Ok(index),
            None => {
                let message = format!("unexpected variant {}", String::from_utf8_lossy(text));
                Err(self.error(message))
            }
        }
    }

    /// `StartBlock` (objistrasn.cpp:1199-1203).
    pub(crate) fn start_block(&mut self) -> Result<bool> {
        self.expect(b'{')?;
        Ok(true)
    }

    /// `NextElement`: whether another element follows (its `,` taken), or the block ends. The
    /// bytes before it are not read again (`AsnSource::release_before`).
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
    pub(crate) fn next_element(&mut self, block_start: &mut bool) -> Result<bool> {
        self.src.release_before(self.pos);
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
    pub(crate) fn end_block(&mut self) -> Result<()> {
        self.expect(b'}')
    }

    pub(crate) fn begin_members(&mut self) -> Result<Members> {
        Ok(Members {
            block_start: self.start_block()?,
            last: None,
        })
    }

    /// The next member of a SEQUENCE whose members are `names` (in module order), as its
    /// position; `None` at the end of the SEQUENCE (its `}` taken). A member out of order or not
    /// in the module is an error.
    pub(crate) fn next_member(
        &mut self,
        members: &mut Members,
        names: &[&str],
    ) -> Result<Option<usize>> {
        if !self.next_element(&mut members.block_start)? {
            self.end_block()?;
            return Ok(None);
        }
        let id = self.member_id_range()?;
        let text = self.bytes(id);
        let Some(index) = names.iter().position(|name| name.as_bytes() == text) else {
            let message = format!("unexpected member {}", String::from_utf8_lossy(text));
            return Err(self.error(message));
        };
        if members.last.is_some_and(|last| index <= last) {
            let message = format!("member {} out of order", String::from_utf8_lossy(text));
            return Err(self.error(message));
        }
        members.last = Some(index);
        Ok(Some(index))
    }

    /// The elements of a SEQUENCE OF or SET OF, each read by `element`.
    pub(crate) fn list<T>(
        &mut self,
        mut element: impl FnMut(&mut Self) -> Result<T>,
    ) -> Result<Vec<T>> {
        let mut out = Vec::new();
        self.each(|reader| {
            out.push(element(reader)?);
            Ok(())
        })?;
        Ok(out)
    }

    /// The elements of a SEQUENCE OF or SET OF, each read by `element` and kept by it.
    pub(crate) fn each(&mut self, mut element: impl FnMut(&mut Self) -> Result<()>) -> Result<()> {
        let mut block_start = self.start_block()?;
        while self.next_element(&mut block_start)? {
            element(self)?;
        }
        self.end_block()
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
    pub(crate) fn integer(&mut self, min: i128, max: i128) -> Result<i128> {
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
        loop {
            let c = self.peek_at(0);
            if !c.is_ascii_digit() {
                break;
            }
            n = n * 10 + i128::from(c - b'0');
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

    pub(crate) fn int4(&mut self) -> Result<i32> {
        Ok(self.integer(i32::MIN.into(), i32::MAX.into())? as i32)
    }

    /// `ReadUint4`: `GetUint4` after `SkipWhiteSpace`, an optional `+` and digits; a `-` is not a
    /// sign there (`'-' - '0'` is more than 9: "bad number", even for `-0`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:783-787
    // ```c++
    // Uint4 CObjectIStreamAsn::ReadUint4(void)
    // {
    //     SkipWhiteSpace();
    //     return m_Input.GetUint4();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:672-700
    // ```c++
    //     char c = SkipWs();
    //     if ( c == '+' )
    //         c = GetChar();
    //
    //     Uint4 n = c - '0';
    //     if ( n > 9 )
    //         BadNumber();
    //
    //     // overflow limits
    //     const Uint4 kMaxBeforeMul = kMax_UI4/10;
    //     const Uint1 kMaxLimitAdd = Uint1(kMax_UI4%10);
    //
    //     for ( ;; ) {
    //         Uint1 d = (Uint1)(PeekCharNoEOF() - '0');
    //         if  ( d > 9 )
    //             break;
    //         SkipChar();
    //
    //         // check multiplication overflow
    //         if ( n > kMaxBeforeMul || (n == kMaxBeforeMul && d > kMaxLimitAdd) )
    //             NumberOverflow();
    //
    //         n = n * 10 + d;
    //     }
    //     return n;
    // ```
    pub(crate) fn uint4(&mut self) -> Result<u32> {
        let mut c = self.skip_white_space()?;
        if c == b'+' {
            self.pos += 1;
            c = self.peek_at(0);
        }
        if !c.is_ascii_digit() {
            return Err(self.error("bad number"));
        }
        let mut n: u32 = 0;
        loop {
            let c = self.peek_at(0);
            if !c.is_ascii_digit() {
                return Ok(n);
            }
            self.pos += 1;
            n = n
                .checked_mul(10)
                .and_then(|n| n.checked_add(u32::from(c - b'0')))
                .ok_or_else(|| self.error("number overflow"))?;
        }
    }

    pub(crate) fn int8(&mut self) -> Result<i64> {
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
    pub(crate) fn boolean(&mut self) -> Result<bool> {
        let c = self.skip_white_space()?;
        if c == b'T'
            && self.peek_at(1) == b'R'
            && self.peek_at(2) == b'U'
            && self.peek_at(3) == b'E'
            && !is_id_char(self.peek_at(4))
        {
            self.pos += 4;
            return Ok(true);
        }
        if c == b'F'
            && self.peek_at(1) == b'A'
            && self.peek_at(2) == b'L'
            && self.peek_at(3) == b'S'
            && self.peek_at(4) == b'E'
            && !is_id_char(self.peek_at(5))
        {
            self.pos += 5;
            return Ok(false);
        }
        Err(self.error("TRUE or FALSE expected"))
    }

    /// `ReadNull`: the word `NULL`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:444-455
    // ```c++
    // void CObjectIStreamAsn::ReadNull(void)
    // {
    //     if ( SkipWhiteSpace() == 'N' &&
    //          m_Input.PeekCharNoEOF(1) == 'U' &&
    //          m_Input.PeekCharNoEOF(2) == 'L' &&
    //          m_Input.PeekCharNoEOF(3) == 'L' &&
    //          !IdChar(m_Input.PeekCharNoEOF(4)) ) {
    //         m_Input.SkipChars(4);
    //     }
    //     else
    //         ThrowError(fFormatError, "'NULL' expected");
    // }
    // ```
    pub(crate) fn null(&mut self) -> Result<()> {
        let c = self.skip_white_space()?;
        if c == b'N'
            && self.peek_at(1) == b'U'
            && self.peek_at(2) == b'L'
            && self.peek_at(3) == b'L'
            && !is_id_char(self.peek_at(4))
        {
            self.pos += 4;
            return Ok(());
        }
        Err(self.error("'NULL' expected"))
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
    pub(crate) fn enumerated(
        &mut self,
        values: &[(&'static str, i32)],
    ) -> Result<(&'static str, i32)> {
        let c = self.skip_white_space()?;
        let id = self.scan_end_of_id(c.is_ascii_lowercase());
        let found = if id.is_empty() {
            let value = self.int4()?;
            values.iter().find(|&&(_, number)| number == value)
        } else {
            let text = self.bytes(id);
            values.iter().find(|&&(name, _)| name.as_bytes() == text)
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
    pub(crate) fn string(&mut self) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.string_into(&mut out)?;
        Ok(out)
    }

    /// A UTF8String: `string` with any byte, as NCBI reads that type (`eFNP_Allow`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:988-991
    // ```c++
    // void CObjectIStreamAsn::ReadString(string& s, EStringType type)
    // {
    //     ReadStringValue(s, type == eStringTypeUTF8? eFNP_Allow: x_FixCharsMethod());
    // }
    // ```
    pub(crate) fn utf8_string(&mut self) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.string_value(&mut out, false)?;
        Ok(out)
    }

    /// `string`, appended to `out`.
    pub(crate) fn string_into(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.string_value(out, self.visible_only)
    }

    /// A string appended to `out`; `visible_only`: a byte outside the visible range is refused.
    fn string_value(&mut self, out: &mut Vec<u8>, visible_only: bool) -> Result<()> {
        self.expect(b'"')?;
        let start_line = self.line;
        loop {
            let Some(c) = self.src.byte_at(self.pos) else {
                return Err(AsnTextError {
                    line: self.line,
                    message: format!("unclosed string starts at line {start_line}"),
                    kind: AsnTextErrorKind::Format,
                });
            };
            self.pos += 1;
            match c {
                b'\r' | b'\n' => self.skip_end_of_line(c),
                b'"' => {
                    if self.peek_at(0) != b'"' {
                        return Ok(());
                    }
                    out.push(b'"');
                    self.pos += 1;
                }
                _ => {
                    // NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:456-459
                    // ```c++
                    // bool GoodVisibleChar(char c)
                    // {
                    //     return c >= ' ' && c <= '~';
                    // }
                    // ```
                    if visible_only && !(b' '..=b'~').contains(&c) {
                        let error = self.not_supported(&format!(
                            "a string with the byte 0x{c:x} (NCBI replaces it by '#' and posts an error)"
                        ));
                        self.refuse(error, ())?;
                    }
                    out.push(c)
                }
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
    pub(crate) fn octets(&mut self) -> Result<Vec<u8>> {
        self.expect(b'\'')?;
        let mut out = Vec::new();
        let mut high: Option<u8> = None;
        loop {
            let Some(c) = self.src.byte_at(self.pos) else {
                return Err(self.end_of_text());
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

    /// `ReadNumber`: a sign or a digit, then digits; returns their positions.
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
    pub(crate) fn number(&mut self) -> Result<Range<usize>> {
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
        Ok(start..self.pos)
    }

    /// Whether the bytes from `start` begin with `name`, ignoring case: `strncasecmp` over the
    /// stream's bytes, whose comparison stops at the first difference, so a byte past the
    /// text's end compares as the terminating 0.
    fn starts_with_ignore_case(&mut self, start: usize, name: &[u8]) -> bool {
        name.iter().enumerate().all(|(i, &expected)| {
            self.src
                .byte_at(start + i)
                .is_some_and(|byte| byte.eq_ignore_ascii_case(&expected))
        })
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
    pub(crate) fn real(&mut self) -> Result<f64> {
        if self.skip_white_space()? != b'{' {
            let start = self.pos;
            let word = self.scan_end_of_id(true);
            if self.starts_with_ignore_case(start, b"PLUS-INFINITY") {
                return Ok(f64::INFINITY);
            }
            if self.starts_with_ignore_case(start, b"MINUS-INFINITY") {
                return Ok(f64::NEG_INFINITY);
            }
            if self.starts_with_ignore_case(start, b"NOT-A-NUMBER") {
                return Ok(f64::NAN);
            }
            return Ok(string_to_double_posix(self.bytes(word), true).value);
        }
        self.pos += 1;
        let is_negative = self.skip_white_space()? == b'-';
        let mantissa_range = self.number()?;
        if mantissa_range.len() >= 128 - 1 {
            return Err(self.error("buffer overflow"));
        }
        let mantissa_length = mantissa_range.len();
        let read = string_to_double_posix(self.bytes(mantissa_range), false);
        let (mut mantissa, used) = (read.value, read.end);
        if used != mantissa_length {
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
    pub(crate) fn skip_member(&mut self, name: &'static str) -> Result<()> {
        self.skipped.push(name);
        self.skip_value()
    }

    /// Skips one value of any type: a block with what it holds, a string, an OCTET STRING or
    /// BIT STRING, or the words and numbers up to the `,` or `}` that ends it (a CHOICE is its
    /// variant's name and value). It checks the blocks, strings and octet strings of the value
    /// but not the grammar of its type: a value that NCBI's typed reader refuses can be skipped.
    pub(crate) fn skip_value(&mut self) -> Result<()> {
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
                        let Some(c) = self.src.byte_at(self.pos) else {
                            return Err(self.end_of_text());
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
}

/// glibc's `pow(10.0, exp)` (the C library NCBI 2.17.0 for x86-64 Linux runs with, version 2.28
/// on): the power of ten rounded to the nearest double, except at the exponents 23 and 210,
/// where glibc's result is the double above the nearest one (its error bound is 0.52 units in
/// the last place, and 10^23 and 10^210 lie almost halfway between two doubles). The table keeps
/// the bytes the same on every platform; `libm_pow10_is_glibcs_pow` checks it against the C
/// library on Linux. The powers from 10^-400 to 10^400 are computed once and kept (a stream of
/// archives reads hundreds of thousands of REALs); the others as they are asked for.
pub(crate) fn libm_pow10(exp: i32) -> f64 {
    static TABLE: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
    if (-400..=400).contains(&exp) {
        let table = TABLE.get_or_init(|| (-400..=400).map(libm_pow10_of).collect());
        return table[(exp + 400) as usize];
    }
    libm_pow10_of(exp)
}

/// `libm_pow10` computed.
fn libm_pow10_of(exp: i32) -> f64 {
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
