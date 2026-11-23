//! The lexical layer of NCBI's ASN.1 text object stream (`CObjectIStreamAsn` over
//! `CIStreamBuffer`), shared by LOSAT's readers of ASN.1 text: the generic reader of NCBI's type
//! information (`report::serial::asn_types`, which reads the PSSM files of `-in_pssm` from a byte
//! slice) and the Blast4-archive reader of `losat rpsbproc` (`rpsbproc::archive`, objects one after
//! another from a stream). The
//! grammar of each type stays with its reader; this module reads what every type is made of:
//! white space and comments, identifiers and `::=`, blocks and their elements, the members of a
//! SEQUENCE in the order of the module, INTEGER, BOOLEAN, ENUMERATED, VisibleString, OCTET
//! STRING and REAL values. For a hand-written grammar (rpsbproc's archive reader) it also gives
//! the frames of NCBI's object stack, which the grammar passes down to the values it reads
//! (`Frame`), and keeps the posts of the bytes that its VisibleStrings read as `#` (`Trail`).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243
//!
//! The text comes from an `AsnSource`: a slice (`SliceSource`, the PSSM file read whole) or a
//! stream read in blocks (`StreamSource`, an archive of any size: only the bytes not yet read
//! and those of the value being read are kept). The readers are generic over the source: both
//! read the kept bytes that `AsnSource::rest` gives and ask `AsnSource::fill` for more at their
//! end, so they share one code path.

use crate::ncbi_exception::ThrowSite;
use crate::ncbi_str::{is_space, string_to_double_posix, NCBI_NAN};
use crate::report::serial::encoder::{object_stack_trace, ReplacedChar, StackStep, StreamPosition};
use crate::rpsbproc::HOST_DEPENDENT_BYTES;
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
    /// The text ends where NCBI's input buffer needed one more byte (`PeekChar`, `GetChar`):
    /// its `CEofException`, which NCBI's object stream turns into a `CSerialException` inside an
    /// object, and passes on at the start of one (a reader of objects one after another takes it
    /// as the end of the input there).
    Eof,
}

/// A document that the reader cannot read: the line (from 1) where it stopped and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AsnTextError {
    pub line: usize,
    pub message: String,
    pub(crate) kind: AsnTextErrorKind,
    /// What NCBI's object stream throws there, where the reader names it (`None`: the reader
    /// does not; a reader that reports NCBI's exceptions, `writedb::mask_info`, builds them from
    /// this).
    pub(crate) ncbi: Option<NcbiThrow>,
}

/// What NCBI throws for an [`AsnTextError`]: the place of the throw and the form of its message.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NcbiThrow {
    /// `ThrowError` of the object stream with a flag and a message (`ThrowError1`,
    /// objistr.cpp:828-853): a `CSerialException` whose message is the stream's position,
    /// `line N: `, then the error's.
    Stream(ThrowSite),
    /// `NCBI_THROW` of a `CSerialException` with the error's message as it is (the values of an
    /// enumerated type: `FindValue`, `FindName`).
    Plain(ThrowSite),
    /// `BadNumber` / `NumberOverflow` of the input buffer: a `CUtilException`, the error's message
    /// and ` in line N`.
    Buffer(ThrowSite),
    /// The `CEofException` that `SkipWhiteSpace` catches (objistrasn.cpp:234-240): it goes on as
    /// it is at a stack depth of 2 or less, and becomes the `CSerialException` of
    /// [`SKIP_WHITE_SPACE_EOF`] deeper. Its kind is [`AsnTextErrorKind::Eof`].
    WhiteSpaceEof,
}

/// A `ThrowError` of `CObjectIStreamAsn` (objistrasn.cpp) at `line` of `function`.
const fn asn_site(line: u32, kind: &'static str, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/serial/objistrasn.cpp",
        line,
        kind,
        function: Some(function),
    }
}

const FORMAT_ERROR: &str = "CSerialException::eFormatError";

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243 (cited at the module and at
// `Reader::skip_white_space`)
/// `SkipWhiteSpace`'s `ThrowError` with `fEOF` and `e.what()` ([`NcbiThrow::WhiteSpaceEof`]).
pub(crate) const SKIP_WHITE_SPACE_EOF: ThrowSite = asn_site(
    239,
    "CSerialException::eEOF",
    "ncbi::CObjectIStreamAsn::SkipWhiteSpace",
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:177-185
const EXPECT: ThrowSite = asn_site(183, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::Expect");
const READ_TYPE_ID: ThrowSite = asn_site(299, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadTypeId");
const READ_MEMBER_ID: ThrowSite =
    asn_site(350, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadMemberId");
const READ_NULL: ThrowSite = asn_site(454, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadNull");
const READ_FILE_HEADER: ThrowSite =
    asn_site(719, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadFileHeader");
const READ_BOOL: ThrowSite = asn_site(762, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadBool");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:866-871
const UNENDED_STRING: ThrowSite =
    asn_site(869, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::UnendedString");
const NEXT_ELEMENT: ThrowSite =
    asn_site(1220, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::NextElement");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:315-329
const READ_NUMBER: ThrowSite = asn_site(319, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadNumber");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857
const READ_DOUBLE_OVERFLOW: ThrowSite = asn_site(
    822,
    "CSerialException::eOverflow",
    "ncbi::CObjectIStreamAsn::ReadDouble",
);
const READ_DOUBLE_BAD: ThrowSite =
    asn_site(828, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadDouble");
const READ_DOUBLE: ThrowSite = asn_site(836, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadDouble");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1482-1508
const GET_HEX_CHAR: ThrowSite = asn_site(1504, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::GetHexChar");
/// `UnexpectedMember` (objistrasn.cpp:1306-1314).
pub(crate) const UNEXPECTED_MEMBER: ThrowSite = asn_site(
    1313,
    FORMAT_ERROR,
    "ncbi::CObjectIStreamAsn::UnexpectedMember",
);
// NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:123-160
// (NCBI's report of these two names no function.)
pub(crate) const FIND_VALUE: ThrowSite = ThrowSite {
    file: "src/serial/enumerated.cpp",
    line: 133,
    kind: "CSerialException::eInvalidData",
    function: None,
};
/// `CEnumeratedTypeValues::FindName` of a number that names no value (enumerated.cpp:156): the
/// binary reader's too.
pub(crate) const FIND_NAME: ThrowSite = ThrowSite {
    file: "src/serial/enumerated.cpp",
    line: 156,
    kind: "CSerialException::eInvalidData",
    function: None,
};
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:570-584
/// `CIStreamBuffer::BadNumber`: the readers of numbers in ASN.1 text (this one and
/// `report::serial::asn_types::read`).
pub(crate) const BAD_NUMBER: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 574,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::BadNumber"),
};
/// `CIStreamBuffer::NumberOverflow`.
pub(crate) const NUMBER_OVERFLOW: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 583,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::NumberOverflow"),
};
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:420-436 (`FillBuffer` at the end of a
// file: `NCBI_THROW` of a `CEofException` with `eEof` and `m_Error`, the error "end of file")
/// The `CEofException` of NCBI's input buffer at the end of the input: what `SkipWhiteSpace`
/// quotes ([`SKIP_WHITE_SPACE_EOF`]).
pub(crate) const FILL_BUFFER_EOF: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 435,
    kind: "CEofException::eEof",
    function: Some("ncbi::CIStreamBuffer::FillBuffer"),
};

impl std::fmt::Display for AsnTextError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for AsnTextError {}

pub(crate) type Result<T> = std::result::Result<T, AsnTextError>;

/// The bytes that an ASN.1 text reader reads, by position from the start of the text.
///
/// A source gives the bytes it keeps as a slice (`rest`), and reads more only when a reader
/// asks for a byte past them (`fill`). The readers hold their position and that slice in
/// locals while they read a value, so the bytes of a value are read without loading the
/// reader's state from memory for each byte; `byte_at` and `scan` are built on these two, one
/// form for every source.
pub(crate) trait AsnSource {
    /// The bytes kept from `pos` on: empty when the byte at `pos` is not kept (yet). A position
    /// before one given to `release_before` is never asked for again.
    fn rest(&self, pos: usize) -> &[u8];
    /// Reads more of the text until the byte at `pos` is kept; `false` when the text ends
    /// before it.
    fn fill(&mut self, pos: usize) -> bool;
    /// The byte at `pos`, `None` past the end of the text.
    #[inline]
    fn byte_at(&mut self, pos: usize) -> Option<u8> {
        match self.rest(pos).first() {
            Some(&c) => Some(c),
            None => self.byte_after_fill(pos),
        }
    }
    /// `byte_at` of a byte not kept yet.
    #[cold]
    #[inline(never)]
    fn byte_after_fill(&mut self, pos: usize) -> Option<u8> {
        if !self.fill(pos) {
            return None;
        }
        self.rest(pos).first().copied()
    }
    /// The bytes from `pos` on for which `more` is true: the position of the first byte for
    /// which it is false and that byte, or the end of the text and `None`. The bytes are read
    /// from the kept slice in locals; the source is asked for more only at its end.
    #[inline]
    fn scan(&mut self, mut pos: usize, mut more: impl FnMut(u8) -> bool) -> (usize, Option<u8>) {
        loop {
            let rest = self.rest(pos);
            if let Some(i) = rest.iter().position(|&c| !more(c)) {
                return (pos + i, Some(rest[i]));
            }
            pos += rest.len();
            if !self.fill(pos) {
                return (pos, None);
            }
        }
    }
    /// The bytes `range`, all of which are kept.
    fn bytes(&self, range: Range<usize>) -> &[u8];
    /// The reader will not ask for a byte before `pos` again (a stream may drop them).
    #[inline]
    fn release_before(&mut self, _pos: usize) {}
    /// Keeps the bytes from `start` on (`None`: no more), whatever `release_before` says, until
    /// they are taken with `bytes` (the bytes of a delayed member, `Input::delay_start`).
    #[inline]
    fn hold(&mut self, _start: Option<usize>) {}
    /// Appends the `count` bytes from `pos` on to `out` (fewer when the text ends first) and
    /// returns how many it appended. `out` grows as the bytes come, so a length that the input
    /// does not hold costs no more than what NCBI reserves for it on every host, and a source may
    /// read a long value straight into `out` instead of keeping it as well.
    fn read_into(&mut self, pos: usize, count: u64, out: &mut Vec<u8>) -> u64 {
        read_kept(self, pos, count, out)
    }
}

/// `AsnSource::read_into` through the bytes the source keeps (`rest`, `fill`).
fn read_kept<S: AsnSource + ?Sized>(
    src: &mut S,
    mut pos: usize,
    count: u64,
    out: &mut Vec<u8>,
) -> u64 {
    let mut done = 0;
    while done < count {
        let rest = src.rest(pos);
        if rest.is_empty() {
            if !src.fill(pos) {
                break;
            }
            continue;
        }
        let take = rest
            .len()
            .min(usize::try_from(count - done).unwrap_or(usize::MAX));
        out.extend_from_slice(&rest[..take]);
        pos += take;
        done += take as u64;
    }
    done
}

/// A text held whole in memory.
pub(crate) struct SliceSource<'a> {
    text: &'a [u8],
}

impl<'a> SliceSource<'a> {
    pub(crate) fn new(text: &'a [u8]) -> Self {
        SliceSource { text }
    }
}

impl AsnSource for SliceSource<'_> {
    #[inline]
    fn rest(&self, pos: usize) -> &[u8] {
        &self.text[pos.min(self.text.len())..]
    }

    /// The whole text is kept: only a byte inside it is there.
    #[inline]
    fn fill(&mut self, pos: usize) -> bool {
        pos < self.text.len()
    }

    #[inline]
    fn bytes(&self, range: Range<usize>) -> &[u8] {
        &self.text[range]
    }
}

/// A text read from a stream in blocks. It keeps the bytes from the last position given to
/// `release_before` on (dropping the bytes before it, when it reads the next block, once they
/// make up most of the buffer), so its memory is about the largest value read between two
/// releases, not the whole text.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:125-139
// NCBI's object stream reads its input in blocks too and gives the bytes it did not use back
// to the stream when an object ends, so the next object is read from the byte after the last
// one. The reader of a stream of objects keeps one `StreamSource` for all of them.
pub(crate) struct StreamSource<R: Read> {
    input: R,
    buf: Vec<u8>,
    /// The position of `buf[0]`.
    base: usize,
    eof: bool,
    /// The first byte that `release_before` keeps (`AsnSource::hold`).
    held: Option<usize>,
    /// The last position given to `release_before`: the bytes before it may be dropped when the
    /// next block is read (`fill_to`), so a release, once per element, is a store.
    released: usize,
    /// The bytes that one fill reads: a block (`STREAM_BLOCK`). The tests also read one byte per
    /// fill (`fill_one_byte_at_a_time`), so that the kept bytes end at every byte of a text.
    fill_size: usize,
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
            held: None,
            released: 0,
            fill_size: STREAM_BLOCK,
        }
    }

    /// Reads one byte per fill from now on (the tests of the readers at the end of the kept
    /// bytes).
    #[cfg(test)]
    pub(crate) fn fill_one_byte_at_a_time(&mut self) {
        self.fill_size = 1;
    }

    /// The bytes kept.
    #[cfg(test)]
    pub(crate) fn buffered(&self) -> usize {
        self.buf.len()
    }

    /// Reads blocks until `pos` is in the buffer or the stream ends, after dropping the released
    /// bytes once they are most of a buffer of more than a few blocks (each byte is moved at
    /// most once on average). A read error ends the text, as a failed read ends an `istream` (its
    /// bad bit, which the object stream reads as the end of its input).
    #[cold]
    fn fill_to(&mut self, pos: usize) -> Option<u8> {
        assert!(pos >= self.base, "a released byte was asked for");
        let keep = self
            .held
            .map_or(self.released, |held| held.min(self.released));
        // (`read_into` may have dropped the bytes up to the end of a long value since the last
        // release: nothing more to drop then.)
        let dead = keep.saturating_sub(self.base);
        if dead > 4 * STREAM_BLOCK && dead * 2 > self.buf.len() {
            self.buf.drain(..dead);
            self.base = keep;
        }
        while pos - self.base >= self.buf.len() {
            if self.eof {
                return None;
            }
            // A block comes straight into the room after the kept bytes (`take` stops at its
            // end), not through a block of its own and a copy. Fewer bytes than `fill_size` is the
            // end of the stream; a read error ends the text, as in `read_into`.
            self.buf.reserve(self.fill_size);
            let before = self.buf.len();
            let failed = std::io::Read::by_ref(&mut self.input)
                .take(self.fill_size as u64)
                .read_to_end(&mut self.buf)
                .is_err();
            if failed || self.buf.len() - before < self.fill_size {
                self.eof = true;
            }
        }
        Some(self.buf[pos - self.base])
    }
}

impl<R: Read> AsnSource for StreamSource<R> {
    #[inline]
    fn rest(&self, pos: usize) -> &[u8] {
        // A released position wraps to past the buffer, so its slice is empty; `fill_to`'s assert
        // then rejects it.
        &self.buf[pos.wrapping_sub(self.base).min(self.buf.len())..]
    }

    #[inline]
    fn fill(&mut self, pos: usize) -> bool {
        self.fill_to(pos).is_some()
    }

    #[inline]
    fn bytes(&self, range: Range<usize>) -> &[u8] {
        &self.buf[range.start - self.base..range.end - self.base]
    }

    /// The bytes before `pos` are dropped when the next block is read (`fill_to`).
    #[inline]
    fn release_before(&mut self, pos: usize) {
        self.released = pos;
    }

    #[inline]
    fn hold(&mut self, start: Option<usize>) {
        self.held = start;
    }

    /// The kept bytes, then (when nothing is held) the rest straight from the stream into `out`:
    /// the bytes before the end of the value are not asked for again, so they are dropped, and a
    /// long value is not kept twice.
    fn read_into(&mut self, pos: usize, count: u64, out: &mut Vec<u8>) -> u64 {
        let rest = self.rest(pos);
        let kept = rest.len().min(usize::try_from(count).unwrap_or(usize::MAX));
        out.extend_from_slice(&rest[..kept]);
        let mut done = kept as u64;
        if done == count || self.held.is_some() || self.eof {
            if done < count {
                done += read_kept(self, pos + kept, count - done, out);
            }
            return done;
        }
        let pos = pos + kept;
        self.buf.clear();
        self.base = pos;
        // NCBI reserves the count before it reads the bytes (`asn_input::reserve`); below the
        // allocations that depend on the host, so does LOSAT, and the bytes come into one block of
        // that size, not into one that doubles as it fills.
        if let Ok(wanted) = usize::try_from(count - done) {
            if (wanted as u64) < HOST_DEPENDENT_BYTES {
                // A reserve that fails leaves `out` to grow as the bytes come.
                let _ = out.try_reserve_exact(wanted);
            }
        }
        let before = out.len();
        // `take` stops at `count`; a read error ends the text, as in `fill_to`.
        let ended = std::io::Read::by_ref(&mut self.input)
            .take(count - done)
            .read_to_end(out)
            .is_err();
        let read = (out.len() - before) as u64;
        done += read;
        self.base = pos + read as usize;
        if ended || done < count {
            self.eof = true;
        }
        done
    }
}

/// `Reader::white_space_from`'s byte at the end of the text: more than any byte.
const END_OF_TEXT: u32 = 0x100;

/// The index of the first byte of `bytes` from `i` (at most its length) on that is not `' '`,
/// or its length when the bytes from `i` on are all spaces.
#[inline]
fn after_spaces(bytes: &[u8], mut i: usize) -> usize {
    while bytes.get(i) == Some(&b' ') {
        i += 1;
    }
    i
}

/// `FirstIdChar` (objistrasn.cpp:118-121).
pub(crate) fn is_first_id_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:123-127
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
    /// What the text is, for the error at its end ("unexpected end of the PSSM").
    document: &'static str,
    /// The reader, for the refusals ("... is not supported by LOSAT's PSSM reader").
    reader: &'static str,
    /// The hand grammar's root and posts (`Trail`).
    pub(crate) trail: Trail,
}

/// What a hand-written grammar over this reader (rpsbproc's archive reader, `rpsbproc::archive`,
/// with its Seq-ids, `seqid::asn_text`) reads beside the values: the name of the document's type,
/// at the root of its frames (`Frame`), and the bytes that its VisibleStrings read as `#`
/// (`string`), each with NCBI's post.
#[derive(Debug, Default)]
pub(crate) struct Trail {
    /// The name of the document's type, the first name of a stack trace (`Blast4-archive`).
    pub(crate) root: &'static str,
    /// The posts of the bytes read as `#`, in order (`ReplaceVisibleChar`'s 803.7).
    pub(crate) replaced: Vec<ReplacedChar>,
}

impl Trail {
    /// The post of a byte that a VisibleString of the hand grammar held by the frame `at` read as
    /// `#` (`Reader::string`): its line, the byte and the piece that holds it, with the stack trace
    /// of `at`.
    #[cold]
    #[inline(never)]
    fn post_replaced(&mut self, at: &Frame<'_>, line: usize, byte: u8, piece: &[u8]) {
        self.replaced.push(ReplacedChar {
            trace: object_stack_trace(self.root, at.steps()),
            byte,
            position: StreamPosition::Line(line as u64),
            string: piece.into(),
        });
    }
}

/// A frame of NCBI's object stack as a hand grammar passes it down to the reader of the value
/// that it holds: the member of a SEQUENCE being read (`Members::at`), the variant of a CHOICE
/// (`Reader::choice`) or the element of a SEQUENCE OF or SET OF (`List::at`), each pointing to the
/// frame that holds it, down to the object's root (`Frame::ROOT`). NCBI reads every value of its
/// object with its typed reader, whose stack trace (`GetStackTrace()`) names the member or variant
/// of each frame and `E` for a container's element; the trace of a frame is made only when a post
/// needs it (`Reader::string`) or the typed reader reads below it (`Frame::steps`,
/// `asn_types::read::ObjectStream::read_value`). So a member costs only the position of the
/// member read, which orders the next (the SEQUENCE's frame takes its parent and names once), an
/// element nothing (its list's frame is made once), and a variant's frame is stored only when a
/// reader below it is given it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Frame<'a> {
    /// The frame that holds this one; `None` at the root.
    parent: Option<&'a Frame<'a>>,
    /// The names of the SEQUENCE's members or of the CHOICE's variants in module order; none for
    /// a container's element and the root.
    names: &'static [&'static str],
    /// The position in `names` of the member or variant read, plus one (0 before the first
    /// member).
    next: usize,
}

impl Frame<'static> {
    /// The object's root, the frame of the document's type (`Trail::root`), which names no step.
    pub(crate) const ROOT: Frame<'static> = Frame {
        parent: None,
        names: &[],
        next: 0,
    };
}

impl<'a> Frame<'a> {
    /// The frame of an element of the list that `at` holds (`E`).
    #[inline]
    fn element(at: &'a Frame<'a>) -> Self {
        Frame {
            parent: Some(at),
            names: &[],
            next: 0,
        }
    }

    /// The position of the variant (`Reader::choice`).
    #[inline]
    pub(crate) fn index(&self) -> usize {
        self.next.wrapping_sub(1)
    }

    /// The name of the member or variant read (`None` for an element, the root, and a SEQUENCE
    /// before its first member).
    pub(crate) fn name(&self) -> Option<&'static str> {
        self.names.get(self.index()).copied()
    }

    /// The steps from the root to this frame: a member or variant name, or `E` for an element.
    #[cold]
    #[inline(never)]
    pub(crate) fn steps(&self) -> Vec<StackStep<'static>> {
        let mut steps = Vec::new();
        let mut frame = self;
        while let Some(parent) = frame.parent {
            if frame.names.is_empty() {
                steps.push(StackStep::Element);
            } else if let Some(name) = frame.name() {
                steps.push(StackStep::Member(name));
            }
            frame = parent;
        }
        steps.reverse();
        steps
    }
}

/// What `ReadMemberId` read: a name, or the digits of a tag in brackets (`[3]`).
pub(crate) enum MemberId {
    Name(Range<usize>),
    Tag(Range<usize>),
}

/// What `ReadEnum` reads ([`Reader::enum_word`]).
#[derive(Debug)]
pub(crate) enum EnumWord {
    /// An identifier: its positions in the text.
    Name(Range<usize>),
    /// A number.
    Number(i32),
}

/// `CObjectIStream::EPointerType`: what `ReadPointerType` finds where NCBI's object stream reads a
/// `CRef<>` that is empty at that moment (`CObjectIStream::ReadPointer`, objistr.cpp:1133-1213).
/// The text stream's is [`Reader::pointer_type`]; the binary one's is the mask reader's
/// (`writedb::mask_info::ber`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PointerType {
    /// `eNullPointer`: `NULL` (text), `05 00` (binary); the `CRef<>` stays empty.
    Null,
    /// `eObjectPointer`: `@` and an object index (text), the application tag 30 (binary).
    Object,
    /// `eOtherPointer`: `:` and a type name (text), the application long tag with the name
    /// (binary); the value of that type follows.
    Other,
    /// `eThisPointer`: anything else, the value of the declared type.
    This,
}

/// The members of a SEQUENCE as they are read: each name in turn, in the order of the module.
pub(crate) struct Members<'a> {
    block_start: bool,
    /// The frame of the member being read: the SEQUENCE's member names and the position of the
    /// last member read, which the next must follow.
    frame: Frame<'a>,
}

impl<'a> Members<'a> {
    /// The frame of the member being read, for the reader of its value.
    #[inline]
    pub(crate) fn at(&self) -> &Frame<'a> {
        &self.frame
    }
}

/// The elements of a SEQUENCE OF or SET OF as a hand grammar reads them (`Reader::begin_list`).
pub(crate) struct List<'a> {
    block_start: bool,
    /// The frame of the elements (`E`).
    frame: Frame<'a>,
}

impl<'a> List<'a> {
    /// The frame of the element being read, for the reader of its value.
    #[inline]
    pub(crate) fn at(&self) -> &Frame<'a> {
        &self.frame
    }
}

impl<S: AsnSource> Reader<S> {
    /// A reader at the start of `src`; `document` and `reader` name the text and the reader in
    /// error messages.
    pub(crate) fn new(src: S, document: &'static str, reader: &'static str) -> Self {
        Reader {
            src,
            pos: 0,
            line: 1,
            document,
            reader,
            trail: Trail::default(),
        }
    }

    #[cold]
    pub(crate) fn error(&self, message: impl Into<String>) -> AsnTextError {
        AsnTextError {
            line: self.line,
            message: message.into(),
            kind: AsnTextErrorKind::Format,
            ncbi: None,
        }
    }

    /// `error`, with what NCBI throws there.
    #[cold]
    pub(crate) fn error_at(&self, ncbi: NcbiThrow, message: impl Into<String>) -> AsnTextError {
        AsnTextError {
            ncbi: Some(ncbi),
            ..self.error(message)
        }
    }

    /// An explicit LOSAT refusal at the current line.
    #[cold]
    pub(crate) fn unsupported(&self, message: impl Into<String>) -> AsnTextError {
        AsnTextError {
            line: self.line,
            message: message.into(),
            kind: AsnTextErrorKind::Unsupported,
            ncbi: None,
        }
    }

    /// "`<what>` is not supported by LOSAT's `<reader>`".
    pub(crate) fn not_supported(&self, what: &str) -> AsnTextError {
        self.unsupported(format!(
            "{what} is not supported by LOSAT's {}",
            self.reader
        ))
    }

    #[cold]
    #[inline(never)]
    pub(crate) fn end_of_text(&self) -> AsnTextError {
        AsnTextError {
            line: self.line,
            message: format!("unexpected end of the {}", self.document),
            kind: AsnTextErrorKind::Eof,
            ncbi: None,
        }
    }

    /// `end_of_text` inside `SkipWhiteSpace` ([`NcbiThrow::WhiteSpaceEof`]).
    #[cold]
    #[inline(never)]
    fn end_of_text_in_white_space(&self) -> AsnTextError {
        AsnTextError {
            ncbi: Some(NcbiThrow::WhiteSpaceEof),
            ..self.end_of_text()
        }
    }

    /// The byte at `pos + offset` (`PeekChar`): the end of the text is [`AsnTextErrorKind::Eof`].
    #[inline]
    pub(crate) fn peek_char(&mut self, offset: usize) -> Result<u8> {
        match self.src.byte_at(self.pos + offset) {
            Some(c) => Ok(c),
            None => Err(self.end_of_text()),
        }
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

    /// `SkipEndOfLine` of `c` after the end-of-line byte `c` was taken.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-537
    pub(crate) fn skip_end_of_line(&mut self, c: u8) {
        self.line += 1;
        if c as u32 + self.peek_at(0) as u32 == b'\r' as u32 + b'\n' as u32 {
            self.pos += 1;
        }
    }

    /// `SkipWhiteSpace`: the next byte that is not white space or in a comment, not taken. The
    /// end of the text is an error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:245-274
    #[inline]
    pub(crate) fn skip_white_space(&mut self) -> Result<u8> {
        let (pos, c) = self.white_space_at(self.pos)?;
        self.pos = pos;
        Ok(c)
    }

    // The readers below that end in `_at` read from the position they are given and return the
    // position after what they read, so that a reader of a value keeps its position in a local
    // from its start to its end and writes `pos` once. On an error they set `pos` to where the
    // reading stopped (the position the error is reported at) themselves.

    /// [`skip_white_space`](Self::skip_white_space) from `pos`: the position of the byte found,
    /// and the byte.
    #[inline]
    pub(crate) fn white_space_at(&mut self, pos: usize) -> Result<(usize, u8)> {
        // A value or a punctuation mark at the position: nothing to skip. A byte not kept yet is
        // read by `white_space_from`.
        if let Some(&c) = self.src.rest(pos).first() {
            if !matches!(c, b' ' | b'\t' | b'\r' | b'\n' | b'-') {
                return Ok((pos, c));
            }
        }
        match self.white_space_from(pos) {
            (pos, END_OF_TEXT) => {
                self.pos = pos;
                Err(self.end_of_text_in_white_space())
            }
            (pos, c) => Ok((pos, c as u8)),
        }
    }

    /// [`white_space_at`](Self::white_space_at) from a byte that may be white space, a comment,
    /// not kept yet or past the end of the text: the position of the byte found and the byte, or
    /// the position where the text ends and [`END_OF_TEXT`]. (A pair of numbers, which a call
    /// returns in registers.)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-274
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:213-253
    #[inline(never)]
    fn white_space_from(&mut self, mut pos: usize) -> (usize, u32) {
        // The line ends found, added to `line` once (when there are any).
        let mut lines = 0;
        loop {
            // Spaces, tabs and line ends are read from the kept bytes in locals; the source is
            // asked for more only at their end.
            let rest = self.src.rest(pos);
            let mut i = 0;
            // What ends the white space in `rest` when it needs another look: the end of `rest`
            // (`None`), a `-` that starts a comment or is the last byte of `rest`, or a line end
            // that is the last byte of `rest` (its complement may be in the next block).
            let c = loop {
                // A run of spaces (the indentation of a text written with line ends) in one
                // tight loop.
                i = after_spaces(rest, i);
                let Some(&c) = rest.get(i) else {
                    break None;
                };
                if !matches!(c, b'\t' | b'\r' | b'\n' | b'-') {
                    // The byte found.
                    self.add_lines(lines);
                    return (pos + i, c.into());
                }
                match c {
                    b'\t' => i += 1,
                    b'-' => match rest.get(i + 1) {
                        // A `-` that does not start a comment (a negative number).
                        Some(&next) if next != b'-' => {
                            self.add_lines(lines);
                            return (pos + i, b'-'.into());
                        }
                        _ => break Some(c),
                    },
                    // A line end.
                    _ => {
                        i += 1;
                        lines += 1;
                        match rest.get(i) {
                            Some(&next) => {
                                if c as u32 + next as u32 == b'\r' as u32 + b'\n' as u32 {
                                    i += 1;
                                }
                            }
                            None => break Some(c),
                        }
                    }
                }
            };
            pos += i;
            match c {
                None => {
                    if !self.src.fill(pos) {
                        break;
                    }
                }
                Some(b'-') => match self.src.byte_at(pos + 1) {
                    // The end of the text after a `-`: `PeekChar` one byte ahead fails at the `-`.
                    None => break,
                    Some(b'-') => match self.comment_at(pos + 2) {
                        Some((end, end_lines)) => (pos, lines) = (end, lines + end_lines),
                        None => {
                            // `comment_at` set `pos` where the text ends.
                            self.add_lines(lines);
                            return (self.pos, END_OF_TEXT);
                        }
                    },
                    Some(_) => {
                        self.add_lines(lines);
                        return (pos, b'-'.into());
                    }
                },
                Some(c) => pos = self.end_of_line_at(pos, c),
            }
        }
        self.add_lines(lines);
        (pos, END_OF_TEXT)
    }

    /// Adds `lines` line ends to `line` (no access to `line` when there are none).
    #[inline]
    fn add_lines(&mut self, lines: usize) {
        if lines != 0 {
            self.line += lines;
        }
    }

    /// The rest of a comment from `pos` (after its `--`): up to the end of the line or the next
    /// `--`. The position after it and the line ends taken (0 or 1), or `None` at the end of the
    /// text (`pos` set there).
    #[cold]
    #[inline(never)]
    fn comment_at(&mut self, mut pos: usize) -> Option<(usize, usize)> {
        loop {
            let (at, c) = self.src.scan(pos, |c| !matches!(c, b'\r' | b'\n' | b'-'));
            pos = at;
            let Some(mut c) = c else {
                self.pos = pos;
                return None;
            };
            pos += 1;
            if c == b'-' {
                let Some(next) = self.src.byte_at(pos) else {
                    self.pos = pos;
                    return None;
                };
                pos += 1;
                if next == b'-' {
                    return Some((pos, 0));
                }
                c = next;
            }
            if matches!(c, b'\r' | b'\n') {
                return Some((self.end_of_line_at(pos, c), 1));
            }
        }
    }

    /// `SkipEndOfLine` of `c` at `pos`, after the end-of-line byte `c`, without its line count:
    /// the position after the complement of `c` when it follows.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-537
    #[inline]
    fn end_of_line_at(&mut self, pos: usize, c: u8) -> usize {
        let next = self.src.byte_at(pos).unwrap_or(0);
        if c as u32 + next as u32 == b'\r' as u32 + b'\n' as u32 {
            return pos + 1;
        }
        pos
    }

    /// `Expect` of `c`, with its flag true.
    pub(crate) fn expect(&mut self, c: u8) -> Result<()> {
        self.pos = self.expect_at(self.pos, c)?;
        Ok(())
    }

    /// [`expect`](Self::expect) from `pos`.
    #[inline]
    fn expect_at(&mut self, pos: usize, c: u8) -> Result<usize> {
        let (pos, found) = self.white_space_at(pos)?;
        if found != c {
            self.pos = pos;
            return Err(self.expected(c));
        }
        Ok(pos + 1)
    }

    /// `Expect`'s error.
    #[cold]
    #[inline(never)]
    fn expected(&self, c: u8) -> AsnTextError {
        self.error_at(
            NcbiThrow::Stream(EXPECT),
            format!("'{}' expected", c as char),
        )
    }

    /// `ScanEndOfId` with `is_id` at the current byte: the byte, then the identifier bytes and the
    /// `-` before one; empty when `is_id` is false. Returns the positions of the identifier.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:276-290
    pub(crate) fn scan_end_of_id(&mut self, is_id: bool) -> Result<Range<usize>> {
        let start = self.pos;
        if !is_id {
            return Ok(start..start);
        }
        self.pos = self.end_of_id_at(start)?;
        Ok(start..self.pos)
    }

    /// `ScanEndOfId` with `true` at `start`: the end of the identifier that starts there.
    #[inline]
    fn end_of_id_at(&mut self, start: usize) -> Result<usize> {
        let mut pos = start + 1;
        loop {
            let (at, c) = self.src.scan(pos, is_id_char);
            if c != Some(b'-') {
                return Ok(at);
            }
            // A `-` belongs to the identifier when an identifier byte follows it.
            match self.src.byte_at(at + 1) {
                Some(next) if is_id_char(next) => pos = at + 1,
                Some(_) => return Ok(at),
                None => {
                    self.pos = start;
                    return Err(self.end_of_text());
                }
            }
        }
    }

    /// `ReadMemberId` after `SkipWhiteSpace`: a member or variant name, bare (from a lower-case
    /// letter) or between `[` and `]` on one line (the positions of the text between them).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:343-363
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:394-397
    // An id in brackets that starts with a digit names the member or variant by its tag
    // (`GetMemberIndex`, and `GetChoiceIndex` at :431-433), which depends on the tags NCBI's
    // generated type information gives each member: `member_or_tag` returns it for a reader that
    // knows them, and `member_id_range`, for the readers that know the members by name only,
    // refuses that form explicitly.
    pub(crate) fn member_id_range(&mut self) -> Result<Range<usize>> {
        match self.member_or_tag()? {
            MemberId::Name(id) => Ok(id),
            MemberId::Tag(_) => Err(self.not_supported("a member or variant named by its tag")),
        }
    }

    /// `ReadMemberId` after `SkipWhiteSpace`, a name or a tag (`member_id_range`).
    pub(crate) fn member_or_tag(&mut self) -> Result<MemberId> {
        let (pos, id) = self.member_or_tag_at(self.pos)?;
        self.pos = pos;
        Ok(id)
    }

    /// [`member_or_tag`](Self::member_or_tag) from `pos`.
    #[inline]
    fn member_or_tag_at(&mut self, pos: usize) -> Result<(usize, MemberId)> {
        let (pos, c) = self.white_space_at(pos)?;
        if c != b'[' {
            let end = if c.is_ascii_lowercase() {
                self.end_of_id_at(pos)?
            } else {
                pos
            };
            return Ok((end, MemberId::Name(pos..end)));
        }
        self.bracketed_member_at(pos)
    }

    /// The member id in brackets at `pos` (its `[`).
    #[inline(never)]
    fn bracketed_member_at(&mut self, pos: usize) -> Result<(usize, MemberId)> {
        let (at, c) = self
            .src
            .scan(pos + 1, |c| !matches!(c, b'\r' | b'\n' | b']'));
        match c {
            Some(b']') => {
                let id = pos + 1..at;
                if self
                    .bytes(id.clone())
                    .first()
                    .is_some_and(u8::is_ascii_digit)
                {
                    return Ok((at + 1, MemberId::Tag(id)));
                }
                Ok((at + 1, MemberId::Name(id)))
            }
            Some(_) => {
                self.pos = pos;
                Err(self.error_at(
                    NcbiThrow::Stream(READ_MEMBER_ID),
                    "end of line: expected ']'",
                ))
            }
            None => {
                self.pos = pos;
                Err(self.end_of_text())
            }
        }
    }

    /// `ReadFileHeader`: the type name before `::=` (empty when there is none), the `::=` taken.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:292-313,709-721
    // A type name in brackets stops before its `]` (`SkipChars` by `i`), so the `::=` test that
    // follows always fails.
    pub(crate) fn file_header(&mut self) -> Result<Range<usize>> {
        let c = self.skip_white_space()?;
        let id = if c == b'[' {
            let mut i = 1;
            loop {
                match self.peek_char(i)? {
                    b'\r' | b'\n' => {
                        return Err(self.error_at(
                            NcbiThrow::Stream(READ_TYPE_ID),
                            "end of line: expected ']'",
                        ))
                    }
                    b']' => {
                        return Err(
                            self.error_at(NcbiThrow::Stream(READ_FILE_HEADER), "'::=' expected")
                        )
                    }
                    _ => i += 1,
                }
            }
        } else {
            self.scan_end_of_id(is_first_id_char(c))?
        };
        if self.skip_white_space()? != b':' || self.peek_at(1) != b':' || self.peek_at(2) != b'=' {
            return Err(self.error_at(NcbiThrow::Stream(READ_FILE_HEADER), "'::=' expected"));
        }
        self.pos += 3;
        Ok(id)
    }

    /// `ReadTypeId` after `SkipWhiteSpace` (cited at `file_header`), as `ReadOtherPointer` reads
    /// the type name after a `:` (objistrasn.cpp:1643-1646): a name from a letter or `_` (empty
    /// otherwise), or the text between `[` and `]` on one line, which NCBI takes without its last
    /// byte (a `CTempString` that starts at `ptr` plus 1 and is `i` less 2 bytes long), leaving
    /// the `]` unread (`SkipChars` by `i`). `None` for `[]`, whose length, `i` less 2, wraps to the
    /// largest `size_t`.
    pub(crate) fn type_id(&mut self) -> Result<Option<Range<usize>>> {
        let c = self.skip_white_space()?;
        if c != b'[' {
            return Ok(Some(self.scan_end_of_id(is_first_id_char(c))?));
        }
        let mut i = 1;
        loop {
            match self.peek_char(i)? {
                b'\r' | b'\n' => {
                    return Err(
                        self.error_at(NcbiThrow::Stream(READ_TYPE_ID), "end of line: expected ']'")
                    )
                }
                b']' => {
                    let start = self.pos;
                    self.pos += i;
                    return Ok(if i >= 2 {
                        Some(start + 1..start + i - 1)
                    } else {
                        None
                    });
                }
                _ => i += 1,
            }
        }
    }

    /// The variant of a CHOICE whose variants are `names` (in module order), held by the frame
    /// `at`: the variant's frame, whose `index` is its position (`variant`). (Inlined, so that a
    /// frame that no reader below it is given is never stored.)
    #[inline]
    pub(crate) fn choice<'a>(
        &mut self,
        at: &'a Frame<'a>,
        names: &'static [&'static str],
    ) -> Result<Frame<'a>> {
        let index = self.variant(names)?;
        Ok(Frame {
            parent: Some(at),
            names,
            next: index + 1,
        })
    }

    /// The position of the variant of a CHOICE whose variants are `names` (in module order); a
    /// name that is not one of them is an error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1475
    // (`NextElement` at the choice's start fails at a `}`, which reads as no id here; the variant
    // names of the modules LOSAT reads are their ASN.1 names, `GetChoiceIndex` (:426-442).)
    fn variant(&mut self, names: &[&str]) -> Result<usize> {
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
    // (Called once per element by the loops over elements, into which it is inlined.)
    #[inline]
    pub(crate) fn next_element(&mut self, block_start: &mut bool) -> Result<bool> {
        let (pos, more) = self.next_element_at(self.pos, block_start)?;
        self.pos = pos;
        Ok(more)
    }

    /// [`next_element`](Self::next_element) from `pos`.
    #[inline]
    fn next_element_at(&mut self, pos: usize, block_start: &mut bool) -> Result<(usize, bool)> {
        self.src.release_before(pos);
        let (pos, c) = self.white_space_at(pos)?;
        if *block_start {
            *block_start = false;
            return Ok((pos, c != b'}'));
        }
        if c == b',' {
            return Ok((pos + 1, true));
        }
        if c != b'}' {
            self.pos = pos;
            return Err(self.error_at(NcbiThrow::Stream(NEXT_ELEMENT), "',' or '}' expected"));
        }
        Ok((pos, false))
    }

    /// `EndBlock`: the `}` that `next_element` stopped at.
    pub(crate) fn end_block(&mut self) -> Result<()> {
        self.expect(b'}')
    }

    /// `StartBlock` of a SEQUENCE of a hand grammar whose members are `names` (in module order),
    /// held by the frame `at`; `next_member` reads the members.
    #[inline]
    pub(crate) fn begin_members<'a>(
        &mut self,
        at: &'a Frame<'a>,
        names: &'static [&'static str],
    ) -> Result<Members<'a>> {
        Ok(Members {
            block_start: self.start_block()?,
            frame: Frame {
                parent: Some(at),
                names,
                next: 0,
            },
        })
    }

    /// The next member of the SEQUENCE, as its position (its frame is `members.at()`); `None` at
    /// the end of the SEQUENCE (its `}` taken). A member out of order or not in the module is an
    /// error.
    pub(crate) fn next_member(&mut self, members: &mut Members) -> Result<Option<usize>> {
        let (pos, more) = self.next_element_at(self.pos, &mut members.block_start)?;
        if !more {
            self.pos = self.expect_at(pos, b'}')?;
            return Ok(None);
        }
        let (pos, id) = self.member_or_tag_at(pos)?;
        self.pos = pos;
        let MemberId::Name(id) = id else {
            return Err(self.not_supported("a member or variant named by its tag"));
        };
        let text = self.bytes(id);
        let Some(index) = members
            .frame
            .names
            .iter()
            .position(|name| name.as_bytes() == text)
        else {
            let message = format!("unexpected member {}", String::from_utf8_lossy(text));
            return Err(self.error(message));
        };
        if index < members.frame.next {
            let message = format!("member {} out of order", String::from_utf8_lossy(text));
            return Err(self.error(message));
        }
        members.frame.next = index + 1;
        Ok(Some(index))
    }

    /// `StartBlock` of a SEQUENCE OF or SET OF of a hand grammar, held by the frame `at`;
    /// `next_item` reads the elements.
    #[inline]
    pub(crate) fn begin_list<'a>(&mut self, at: &'a Frame<'a>) -> Result<List<'a>> {
        Ok(List {
            block_start: self.start_block()?,
            frame: Frame::element(at),
        })
    }

    /// `NextElement` of the list (`next_element`; its frame is `list.at()`).
    #[inline]
    pub(crate) fn next_item(&mut self, list: &mut List) -> Result<bool> {
        self.next_element(&mut list.block_start)
    }

    /// The elements of a SEQUENCE OF or SET OF held by the frame `at`, each read by `element`
    /// (given its frame) and kept by it.
    pub(crate) fn each(
        &mut self,
        at: &Frame<'_>,
        mut element: impl FnMut(&mut Self, &Frame<'_>) -> Result<()>,
    ) -> Result<()> {
        let mut list = self.begin_list(at)?;
        while self.next_item(&mut list)? {
            element(self, list.at())?;
        }
        self.end_block()
    }

    /// `GetInt4`, `GetUint4` (`signed` false: only `+`, a `-` is "bad number" even for `-0`) and
    /// `GetInt8` after `SkipWhiteSpace` (`ReadInt4`, `ReadUint4`, `ReadInt8`): the digits checked
    /// against `max` as NCBI checks them (the digit that overflows is taken first). The one port
    /// of NCBI's text numbers, for the archive reader and the object streams (`asn_input`).
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-668
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:783-787
    pub(crate) fn get_number(&mut self, max: u64, signed: bool) -> Result<i64> {
        // The position is kept in a local and written once (`white_space_at`).
        let (pos, c) = self.white_space_at(self.pos)?;
        // `SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
        let (pos, mut c) = self.ws_at(pos, c)?;
        let mut pos = pos + 1;
        let mut negative = false;
        if (signed && c == b'-') || c == b'+' {
            negative = c == b'-';
            c = self.char_at(pos)?;
            pos += 1;
        }
        if !c.is_ascii_digit() {
            self.pos = pos;
            return Err(self.bad_number());
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        let mut overflow = false;
        let (end, _) = self.src.scan(pos, |d| {
            if !d.is_ascii_digit() {
                return false;
            }
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                overflow = true;
                return false;
            }
            n = n * 10 + d;
            true
        });
        if overflow {
            // The digit that overflows is taken before the check.
            self.pos = end + 1;
            return Err(self.number_overflow());
        }
        self.pos = end;
        // `-Int4(n)` of n = 2^31 is the least Int4 (and likewise in 64 bits).
        Ok(if negative {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    /// `PeekChar()` at `pos`: the end of the text is [`AsnTextErrorKind::Eof`] at `pos`.
    #[inline]
    pub(crate) fn char_at(&mut self, pos: usize) -> Result<u8> {
        match self.src.byte_at(pos) {
            Some(c) => Ok(c),
            None => {
                self.pos = pos;
                Err(self.end_of_text())
            }
        }
    }

    /// `BadNumber` of the input buffer.
    #[cold]
    #[inline(never)]
    fn bad_number(&self) -> AsnTextError {
        self.error_at(NcbiThrow::Buffer(BAD_NUMBER), "bad number")
    }

    /// `NumberOverflow` of the input buffer.
    #[cold]
    #[inline(never)]
    fn number_overflow(&self) -> AsnTextError {
        self.error_at(NcbiThrow::Buffer(NUMBER_OVERFLOW), "number overflow")
    }

    /// `CIStreamBuffer::SkipWs` before the first byte of a number, which `GetInt4` and the other
    /// readers of numbers call after `SkipWhiteSpace`: any C `isspace` byte (also `\v` and `\f`,
    /// which `SkipWhiteSpace` stops at) is skipped, line ends uncounted, comments not seen. Returns
    /// the next byte, not taken.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:617-624
    // `c` is the byte at `pos`, which `SkipWhiteSpace` found.
    #[inline]
    pub(crate) fn ws_at(&mut self, pos: usize, c: u8) -> Result<(usize, u8)> {
        if !is_space(c) {
            return Ok((pos, c));
        }
        self.spaces_from(pos + 1)
    }

    /// `ws_at` past a first `isspace` byte.
    #[inline(never)]
    fn spaces_from(&mut self, pos: usize) -> Result<(usize, u8)> {
        match self.src.scan(pos, is_space) {
            (pos, Some(c)) => Ok((pos, c)),
            (pos, None) => {
                self.pos = pos;
                Err(self.end_of_text())
            }
        }
    }

    /// `ReadInt4` ([`get_number`](Self::get_number)).
    pub(crate) fn int4(&mut self) -> Result<i32> {
        Ok(self.get_number(i32::MAX as u64, true)? as i32)
    }

    /// `ReadUint4` ([`get_number`](Self::get_number)).
    pub(crate) fn uint4(&mut self) -> Result<u32> {
        Ok(self.get_number(u64::from(u32::MAX), false)? as u32)
    }

    /// `ReadInt8` ([`get_number`](Self::get_number)).
    pub(crate) fn int8(&mut self) -> Result<i64> {
        self.get_number(i64::MAX as u64, true)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:739-764
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
        Err(self.error_at(NcbiThrow::Stream(READ_BOOL), "TRUE or FALSE expected"))
    }

    /// `ReadNull`: the word `NULL`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:444-455
    pub(crate) fn null(&mut self) -> Result<()> {
        let c = self.skip_white_space()?;
        if self.take_null_word(c) {
            return Ok(());
        }
        Err(self.error_at(NcbiThrow::Stream(READ_NULL), "'NULL' expected"))
    }

    /// Whether the word `NULL`, not followed by an identifier byte, starts at the current byte `c`
    /// (`ReadNull`, `ReadPointerType`); it is taken when it does.
    fn take_null_word(&mut self, c: u8) -> bool {
        if c == b'N'
            && self.peek_at(1) == b'U'
            && self.peek_at(2) == b'L'
            && self.peek_at(3) == b'L'
            && !is_id_char(self.peek_at(4))
        {
            self.pos += 4;
            return true;
        }
        false
    }

    /// `ReadPointerType`: `NULL`, `@` or `:` (taken), or the value of the declared type (nothing
    /// taken), after `SkipWhiteSpace` (`PeekChar` with `true`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1606-1630
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:161-164
    pub(crate) fn pointer_type(&mut self) -> Result<PointerType> {
        let c = self.skip_white_space()?;
        if self.take_null_word(c) {
            return Ok(PointerType::Null);
        }
        match c {
            b'@' => {
                self.pos += 1;
                Ok(PointerType::Object)
            }
            b':' => {
                self.pos += 1;
                Ok(PointerType::Other)
            }
            _ => Ok(PointerType::This),
        }
    }

    /// The name or number that `ReadEnum` reads after `SkipWhiteSpace`: an identifier when the
    /// byte there is a lower-case letter, else a number (`ReadInt4`). The one reader of both for
    /// the archive reader and the object streams (`asn_input`), which look it up in their values.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    #[inline(always)]
    pub(crate) fn enum_word(&mut self) -> Result<EnumWord> {
        let (pos, c) = self.white_space_at(self.pos)?;
        if !c.is_ascii_lowercase() {
            self.pos = pos;
            return Ok(EnumWord::Number(self.int4()?));
        }
        let end = self.end_of_id_at(pos)?;
        self.pos = end;
        Ok(EnumWord::Name(pos..end))
    }

    /// `ReadEnum` of an ENUMERATED type with the values `values`: a name or a number.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    pub(crate) fn enumerated(
        &mut self,
        values: &[(&'static str, i32)],
    ) -> Result<(&'static str, i32)> {
        match self.enum_item(values)? {
            Ok(item) => Ok(item),
            Err(value) => Err(self.error_at(
                NcbiThrow::Plain(FIND_NAME),
                format!("invalid value of enumerated type: {value}"),
            )),
        }
    }

    /// The item of `values` named or numbered by the text, or the number that is none of them
    /// (a name that is none of them is an error).
    fn enum_item(
        &mut self,
        values: &[(&'static str, i32)],
    ) -> Result<std::result::Result<(&'static str, i32), i32>> {
        let id = match self.enum_word()? {
            EnumWord::Name(id) => id,
            EnumWord::Number(value) => {
                return Ok(values
                    .iter()
                    .find(|&&(_, number)| number == value)
                    .copied()
                    .ok_or(value));
            }
        };
        let text = self.bytes(id);
        match values.iter().find(|&&(name, _)| name.as_bytes() == text) {
            Some(&item) => Ok(Ok(item)),
            None => {
                let message = format!(
                    "invalid value of enumerated type: {}",
                    String::from_utf8_lossy(text)
                );
                Err(self.error_at(NcbiThrow::Plain(FIND_VALUE), message))
            }
        }
    }

    /// A VisibleString of a hand grammar (and a StringStore, which the text stream reads as one)
    /// held by the frame `at`: `visible_string`, each byte outside 0x20..0x7E read as `#` and
    /// posted with the stack trace of `at` (`Trail::replaced`), as NCBI's `ReadString` of the type
    /// posts it.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:988-991
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (cited at
    // `report::serial::encoder::ReplacedChar`)
    pub(crate) fn string(&mut self, at: &Frame<'_>) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.visible_string(
            &mut out,
            true,
            false,
            &mut move |trail, line, byte, piece| trail.post_replaced(at, line, byte, piece),
        )?;
        Ok(out)
    }

    /// A UTF8String: its bytes as they are, as NCBI reads that type (`eFNP_Allow`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:988-991
    pub(crate) fn utf8_string(&mut self) -> Result<Vec<u8>> {
        let mut out = Vec::new();
        self.visible_string(&mut out, false, false, &mut |_, _, _, _| {})?;
        Ok(out)
    }

    /// The `CEofException` inside a string, which `ReadStringValue` turns into its error
    /// (`UnendedString`, cited at [`UNENDED_STRING`]) at the line where the text ends.
    fn unended_string(&self, start_line: usize) -> AsnTextError {
        self.error_at(
            NcbiThrow::Stream(UNENDED_STRING),
            format!("unclosed string starts at line {start_line}"),
        )
    }

    /// `ReadStringValue` (`SkipString` when `skip`) of a VisibleString (`visible`) or a
    /// UTF8String: the string's lines joined and `""` read as one quote, appended to `out`
    /// (nothing when skipping). In a VisibleString each byte outside 0x20..0x7E is read as `#`,
    /// and `bad_char` gets the reader's trail (where a hand grammar posts it, `string`), the line,
    /// the byte and the piece of the string that NCBI's message
    /// quotes (`ReplaceVisibleChar`'s post, `report::serial::encoder::ReplacedChar`): the piece
    /// being appended (`AppendStringData`: up to the end of the line, the quote, or 128 bytes;
    /// after `""` the next piece starts with the second quote), none when skipping. The one
    /// reader of a VisibleString in ASN.1 text: the typed reader's (`asn_types::read`, the PSSM
    /// reader and rpsbproc's), the mask reader's (`writedb::mask_info`) and a hand grammar's
    /// (`string`, rpsbproc's archive reader).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:877-912 (`AppendStringData`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-986 (`ReadStringValue`, quoted
    // at `string`; its `catch` of a `CEofException` calls `UnendedString`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1093-1155
    // (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0, objistr.hpp:1163: `valid` stays 0;
    // the fix is `eFNP_ReplaceAndWarn`, the default of `SERIAL_WRONG_CHARS_READ`, objistr.cpp:302.)
    pub(crate) fn visible_string(
        &mut self,
        out: &mut Vec<u8>,
        visible: bool,
        skip: bool,
        bad_char: &mut impl FnMut(&mut Trail, usize, u8, &[u8]),
    ) -> Result<()> {
        // `pos` is the start of the run being read (the bytes up to the next line end or quote,
        // which NCBI appends in pieces of 128 bytes from its start), `i` its bytes so far (after
        // `""` the run starts with the second quote).
        let mut pos = self.expect_at(self.pos, b'"')?;
        let mut line = self.line;
        let start_line = line;
        let mut i = 0usize;
        loop {
            // Whether the run holds a byte that a VisibleString reads as `#`.
            let mut clean = true;
            let (end, c) = self.src.scan(pos + i, |c| {
                if (b' '..=b'~').contains(&c) {
                    return c != b'"';
                }
                if matches!(c, b'\r' | b'\n') {
                    return false;
                }
                clean &= !visible;
                true
            });
            let Some(c) = c else {
                // The pieces of 128 bytes before the end of the text were appended; the rest was
                // not (`PeekChar` throws first).
                let whole = (end - pos) / 128 * 128;
                self.append_run(pos..end, pos + whole, line, out, clean, skip, bad_char);
                self.pos = pos;
                self.line = line;
                return Err(self.unended_string(start_line));
            };
            // A run without such a byte is appended whole (nothing when skipping).
            if !clean {
                self.append_run(pos..end, end, line, out, clean, skip, bad_char);
            } else if !skip {
                out.extend_from_slice(self.src.bytes(pos..end));
            }
            pos = end + 1;
            i = 0;
            if c == b'"' {
                if self.src.byte_at(pos) != Some(b'"') {
                    self.pos = pos;
                    self.line = line;
                    return Ok(());
                }
                // `""`: the second quote starts the next run (`SkipString` takes both).
                if skip {
                    pos += 1;
                } else {
                    i = 1;
                }
            } else {
                line += 1;
                pos = self.end_of_line_at(pos, c);
            }
        }
    }

    /// `AppendStringData` of the pieces of the run `run` (on line `line`) up to `appended`: 128
    /// bytes each from the run's start, the last up to the run's end. `clean`: no byte of the run
    /// is one that the string reads as `#` (then the run is appended whole). A skipped string
    /// appends nothing and posts each such byte of the whole run without a piece (`SkipString`).
    #[allow(clippy::too_many_arguments)]
    #[cold]
    #[inline(never)]
    fn append_run(
        &mut self,
        run: Range<usize>,
        appended: usize,
        line: usize,
        out: &mut Vec<u8>,
        clean: bool,
        skip: bool,
        bad_char: &mut dyn FnMut(&mut Trail, usize, u8, &[u8]),
    ) {
        let (src, trail) = (&self.src, &mut self.trail);
        if skip {
            if !clean {
                for &c in src.bytes(run) {
                    if !(b' '..=b'~').contains(&c) {
                        bad_char(trail, line, c, &[]);
                    }
                }
            }
            return;
        }
        let data = src.bytes(run.start..appended);
        if clean {
            out.extend_from_slice(data);
            return;
        }
        for piece in data.chunks(128) {
            let mut done = 0;
            for (i, &c) in piece.iter().enumerate() {
                if !(b' '..=b'~').contains(&c) {
                    out.extend_from_slice(&piece[done..i]);
                    bad_char(trail, line, c, piece);
                    out.push(b'#');
                    done = i + 1;
                }
            }
            out.extend_from_slice(&piece[done..]);
        }
    }

    /// An OCTET STRING: `'`, hex digits across lines, `'H`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1477-1538
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
                    return Err(self.error_at(
                        NcbiThrow::Stream(GET_HEX_CHAR),
                        format!("bad char in octet string: #{}", c as i8),
                    ));
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
            return Err(self.expected(b'H'));
        }
        self.pos += 1;
        Ok(out)
    }

    /// `ReadNumber`: a sign or a digit, then digits; returns their positions.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:315-329
    pub(crate) fn number(&mut self) -> Result<Range<usize>> {
        let (start, c) = self.white_space_at(self.pos)?;
        if c != b'-' && c != b'+' && !c.is_ascii_digit() {
            self.pos = start;
            return Err(self.error_at(NcbiThrow::Stream(READ_NUMBER), "invalid number"));
        }
        let (end, _) = self.src.scan(start + 1, |c| c.is_ascii_digit());
        self.pos = end;
        Ok(start..end)
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

    /// `CObjectIStreamAsn::ReadDouble`: a word, or `{ mantissa, base, exponent }` whose base is
    /// read with `ReadUint4` (a `CUtilException` for a `-`) and exponent with `ReadInt4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857
    pub(crate) fn real(&mut self) -> Result<f64> {
        if let Some(word) = self.real_word()? {
            return Ok(word);
        }
        let (mantissa, is_negative) = self.real_mantissa()?;
        self.expect(b',')?;
        let base = self.uint4()?;
        self.expect(b',')?;
        let exp = self.int4()?;
        self.expect(b'}')?;
        real_of_parts(mantissa, is_negative, base, exp).ok_or_else(|| {
            self.error_at(
                NcbiThrow::Stream(READ_DOUBLE),
                "illegal REAL base (must be 2 or 10)",
            )
        })
    }

    /// The word form of `ReadDouble` (`PLUS-INFINITY`, `MINUS-INFINITY`, `NOT-A-NUMBER` or a
    /// decimal number), or `None` at a `{` (not taken).
    fn real_word(&mut self) -> Result<Option<f64>> {
        if self.skip_white_space()? == b'{' {
            return Ok(None);
        }
        let start = self.pos;
        let word = self.scan_end_of_id(true)?;
        if self.starts_with_ignore_case(start, b"PLUS-INFINITY") {
            return Ok(Some(f64::INFINITY));
        }
        if self.starts_with_ignore_case(start, b"MINUS-INFINITY") {
            return Ok(Some(f64::NEG_INFINITY));
        }
        if self.starts_with_ignore_case(start, b"NOT-A-NUMBER") {
            return Ok(Some(NCBI_NAN));
        }
        Ok(Some(string_to_double_posix(self.bytes(word), true).value))
    }

    /// The `{` (at the reader's position) and the mantissa of `ReadDouble`: its value and whether
    /// it starts with `-`.
    fn real_mantissa(&mut self) -> Result<(f64, bool)> {
        self.pos += 1;
        let is_negative = self.skip_white_space()? == b'-';
        let mantissa_range = self.number()?;
        if mantissa_range.len() >= 128 - 1 {
            return Err(self.error_at(NcbiThrow::Stream(READ_DOUBLE_OVERFLOW), "buffer overflow"));
        }
        let mantissa_length = mantissa_range.len();
        let read = string_to_double_posix(self.bytes(mantissa_range), false);
        if read.end != mantissa_length {
            return Err(self.error_at(
                NcbiThrow::Stream(READ_DOUBLE_BAD),
                format!("bad double in line {}", self.line),
            ));
        }
        Ok((read.value, is_negative))
    }
}

/// The end of `ReadDouble` (cited at `Reader::real`): the value of the mantissa `mantissa` (read
/// with its sign; `is_negative` whether it had a `-`), `base` and `exp`; `None` for a base other
/// than 2 and 10 ("illegal REAL base").
fn real_of_parts(mut mantissa: f64, is_negative: bool, base: u32, exp: i32) -> Option<f64> {
    if base != 2 && base != 10 {
        return None;
    }
    if mantissa == 0.0 {
        return Some(mantissa);
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
    Some(if is_negative { -result } else { result })
}

/// glibc's `pow(10.0, exp)` (the C library NCBI 2.17.0 for x86-64 Linux runs with, version 2.28
/// on): the power of ten rounded to the nearest double, except at the exponents 23 and 210,
/// where glibc's result is the double above the nearest one (its error bound is 0.52 units in
/// the last place, and 10^23 and 10^210 lie almost halfway between two doubles). The table keeps
/// the bytes the same on every platform; `libm_pow10_is_glibcs_pow` checks it against the C
/// library on Linux. The powers from 10^-400 to 10^400 are kept once computed, each the first time
/// it is asked for (a stream of archives reads hundreds of thousands of REALs, a PSSM a few hundred,
/// of a handful of exponents: computing all 801 at the first REAL cost 735,000 instructions); the
/// others as they are asked for.
pub(crate) fn libm_pow10(exp: i32) -> f64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    // A power not yet computed: a NaN's bits, which no power of ten has (10^-400 is +0.0).
    const NOT_YET: u64 = u64::MAX;
    static TABLE: [AtomicU64; 801] = [const { AtomicU64::new(NOT_YET) }; 801];
    if !(-400..=400).contains(&exp) {
        return libm_pow10_of(exp);
    }
    let slot = &TABLE[(exp + 400) as usize];
    let bits = slot.load(Ordering::Relaxed);
    if bits != NOT_YET {
        return f64::from_bits(bits);
    }
    let value = libm_pow10_of(exp);
    slot.store(value.to_bits(), Ordering::Relaxed);
    value
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

#[cfg(test)]
#[path = "asn_text_reader_tests.rs"]
pub(crate) mod tests;
