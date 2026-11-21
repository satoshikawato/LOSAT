//! The primitives of NCBI's two ASN.1 object input streams as the object stacks of LOSAT's readers
//! call them (`Input`): the ASN.1 text stream (`CObjectIStreamAsn`, [`TextInput`], over the
//! lexical layer of `asn_text_reader`) and the binary one (`CObjectIStreamAsnBinary`,
//! `asn_ber_reader::BerReader`). The object stacks are NCBI's `CObjectIStream` (the frames,
//! `ReadClassSequential`, `ReadChoiceSimple`, `ReadContainer`): the reader of any type of the type
//! table (`asn_types::read`) and the reader of makeblastdb's mask files (`writedb::mask_info`).
//! Each is generic over `Input`, so a file of either encoding is read by the same frames, with
//! the primitives monomorphised (no dynamic call on a hot path).
//!
//! What a primitive throws is NCBI's exception class (`Fail`): the end of the input (NCBI's
//! `CEofException`, which the frames turn into a `CSerialException` unless fewer than two frames
//! stay, `HandleEOF`), a `CSerialException` at its place in NCBI's source with the stream's
//! position (`throw_error`), another `CException` (the numbers of the text stream), a
//! `std::bad_alloc` or another `std::exception` of a string's reserve, or LOSAT's refusal of what
//! NCBI's process may or may not read by the host's memory. What a primitive posts while it reads
//! (a byte of a VisibleString read as `#`, a UTF8String tag for a VisibleString) is handed to the
//! caller in order (`StringPost`), which adds the stack trace that NCBI's message names.

use super::asn_text_reader::{
    AsnSource, AsnTextError, AsnTextErrorKind, MemberId, NcbiThrow, PointerType, Reader,
    BAD_NUMBER, FILL_BUFFER_EOF, NUMBER_OVERFLOW, SKIP_WHITE_SPACE_EOF, UNEXPECTED_MEMBER,
};
use super::encoder::StreamPosition;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::string_to_int;
use crate::rpsbproc::{HOST_DEPENDENT_BYTES, USER_SPACE_BYTES};

/// What a primitive or a frame reports: the class of what NCBI's reader throws. What a variant
/// holds is in a box of one word (`Fail::serial`, `Fail::exception`, `Fail::unsupported`), so
/// `Fail` is a tag and a word, and the result of a primitive (`R<u32>`, `R<()>`, `R<bool>`) is
/// returned in two registers (an unboxed `String` made every result a niche of its capacity, a
/// 64-bit constant at each `Ok` and `?`, returned through memory); no error is built on the
/// paths that read a valid file.
#[derive(Debug)]
pub(crate) enum Fail {
    /// NCBI's `CEofException` (see the module comment).
    Eof,
    /// A `CSerialException`, with the frame stack that `END_OBJECT_FRAME` adds to it.
    Serial(Box<NcbiException>),
    /// Another `CException` (the `CUtilException` of a number in text, the `CStringException` of
    /// a member's tag), with a line for each frame it leaves (`NCBI_RETHROW_SAME`).
    Exception(Box<NcbiException>),
    /// `std::bad_alloc` (the string that a length announces, `GetChars`' `reserve`).
    OutOfMemory,
    /// Another `std::exception`, by its `what()`: `std::length_error` of a reserve
    /// ([`STRING_CREATE_LENGTH_ERROR`], `vector::reserve` of a reserve hook).
    Std(&'static str),
    /// NCBI reads it, or not, by what LOSAT does not know (the host's memory); LOSAT refuses it.
    Unsupported(Box<Refusal>),
}

/// What [`Fail::Unsupported`] says of what LOSAT does not read.
#[derive(Debug)]
pub(crate) struct Refusal(pub(crate) String);

pub(crate) type R<T> = Result<T, Fail>;

impl Fail {
    /// [`Fail::Serial`] of `e`.
    #[cold]
    pub(crate) fn serial(e: NcbiException) -> Fail {
        Fail::Serial(Box::new(e))
    }

    /// [`Fail::Exception`] of `e`.
    #[cold]
    pub(crate) fn exception(e: NcbiException) -> Fail {
        Fail::Exception(Box::new(e))
    }

    /// [`Fail::Unsupported`] of `what`.
    #[cold]
    pub(crate) fn unsupported(what: String) -> Fail {
        Fail::Unsupported(Box::new(Refusal(what)))
    }
}

/// `ThrowError(flag, message)` of NCBI's object stream at `site` (`ThrowError1`,
/// objistr.cpp:828-853): a `CSerialException` whose message is the stream's position
/// (`GetPosition()`: `line N` in text, `byte N` in binary) and `message`. It returns the
/// [`Fail`], not an `R`: the path that calls it ends there (`return Err(throw_error(..))`), so
/// the reading path around it keeps nothing for its return.
#[cold]
#[inline(never)]
pub(crate) fn throw_error(
    site: ThrowSite,
    position: StreamPosition,
    message: impl AsRef<[u8]>,
) -> Fail {
    Fail::serial(NcbiException::new(
        site,
        [position.to_string().as_bytes(), b": ", message.as_ref()].concat(),
    ))
}

/// `what()` of the `std::length_error` that libstdc++'s `_M_create` throws for a capacity over
/// `max_size()` (`rpsbproc::prosite::CxxString`).
pub(crate) const STRING_CREATE_LENGTH_ERROR: &str = "basic_string::_M_create";

/// The bytes of the input that one string of `length` bytes asks NCBI for (`GetChars`:
/// `str.reserve(count)` before the bytes are read, when they are not all in the input buffer):
/// a length over `max_size()` (2^63 - 1, `rpsbproc::prosite::CxxString`) is `std::length_error`;
/// an allocation of [`USER_SPACE_BYTES`] or more `std::bad_alloc` on every host; one of
/// [`HOST_DEPENDENT_BYTES`] or more succeeds or fails by the host's memory, which LOSAT refuses
/// (DECISIONS 40 of line PS: memory that NCBI builds from a length field of the input, as
/// `losat rpsbproc` does). Below that the allocation succeeds and the read goes on to the end of
/// the input. The readers ask it when the input ends inside the string: a string whose bytes are
/// there costs what they cost, whatever NCBI reserved first.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:496-523
// (checked with makeblastdb 2.17.0 on lengths 2^62 - 1, 2^62, 2^63 - 2 and 2^63 - 1: "BLAST ran
// out of memory", status 4; 2^63: "Error: basic_string::_M_create", status 255.)
pub(crate) fn reserve_string(length: u64) -> R<()> {
    if length > (1 << 63) - 1 {
        return Err(Fail::Std(STRING_CREATE_LENGTH_ERROR));
    }
    let bytes = length + 1;
    if bytes >= USER_SPACE_BYTES {
        return Err(Fail::OutOfMemory);
    }
    if bytes >= HOST_DEPENDENT_BYTES {
        return Err(host_dependent_string(length));
    }
    Ok(())
}

/// The refusal of a string of `length` bytes that `reserve_string` leaves to the host.
#[cold]
#[inline(never)]
fn host_dependent_string(length: u64) -> Fail {
    let bytes = length + 1;
    Fail::unsupported(format!(
        "a string whose length ({length} bytes) passes the end of the input: NCBI's reader \
         allocates {bytes} bytes for it, which succeeds or fails by the host's memory"
    ))
}

/// What reading a string posts, in the order NCBI posts it.
#[derive(Debug)]
pub(crate) enum StringPost {
    /// The warning of a UTF8String tag for a VisibleString (binary, `ExpectStringTag`,
    /// objistrasnb.cpp:296-305; [`utf8_tag_message`]).
    Utf8Tag,
    /// `ReplaceVisibleChar`'s post of a byte read as `#`: where, the byte, and the string (or
    /// piece of it) that the message quotes.
    BadChar {
        position: StreamPosition,
        byte: u8,
        string: Vec<u8>,
    },
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-310
/// The message of `ERR_POST_X_ONCE(10, Warning << ...)` of a UTF8String tag for the VisibleString
/// member at the stack trace `trace` (error 803.10, once per process): the one owner of its text,
/// for the mask reader's and the table reader's posts.
pub(crate) fn utf8_tag_message(trace: &str) -> String {
    format!(
        "CObjectIStreamAsnBinary: UTF8String data for VisibleString member {trace}, ASN.1 \
         specification may need an update"
    )
}

/// The kind of a string, which decides its tag in binary and what is done with its bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StringKind {
    /// VisibleString: a byte outside 0x20..0x7E is read as `#` with NCBI's error message.
    Visible,
    /// UTF8String: the bytes as they are (`eFNP_Allow`).
    Utf8,
    /// StringStore: in text as a VisibleString (`CObjectIStream::ReadStringStore` is
    /// `ReadString`), in binary its bytes as they are (`CObjectIStreamAsnBinary::ReadStringStore`,
    /// `eFNP_Allow`).
    Store,
}

/// The encoding of an input (`ESerialDataFormat`): what a delayed member's bytes are
/// (`Value::Delayed`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Encoding {
    /// `eSerial_AsnText`
    Text,
    /// `eSerial_AsnBinary`
    Binary,
}

/// The members of a SEQUENCE or the variants of a CHOICE, as an encoding looks them up: their
/// names in module order. NCBI's modules have automatic tags, so the tag of each is its position
/// (`[n]` in text, `a0 | n` in binary; `scripts/asn_spec_gen.py` checks it for the table).
pub(crate) trait Items {
    fn count(&self) -> usize;
    fn name(&self, index: usize) -> &'static str;
}

/// The named values of an ENUMERATED (`integer` false) or of an INTEGER with named values.
pub(crate) trait EnumItems {
    /// `IsInteger()`: any number is a value.
    fn integer(&self) -> bool;
    /// The value named `name`.
    fn value_of(&self, name: &[u8]) -> Option<i32>;
    /// Whether `value` is one of the named values.
    fn has_value(&self, value: i32) -> bool;
}

/// The primitives of one encoding (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`) that NCBI's
/// object stack calls (the names are NCBI's functions).
pub(crate) trait Input {
    /// The encoding.
    const ENCODING: Encoding;
    /// `GetPosition()`: `line N` in text, `byte N` (from the start of the object) in binary.
    fn position(&self) -> StreamPosition;
    /// The state of a new object stream at the next byte (each object of a file is read by a
    /// stream of its own in NCBI: its line count, or its byte count and data limits, start
    /// again).
    fn new_object(&mut self);
    /// The depth of NCBI's object stack (`GetStackDepth()`), which the text reader's
    /// `SkipWhiteSpace` asks at the end of the text.
    #[inline]
    fn set_stack_depth(&mut self, _depth: usize) {}
    /// `ReadFileHeader`: the type name that precedes an object in text; binary files have none
    /// (the first byte peeked at: an empty input is NCBI's `CEofException`).
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>>;
    fn begin_class(&mut self) -> R<()>;
    /// `BeginClassMember(classType, pos)`: the next member of `items`, searched from `pos` on,
    /// or `None` at the end of the class.
    fn begin_member(&mut self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>>;
    fn end_member(&mut self) -> R<()>;
    fn end_class(&mut self) -> R<()>;
    fn begin_choice(&mut self) -> R<()>;
    /// `BeginChoiceVariant`: the variant of `items`.
    fn begin_variant(&mut self, items: &(impl Items + ?Sized)) -> R<usize>;
    fn end_variant(&mut self) -> R<()>;
    fn end_choice(&mut self) -> R<()>;
    /// `BeginContainer` of a SEQUENCE OF or (`set`) SET OF.
    fn begin_container(&mut self, set: bool) -> R<()>;
    /// `BeginContainerElement`: whether an element follows.
    fn begin_element(&mut self) -> R<bool>;
    fn end_container(&mut self) -> R<()>;
    /// `ReadInt4`, `ReadInt8`, `ReadUint4`.
    fn read_i32(&mut self) -> R<i32>;
    fn read_i64(&mut self) -> R<i64>;
    fn read_u32(&mut self) -> R<u32>;
    /// `SkipSNumber` (`signed`) / `SkipUNumber`.
    fn skip_number(&mut self, signed: bool) -> R<()>;
    fn read_bool(&mut self) -> R<bool>;
    fn read_null(&mut self) -> R<()>;
    /// `ReadDouble`, `SkipFNumber`.
    fn read_real(&mut self) -> R<f64>;
    fn skip_real(&mut self) -> R<()>;
    /// `ReadString` (`ReadStringStore`) of a string of `kind`, appended to `out`, with what NCBI
    /// posts while it reads it in `posts` (also when it then fails).
    fn read_string(
        &mut self,
        kind: StringKind,
        out: &mut Vec<u8>,
        posts: &mut Vec<StringPost>,
    ) -> R<()>;
    /// `SkipString` (`SkipStringStore`).
    fn skip_string(&mut self, kind: StringKind, posts: &mut Vec<StringPost>) -> R<()>;
    /// An OCTET STRING (`ReadBytes` of a `vector<char>`), `SkipByteBlock`.
    fn read_octets(&mut self) -> R<Vec<u8>>;
    fn skip_octets(&mut self) -> R<()>;
    /// `ReadEnum(values)` (which `SkipEnum` calls too): the value (`CEnumeratedTypeInfo::ReadEnum`
    /// adds its line to an error).
    fn read_enum(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32>;
    /// `ReadPointerType()`, where NCBI reads a `CRef<>` that is empty at that moment: a `NULL` is
    /// taken, the `@` or `:` that starts the other two forms too.
    fn read_pointer_type(&mut self) -> R<PointerType>;
    /// `ReadObjectPointer()`: the object index after `@` (text: `ReadInt8`, `TObjectIndex` being
    /// 64 bits) or of the object reference tag (binary: `ReadStdSigned`). Its value does not
    /// matter: `GetRegisteredObject` throws for every index.
    fn read_object_pointer(&mut self) -> R<()>;
    /// `ReadOtherPointer()`: the name of the type whose value follows.
    fn read_other_pointer(&mut self) -> R<Vec<u8>>;
    /// `ReadOtherPointerEnd()`, after that value.
    fn read_other_pointer_end(&mut self) -> R<()>;
    /// `StartDelayBuffer`: the position in the source where a delayed member's bytes start (text:
    /// after white space). The source keeps the bytes from there until [`Input::delayed`].
    fn delay_start(&mut self) -> R<usize>;
    /// `EndDelayBuffer`: the bytes of the delayed member from `start` (`delay_start`) to the
    /// current position.
    fn delayed(&mut self, start: usize) -> Vec<u8>;
}

/// ASN.1 text (`MSerial_AsnText`, `CObjectIStreamAsn`): the members and variants found by name or
/// by tag, `m_BlockStart`, the numbers, and what NCBI's text stream throws for an error of the
/// lexical reader (`TextInput::fail`: the end of the text is a `CEofException`, except inside
/// `SkipWhiteSpace` at a stack depth over 2).
///
/// A text read from a stream in blocks (`StreamSource`) is read by one `TextInput` for all objects
/// of the stream, as NCBI's object stream gives back the bytes it did not use when an object ends
/// (`CIStreamBuffer::Close`, strbuffer.cpp:125-139), so that the next object starts at the byte
/// after the last one read. Each object gets a new `CObjectIStream` in NCBI (`s_ReadObject`), so
/// its line count starts at 1 again (checked with makeblastdb 2.17.0: a string on the first line
/// of a second object is "at line 1").
pub(crate) struct TextInput<S: AsnSource> {
    pub(crate) r: Reader<S>,
    /// `m_BlockStart`
    block_start: bool,
    /// `GetStackDepth()` (`Input::set_stack_depth`).
    depth: usize,
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
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1455-1463 (`BeginChoiceVariant`,
// cited at `begin_variant`)
const BEGIN_CHOICE_VARIANT: ThrowSite = asn_site(
    1463,
    FORMAT_ERROR,
    "ncbi::CObjectIStreamAsn::BeginChoiceVariant",
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1488-1491 (`ReadChoiceSimple`, cited at
// `begin_variant`)
const READ_CHOICE_SIMPLE: ThrowSite = ThrowSite {
    file: "src/serial/objistr.cpp",
    line: 1490,
    kind: FORMAT_ERROR,
    function: Some("ncbi::CObjectIStream::ReadChoiceSimple"),
};
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1074
const SKIP_S_NUMBER: ThrowSite =
    asn_site(1042, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::SkipSNumber");
const SKIP_U_NUMBER: ThrowSite =
    asn_site(1067, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::SkipUNumber");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
const SKIP_F_NUMBER: ThrowSite =
    asn_site(1090, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::SkipFNumber");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (cited at `Reader::real`)
const READ_DOUBLE: ThrowSite = asn_site(836, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::ReadDouble");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
const SKIP_BYTE_BLOCK: ThrowSite =
    asn_site(1192, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::SkipByteBlock");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:187-205
const EXPECT_CHOICE: ThrowSite = asn_site(203, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::Expect");

impl<S: AsnSource> TextInput<S> {
    /// The text input over `src`; `document` and `reader` name the text and the reader in
    /// LOSAT's own messages (`Reader::new`).
    pub(crate) fn new(src: S, document: &'static str, reader: &'static str) -> Self {
        TextInput {
            r: Reader::new(src, document, reader),
            block_start: false,
            depth: 0,
        }
    }

    /// What NCBI's text stream throws for `error` of the lexical reader (`AsnTextError::ncbi`):
    /// the end of the text is a `CEofException`, except inside `SkipWhiteSpace` at a stack depth
    /// over 2, where it becomes the `CSerialException` that quotes it
    /// (`asn_text_reader::SKIP_WHITE_SPACE_EOF`).
    #[cold]
    fn fail(&self, error: AsnTextError) -> Fail {
        let line = error.line;
        match (error.kind, error.ncbi) {
            (AsnTextErrorKind::Eof, Some(NcbiThrow::WhiteSpaceEof)) if self.depth > 2 => {
                let eof = NcbiException::new(FILL_BUFFER_EOF, "end of file").what();
                Fail::serial(NcbiException::new(
                    SKIP_WHITE_SPACE_EOF,
                    [format!("line {line}: ").as_bytes(), &eof].concat(),
                ))
            }
            (AsnTextErrorKind::Eof, _) => Fail::Eof,
            (AsnTextErrorKind::Format, Some(NcbiThrow::Stream(site))) => Fail::serial(
                NcbiException::new(site, format!("line {line}: {}", error.message)),
            ),
            (AsnTextErrorKind::Format, Some(NcbiThrow::Plain(site))) => {
                Fail::serial(NcbiException::new(site, error.message))
            }
            (AsnTextErrorKind::Format, Some(NcbiThrow::Buffer(site))) => Fail::exception(
                NcbiException::new(site, format!("{} in line {line}", error.message)),
            ),
            // Every error of the lexical reader that this reader meets names NCBI's throw; one
            // that did not would be refused, not reported in LOSAT's words.
            _ => Fail::unsupported(error.message),
        }
    }

    /// `ThrowError(fFormatError, message)` of the text stream at `site`.
    #[cold]
    #[inline(never)]
    fn throw(&self, site: ThrowSite, message: impl AsRef<[u8]>) -> Fail {
        throw_error(site, self.position(), message)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:389-424
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:432-439,462-473
    // (The members and variants have the automatic tags 0, 1, ... in module order. Their names
    // start with a lower case letter, so `GetAltItemIndex`, which looks up the name with an upper
    // case first letter, finds nothing: `Biomol` is no member, c21 of SURVEY_D3.)
    /// `ReadMemberId(SkipWhiteSpace())` and the member or variant it names among `items`,
    /// searched from `pos` on; one that is not there is an error (`UnexpectedMember`, which lists
    /// them all). `choice`: an empty name is `BeginChoiceVariant`'s error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1305-1314
    fn find_member(
        &mut self,
        items: &(impl Items + ?Sized),
        pos: usize,
        choice: bool,
    ) -> R<usize> {
        let count = items.count();
        let member = self.r.member_or_tag().map_err(|e| self.fail(e))?;
        let (id, found) = match member {
            MemberId::Tag(digits) => {
                // `NStr::StringToInt(id)` of the tag's digits: a number that does not fit is its
                // `CStringException` (NCBI 2.17.0: `[3x]`, `[99999999999]`, probe e03, e07).
                let id = self.r.bytes(digits);
                let tag = string_to_int(id, 0).map_err(Fail::exception)?;
                let found = usize::try_from(tag)
                    .ok()
                    .filter(|tag| (pos..count).contains(tag));
                (id, found)
            }
            MemberId::Name(id) => {
                let id = self.r.bytes(id);
                if choice && id.is_empty() {
                    return Err(self.throw(BEGIN_CHOICE_VARIANT, "choice variant id expected"));
                }
                let found = (pos..count).find(|&i| items.name(i).as_bytes() == id);
                (id, found)
            }
        };
        match found {
            Some(index) => Ok(index),
            None => Err(self.unexpected_member(id, items)),
        }
    }

    /// `UnexpectedMember(id, items)` of the members or variants `items`.
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, id: &[u8], items: &(impl Items + ?Sized)) -> Fail {
        let mut message = [b"\"", id, b"\": unexpected member, should be one of: "].concat();
        for i in 0..items.count() {
            message.extend_from_slice(format!("\"{}\" ", items.name(i)).as_bytes());
        }
        self.throw(UNEXPECTED_MEMBER, message)
    }

    /// `CIStreamBuffer::BadNumber` / `NumberOverflow`: a `CUtilException`, ` in line N`.
    #[cold]
    #[inline(never)]
    fn util_error(&self, site: ThrowSite, what: &str) -> Fail {
        Fail::exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.r.line),
        ))
    }

    /// `GetInt4`, `GetUint4` (`signed` false: no `-`) and `GetInt8` after `SkipWhiteSpace`, the
    /// digits checked against `max` as NCBI checks them.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-668
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let (pos, c) = match self.r.white_space_at(self.r.pos) {
            Ok(found) => found,
            Err(e) => return Err(self.fail(e)),
        };
        // `SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
        let (pos, mut c) = match self.r.ws_at(pos, c) {
            Ok(found) => found,
            Err(e) => return Err(self.fail(e)),
        };
        let mut pos = pos + 1;
        let mut negative = false;
        if (signed && c == b'-') || c == b'+' {
            negative = c == b'-';
            c = match self.r.char_at(pos) {
                Ok(c) => c,
                Err(e) => return Err(self.fail(e)),
            };
            pos += 1;
        }
        let r = &mut self.r;
        if !c.is_ascii_digit() {
            r.pos = pos;
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        let mut overflow = false;
        let (end, _) = r.src.scan(pos, |d| {
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
            // `GetInt4` takes the digit before it checks the overflow.
            r.pos = end + 1;
            return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
        }
        r.pos = end;
        // `-Int4(n)` of n = 2^31 is the least Int4 (and likewise in 64 bits).
        Ok(if negative {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    /// The end of the text inside a number that starts at `start`: `PeekChar`'s, at `start`.
    #[cold]
    fn end_of_number(&mut self, start: usize) -> Fail {
        self.r.pos = start;
        let error = self.r.end_of_text();
        self.fail(error)
    }

    /// `SkipWhiteSpace` with NCBI's error for the end of the text (`fail`).
    fn skip_white_space(&mut self) -> R<u8> {
        let c = self.r.skip_white_space();
        c.map_err(|e| self.fail(e))
    }

    /// `Expect(c, true)`.
    fn expect(&mut self, c: u8) -> R<()> {
        let result = self.r.expect(c);
        result.map_err(|e| self.fail(e))
    }
}

#[cfg(test)]
impl<R: std::io::Read> TextInput<super::asn_text_reader::StreamSource<R>> {
    /// The bytes the stream keeps.
    pub(crate) fn buffered(&self) -> usize {
        self.r.src.buffered()
    }
}

impl<S: AsnSource> Input for TextInput<S> {
    const ENCODING: Encoding = Encoding::Text;

    /// `GetPosition()` of the text stream (objistr.cpp:807-810): `line N`.
    fn position(&self) -> StreamPosition {
        StreamPosition::Line(self.r.line as u64)
    }

    /// The state of a new `CObjectIStreamAsn` at the next byte: line 1, no block started.
    fn new_object(&mut self) {
        self.r.line = 1;
        self.block_start = false;
    }

    #[inline]
    fn set_stack_depth(&mut self, depth: usize) {
        self.depth = depth;
    }

    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        let name = self.r.file_header().map_err(|e| self.fail(e))?;
        Ok(Some(self.r.bytes(name).to_vec()))
    }

    fn begin_class(&mut self) -> R<()> {
        self.block_start = self.r.start_block().map_err(|e| self.fail(e))?;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1329-1362
    // (`CanSkipUnknownMembers()` is false by default, `m_TypeAlias` is null.)
    fn begin_member(&mut self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>> {
        let more = self.r.next_element(&mut self.block_start);
        if !more.map_err(|e| self.fail(e))? {
            return Ok(None);
        }
        Ok(Some(self.find_member(items, pos, false)?))
    }

    fn end_member(&mut self) -> R<()> {
        Ok(())
    }

    fn end_class(&mut self) -> R<()> {
        let result = self.r.end_block();
        result.map_err(|e| self.fail(e))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1478
    // (BeginChoice has just set `m_BlockStart`, so `skipname` is false. A `}` there is no
    // variant: `ReadChoiceSimple`'s "choice variant id expected".)
    fn begin_choice(&mut self) -> R<()> {
        self.block_start = true;
        Ok(())
    }

    fn begin_variant(&mut self, items: &(impl Items + ?Sized)) -> R<usize> {
        let more = self.r.next_element(&mut self.block_start);
        if !more.map_err(|e| self.fail(e))? {
            return Err(self.throw(READ_CHOICE_SIMPLE, "choice variant id expected"));
        }
        self.find_member(items, 0, true)
    }

    fn end_variant(&mut self) -> R<()> {
        Ok(())
    }

    fn end_choice(&mut self) -> R<()> {
        self.block_start = false;
        Ok(())
    }

    fn begin_container(&mut self, _set: bool) -> R<()> {
        self.begin_class()
    }

    fn begin_element(&mut self) -> R<bool> {
        let more = self.r.next_element(&mut self.block_start);
        more.map_err(|e| self.fail(e))
    }

    fn end_container(&mut self) -> R<()> {
        self.end_class()
    }

    #[inline]
    fn read_i32(&mut self) -> R<i32> {
        Ok(self.get_number(i32::MAX as u64, true)? as i32)
    }

    #[inline]
    fn read_i64(&mut self) -> R<i64> {
        self.get_number(i64::MAX as u64, true)
    }

    #[inline]
    fn read_u32(&mut self) -> R<u32> {
        Ok(self.get_number(u64::from(u32::MAX), false)? as u32)
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1048
    // (`SkipUNumber`, 1050-1074, takes only `+` and peeks past the digits with `PeekCharNoEOF`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let found = self.r.white_space_at(self.r.pos);
        let (start, mut c) = found.map_err(|e| self.fail(e))?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = match self.r.src.byte_at(start + 1) {
                Some(c) => c,
                None => return Err(self.end_of_number(start)),
            };
            i = 2;
        }
        if !c.is_ascii_digit() {
            self.r.pos = start;
            let (site, kind) = if signed {
                (SKIP_S_NUMBER, "signed")
            } else {
                (SKIP_U_NUMBER, "unsigned")
            };
            let message = format!("bad {kind} integer in line {}", self.r.line);
            return Err(self.throw(site, message));
        }
        let (end, c) = self.r.src.scan(start + i, |c| c.is_ascii_digit());
        // `SkipSNumber` peeks past the digits with `PeekChar`, which fails at the end of the
        // text; `SkipUNumber` with `PeekCharNoEOF`.
        if signed && c.is_none() {
            return Err(self.end_of_number(start));
        }
        self.r.pos = end;
        Ok(())
    }

    fn read_bool(&mut self) -> R<bool> {
        let value = self.r.boolean();
        value.map_err(|e| self.fail(e))
    }

    fn read_null(&mut self) -> R<()> {
        let value = self.r.null();
        value.map_err(|e| self.fail(e))
    }

    /// `ReadDouble`: a word, or `{ mantissa, base, exponent }` whose base is read with
    /// `ReadUint4` (a `CUtilException` for a `-`) and exponent with `ReadInt4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (cited at `Reader::real`)
    fn read_real(&mut self) -> R<f64> {
        let word = self.r.real_word();
        if let Some(word) = word.map_err(|e| self.fail(e))? {
            return Ok(word);
        }
        let mantissa = self.r.real_mantissa();
        let (mantissa, is_negative) = mantissa.map_err(|e| self.fail(e))?;
        self.expect(b',')?;
        let base = self.read_u32()?;
        self.expect(b',')?;
        let exp = self.read_i32()?;
        self.expect(b'}')?;
        super::asn_text_reader::real_of_parts(mantissa, is_negative, base, exp)
            .ok_or_else(|| self.throw(READ_DOUBLE, "illegal REAL base (must be 2 or 10)"))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
    fn skip_real(&mut self) -> R<()> {
        if self.skip_white_space()? != b'{' {
            let id = self.r.scan_end_of_id(true);
            id.map_err(|e| self.fail(e))?;
            return Ok(());
        }
        self.r.pos += 1;
        self.skip_number(true)?;
        self.expect(b',')?;
        let base = self.read_u32()?;
        self.expect(b',')?;
        self.skip_number(true)?;
        self.expect(b'}')?;
        if base != 2 && base != 10 {
            return Err(self.throw(SKIP_F_NUMBER, "illegal REAL base (must be 2 or 10)"));
        }
        Ok(())
    }

    /// `ReadStringValue` with the fix of `kind` (`Reader::visible_string`): for a VisibleString
    /// and a StringStore each byte outside 0x20..0x7E read as `#`, posted at its line with the
    /// piece being appended.
    fn read_string(
        &mut self,
        kind: StringKind,
        out: &mut Vec<u8>,
        posts: &mut Vec<StringPost>,
    ) -> R<()> {
        self.string(kind, Some(out), posts)
    }

    fn skip_string(&mut self, kind: StringKind, posts: &mut Vec<StringPost>) -> R<()> {
        self.string(kind, None, posts)
    }

    /// `ReadBytes` of an OCTET STRING (`Reader::octets`).
    fn read_octets(&mut self) -> R<Vec<u8>> {
        let value = self.r.octets();
        value.map_err(|e| self.fail(e))
    }

    /// `SkipByteBlock` (the skip takes `H` or `B` after white space).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
    fn skip_octets(&mut self) -> R<()> {
        self.expect(b'\'')?;
        loop {
            let c = self.r.peek_char(0);
            let c = c.map_err(|e| self.fail(e))?;
            self.r.pos += 1;
            match c {
                b'0'..=b'9' | b'A'..=b'F' | b'a'..=b'f' => {}
                b'\'' => break,
                b'\r' | b'\n' => self.r.skip_end_of_line(c),
                _ => {
                    self.r.pos -= 1;
                    let message = format!("bad char in octet string: #{}", c as i8);
                    return Err(self.throw(SKIP_BYTE_BLOCK, message));
                }
            }
        }
        match self.skip_white_space()? {
            b'H' | b'B' => {
                self.r.pos += 1;
                Ok(())
            }
            _ => Err(self.throw(EXPECT_CHOICE, "'H' or 'B' expected")),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    /// `ReadEnum(values)`: a lower-case name of the list, or a number (one of the list for an
    /// ENUMERATED).
    fn read_enum(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32> {
        let name = self.r.enum_name();
        if let Some(id) = name.map_err(|e| self.fail(e))? {
            return match values.value_of(self.r.bytes(id.clone())) {
                Some(value) => Ok(value),
                None => Err(self.fail(self.r.unknown_enum_name(id))),
            };
        }
        let value = self.read_i32()?;
        if !values.integer() && !values.has_value(value) {
            return Err(self.fail(self.r.unknown_enum_value(value)));
        }
        Ok(value)
    }

    /// `ReadPointerType` (`Reader::pointer_type`).
    fn read_pointer_type(&mut self) -> R<PointerType> {
        let value = self.r.pointer_type();
        value.map_err(|e| self.fail(e))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1631-1640
    // (`TObjectIndex` is `size_t`, objistr.hpp:1029.)
    fn read_object_pointer(&mut self) -> R<()> {
        self.read_i64().map(|_| ())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1643-1646
    /// `ReadOtherPointer` (`Reader::type_id`). The name of `[]` is the `std::string` of
    /// `(size_t)-1` bytes that `CTempString::operator string` asks for: `std::length_error`.
    fn read_other_pointer(&mut self) -> R<Vec<u8>> {
        let id = self.r.type_id();
        match id.map_err(|e| self.fail(e))? {
            Some(id) => Ok(self.r.bytes(id).to_vec()),
            None => Err(Fail::Std(STRING_CREATE_LENGTH_ERROR)),
        }
    }

    /// `ReadOtherPointerEnd` of the text stream: `CObjectIStream`'s, which reads nothing
    /// (objistr.cpp:1223-1225).
    fn read_other_pointer_end(&mut self) -> R<()> {
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1648-1652
    // (psiblast and tblastn set no hook on their stream: `ShouldParseDelayBuffer` is false,
    // objistr.cpp:745-760.)
    /// `StartDelayBuffer`: white space skipped, the text kept from there.
    fn delay_start(&mut self) -> R<usize> {
        self.skip_white_space()?;
        self.r.src.hold(Some(self.r.pos));
        Ok(self.r.pos)
    }

    fn delayed(&mut self, start: usize) -> Vec<u8> {
        let bytes = self.r.bytes(start..self.r.pos).to_vec();
        self.r.src.hold(None);
        bytes
    }
}

impl<S: AsnSource> TextInput<S> {
    /// A string of `kind`, appended to `out` (`None`: skipped), with NCBI's posts.
    fn string(
        &mut self,
        kind: StringKind,
        out: Option<&mut Vec<u8>>,
        posts: &mut Vec<StringPost>,
    ) -> R<()> {
        let skip = out.is_none();
        let mut none = Vec::new();
        let out = out.unwrap_or(&mut none);
        let result = self.r.visible_string(
            out,
            kind != StringKind::Utf8,
            skip,
            &mut |line, byte, piece| {
                posts.push(StringPost::BadChar {
                    position: StreamPosition::Line(line as u64),
                    byte,
                    string: piece.to_vec(),
                })
            },
        );
        result.map_err(|e| self.fail(e))
    }
}
