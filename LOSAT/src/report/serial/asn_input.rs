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
//! caller as it is made (`StringPost`), which adds the stack trace that NCBI's message names.

use super::asn_text_reader::{
    AsnSource, AsnTextError, AsnTextErrorKind, EnumWord, MemberId, NcbiThrow, PointerType, Reader,
    FILL_BUFFER_EOF, FIND_NAME, FIND_VALUE, SKIP_WHITE_SPACE_EOF, UNEXPECTED_MEMBER,
};
use super::asn_types::read::PostSink;
use super::asn_types::Member;
use super::encoder::StreamPosition;
use super::object_frame::{text_frame_site, FrameKind};
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
    /// A `CSerialException`, with the frame stack that `END_OBJECT_FRAME` adds to it
    /// ([`SerialFail`]).
    Serial(Box<SerialFail>),
    /// Another `CException` (the `CUtilException` of a number in text, the `CStringException` of
    /// a member's tag), with a line for each frame it leaves (`NCBI_RETHROW_SAME`).
    Exception(Box<NcbiException>),
    /// `std::bad_alloc` (the string that a length announces, `GetChars`' `reserve`).
    OutOfMemory,
    /// Another `std::exception`, by its `what()`: `std::length_error` of a reserve
    /// ([`STRING_CREATE_LENGTH_ERROR`], `vector::reserve` of a reserve hook). The `what()` is
    /// held by a thin reference to its constant, so `Fail` stays a tag and a word (see above).
    Std(&'static &'static str),
    /// NCBI reads it, or not, by what LOSAT does not know (the host's memory); LOSAT refuses it.
    Unsupported(Box<Refusal>),
}

/// A `CSerialException` on its way out through the object frames: the exception, and the names
/// (`GetFrameName()`) of the frames it has left since it was last used whole, innermost first,
/// which `END_OBJECT_FRAME` adds to its frame stack (`AddFrameInfo`). They are added at once where
/// the exception is used whole ([`SerialFail::exception`]): NCBI puts each name in front of the
/// stack, so an error deep inside a nest of any depth built its stack in time quadratic in the
/// depth; this builds it in linear time.
#[derive(Debug)]
pub(crate) struct SerialFail {
    exception: NcbiException,
    frames: Vec<&'static str>,
}

impl SerialFail {
    /// `AddFrameInfo` of the frame's name `name` as the exception leaves the frame (an empty name
    /// adds nothing).
    pub(crate) fn add_frame(&mut self, name: &'static str) {
        if !name.is_empty() {
            self.frames.push(name);
        }
    }

    /// `GetErrCode()` of the exception, by its name (`NcbiException::kind`).
    pub(crate) fn kind(&self) -> &'static str {
        self.exception.kind()
    }

    /// The exception with the names of the frames it left added to its frame stack.
    pub(crate) fn exception(self) -> NcbiException {
        let SerialFail {
            mut exception,
            frames,
        } = self;
        exception.add_frame_infos(&frames);
        exception
    }
}

/// What [`Fail::Unsupported`] says of what LOSAT does not read.
#[derive(Debug)]
pub(crate) struct Refusal(pub(crate) String);

pub(crate) type R<T> = Result<T, Fail>;

// `Fail` is a tag and a word: a primitive's `R<()>` comes back in two registers, its `R<u32>` or
// `R<bool>` in 16 bytes. A third word (a `&'static str` in a variant) returned every result of every
// primitive through 24 bytes of memory (part of +16.6 % instructions on a binary mask file, line RG
// item D3c-3a).
const _: () = assert!(std::mem::size_of::<Fail>() == 2 * std::mem::size_of::<usize>());
const _: () = assert!(std::mem::size_of::<R<u32>>() == std::mem::size_of::<Fail>());

impl Fail {
    /// [`Fail::Serial`] of `e`.
    #[cold]
    pub(crate) fn serial(e: NcbiException) -> Fail {
        Fail::Serial(Box::new(SerialFail {
            exception: e,
            frames: Vec::new(),
        }))
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

/// `ThrowError` of NCBI's object stream with a flag and a message, at `site` (`ThrowError1`,
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

/// The `CSerialException` (`eInvalidData`) of `CEnumeratedTypeValues::FindValue` (a name) or
/// `FindName` (a number) at `site` for `value`, which names no value of the type (the binary
/// stream's `ReadEnum` too).
#[cold]
#[inline(never)]
pub(crate) fn invalid_enum(site: ThrowSite, value: String) -> Fail {
    Fail::serial(NcbiException::new(
        site,
        format!("invalid value of enumerated type: {value}"),
    ))
}

/// `what()` of the `std::length_error` that libstdc++'s `_M_create` throws for a capacity over
/// `max_size()` (`rpsbproc::prosite::CxxString`).
pub(crate) const STRING_CREATE_LENGTH_ERROR: &str = "basic_string::_M_create";

/// `what()` of the `std::length_error` that libstdc++'s `vector::reserve` throws for a capacity
/// over `max_size()` (2^63 - 1 bytes for a `vector<char>`): a byte block's reserve, a reserve hook
/// (`asn_types::hooks`), `losat rpsbproc`'s tables.
pub(crate) const VECTOR_RESERVE_LENGTH_ERROR: &str = "vector::reserve";

/// What LOSAT's own messages call the document that a reader reads, the reader, the input that a
/// length can pass the end of, and the part of NCBI's program that allocates what a length asks
/// for (the lexical layer's refusals, `asn_text_reader::Reader::new`, and [`reserve`]'s): the one
/// owner of those words for each caller.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Names {
    pub(crate) document: &'static str,
    pub(crate) reader: &'static str,
    pub(crate) input: &'static str,
    pub(crate) allocator: &'static str,
    /// Whether the program reports every `CSerialException` of the read with a message of its own
    /// (psiblast's and tblastn's "Unrecognized format for PSSM", `CPsiBlastArgs`), so that what
    /// the exception says does not show (`read::TypedReader::other_type`).
    pub(crate) own_serial_message: bool,
}

impl Names {
    /// The PSSM file of `-in_pssm` (psiblast and tblastn).
    pub(crate) const PSSM: Names = Names {
        document: "PSSM",
        reader: "PSSM reader",
        input: "input",
        allocator: "NCBI's reader",
        own_serial_message: true,
    };
    /// A profile file of makeprofiledb (`-in`, ASN.1 text or `-binary`). makeprofiledb prints the
    /// exception's own text, so `:OtherType` depends on the types its process has met before.
    pub(crate) const PROFILE: Names = Names {
        document: "profile",
        reader: "PSSM reader",
        input: "file",
        allocator: "makeprofiledb",
        own_serial_message: false,
    };
    /// Binary definition lines (`ASN1_BlastDefLine` user object) of makeprofiledb.
    pub(crate) const BINARY_DEFLINES: Names = Names {
        document: "binary definition lines",
        reader: "binary definition lines reader",
        input: "binary header",
        allocator: "makeprofiledb",
        own_serial_message: false,
    };
    /// A mask data file of makeblastdb's `-mask_data`.
    pub(crate) const MASK: Names = Names {
        document: "mask data",
        reader: "mask data reader",
        input: "file",
        allocator: "makeblastdb",
        own_serial_message: false,
    };
    /// The input of `losat rpsbproc`: `Blast4-archive` objects in ASN.1 text (`rpsbproc::archive`).
    pub(crate) const ARCHIVE: Names = Names {
        document: "archive",
        reader: "rpsbproc",
        input: "input",
        allocator: "NCBI's rpsbproc",
        own_serial_message: false,
    };
    /// The bytes of a delayed member that a writer reads again (`CDelayBuffer::Update`). The read
    /// that kept them took them whole, so no length passes their end.
    pub(crate) const DELAYED: Names = Names {
        document: "delayed member",
        reader: "delayed member reader",
        input: "member's bytes",
        allocator: "NCBI's reader",
        own_serial_message: false,
    };
    /// The documents of the tests of the readers.
    #[cfg(test)]
    pub(crate) const TEST: Names = Names {
        document: "test document",
        reader: "test reader",
        input: "input",
        allocator: "NCBI's reader",
        own_serial_message: true,
    };
}

/// What NCBI reserves for a length before it reads the bytes ([`reserve`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Reserve {
    /// A string: `GetChars` reserves the count in the string when the bytes are not all in the
    /// input buffer (a capacity of the length and its terminating 0).
    String,
    /// A byte block of known length (a binary OCTET STRING or BIT STRING) read into a
    /// `vector<char>`: `reserve(length)` before its bytes.
    Bytes,
}

/// The reserve that NCBI makes for `length` bytes of `what` before it reads them, the one owner
/// of its outcome for strings and byte blocks: a capacity over `max_size()` (2^63 - 1) is
/// `std::length_error`; an allocation of [`USER_SPACE_BYTES`] or more `std::bad_alloc` on every
/// host; one of [`HOST_DEPENDENT_BYTES`] or more succeeds or fails by the host's memory, which
/// LOSAT refuses in the words of `names` (DECISIONS 40 of line PS: memory that NCBI builds from a
/// length field of the input, as `losat rpsbproc` does). Below that the allocation succeeds and
/// the read goes on to the end of the input. The readers ask it when the input ends inside the
/// bytes: bytes that are there cost what they cost, whatever NCBI reserved first.
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:496-523
// NCBI reference (598d8ae6): c++/src/serial/stdtypes.cpp:1538-1548
// (checked with makeblastdb 2.17.0 on string lengths 2^62 - 1, 2^62, 2^63 - 2 and 2^63 - 1: "BLAST
// ran out of memory", status 4; 2^63: "Error: basic_string::_M_create", status 255; with psiblast
// 2.17.0 on byte blocks of 2^62 and 2^63 - 1 bytes: status 4; of 2^63 and 2^64 - 1 bytes: "Error:
// vector::reserve", status 255.)
pub(crate) fn reserve(length: u64, what: Reserve, names: Names) -> R<()> {
    if length > (1 << 63) - 1 {
        return Err(Fail::Std(match what {
            Reserve::String => &STRING_CREATE_LENGTH_ERROR,
            Reserve::Bytes => &VECTOR_RESERVE_LENGTH_ERROR,
        }));
    }
    let bytes = match what {
        Reserve::String => length + 1,
        Reserve::Bytes => length,
    };
    if bytes >= USER_SPACE_BYTES {
        return Err(Fail::OutOfMemory);
    }
    if bytes >= HOST_DEPENDENT_BYTES {
        return Err(host_dependent(length, bytes, what, names));
    }
    Ok(())
}

/// The refusal of `length` bytes of `what` (an allocation of `bytes`) that [`reserve`] leaves to
/// the host.
#[cold]
#[inline(never)]
fn host_dependent(length: u64, bytes: u64, what: Reserve, names: Names) -> Fail {
    let what = match what {
        Reserve::String => "a string",
        Reserve::Bytes => "an OCTET STRING or BIT STRING",
    };
    let Names {
        input, allocator, ..
    } = names;
    Fail::unsupported(format!(
        "{what} whose length ({length} bytes) passes the end of the {input}: {allocator} \
         allocates {bytes} bytes for it, which succeeds or fails by the host's memory"
    ))
}

/// What reading a string posts, given to the caller when NCBI posts it.
#[derive(Debug)]
pub(crate) enum StringPost<'a> {
    /// The warning of a UTF8String tag for a VisibleString (binary, `ExpectStringTag`,
    /// objistrasnb.cpp:296-305; [`utf8_tag_message`]).
    Utf8Tag,
    /// `ReplaceVisibleChar`'s post of a byte read as `#`: where, the byte, and the string (or
    /// piece of it) that the message quotes.
    BadChar {
        position: StreamPosition,
        byte: u8,
        string: &'a [u8],
    },
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:296-310
/// The message of the warning 803.10 of a UTF8String tag for the VisibleString member at the
/// stack trace `trace` (once per run, [`Utf8TagOnce`]): the one owner of its text, for the mask
/// reader's and the table reader's posts.
pub(crate) fn utf8_tag_message(trace: &str) -> String {
    format!(
        "CObjectIStreamAsnBinary: UTF8String data for VisibleString member {trace}, ASN.1 \
         specification may need an update"
    )
}

/// `ERR_POST_X_ONCE` of the UTF8String tag warning ([`utf8_tag_message`]): NCBI posts it at the
/// first such tag its process reads and never again, whichever file or object holds the tag. The
/// one owner of that rule, for the table reader (`asn_types::read`) and makeblastdb's mask files.
/// One LOSAT run is one NCBI process, so the flag lives as long as the run's reader, not in a
/// static: a library that runs LOSAT twice, and the tests, run in one process.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Utf8TagOnce {
    posted: bool,
}

impl Utf8TagOnce {
    /// The flag of a run that has posted the warning. `CDelayBuffer::Update` re-reads bytes that
    /// the read before it passed, so the tags it meets were posted for (or after an earlier one).
    pub(crate) const POSTED: Self = Utf8TagOnce { posted: true };

    /// Whether this tag is the run's first, whose warning is posted; the next ones are not.
    pub(crate) fn first(&mut self) -> bool {
        !std::mem::replace(&mut self.posted, true)
    }
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
pub enum Encoding {
    /// `eSerial_AsnText`
    Text,
    /// `eSerial_AsnBinary`
    Binary,
    /// `eSerial_Xml` (NCBI's ASN.1-style XML, `asn_types::read_xml`)
    Xml,
}

/// The members of a SEQUENCE or the variants of a CHOICE, as an encoding looks them up: their
/// names in module order. NCBI's modules have automatic tags, so the tag of each is its position
/// (`[n]` in text, `a0 | n` in binary; `scripts/asn_spec_gen.py` checks it for the table).
pub(crate) trait Items {
    fn count(&self) -> usize;
    fn name(&self, index: usize) -> &'static str;
    /// Whether the member or variant at `index` is named `id`.
    fn is_named(&self, index: usize, id: &[u8]) -> bool {
        self.name(index).as_bytes() == id
    }
}

/// The named values of an ENUMERATED (`integer` false) or of an INTEGER with named values.
pub(crate) trait EnumItems {
    /// `IsInteger()`: any number is a value.
    fn integer(&self) -> bool;
    /// The value named `name`.
    fn value_of(&self, name: &[u8]) -> Option<i32>;
    /// Whether `value` is one of the named values.
    fn has_value(&self, value: i32) -> bool;
    /// `GetName()`: the name of the ENUMERATED's type, empty for an anonymous one (the XML
    /// stream's `ReadEnum` reads the value inside a tag of that name).
    fn type_name(&self) -> &'static str {
        ""
    }
}

/// The primitives of one encoding (`CObjectIStreamAsn`, `CObjectIStreamAsnBinary`,
/// `CObjectIStreamXml`) that NCBI's object stack calls (the names are NCBI's functions).
///
/// The XML stream's object stack has frames of its own and names its tags by them
/// (`SkipStackTagName`): the consts and the default methods below (`NAMED_TYPES` to
/// `drain_posts`) are what the object stack asks of an encoding for them. Their defaults are the
/// text and binary streams', which inline to nothing (no runtime flag: the stack is monomorphised
/// for each encoding).
pub(crate) trait Input {
    /// The encoding.
    const ENCODING: Encoding;
    /// Whether a named type that is neither a SEQUENCE, a CHOICE nor an ENUMERATED (an alias, a
    /// named SEQUENCE OF such as `Seq-descr`) is read in a frame of its own with a tag of its own
    /// (`ReadNamedType`: `BeginNamedType`, `EndNamedType`; the text and binary streams read it
    /// without either, `VIRTUAL_MID_LEVEL_IO`).
    const NAMED_TYPES: bool = false;
    /// Whether a SEQUENCE OF or SET OF of the table is read in an array frame (`eFrameArray`; the
    /// XML stream's `ReadContainer` pushes none for a container without a name, and every list
    /// type of the table is read as one: a named one through its named frame).
    const ARRAY_FRAMES: bool = true;
    /// Whether the stream throws `eMissingValue` (the XML stream's ENUMERATED without its
    /// attribute), which `ReadWithSetFlagMember` takes for an OPTIONAL or DEFAULT member.
    const MISSING_VALUES: bool = false;
    /// Whether the stream holds what it posts outside its strings until `drain_posts` (the XML
    /// stream), which the object stack drains at the end of each value it reads.
    const HELD_POSTS: bool = false;
    /// Whether the elements of a list of the type `element` are read in an element frame
    /// (`eFrameArrayElement`; the XML stream's only when the type will have no name,
    /// `WillHaveName`).
    #[inline(always)]
    fn element_frame(_element: u16) -> bool {
        true
    }
    /// The source file, line and function of the `END_OBJECT_FRAME` that leaves a frame of `kind`
    /// (of a skip function: `skip`), whose line `NCBI_RETHROW_SAME` adds to a `CException` other
    /// than a `CSerialException` passing through it.
    #[inline(always)]
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
        text_frame_site(kind, skip)
    }
    /// A frame of `kind` (of the table's type `ty`: none, `u16::MAX`, for a member frame, the
    /// element type for an element frame) was pushed; `depth` frames are on the stack.
    #[inline(always)]
    fn frame_pushed(&mut self, _kind: FrameKind, _ty: u16, depth: usize) {
        self.set_stack_depth(depth);
    }
    /// The member or variant of the top frame (`SetTopMemberId`).
    #[inline(always)]
    fn frame_member(&mut self, _member: &'static Member) {}
    /// The top frame was left; `depth` frames stay.
    #[inline(always)]
    fn frame_popped(&mut self, depth: usize) {
        self.set_stack_depth(depth);
    }
    /// `BeginNamedType` of the type of the named frame on top (`NAMED_TYPES`).
    #[inline(always)]
    fn begin_named(&mut self) -> R<()> {
        Ok(())
    }
    /// `EndNamedType`, before the named frame is left.
    #[inline(always)]
    fn end_named(&mut self) -> R<()> {
        Ok(())
    }
    /// `EndArrayElement` after each element of a list (the XML stream's closing tag of an
    /// element read in an element frame).
    #[inline(always)]
    fn end_element(&mut self) -> R<()> {
        Ok(())
    }
    /// `SkipNull` (the text and binary streams read what `ReadNull` reads).
    #[inline(always)]
    fn skip_null(&mut self) -> R<()> {
        self.read_null()
    }
    /// Moves what the stream posted outside the strings it read (the XML stream posts each byte
    /// it replaces at once, also in the attributes it skips) to `posts`, in the order posted.
    #[inline(always)]
    fn drain_posts(&mut self, _posts: &mut dyn PostSink) {}
    /// `GetPosition()`: `line N` in text, `byte N` (from the start of the object) in binary.
    fn position(&self) -> StreamPosition;
    /// The state of a new object stream at the next byte (each object of a file is read by a
    /// stream of its own in NCBI: its line count, or its byte count and data limits, start
    /// again).
    fn new_object(&mut self);
    /// `EndOfData()`: whether no object follows (binary: no byte left; text: nothing but white
    /// space and comments, or what `SkipWhiteSpace` stops at). Only the tests read several
    /// objects of one stream so far; makeblastdb's Seq-entry input will (line RG item D3c-3b).
    #[cfg(test)]
    fn end_of_data(&mut self) -> bool;
    /// The depth of NCBI's object stack (`GetStackDepth()`), which the text reader's
    /// `SkipWhiteSpace` asks at the end of the text.
    #[inline]
    fn set_stack_depth(&mut self, _depth: usize) {}
    /// `ReadFileHeader`: the type name that precedes an object in text; binary files have none
    /// (the first byte peeked at: an empty input is NCBI's `CEofException`).
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>>;
    fn begin_class(&mut self) -> R<()>;
    /// `BeginClassMember` of the class type and `pos`: the next member of `items`, searched from `pos` on,
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
    /// posts while it reads it given to `posts` when it is posted (also when the read then
    /// fails).
    fn read_string(
        &mut self,
        kind: StringKind,
        out: &mut Vec<u8>,
        posts: &mut dyn FnMut(StringPost<'_>),
    ) -> R<()>;
    /// `SkipString` (`SkipStringStore`).
    fn skip_string(&mut self, kind: StringKind, posts: &mut dyn FnMut(StringPost<'_>)) -> R<()>;
    /// An OCTET STRING (`ReadBytes` of a `vector<char>`), `SkipByteBlock`.
    fn read_octets(&mut self) -> R<Vec<u8>>;
    fn skip_octets(&mut self) -> R<()>;
    /// `ReadEnum` of the enum's `values` (which `SkipEnum` calls too): the value (`CEnumeratedTypeInfo::ReadEnum`
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
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
const SKIP_BYTE_BLOCK: ThrowSite =
    asn_site(1192, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::SkipByteBlock");
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:187-205
const EXPECT_CHOICE: ThrowSite = asn_site(203, FORMAT_ERROR, "ncbi::CObjectIStreamAsn::Expect");

impl<S: AsnSource> TextInput<S> {
    /// The text input over `src`, which LOSAT's own messages call by `names` (`Reader::new`).
    pub(crate) fn new(src: S, names: Names) -> Self {
        TextInput {
            r: Reader::new(src, names.document, names.reader),
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

    /// `ThrowError` with `fFormatError` and a message, of the text stream at `site`.
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
    /// `ReadMemberId` of what `SkipWhiteSpace` returns, and the member or variant it names among `items`,
    /// searched from `pos` on; one that is not there is an error (`UnexpectedMember`, which lists
    /// them all). `choice`: an empty name is `BeginChoiceVariant`'s error.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1305-1314
    #[inline(always)]
    fn find_member(&mut self, items: &(impl Items + ?Sized), pos: usize, choice: bool) -> R<usize> {
        let count = items.count();
        let member = self.r.member_or_tag().map_err(|e| self.fail(e))?;
        let (id, found) = match member {
            MemberId::Tag(digits) => {
                // `NStr::StringToInt` of `id`, the tag's digits: a number that does not fit is its
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
                // A plain loop: the lookup, with the compare it inlines, stays in the caller (an
                // iterator's `find` is compiled as a function of its own).
                let mut found = None;
                let mut i = pos;
                while i < count {
                    if items.is_named(i, id) {
                        found = Some(i);
                        break;
                    }
                    i += 1;
                }
                (id, found)
            }
        };
        match found {
            Some(index) => Ok(index),
            None => Err(self.unexpected_member(id, items)),
        }
    }

    /// `UnexpectedMember` of `id` and the items, of the members or variants `items`.
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, id: &[u8], items: &(impl Items + ?Sized)) -> Fail {
        let mut message = [b"\"", id, b"\": unexpected member, should be one of: "].concat();
        for i in 0..items.count() {
            message.extend_from_slice(format!("\"{}\" ", items.name(i)).as_bytes());
        }
        self.throw(UNEXPECTED_MEMBER, message)
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

    /// `Expect` of `c` with `true`.
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

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:98-109
    #[cfg(test)]
    fn end_of_data(&mut self) -> bool {
        self.r.skip_white_space().is_err()
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

    #[inline(always)]
    fn begin_element(&mut self) -> R<bool> {
        let more = self.r.next_element(&mut self.block_start);
        more.map_err(|e| self.fail(e))
    }

    fn end_container(&mut self) -> R<()> {
        self.end_class()
    }

    #[inline]
    fn read_i32(&mut self) -> R<i32> {
        self.r.int4().map_err(|e| self.fail(e))
    }

    #[inline]
    fn read_i64(&mut self) -> R<i64> {
        self.r.int8().map_err(|e| self.fail(e))
    }

    #[inline]
    fn read_u32(&mut self) -> R<u32> {
        self.r.uint4().map_err(|e| self.fail(e))
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

    /// `ReadDouble` (`Reader::real`).
    fn read_real(&mut self) -> R<f64> {
        let value = self.r.real();
        value.map_err(|e| self.fail(e))
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
        posts: &mut dyn FnMut(StringPost<'_>),
    ) -> R<()> {
        self.string(kind, Some(out), posts)
    }

    fn skip_string(&mut self, kind: StringKind, posts: &mut dyn FnMut(StringPost<'_>)) -> R<()> {
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
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:123-160 (`FindValue`, `FindName`;
    // their reports name no function)
    /// `ReadEnum` of the enum's `values`: a lower-case name of the list, or a number (one of the list for an
    /// ENUMERATED).
    fn read_enum(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32> {
        let word = self.r.enum_word();
        let value = match word.map_err(|e| self.fail(e))? {
            EnumWord::Name(id) => {
                let name = self.r.bytes(id);
                return match values.value_of(name) {
                    Some(value) => Ok(value),
                    None => Err(invalid_enum(
                        FIND_VALUE,
                        String::from_utf8_lossy(name).into_owned(),
                    )),
                };
            }
            EnumWord::Number(value) => value,
        };
        if !values.integer() && !values.has_value(value) {
            return Err(invalid_enum(FIND_NAME, value.to_string()));
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
    /// the largest `size_t` count of bytes (-1 cast to `size_t`) that `CTempString::operator string` asks for: `std::length_error`.
    fn read_other_pointer(&mut self) -> R<Vec<u8>> {
        let id = self.r.type_id();
        match id.map_err(|e| self.fail(e))? {
            Some(id) => Ok(self.r.bytes(id).to_vec()),
            None => Err(Fail::Std(&STRING_CREATE_LENGTH_ERROR)),
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
        posts: &mut dyn FnMut(StringPost<'_>),
    ) -> R<()> {
        let skip = out.is_none();
        let mut none = Vec::new();
        let out = out.unwrap_or(&mut none);
        let result = self.r.visible_string(
            out,
            kind != StringKind::Utf8,
            skip,
            &mut |_, line, byte, piece| {
                posts(StringPost::BadChar {
                    position: StreamPosition::Line(line as u64),
                    byte,
                    string: piece,
                })
            },
        );
        result.map_err(|e| self.fail(e))
    }
}
