//! BER decoding of the header record of a BLAST database (`.phr` / `.nhr`): the reader of
//! the encoding (`Reader`, shared with the ordinal-id decoder of `seqdb::header`) and the
//! decoder of a whole `Blast-def-line-set` with every member and every Seq-id alternative
//! (`decode_full_set`). An INTEGER is read as NCBI's `ReadInt4` / `ReadInt8` read it, at the
//! width of the member's NCBI type (`read_std_signed` of the binary object stream,
//! `report::serial::asn_ber_reader`, whose tag bytes it shares too).
//!
//! `decode_full_set` reads the record as NCBI's object stream does (`CObjectIStreamAsnBinary`
//! with the type information of `Blast-def-line-set`): where that read fails, its error is
//! NCBI's `CSerialException` ([`BerError::exception`]), made where NCBI throws, with NCBI's
//! position, message and code, the frames of NCBI's object stack that the exception leaves, and
//! the backlog that each member with a set flag adds (`error while reading <member>`). The
//! ordinal-id decoder of `seqdb::header` fails wherever this read fails, and takes its error.
//!
//! NCBI writes the record with `MSerial_AsnBinary`: EXPLICIT context tags, indefinite length
//! (`80 ... 00 00`) on every constructed value, definite length on INTEGER and VisibleString,
//! minimal two's-complement integers. NCBI's reader (`CObjectIStreamAsnBinary`) reads both
//! length forms, a constructed tag in the long form (`Reader::tag`) and a UTF8String where a
//! VisibleString is expected (`Reader::visible_string`, with NCBI's warning once per process,
//! [`utf8_string_path`]); so does this one.
//!
//! ASN.1 definitions (NCBI 598d8ae6): `c++/src/objects/blastdb/blastdb.asn:31-47`
//! (`Blast-def-line-set`, `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-74`
//! (`Seq-id` and its member types), `c++/src/objects/general/general.asn:54-80` (`Date`,
//! `Dbtag`, `Object-id`), `c++/src/objects/biblio/biblio.asn:119-124` (`Id-pat`).
//!
//! NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:32-47

use super::defline_set::{BlastDefline, DeflineSet};
use super::model::{
    Date, DateStd, Dbtag, GiimportId, IdPat, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId,
    PdbSeqId, TextseqId,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_ber_reader::{
    read_std_signed, StdSignedError, TAG_APPLICATION_INTEGER, TAG_INTEGER, TAG_UTF8_STRING,
    TAG_VISIBLE_STRING, UNIVERSAL_CONSTRUCTED as TAG_CONSTRUCTED,
};
pub(crate) use crate::report::serial::asn_ber_reader::{
    CONTEXT_CONSTRUCTED as TAG_CONTEXT_CONSTRUCTED, TAG_SEQUENCE,
};
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1032-1044
// NCBI reference (598d8ae6): c++/src/serial/member.cpp:464-477
// `ExpectedMember`'s throw: the members' read functions (`ReadMissingSimpleMember`) call it for a
// mandatory member that the record skips; it is the read check that `SERIAL_VERIFY_DATA_READ`
// turns off (`object_frame::DataCheck`, the database reader's refusal:
// `seqdb::SeqDb::refuse_unchecked_reads`).
use crate::report::serial::object_frame::{frame_info, frame_name, FrameKind, EXPECTED_MEMBER};
pub(crate) const TAG_NUMBER_MASK: u8 = 31;
const INDEFINITE_LENGTH_BYTE: u8 = 0x80;
/// Deepest nesting accepted inside a skipped member.
const MAX_SKIP_DEPTH: u32 = 8;

/// Member tags of `Blast-def-line` (blastdb.asn:34-46, in order).
pub(crate) const DEFLINE_TITLE: u8 = 0;
pub(crate) const DEFLINE_SEQID: u8 = 1;
pub(crate) const DEFLINE_TAXID: u8 = 2;
const DEFLINE_MEMBERSHIPS: u8 = 3;
const DEFLINE_LINKS: u8 = 4;
const DEFLINE_OTHER_INFO: u8 = 5;
pub(crate) const DEFLINE_LAST_MEMBER: u8 = DEFLINE_OTHER_INFO;
/// The names of the members of `Blast-def-line`, by tag.
pub(crate) const DEFLINE_MEMBERS: [&str; 6] = [
    "title",
    "seqid",
    "taxid",
    "memberships",
    "links",
    "other-info",
];

// The places where NCBI's object stream throws while it reads a header record (NCBI 598d8ae6,
// `CObjectIStreamAsnBinary` and `CObjectIStream`), each cited below where this reader throws
// it. `ThrowError` (objistr.cpp:829-852) makes the message `byte <N>: <text>` (`GetPosition`,
// objistr.cpp:807-819), N the position of NCBI's stream in the record; its flag gives the code
// (`fFormatError`: eFormatError, `fOverflow`: eOverflow, `fEOF`: eEOF). The functions are what
// NCBI's report names (`ReadStdSigned`, a template, names none); every place was checked with
// NCBI 2.17.0 (segmasker, dustmasker and blastdbcmd on edited records, RG-T5d probes
// `T/t5d/survey/probe` and `T/t5d/impl/probe1`).
const OBJISTRASNB: &str = "src/serial/objistrasnb.cpp";
const OBJISTR: &str = "src/serial/objistr.cpp";
const FORMAT_ERROR: &str = "CSerialException::eFormatError";
const OVERFLOW_ERROR: &str = "CSerialException::eOverflow";

const fn site(
    file: &'static str,
    line: u32,
    kind: &'static str,
    function: Option<&'static str>,
) -> ThrowSite {
    ThrowSite {
        file,
        line,
        kind,
        function,
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:174-179
const UNEXPECTED_TAG_CLASS_BYTE: ThrowSite = site(
    OBJISTRASNB,
    177,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedTagClassByte"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:151-172
// NCBI 2.17.0 reports the throw at line 161 (probe RG-T5d review l07).
const PEEK_LONG_TAG: ThrowSite = site(
    OBJISTRASNB,
    161,
    OVERFLOW_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::PeekLongTag"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:195-198
const UNEXPECTED_CONTINUATION: ThrowSite = site(
    OBJISTRASNB,
    197,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedContinuation"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:200-203
const UNEXPECTED_LONG_LENGTH: ThrowSite = site(
    OBJISTRASNB,
    202,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedLongLength"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:205-217
const UNEXPECTED_TAG_VALUE: ThrowSite = site(
    OBJISTRASNB,
    216,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedTagValue"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:374-379
const UNEXPECTED_SYS_TAG_BYTE: ThrowSite = site(
    OBJISTRASNB,
    377,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedSysTagByte"),
);
/// `ReadLengthLong` (cited at `Reader::read_length_long`).
const READ_LENGTH_LONG: Option<&str> = Some("ncbi::CObjectIStreamAsnBinary::ReadLengthLong");
const INDEFINITE_LENGTH: ThrowSite = site(OBJISTRASNB, 398, FORMAT_ERROR, READ_LENGTH_LONG);
const LENGTH_OVERFLOW: ThrowSite = site(OBJISTRASNB, 401, OVERFLOW_ERROR, READ_LENGTH_LONG);
const ILLEGAL_LENGTH_START: ThrowSite = site(OBJISTRASNB, 405, FORMAT_ERROR, READ_LENGTH_LONG);
/// `ReadStdSigned` (cited at `read_std_signed`): its report names no function.
const ZERO_LENGTH_OF_NUMBER: ThrowSite = site(OBJISTRASNB, 480, FORMAT_ERROR, None);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1057-1067
const UNEXPECTED_MEMBER: ThrowSite = site(
    OBJISTRASNB,
    1066,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedMember"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608
const HANDLE_EOF: ThrowSite = site(
    OBJISTR,
    606,
    "CSerialException::eEOF",
    Some("ncbi::CObjectIStream::HandleEOF"),
);
/// `HandleEOF`'s `NCBI_RETHROW_SAME` (the code stays `CEofException::eEof`).
const HANDLE_EOF_RETHROW: (&str, u32, &str) = (OBJISTR, 604, "ncbi::CObjectIStream::HandleEOF");
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:352-358
const END_OF_FILE: ThrowSite = site(
    "src/util/strbuffer.cpp",
    358,
    "CEofException::eEof",
    Some("ncbi::CIStreamBuffer::FillBuffer"),
);

// NCBI reference (598d8ae6): c++/include/serial/impl/objstack.hpp:201-214
// NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:399-451
/// A frame of NCBI's object stack that the read of a header record enters (`BEGIN_OBJECT_FRAME`):
/// the root, a container (`ReadContainer`) and its elements, a class (`ReadClassSequential`)
/// and its members, a CHOICE's variant (`ReadChoice`; the frame of the CHOICE itself has no
/// name, and no read happens in it). An exception that leaves a frame gets its name
/// (`AddFrameInfo`); the end of the record becomes NCBI's exception in the innermost frame
/// (`HandleEOF`, whose text is the frame's `GetFrameInfo`).
#[derive(Clone, Copy, Debug)]
enum Frame {
    /// `eFrameNamed`: the root, `Blast-def-line-set`.
    Root,
    /// `eFrameArray` of the set, the second frame of the stack: there `HandleEOF` rethrows the
    /// `CEofException` (a `GetStackDepth` below 2 after the frame is left).
    SetArray,
    /// `eFrameArray` of a member's container (its type has no name).
    Array,
    /// `eFrameArrayElement`.
    ArrayElement,
    /// `eFrameClass` of the class type named.
    Class(&'static str),
    /// `eFrameClassMember`: the member last begun (`SetTopMemberId`), none before the first.
    Member(Option<&'static str>),
    /// `eFrameChoiceVariant`: the variant, none before it is known.
    Variant(Option<&'static str>),
}

/// The type name of the root (`Blast-def-line-set ::= SEQUENCE OF Blast-def-line`).
const SET_TYPE: &str = "Blast-def-line-set";

impl Frame {
    /// The frame as NCBI's object stack holds it: its type (`EFrameType`), the name of its type
    /// (empty for the containers' anonymous types) and its member or variant
    /// (`report::serial::object_frame`, the one owner of the frames' texts).
    fn parts(self) -> (FrameKind, Option<&'static str>, Option<&'static str>) {
        match self {
            Frame::Root => (FrameKind::Named, Some(SET_TYPE), None),
            Frame::SetArray | Frame::Array => (FrameKind::Array, Some(""), None),
            Frame::ArrayElement => (FrameKind::ArrayElement, None, None),
            Frame::Class(name) => (FrameKind::Class, Some(name), None),
            Frame::Member(member) => (FrameKind::ClassMember, None, member),
            Frame::Variant(member) => (FrameKind::ChoiceVariant, None, member),
        }
    }

    /// `GetFrameName`, which an exception that leaves the frame adds to its frames
    /// (`object_frame::frame_name`).
    fn name(self) -> &'static str {
        let (kind, type_name, member) = self.parts();
        frame_name(kind, type_name, member)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objstack.cpp:79-119
    /// What the frame adds to `GetStackTraceASN` (after a `.`; the root's is the start).
    fn trace_name(self) -> Option<&'static str> {
        match self {
            Frame::Root => Some(SET_TYPE),
            Frame::ArrayElement => Some("E"),
            Frame::Member(Some(name)) | Frame::Variant(Some(name)) => Some(name),
            Frame::SetArray
            | Frame::Array
            | Frame::Class(_)
            | Frame::Member(None)
            | Frame::Variant(None) => None,
        }
    }

    /// `GetFrameInfo` (`object_frame::frame_info`).
    fn info(self) -> String {
        let (kind, type_name, member) = self.parts();
        frame_info(kind, type_name, member)
    }
}

/// Why the decode of a header record stops.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Failure {
    /// NCBI's `CSerialException` as it leaves the frames of NCBI's object stack: the exception
    /// with its backlogs (`NCBI_RETHROW_SAME`) and its frame stack (`m_FrameStack`, which
    /// `AddFrameInfo` extends, `NcbiException::add_frame_infos`, and `ReportExtra` prints as
    /// ` at <frames>`).
    Serial { exception: NcbiException },
    /// The record ends where NCBI's stream needs a byte (the `CEofException` of
    /// `CIStreamBuffer::FillBuffer`), not yet in a frame: the innermost frame makes it NCBI's
    /// exception (`HandleEOF`), with the position of NCBI's stream then: the record's length
    /// where a byte is read, the tag's first byte where a byte after it is peeked
    /// (`PeekTagByte` of a later byte of a tag in the long form leaves the stream where it is).
    Eof { pos: usize },
    /// The `CEofException` that `HandleEOF` rethrows in the outer two frames (the end of the
    /// record where the set's container is read, outside its elements): no frames, and a
    /// backlog for each of the two.
    Rethrown(NcbiException),
    /// A failure whose report by NCBI's reader is not ported, with LOSAT's text: the checks of
    /// the ordinal-id decoder (`seqdb::header`), which NCBI's read replaces where it fails too.
    Unported { pos: usize, message: String },
    /// Not a failure of NCBI's read: the read that looks for the first UTF8String where a
    /// VisibleString is expected (`utf8_string_path`) stops there, and each frame it leaves
    /// adds its part of NCBI's `GetStackTraceASN` to `path` (`Frame::trace_name`).
    Utf8 { path: String },
}

/// An error of the BER decoding of a header record: NCBI's exception
/// ([`BerError::exception`]), or LOSAT's text where NCBI's report is not ported (the
/// `Display` text).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BerError(Box<Failure>);

impl BerError {
    /// NCBI's exception of this failure (`CSerialException`, or a `CEofException` of the
    /// stream): what NCBI's application reports where its read of the record fails. `None`:
    /// NCBI's report is not ported (LOSAT's refusal, in the `Display` text).
    pub fn exception(&self) -> Option<NcbiException> {
        match &*self.0 {
            Failure::Serial { exception } => Some(exception.clone()),
            Failure::Rethrown(exception) => Some(exception.clone()),
            Failure::Eof { .. } | Failure::Unported { .. } | Failure::Utf8 { .. } => None,
        }
    }

    /// The exception leaves `frame` (`END_OBJECT_FRAME`): a `CSerialException` gets its name;
    /// the end of the record becomes NCBI's exception there (`HandleEOF`).
    #[cold]
    #[inline(never)]
    fn left(self, frame: Frame) -> Self {
        let failure = match *self.0 {
            Failure::Serial { mut exception } => {
                exception.add_frame_infos(&[frame.name()]);
                Failure::Serial { exception }
            }
            Failure::Eof { pos } => match frame {
                Frame::SetArray | Frame::Root => {
                    let (file, line, function) = HANDLE_EOF_RETHROW;
                    let eof = NcbiException::new(END_OF_FILE, "end of file").rethrow_same(
                        file,
                        line,
                        function,
                        Frame::SetArray.info(),
                    );
                    match frame {
                        Frame::Root => Failure::Rethrown(eof.rethrow_same(
                            file,
                            line,
                            function,
                            Frame::Root.info(),
                        )),
                        _ => Failure::Rethrown(eof),
                    }
                }
                _ => Failure::Serial {
                    exception: NcbiException::new(
                        HANDLE_EOF,
                        format!("byte {pos}: {}", frame.info()),
                    ),
                },
            },
            Failure::Rethrown(exception) => match frame {
                Frame::Root => {
                    let (file, line, function) = HANDLE_EOF_RETHROW;
                    Failure::Rethrown(exception.rethrow_same(
                        file,
                        line,
                        function,
                        Frame::Root.info(),
                    ))
                }
                _ => Failure::Rethrown(exception),
            },
            unported @ Failure::Unported { .. } => unported,
            Failure::Utf8 { mut path } => {
                if let Some(name) = frame.trace_name() {
                    if !path.is_empty() {
                        path.insert(0, '.');
                    }
                    path.insert_str(0, name);
                }
                Failure::Utf8 { path }
            }
        };
        BerError(Box::new(failure))
    }

    /// The read of a member with a set flag fails (`ReadWithSetFlagMember`): a
    /// `CSerialException` gets the backlog `error while reading <member>`, whose predecessor
    /// keeps the frames so far. The end of the record (a `CEofException`) passes.
    #[cold]
    #[inline(never)]
    fn member_read(self, member: &'static str) -> Self {
        match *self.0 {
            Failure::Serial { exception } => BerError(Box::new(Failure::Serial {
                exception: crate::report::serial::asn_types::read::rethrow_set_flag_member(
                    exception, member,
                ),
            })),
            failure => BerError(Box::new(failure)),
        }
    }
}

impl std::fmt::Display for BerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &*self.0 {
            Failure::Unported { pos, message } => {
                write!(f, "header record, byte {pos}: {message}")
            }
            Failure::Eof { pos } => {
                write!(
                    f,
                    "header record, byte {pos}: the record ends inside a value"
                )
            }
            Failure::Utf8 { path } => {
                write!(
                    f,
                    "header record: a UTF8String for VisibleString member {path}"
                )
            }
            Failure::Serial { .. } | Failure::Rethrown(_) => {
                let exception = self.exception().expect("NCBI's exception");
                write!(
                    f,
                    "header record: {}",
                    String::from_utf8_lossy(&exception.what()).trim_end()
                )
            }
        }
    }
}

impl std::error::Error for BerError {}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:318-372
/// `TagToString`: how NCBI's reader names a tag byte in its messages.
fn tag_to_string(byte: u8) -> String {
    let class = match byte & 0xc0 {
        0x40 => "application/",
        0x80 => "contextspecific/",
        0xc0 => "private/",
        _ => "",
    };
    let constructed = if byte & 0x20 != 0 { "constructed/" } else { "" };
    let value = byte & TAG_NUMBER_MASK;
    let value = if byte & 0xc0 == 0 {
        // The universal tags of objstrasnb.hpp:81-124 (`ETagValue`) that `TagToString` names.
        match value {
            0 => "None",
            1 => "Boolean",
            2 => "Integer",
            3 => "BitString",
            4 => "OctetString",
            5 => "Null",
            6 => "ObjectIdentifier",
            7 => "ObjectDescriptor",
            8 => "External",
            9 => "Real",
            10 => "Enumerated",
            12 => "UTF8String",
            16 => "Sequence",
            17 => "Set",
            18 => "NumericString",
            19 => "PrintableString",
            20 => "TeletextString",
            21 => "VideotextString",
            22 => "IA5String",
            23 => "UTCTime",
            24 => "GeneralizedTime",
            25 => "GraphicString",
            26 => "VisibleString",
            27 => "GeneralString",
            29 => "MemberReference",
            30 => "ObjectReference",
            _ => "unknown",
        }
        .to_string()
    } else {
        value.to_string()
    };
    format!("{class}{constructed}{value} ({byte})")
}

/// Where a constructed value ends.
#[derive(Clone, Copy)]
pub(crate) enum End {
    /// Indefinite length: two zero bytes close the value.
    Marker,
    /// Definite length: the value ends at this position.
    At(usize),
}

/// The identifier of a value as NCBI's `PeekTag` reads it: the class and form bits of its first
/// byte, its number (in the long form, from the bytes after: `PeekLongTag`), and how many bytes
/// it takes (`m_CurrentTagLength`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Tag {
    /// The first byte without the number (`first & !TAG_NUMBER_MASK`).
    pub(crate) class: u8,
    /// `TLongTag` (an `Int4`; `PeekLongTag` keeps it below 2^31).
    pub(crate) number: u32,
    len: u8,
}

impl Tag {
    /// The identifier whose first byte is `first` in the short form (`PeekTag`'s one byte,
    /// `m_CurrentTagLength` set to 1); `None` in the long form (`PeekLongTag`).
    #[inline(always)]
    fn short(first: u8) -> Option<Tag> {
        (first & TAG_NUMBER_MASK != TAG_NUMBER_MASK).then_some(Tag {
            class: first & !TAG_NUMBER_MASK,
            number: u32::from(first & TAG_NUMBER_MASK),
            len: 1,
        })
    }

    /// Whether this is the tag of the one-byte identifier `byte` (whose number is below 31):
    /// the same class and form, and the same number in either form.
    pub(crate) fn is(self, byte: u8) -> bool {
        self.class == byte & !TAG_NUMBER_MASK && self.number == u32::from(byte & TAG_NUMBER_MASK)
    }

    /// Whether the value is constructed (`eConstructed`).
    pub(crate) fn constructed(self) -> bool {
        self.class & TAG_CONSTRUCTED != 0
    }
}

/// A reader of one BER-encoded record.
pub(crate) struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    /// Whether a UTF8String was read where a VisibleString is expected (`visible_string`).
    utf8: bool,
    /// Whether such a string stops the read (`Failure::Utf8`; `utf8_string_path`).
    find_utf8: bool,
    /// Where the read keeps the titles that hold a byte outside ' ' ..= '~' (the range of each
    /// stored string in the record; `title_replaced_chars`), else `None`.
    replaced_titles: Option<Vec<std::ops::Range<usize>>>,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Reader {
            data,
            pos: 0,
            utf8: false,
            find_utf8: false,
            replaced_titles: None,
        }
    }

    /// The title just read (`length` bytes before the position) holds a byte that its read
    /// replaced: kept where the read keeps such titles (`title_replaced_chars`).
    #[cold]
    #[inline(never)]
    fn note_replaced_title(&mut self, length: usize) {
        let end = self.pos;
        if let Some(titles) = self.replaced_titles.as_mut() {
            titles.push(end - length..end);
        }
    }

    /// Whether the reader read a UTF8String where a VisibleString is expected, which NCBI's
    /// read takes with a warning (`visible_string`; `utf8_string_path` names the member).
    pub(crate) fn read_utf8(&self) -> bool {
        self.utf8
    }

    /// Whether every byte of the record was read.
    #[cfg(test)]
    pub(crate) fn at_record_end(&self) -> bool {
        self.pos == self.data.len()
    }

    /// LOSAT's own failure at the reader's position (`Failure::Unported`): a check of the
    /// ordinal-id decoder.
    #[cold]
    #[inline(never)]
    pub(crate) fn error(&self, message: impl Into<String>) -> BerError {
        BerError(Box::new(Failure::Unported {
            pos: self.pos,
            message: message.into(),
        }))
    }

    /// NCBI's `ThrowError` at `site` with the stream at `pos`.
    #[cold]
    #[inline(never)]
    fn throw(&self, site: ThrowSite, pos: usize, message: &str) -> BerError {
        BerError(Box::new(Failure::Serial {
            exception: NcbiException::new(site, format!("byte {pos}: {message}")),
        }))
    }

    /// The end of the record where a byte is needed.
    #[cold]
    #[inline(never)]
    fn eof(&self) -> BerError {
        BerError(Box::new(Failure::Eof {
            pos: self.data.len(),
        }))
    }

    fn byte(&mut self) -> Result<u8, BerError> {
        let Some(&byte) = self.data.get(self.pos) else {
            return Err(self.eof());
        };
        self.pos += 1;
        Ok(byte)
    }

    /// `PeekTagByte`: the byte at the position, which stays.
    fn peek(&self) -> Result<u8, BerError> {
        match self.data.get(self.pos) {
            Some(&byte) => Ok(byte),
            None => Err(self.eof()),
        }
    }

    fn take(&mut self, len: u64) -> Result<&'a [u8], BerError> {
        let end = usize::try_from(len)
            .ok()
            .and_then(|len| self.pos.checked_add(len))
            .filter(|&end| end <= self.data.len());
        let Some(end) = end else {
            return Err(self.eof());
        };
        let slice = &self.data[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:175-189
    /// `PeekTag` of the identifier whose first byte `first` is at the position, which stays
    /// (`PeekTagByte` does not move the stream): the one tag reader of both decoders.
    #[inline]
    fn peek_tag(&self, first: u8) -> Result<Tag, BerError> {
        match Tag::short(first) {
            Some(tag) => Ok(tag),
            None => self.peek_long_tag(first),
        }
    }

    /// `PeekLongTag` (cited at `PEEK_LONG_TAG`): the number in the bytes after `first`, seven
    /// bits each, the last without the high bit; NCBI's overflow before a byte that could pass
    /// 2^31, and the end of the record at the tag's first byte.
    #[cold]
    #[inline(never)]
    fn peek_long_tag(&self, first: u8) -> Result<Tag, BerError> {
        let mut tag: u32 = 0;
        let mut i = 1;
        loop {
            if tag >= 1 << 24 {
                let message = format!("tag number is too big: {tag}");
                return Err(self.throw(PEEK_LONG_TAG, self.pos, &message));
            }
            let Some(&byte) = self.data.get(self.pos + i) else {
                return Err(BerError(Box::new(Failure::Eof { pos: self.pos })));
            };
            i += 1;
            tag = (tag << 7) | u32::from(byte & 0x7f);
            if byte & 0x80 == 0 {
                break;
            }
        }
        Ok(Tag {
            class: first & !TAG_NUMBER_MASK,
            number: tag,
            len: i as u8,
        })
    }

    /// Reads one identifier (`peek_tag`), short or long.
    #[inline]
    pub(crate) fn tag(&mut self) -> Result<Tag, BerError> {
        let first = self.peek()?;
        if let Some(tag) = Tag::short(first) {
            self.pos += 1;
            return Ok(tag);
        }
        self.long_tag(first)
    }

    /// `tag` of an identifier in the long form.
    #[cold]
    #[inline(never)]
    fn long_tag(&mut self, first: u8) -> Result<Tag, BerError> {
        let tag = self.peek_long_tag(first)?;
        self.pos += usize::from(tag.len);
        Ok(tag)
    }

    /// Reads the next identifier where it is that of the one-byte tag `expected` in either form
    /// (`Tag::is`): whether it was.
    #[inline]
    pub(crate) fn take_tag(&mut self, expected: u8) -> Result<bool, BerError> {
        if self.data.get(self.pos) == Some(&expected) {
            self.pos += 1;
            return Ok(true);
        }
        self.take_other_tag(expected)
    }

    /// `take_tag` of another first byte: the long form of `expected`, or another tag (the
    /// position stays), or the end of the record.
    #[cold]
    #[inline(never)]
    fn take_other_tag(&mut self, expected: u8) -> Result<bool, BerError> {
        let at = self.pos;
        let tag = self.tag()?;
        if tag.is(expected) {
            return Ok(true);
        }
        self.pos = at;
        Ok(false)
    }

    /// Reads the identifier of the one-byte tag `expected` (`Tag::is`), in either form.
    pub(crate) fn expect_tag(&mut self, expected: u8, what: &str) -> Result<(), BerError> {
        if self.take_tag(expected)? {
            return Ok(());
        }
        let tag = self.tag()?;
        Err(self.error(format!(
            "tag {}/{} where {what} (0x{expected:02x}) is expected",
            tag.class, tag.number
        )))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:394-411
    /// `ReadLengthLong` of the first length byte `first` (0x80 or above, read): the length in
    /// the long form, as NCBI's 64-bit `size_t` holds it.
    fn read_length_long(&mut self, first: u8) -> Result<u64, BerError> {
        let mut length_length = first - INDEFINITE_LENGTH_BYTE;
        if length_length == 0 {
            return Err(self.throw(INDEFINITE_LENGTH, self.pos, "unexpected indefinite length"));
        }
        if usize::from(length_length) > std::mem::size_of::<u64>() {
            return Err(self.throw(LENGTH_OVERFLOW, self.pos, "length overflow"));
        }
        let mut length = u64::from(self.byte()?);
        if length == 0 {
            return Err(self.throw(ILLEGAL_LENGTH_START, self.pos, "illegal length start"));
        }
        length_length -= 1;
        while length_length > 0 {
            length = (length << 8) | u64::from(self.byte()?);
            length_length -= 1;
        }
        Ok(length)
    }

    /// `ReadLength`: a definite length (short or long form) of a primitive value.
    fn read_length(&mut self) -> Result<u64, BerError> {
        let first = self.byte()?;
        if first < INDEFINITE_LENGTH_BYTE {
            return Ok(u64::from(first));
        }
        self.read_length_long(first)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:277-312
    // `USE_DEF_LEN` is 1 (objistrasnb.hpp:46). Of the long form NCBI 2.17.0 takes the position
    // before the length bytes after the first (the limit of `a0 81 18` and 24 content bytes ends
    // one byte before them: `UnexpectedContinuation` at the member's end, probe c30).
    /// `ExpectIndefiniteLength` after a constructed tag: the end of the value, an indefinite
    /// length or a definite one (NCBI's reader takes both).
    fn expect_indefinite_length(&mut self) -> Result<End, BerError> {
        let first = self.byte()?;
        if first == INDEFINITE_LENGTH_BYTE {
            return Ok(End::Marker);
        }
        if first < INDEFINITE_LENGTH_BYTE {
            return Ok(End::At(self.pos + usize::from(first)));
        }
        let at = self.pos;
        let length = self.read_length_long(first)?;
        let length = usize::try_from(length).unwrap_or(usize::MAX);
        Ok(End::At(at.saturating_add(length)))
    }

    /// Reads the length of a constructed value as NCBI does
    /// (`expect_indefinite_length`); a definite length must stay inside the record.
    pub(crate) fn open(&mut self) -> Result<End, BerError> {
        let end = self.expect_indefinite_length()?;
        match end {
            End::At(end) if end > self.data.len() => {
                Err(self.error("a value is longer than the record"))
            }
            end => Ok(end),
        }
    }

    /// Whether the constructed value opened with `end` has no further
    /// element; consumes the end marker of an indefinite length.
    pub(crate) fn at_end(&mut self, end: End) -> Result<bool, BerError> {
        match end {
            End::At(limit) => {
                if self.pos > limit {
                    return Err(self.error("a value passes the end of its container"));
                }
                Ok(self.pos == limit)
            }
            End::Marker => {
                if self.pos >= self.data.len() {
                    return Err(self.eof());
                }
                if self.data.get(self.pos) == Some(&0) {
                    if self.data.get(self.pos + 1) != Some(&0) {
                        return Err(self.error("end-of-contents marker cut short"));
                    }
                    self.pos += 2;
                    return Ok(true);
                }
                Ok(false)
            }
        }
    }

    pub(crate) fn close(&mut self, end: End, what: &str) -> Result<(), BerError> {
        if self.at_end(end)? {
            Ok(())
        } else {
            Err(self.error(format!("unexpected data after {what}")))
        }
    }

    /// Skips one value of any type. `depth` bounds the nesting: the skipped
    /// members of `Blast-def-line` are an INTEGER or a SEQUENCE OF INTEGER.
    pub(crate) fn skip_value(&mut self, depth: u32) -> Result<(), BerError> {
        let tag = self.tag()?;
        if tag.constructed() {
            if depth >= MAX_SKIP_DEPTH {
                return Err(self.error("values nested too deeply"));
            }
            let end = self.open()?;
            while !self.at_end(end)? {
                self.skip_value(depth + 1)?;
            }
        } else {
            let length = self.read_length()?;
            self.take(length)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:281-312
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:785-790
    // NCBI's reader takes a UTF8String (tag 0x0c) where a VisibleString is expected: the
    // parameter `SERIAL_READ_ANY_VISIBLESTRING_TAG` is 1 by default (objistrasnb.cpp:273-279),
    // and its bytes are read as a VisibleString's (`x_FixCharsMethod`). The warning is posted
    // once per process (`ERR_POST_X_ONCE`, error code 803 `Serial_IStream`, error_codes.hpp:47),
    // at NCBI's first read of such a string: the programs post it where NCBI first reads the
    // record (`utf8_string_warning`, `utf8_string_path`; NCBI 2.17.0, probe RG-T5d review u01).
    /// `ReadString` of a VisibleString: the stored bytes. A UTF8String is read the same way
    /// (`read_utf8` tells it).
    pub(crate) fn visible_string(&mut self) -> Result<&'a [u8], BerError> {
        let tag = self.peek()?;
        if tag != TAG_VISIBLE_STRING {
            self.other_string_tag(tag)?;
        }
        self.pos += 1;
        let length = self.read_length()?;
        self.take(length)
    }

    /// The tag `tag` at the position where a VisibleString is expected (`ExpectStringTag`): a
    /// UTF8String, noted (`read_utf8`; with `find_utf8` the read stops there), or NCBI's
    /// `UnexpectedSysTagByte`.
    #[cold]
    #[inline(never)]
    fn other_string_tag(&mut self, tag: u8) -> Result<(), BerError> {
        if tag != TAG_UTF8_STRING {
            return Err(self.unexpected_sys_tag_byte(tag, TAG_VISIBLE_STRING));
        }
        if self.find_utf8 {
            return Err(BerError(Box::new(Failure::Utf8 {
                path: String::new(),
            })));
        }
        self.utf8 = true;
        Ok(())
    }

    /// `UnexpectedSysTagByte` of the expected byte `expected`, with the tag byte `got` at the
    /// position.
    #[cold]
    #[inline(never)]
    fn unexpected_sys_tag_byte(&self, got: u8, expected: u8) -> BerError {
        let message = format!(
            "unexpected tag: {}, should be: {}",
            tag_to_string(got),
            tag_to_string(expected)
        );
        self.throw(UNEXPECTED_SYS_TAG_BYTE, self.pos, &message)
    }

    /// `UnexpectedTagClassByte(got, expected_class)` with the tag byte at the position.
    #[cold]
    #[inline(never)]
    fn unexpected_tag_class_byte(&self, got: u8, expected_class: u8) -> BerError {
        let message = format!(
            "unexpected tag: {}, should be {}",
            tag_to_string(got),
            tag_to_string(expected_class)
        );
        self.throw(UNEXPECTED_TAG_CLASS_BYTE, self.pos, &message)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:145-160
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:337-345
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:569-575
    /// An INTEGER as NCBI's `ReadInt4` (`width` 4) or `ReadInt8` (`width` 8) reads it: the tag
    /// of an INTEGER or of an application INTEGER (`eReadAsBigInt` changes nothing in the value),
    /// a short length, the value by `read_std_signed`, with NCBI's errors.
    fn std_signed(&mut self, width: usize) -> Result<i64, BerError> {
        let tag = self.peek()?;
        if tag != TAG_INTEGER && tag != TAG_APPLICATION_INTEGER {
            return Err(self.unexpected_sys_tag_byte(tag, TAG_INTEGER));
        }
        self.pos += 1;
        let first = self.byte()?;
        if first >= INDEFINITE_LENGTH_BYTE {
            return Err(self.throw(UNEXPECTED_LONG_LENGTH, self.pos, "ShortLength expected"));
        }
        match read_std_signed(usize::from(first), width, || self.byte()) {
            Ok(value) => Ok(value),
            Err(StdSignedError::ZeroLength) => {
                Err(self.throw(ZERO_LENGTH_OF_NUMBER, self.pos, "zero length of number"))
            }
            Err(StdSignedError::Overflow { line }) => {
                let site = site(OBJISTRASNB, line, OVERFLOW_ERROR, None);
                Err(self.throw(site, self.pos, "overflow error"))
            }
            Err(StdSignedError::Read(error)) => Err(error),
        }
    }

    /// An INTEGER member that NCBI reads into an `int` (`ReadInt4`): the taxid (`TTaxId`), the
    /// lists of `Blast-def-line`, and every INTEGER of a Seq-id but the GI. A value outside
    /// `int` fails as in NCBI ("overflow error").
    pub(crate) fn int(&mut self) -> Result<i32, BerError> {
        let value = self.std_signed(std::mem::size_of::<i32>())?;
        // `read_std_signed` of width 4 keeps at most four content bytes: the value is an `int`.
        Ok(value as i32)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:585-591
    /// An INTEGER member that NCBI reads into an `Int8` (`ReadInt8`): the GI (`TGi`, an `Int8`
    /// with `NCBI_INT8_GI`, corelib/ncbimisc.hpp:975-1003).
    fn int8(&mut self) -> Result<i64, BerError> {
        self.std_signed(std::mem::size_of::<i64>())
    }

    pub(crate) fn open_context(&mut self, number: u8, what: &str) -> Result<End, BerError> {
        self.expect_tag(TAG_CONTEXT_CONSTRUCTED | number, what)?;
        self.open()
    }

    /// A VisibleString as NCBI's reader gives it (`fix_visible_chars`), written over `value`
    /// (whose buffer is kept), and the number of bytes it replaced.
    fn fixed_string_in(&mut self, value: &mut Vec<u8>) -> Result<usize, BerError> {
        let stored = self.visible_string()?;
        Ok(fix_visible_chars(stored, value))
    }

    /// Opens a SEQUENCE value of `type_name` (the ordinal-id decoder's lists).
    fn open_sequence(&mut self, type_name: &str) -> Result<End, BerError> {
        self.expect_tag(TAG_SEQUENCE, type_name)?;
        self.open()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1029-1044
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:210-219
    /// The start of a class or of a container (`BeginClass`, `ReadContainer`): the tag of a
    /// universal constructed SEQUENCE (`ExpectTag`, which also takes a SET's class and form),
    /// and the length (`ExpectIndefiniteLength`).
    fn begin_sequence(&mut self) -> Result<End, BerError> {
        let first = self.peek()?;
        if first == TAG_SEQUENCE {
            self.pos += 1;
        } else {
            self.other_sequence_tag(first)?;
        }
        self.expect_indefinite_length()
    }

    /// `begin_sequence`'s tag other than SEQUENCE's one byte: NCBI's errors, or SEQUENCE's
    /// number in the long form (read).
    #[cold]
    #[inline(never)]
    fn other_sequence_tag(&mut self, first: u8) -> Result<(), BerError> {
        if first & !TAG_NUMBER_MASK != TAG_CONSTRUCTED {
            return Err(self.unexpected_tag_class_byte(first, TAG_CONSTRUCTED));
        }
        let tag = self.peek_tag(first)?;
        let expected = u32::from(TAG_SEQUENCE & TAG_NUMBER_MASK);
        if tag.number != expected {
            let message = format!("Unexpected tag: {}, expected: {expected}", tag.number);
            return Err(self.throw(UNEXPECTED_TAG_VALUE, self.pos, &message));
        }
        self.pos += usize::from(tag.len);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1122-1160
    // The types of a header record have automatic tags (the cited path); a member's tag must be
    // that of a member after the last one read (`Find` from `pos`, memberlist.cpp:462-472), and
    // `SERIAL_SKIP_UNKNOWN_MEMBERS` (and `..._VARIANTS` of a CHOICE) is no by default.
    /// The tag of a member or variant (`[n]` EXPLICIT, context-specific constructed, in either
    /// form) and its length: the number `n` and the end of the value.
    fn context_tag(&mut self) -> Result<(u32, End), BerError> {
        let first = self.peek()?;
        // A context-specific constructed tag in the short form (`[0]` to `[30]`).
        let number = if first.wrapping_sub(TAG_CONTEXT_CONSTRUCTED) < TAG_NUMBER_MASK {
            self.pos += 1;
            u32::from(first & TAG_NUMBER_MASK)
        } else {
            self.other_context_tag(first)?
        };
        Ok((number, self.expect_indefinite_length()?))
    }

    /// `context_tag`'s other first bytes: NCBI's error of another class or form, or the number in
    /// the long form (read).
    #[cold]
    #[inline(never)]
    fn other_context_tag(&mut self, first: u8) -> Result<u32, BerError> {
        if first & !TAG_NUMBER_MASK != TAG_CONTEXT_CONSTRUCTED {
            return Err(self.unexpected_tag_class_byte(first, TAG_CONTEXT_CONSTRUCTED));
        }
        let tag = self.peek_tag(first)?;
        self.pos += usize::from(tag.len);
        Ok(tag.number)
    }

    /// `UnexpectedMember` of `tag` in the type whose members (or variants) are `names`.
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, tag: u32, names: &[&str]) -> BerError {
        let mut message = format!("unexpected member: [{tag}], should be one of: ");
        for (at, name) in names.iter().enumerate() {
            message.push_str(&format!("{name}[{at}] "));
        }
        self.throw(UNEXPECTED_MEMBER, self.pos, &message)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:243-251
    /// `HaveMoreElements` of the constructed value that ends at `end`.
    fn more(&self, end: End) -> Result<bool, BerError> {
        match end {
            End::At(limit) => Ok(limit > self.pos),
            End::Marker => Ok(self.peek()? != 0),
        }
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:374-405
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:73-88
    /// `ExpectEndOfContent` of the constructed value that ends at `end`: the first zero byte of
    /// the end marker is taken before the second is checked.
    fn end_of_content(&mut self, end: End) -> Result<(), BerError> {
        match end {
            End::At(limit) => {
                if limit != self.pos {
                    return Err(self.unexpected_continuation());
                }
            }
            End::Marker => {
                for _ in 0..2 {
                    if self.peek()? != 0 {
                        return Err(self.unexpected_continuation());
                    }
                    self.pos += 1;
                }
            }
        }
        Ok(())
    }

    #[cold]
    #[inline(never)]
    fn unexpected_continuation(&self) -> BerError {
        self.throw(
            UNEXPECTED_CONTINUATION,
            self.pos,
            "EndOfContentsByte expected",
        )
    }

    /// `ExpectedMember` of the mandatory member `name` that the record skips.
    #[cold]
    #[inline(never)]
    fn expected_member(&self, name: &str) -> BerError {
        self.throw(
            EXPECTED_MEMBER,
            self.pos,
            &format!("member {name} expected"),
        )
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:724-738
// A VisibleString is read with the method `eFNP_ReplaceAndWarn` unless the
// parameter SERIAL / WRONG_CHARS_READ says otherwise (objistr.cpp:302-303):
// every byte outside ' ' ..= '~' (`GoodVisibleChar`, objistr.inl:456-459)
// becomes the substitute '#', UTF-8 included
// (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0, objistr.hpp:1163).
// NCBI also writes a warning for each such string; this reader reports none.
// makeblastdb stores titles already repaired this way.
/// The stored VisibleString `value` as NCBI's reader gives it, written over `out`, and the number
/// of bytes it replaced. One pass makes the bytes and counts the replaced ones, `FIX_LANES` bytes
/// a step: each lane counts in a `u8`, added up every 255 steps (before it could pass 255), so
/// the step takes as many bytes as the vector registers hold. An unallocated `out` is allocated
/// once with `Vec::with_capacity` so fresh buffers skip `Vec::resize`'s amortized grow path; a
/// kept buffer is cut or lengthened in place without clearing bytes it already holds.
pub(crate) fn fix_visible_chars(value: &[u8], out: &mut Vec<u8>) -> usize {
    if out.capacity() == 0 {
        *out = Vec::with_capacity(value.len());
    }
    out.truncate(value.len());
    out.resize(value.len(), 0);
    let (steps, rest) = value.as_chunks::<FIX_LANES>();
    let (out_steps, out_rest) = out.split_at_mut(steps.len() * FIX_LANES);
    let (out_steps, _) = out_steps.as_chunks_mut::<FIX_LANES>();
    let block = usize::from(u8::MAX);
    let mut replaced = 0;
    for (steps, out_steps) in steps.chunks(block).zip(out_steps.chunks_mut(block)) {
        let mut lanes = [0u8; FIX_LANES];
        for (step, fixed) in steps.iter().zip(out_steps) {
            for lane in 0..FIX_LANES {
                let (byte, kept) = visible_byte(step[lane]);
                fixed[lane] = byte;
                lanes[lane] += u8::from(!kept);
            }
        }
        replaced += lanes.iter().map(|&n| usize::from(n)).sum::<usize>();
    }
    for (&stored, fixed) in rest.iter().zip(out_rest) {
        let (byte, kept) = visible_byte(stored);
        *fixed = byte;
        replaced += usize::from(!kept);
    }
    replaced
}

/// The bytes of one step of `fix_visible_chars`.
const FIX_LANES: usize = 16;

/// The byte `c` as `fix_visible_chars` gives it, and whether it is a byte of ' ' ..= '~', which
/// it keeps.
fn visible_byte(c: u8) -> (u8, bool) {
    let kept = (b' '..=b'~').contains(&c);
    (if kept { c } else { b'#' }, kept)
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:955-1000
/// `ReadContainer` of a SEQUENCE OF whose array frame is `array`: each element read by
/// `element` in the element frame.
fn read_container<'a>(
    ber: &mut Reader<'a>,
    array: Frame,
    mut element: impl FnMut(&mut Reader<'a>) -> Result<(), BerError>,
) -> Result<(), BerError> {
    let end = ber.begin_sequence().map_err(|e| e.left(array))?;
    loop {
        let more = ber
            .more(end)
            .map_err(|e| e.left(Frame::ArrayElement).left(array))?;
        if !more {
            break;
        }
        element(ber).map_err(|e| e.left(Frame::ArrayElement).left(array))?;
    }
    ber.end_of_content(end).map_err(|e| e.left(array))
}

/// The members of a class type that its read checks (generated `*_.cpp` of NCBI 598d8ae6, each
/// `ADD_NAMED_*_MEMBER` in order).
struct Class {
    /// The type's name (`eFrameClass`'s `Object type`).
    name: &'static str,
    /// The members' names, by tag (automatic tags: 0, 1, ...).
    members: &'static [&'static str],
    /// The members that are neither `SetOptional()` nor with a default (bit n: member n).
    mandatory: u8,
    /// The members with a set flag (`SetSetFlag`), whose read failure gets a backlog
    /// (`ReadWithSetFlagMember`); the others are read by `ReadSimpleMember`.
    set_flag: u8,
}

// NCBI reference (598d8ae6): c++/include/serial/impl/objistrimpl.hpp:86-111
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1209-1232
// The generated classes call `info->RandomOrder()`, a getter (classinfo.hpp:100) that leaves them
// sequential.
/// `ReadClassSequential` of `class`: each member present read by `member` (with its tag) in the
/// member frame, the mandatory members checked, the frames and backlogs of NCBI's read.
fn read_class<'a>(
    ber: &mut Reader<'a>,
    class: &Class,
    mut member: impl FnMut(&mut Reader<'a>, u8) -> Result<(), BerError>,
) -> Result<(), BerError> {
    let frame = Frame::Class(class.name);
    let end = ber.begin_sequence().map_err(|e| e.left(frame))?;
    let mut next = 0u8;
    let mut current = None;
    loop {
        let begun = begin_class_member(ber, end, class.members, next)
            .map_err(|e| e.left(Frame::Member(current)).left(frame))?;
        let Some((number, member_end)) = begun else {
            break;
        };
        let name = class.members[usize::from(number)];
        current = Some(name);
        let in_member = |e: BerError| e.left(Frame::Member(Some(name))).left(frame);
        if let Some(skipped) = (next..number).find(|&at| class.mandatory & (1 << at) != 0) {
            return Err(in_member(
                ber.expected_member(class.members[usize::from(skipped)]),
            ));
        }
        member(ber, number).map_err(|e| {
            in_member(if class.set_flag & (1 << number) != 0 {
                e.member_read(name)
            } else {
                e
            })
        })?;
        ber.end_of_content(member_end).map_err(in_member)?;
        next = number + 1;
    }
    let count = class.members.len() as u8;
    if let Some(skipped) = (next..count).find(|&at| class.mandatory & (1 << at) != 0) {
        let name = class.members[usize::from(skipped)];
        return Err(ber.expected_member(name).left(frame));
    }
    ber.end_of_content(end).map_err(|e| e.left(frame))
}

/// `BeginClassMember` of a class with the members `names` whose next member is `next` (NCBI's
/// `pos`): the member's tag and end, `None` at the end of the members.
fn begin_class_member(
    ber: &mut Reader<'_>,
    end: End,
    names: &[&str],
    next: u8,
) -> Result<Option<(u8, End)>, BerError> {
    if !ber.more(end)? {
        return Ok(None);
    }
    let (number, member_end) = ber.context_tag()?;
    match member_index(number, names) {
        Some(member) if member >= next => Ok(Some((member, member_end))),
        _ => Err(ber.unexpected_member(number, names)),
    }
}

/// The index of the member or variant whose tag is `number` among `names` (automatic tags: 0,
/// 1, ...), `None` for a tag of none.
fn member_index(number: u32, names: &[&str]) -> Option<u8> {
    (number < names.len() as u32).then_some(number as u8)
}

// NCBI reference (598d8ae6): c++/src/serial/choice.cpp:371-411
// `BeginChoiceVariant` of a CHOICE with automatic tags (objistrasnb.cpp:1323-1341) reads the
// tag and length as `BeginClassMember` does and finds the variant by its tag alone.
/// `ReadChoice` of the CHOICE whose variants are `names`: the variant read by `variant` (with
/// its tag) in the variant frame.
fn read_choice<'a, T>(
    ber: &mut Reader<'a>,
    names: &'static [&'static str],
    variant: impl FnOnce(&mut Reader<'a>, u8) -> Result<T, BerError>,
) -> Result<T, BerError> {
    let (number, end) =
        begin_choice_variant(ber, names).map_err(|e| e.left(Frame::Variant(None)))?;
    let frame = Frame::Variant(Some(names[usize::from(number)]));
    let value = variant(ber, number).map_err(|e| e.left(frame))?;
    ber.end_of_content(end).map_err(|e| e.left(frame))?;
    Ok(value)
}

/// `BeginChoiceVariant` of a CHOICE whose variants are `names`.
fn begin_choice_variant(ber: &mut Reader<'_>, names: &[&str]) -> Result<(u8, End), BerError> {
    let (number, end) = ber.context_tag()?;
    match member_index(number, names) {
        Some(variant) => Ok((variant, end)),
        None => Err(ber.unexpected_member(number, names)),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2208-2218
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:959-969
// `Blast-def-line-set` is an implicit class (`BEGIN_NAMED_BASE_IMPLICIT_CLASS_INFO`,
// Blast_def_line_set_.cpp:56-63) whose one member is the container. `EndOfRead` does not look at
// what follows the object: NCBI reads a record with bytes after the set (probe c25).
/// Decodes the `Blast-def-line-set` of one sequence with every member: the titles, the
/// Seq-ids of every alternative, the taxids, memberships, links and other-info, as NCBI's
/// object stream reads it, failing where it fails with its exception
/// ([`BerError::exception`]). An empty record gives no definition line, as in NCBI. The
/// `general` ids `BL_ORD_ID` keep the volume-relative tag that the record stores
/// (`SeqDb::defline_set` adds the volume start).
pub fn decode_full_set(record: &[u8]) -> Result<DeflineSet, BerError> {
    let mut set = DeflineSet::default();
    decode_full_set_in(record, &mut set)?;
    Ok(set)
}

/// [`decode_full_set`] written over `set`, whose buffers (its lines, their titles, Seq-id lists
/// and Seq-id strings) are kept: the decode of a reader that reads the header of one sequence after another
/// (blastdbcmd, `SeqDb::header_set_in`) allocates them once, not for each sequence. Where the
/// decode fails, `set` holds part of the record.
pub fn decode_full_set_in(record: &[u8], set: &mut DeflineSet) -> Result<(), BerError> {
    set.utf8 = false;
    if record.is_empty() {
        set.deflines.clear();
        return Ok(());
    }
    let mut ber = Reader::new(record);
    read_set(&mut ber, &mut set.deflines)?;
    set.utf8 = ber.read_utf8();
    Ok(())
}

/// `CObjectIStream::Read` of the root `Blast-def-line-set` (cited at `decode_full_set`), written
/// over `deflines`.
fn read_set(ber: &mut Reader<'_>, deflines: &mut Vec<BlastDefline>) -> Result<(), BerError> {
    let mut count = 0;
    read_container(ber, Frame::SetArray, |ber| {
        if count == deflines.len() {
            deflines.push(BlastDefline::default());
        }
        decode_defline(ber, &mut deflines[count])?;
        count += 1;
        Ok(())
    })
    .map_err(|e| e.left(Frame::Root))?;
    deflines.truncate(count);
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:298-309
/// The member that NCBI's warning names where its read of `record` reaches a UTF8String where a
/// VisibleString is expected (`GetStackTraceASN` there: `Blast-def-line-set.E.title`,
/// `Blast-def-line-set.E.seqid.E.other.accession`): the first in the order of the read, `None`
/// where the read meets none or fails before one (NCBI then posts nothing). The record is read
/// again, as `decode_full_set` reads it, up to that string: asked only of a record whose read
/// noted one (`DeflineSet::utf8`, `seqdb::header::Defline::utf8`) or failed.
#[cold]
#[inline(never)]
pub fn utf8_string_path(record: &[u8]) -> Option<String> {
    if record.is_empty() {
        return None;
    }
    let mut ber = Reader::new(record);
    ber.find_utf8 = true;
    match read_set(&mut ber, &mut Vec::new()) {
        Err(BerError(failure)) => match *failure {
            Failure::Utf8 { path } => Some(path),
            _ => None,
        },
        Ok(_) => None,
    }
}

/// The member path of NCBI's post of a byte that its read of a title replaces
/// (`title_replaced_chars`: `GetStackTrace` of the read, `Blast-def-line-set.E.title`).
const TITLE_TRACE: &str = "Blast-def-line-set.E.title";

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1865-1886
/// What NCBI's read of `record` (a `Blast-def-line-set`, as `decode_full_set` reads it) posts for
/// the bytes outside ' ' ..= '~' of its titles: for each such byte, in the order of the read, the
/// byte, the position of the stream after the whole string (the same for every byte of one
/// string) and the string as it is when NCBI posts it, the bytes before replaced by '#' and this
/// one and those after still stored (`FixVisibleChars`, whose one port is
/// `report::serial::asn_ber_reader::fix_visible_chars`; `report::serial::encoder::ReplacedChar`,
/// whose `post` is the text). The record is read again up to a failure (NCBI posts what it read
/// before it). Asked at each read of a record that the model of NCBI's cache of header records
/// does not answer, once a title with such a byte was read (`SeqDb::title_posts_at`), and of a
/// read that failed (`blastinput::app::ErrorAfterBatchReports::title_posts`): a record whose
/// titles hold none returns nothing.
#[cold]
#[inline(never)]
pub fn title_replaced_chars(record: &[u8]) -> Vec<crate::report::serial::encoder::ReplacedChar> {
    use crate::report::serial::asn_input::StringPost;
    use crate::report::serial::encoder::{ReplacedChar, StreamPosition};
    let mut ber = Reader::new(record);
    ber.replaced_titles = Some(Vec::new());
    let _ = read_set(&mut ber);
    let mut posts = Vec::new();
    for range in ber.replaced_titles.take().unwrap_or_default() {
        let end = StreamPosition::Byte(range.end as u64);
        crate::report::serial::asn_ber_reader::fix_visible_chars(
            &mut record[range].to_vec(),
            end,
            &mut |post| {
                if let StringPost::BadChar {
                    position,
                    byte,
                    string,
                } = post
                {
                    posts.push(ReplacedChar {
                        trace: TITLE_TRACE.to_string(),
                        byte,
                        position,
                        string: string.into(),
                    });
                }
            },
        );
    }
    posts
}

/// NCBI's warning of a UTF8String read where a VisibleString is expected (cited at
/// `utf8_string_path`), whose member is `path`, as the diagnostics of the program with the post
/// prefix `prefix` (`SetDiagPostPrefix`: `[blastp]`, `[blastdbcmd]`) print it: `Warning:
/// (803.10) [<prefix>] CObjectIStreamAsnBinary: UTF8String data for VisibleString member <path>,
/// ASN.1 specification may need an update` and a line feed. NCBI posts it once per process. The
/// text is the serial reader's (`asn_input::utf8_tag_message`, also of makeblastdb's mask files).
pub fn utf8_string_warning(prefix: Option<&str>, path: &str) -> Vec<u8> {
    crate::ncbi_exception::post_line_in(
        crate::ncbi_exception::Severity::Warning,
        Some("803.10"),
        prefix,
        crate::report::serial::asn_input::utf8_tag_message(path).as_bytes(),
        false,
    )
}

/// Blast_def_line_.cpp:91-103: every member has a set flag; `seqid` is the one mandatory.
const BLAST_DEF_LINE: Class = Class {
    name: "Blast-def-line",
    members: &DEFLINE_MEMBERS,
    mandatory: 1 << DEFLINE_SEQID,
    set_flag: 0b11_1111,
};

/// Decodes one `Blast-def-line` over `defline`, whose title, Seq-id list and Seq-id strings keep
/// their buffers (`IdStrings`).
fn decode_defline(ber: &mut Reader<'_>, defline: &mut BlastDefline) -> Result<(), BerError> {
    let mut title = defline.title.take().unwrap_or_default();
    let mut has_title = false;
    defline.replaced_title_bytes = 0;
    defline.replaced_id_bytes = 0;
    defline.taxid = None;
    defline.memberships = None;
    defline.links = None;
    defline.other_info = None;
    read_class(ber, &BLAST_DEF_LINE, |ber, member| {
        match member {
            DEFLINE_TITLE => {
                let replaced = ber.fixed_string_in(&mut title)?;
                if replaced > 0 {
                    ber.note_replaced_title(title.len());
                }
                defline.replaced_title_bytes = replaced;
                has_title = true;
            }
            DEFLINE_SEQID => {
                let ids = &mut defline.ids;
                let mut strings = IdStrings::default();
                let mut count = 0;
                read_container(ber, Frame::Array, |ber| {
                    match ids.get_mut(count) {
                        Some(place) => {
                            strings.keep(std::mem::replace(place, ParsedSeqId::Gi(0)));
                            *place = decode_seq_id(ber, &mut strings)?;
                        }
                        None => ids.push(decode_seq_id(ber, &mut strings)?),
                    }
                    count += 1;
                    Ok(())
                })?;
                ids.truncate(count);
                defline.replaced_id_bytes = strings.replaced;
            }
            DEFLINE_TAXID => defline.taxid = Some(ber.int()?),
            DEFLINE_MEMBERSHIPS => defline.memberships = Some(read_int_list(ber)?),
            DEFLINE_LINKS => defline.links = Some(read_int_list(ber)?),
            _ => defline.other_info = Some(read_int_list(ber)?),
        }
        Ok(())
    })?;
    defline.title = has_title.then_some(title);
    Ok(())
}

/// A `SEQUENCE OF INTEGER` member of `Blast-def-line` (NCBI's `list<int>`).
fn read_int_list(ber: &mut Reader<'_>) -> Result<Vec<i32>, BerError> {
    let mut values = Vec::new();
    read_container(ber, Frame::Array, |ber| {
        values.push(ber.int()?);
        Ok(())
    })?;
    Ok(values)
}

/// A `SEQUENCE OF INTEGER` member of `Blast-def-line` as the ordinal-id decoder reads it
/// (`seqdb::header`).
pub(crate) fn int_list(ber: &mut Reader<'_>, what: &str) -> Result<Vec<i32>, BerError> {
    let end = ber.open_sequence(what)?;
    let mut values = Vec::new();
    while !ber.at_end(end)? {
        values.push(ber.int()?);
    }
    Ok(values)
}

/// The variants of `Seq-id` (Seq_id_.cpp:570-595), by tag.
const SEQ_ID: [&str; 20] = [
    "local",
    "gibbsq",
    "gibbmt",
    "giim",
    "genbank",
    "embl",
    "pir",
    "swissprot",
    "patent",
    "other",
    "general",
    "gi",
    "ddbj",
    "prf",
    "pdb",
    "tpg",
    "tpe",
    "tpd",
    "gpipe",
    "named-annot-track",
];

/// Decodes one `Seq-id`, its strings read by `strings` (which counts the bytes NCBI's reader
/// replaces, `fix_visible_chars`).
///
/// Not inlined: inlined into the member closure of `decode_defline` (the release build does
/// since the INTEGER reads became small), it keeps that closure, which runs for every member,
/// out of `decode_full_set`, with a 1.3 KiB frame: about 100 instructions more per record
/// (callgrind of `blastdbcmd -entry all`, RG-T5c fix 2).
#[inline(never)]
pub(crate) fn decode_seq_id(
    ber: &mut Reader<'_>,
    strings: &mut IdStrings,
) -> Result<ParsedSeqId, BerError> {
    read_choice(ber, &SEQ_ID, |ber, alternative| {
        Ok(match alternative {
            0 => ParsedSeqId::Local(decode_object_id(ber, strings)?),
            1 => ParsedSeqId::Gibbsq(ber.int()?),
            2 => ParsedSeqId::Gibbmt(ber.int()?),
            3 => ParsedSeqId::Giim(decode_giimport_id(ber, strings)?),
            4 => ParsedSeqId::Genbank(decode_textseq_id(ber, strings)?),
            5 => ParsedSeqId::Embl(decode_textseq_id(ber, strings)?),
            6 => ParsedSeqId::Pir(decode_textseq_id(ber, strings)?),
            7 => ParsedSeqId::Swissprot(decode_textseq_id(ber, strings)?),
            8 => ParsedSeqId::Patent(decode_patent_seq_id(ber, strings)?),
            9 => ParsedSeqId::Other(decode_textseq_id(ber, strings)?),
            10 => ParsedSeqId::General(decode_dbtag(ber, strings)?),
            11 => ParsedSeqId::Gi(ber.int8()?),
            12 => ParsedSeqId::Ddbj(decode_textseq_id(ber, strings)?),
            13 => ParsedSeqId::Prf(decode_textseq_id(ber, strings)?),
            14 => ParsedSeqId::Pdb(decode_pdb_seq_id(ber, strings)?),
            15 => ParsedSeqId::Tpg(decode_textseq_id(ber, strings)?),
            16 => ParsedSeqId::Tpe(decode_textseq_id(ber, strings)?),
            17 => ParsedSeqId::Tpd(decode_textseq_id(ber, strings)?),
            18 => ParsedSeqId::Gpipe(decode_textseq_id(ber, strings)?),
            _ => ParsedSeqId::NamedAnnotTrack(decode_textseq_id(ber, strings)?),
        })
    })
}

/// The strings that the decode of the Seq-ids of a definition line reads: the number of their
/// bytes that NCBI's reader replaces (`fix_visible_chars`), and the buffers of the Seq-id that
/// held the place before (`decode_defline`, which writes over the lines of the record read
/// before), which the strings are read into in the order of their members: a line of the shape
/// of the one before allocates no string.
#[derive(Default)]
pub(crate) struct IdStrings {
    /// The number of the strings' bytes that read as '#'.
    pub replaced: usize,
    spare: [Vec<u8>; ID_STRINGS],
    /// `spare[next..count]` are the buffers not yet read into.
    next: usize,
    count: usize,
}

/// The most strings of one `Seq-id` (a `pdb` id's `mol`, release date and `chain-id`; a patent's
/// country, number and document type; a `Textseq-id`'s name, accession and release).
const ID_STRINGS: usize = 3;

impl IdStrings {
    /// A VisibleString member of a Seq-id (`Reader::fixed_string_in`), read into a spare buffer.
    #[inline]
    fn read(&mut self, ber: &mut Reader<'_>) -> Result<Vec<u8>, BerError> {
        let mut value = if self.next < self.count {
            self.next += 1;
            std::mem::take(&mut self.spare[self.next - 1])
        } else {
            Vec::new()
        };
        self.replaced += ber.fixed_string_in(&mut value)?;
        Ok(value)
    }

    /// Makes the string buffers of `id`, the Seq-id at the place the next one is read into, the
    /// spare buffers, in the order of its members. Kept out of line so its 20-variant match does
    /// not bloat `read_set`'s stack frame or spill registers on fresh sets where it never runs.
    #[inline(never)]
    fn keep(&mut self, id: ParsedSeqId) {
        (self.next, self.count) = (0, 0);
        let mut put = |buffer: Vec<u8>| {
            if self.count < ID_STRINGS {
                self.spare[self.count] = buffer;
                self.count += 1;
            }
        };
        let text = |id: TextseqId, put: &mut dyn FnMut(Vec<u8>)| {
            for buffer in [id.name, id.accession, id.release].into_iter().flatten() {
                put(buffer);
            }
        };
        match id {
            ParsedSeqId::Local(ObjectId::Str(value)) => put(value),
            ParsedSeqId::Local(ObjectId::Id(_))
            | ParsedSeqId::Gibbsq(_)
            | ParsedSeqId::Gibbmt(_)
            | ParsedSeqId::Gi(_) => {}
            ParsedSeqId::Giim(GiimportId { db, release, .. }) => {
                for buffer in [db, release].into_iter().flatten() {
                    put(buffer);
                }
            }
            ParsedSeqId::Genbank(id)
            | ParsedSeqId::Embl(id)
            | ParsedSeqId::Pir(id)
            | ParsedSeqId::Swissprot(id)
            | ParsedSeqId::Other(id)
            | ParsedSeqId::Ddbj(id)
            | ParsedSeqId::Prf(id)
            | ParsedSeqId::Tpg(id)
            | ParsedSeqId::Tpe(id)
            | ParsedSeqId::Tpd(id)
            | ParsedSeqId::Gpipe(id)
            | ParsedSeqId::NamedAnnotTrack(id) => text(id, &mut put),
            ParsedSeqId::Patent(PatentSeqId { cit, .. }) => {
                put(cit.country);
                match cit.id {
                    PatentNumber::Number(number) | PatentNumber::AppNumber(number) => put(number),
                }
                if let Some(doc_type) = cit.doc_type {
                    put(doc_type);
                }
            }
            ParsedSeqId::General(Dbtag { db, tag }) => {
                put(db);
                if let ObjectId::Str(value) = tag {
                    put(value);
                }
            }
            ParsedSeqId::Pdb(PdbSeqId {
                mol, rel, chain_id, ..
            }) => {
                put(mol);
                match rel {
                    Some(Date::Str(value)) => put(value),
                    Some(Date::Std(date)) => {
                        if let Some(season) = date.season {
                            put(season);
                        }
                    }
                    None => {}
                }
                if let Some(chain_id) = chain_id {
                    put(chain_id);
                }
            }
        }
    }
}

/// Object_id_.cpp:108-115.
const OBJECT_ID: [&str; 2] = ["id", "str"];

fn decode_object_id(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<ObjectId, BerError> {
    read_choice(ber, &OBJECT_ID, |ber, alternative| {
        Ok(match alternative {
            0 => ObjectId::Id(ber.int()?),
            _ => ObjectId::Str(strings.read(ber)?),
        })
    })
}

/// Dbtag_.cpp:77-85: `db` with a set flag, `tag` a reference member without one; both
/// mandatory.
const DBTAG: Class = Class {
    name: "Dbtag",
    members: &["db", "tag"],
    mandatory: 0b11,
    set_flag: 0b01,
};

fn decode_dbtag(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<Dbtag, BerError> {
    let (mut db, mut tag) = (Vec::new(), ObjectId::Id(0));
    read_class(ber, &DBTAG, |ber, member| {
        match member {
            0 => db = strings.read(ber)?,
            _ => tag = decode_object_id(ber, strings)?,
        }
        Ok(())
    })?;
    Ok(Dbtag { db, tag })
}

/// Textseq_id_.cpp:76-86: every member optional, with a set flag.
const TEXTSEQ_ID: Class = Class {
    name: "Textseq-id",
    members: &["name", "accession", "release", "version"],
    mandatory: 0,
    set_flag: 0b1111,
};

fn decode_textseq_id(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<TextseqId, BerError> {
    let mut id = TextseqId::default();
    read_class(ber, &TEXTSEQ_ID, |ber, member| {
        match member {
            0 => id.name = Some(strings.read(ber)?),
            1 => id.accession = Some(strings.read(ber)?),
            2 => id.release = Some(strings.read(ber)?),
            _ => id.version = Some(ber.int()?),
        }
        Ok(())
    })?;
    Ok(id)
}

/// Giimport_id_.cpp:69-78: `id` mandatory; every member with a set flag.
const GIIMPORT_ID: Class = Class {
    name: "Giimport-id",
    members: &["id", "db", "release"],
    mandatory: 0b001,
    set_flag: 0b111,
};

fn decode_giimport_id(
    ber: &mut Reader<'_>,
    strings: &mut IdStrings,
) -> Result<GiimportId, BerError> {
    let (mut id, mut db, mut release) = (0, None, None);
    read_class(ber, &GIIMPORT_ID, |ber, member| {
        match member {
            0 => id = ber.int()?,
            1 => db = Some(strings.read(ber)?),
            _ => release = Some(strings.read(ber)?),
        }
        Ok(())
    })?;
    Ok(GiimportId { id, db, release })
}

/// Patent_seq_id_.cpp:71-79: both mandatory; `cit` a reference member without a set flag.
const PATENT_SEQ_ID: Class = Class {
    name: "Patent-seq-id",
    members: &["seqid", "cit"],
    mandatory: 0b11,
    set_flag: 0b01,
};

fn decode_patent_seq_id(
    ber: &mut Reader<'_>,
    strings: &mut IdStrings,
) -> Result<PatentSeqId, BerError> {
    let (mut seqid, mut cit) = (0, None);
    read_class(ber, &PATENT_SEQ_ID, |ber, member| {
        match member {
            0 => seqid = ber.int()?,
            _ => cit = Some(decode_id_pat(ber, strings)?),
        }
        Ok(())
    })?;
    let cit = cit.expect("the mandatory member cit was read");
    Ok(PatentSeqId { seqid, cit })
}

/// Id_pat_.cpp:170-179: `country` and `id` mandatory; `id` a reference member (the unnamed
/// CHOICE `C_Id`, :113-121) without a set flag.
const ID_PAT: Class = Class {
    name: "Id-pat",
    members: &["country", "id", "doc-type"],
    mandatory: 0b011,
    set_flag: 0b101,
};

/// The variants of `Id-pat.id` (Id_pat_.cpp:113-121).
const ID_PAT_ID: [&str; 2] = ["number", "app-number"];

fn decode_id_pat(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<IdPat, BerError> {
    let (mut country, mut id, mut doc_type) = (Vec::new(), None, None);
    read_class(ber, &ID_PAT, |ber, member| {
        match member {
            0 => country = strings.read(ber)?,
            1 => {
                id = Some(read_choice(ber, &ID_PAT_ID, |ber, alternative| {
                    let number = strings.read(ber)?;
                    Ok(match alternative {
                        0 => PatentNumber::Number(number),
                        _ => PatentNumber::AppNumber(number),
                    })
                })?)
            }
            _ => doc_type = Some(strings.read(ber)?),
        }
        Ok(())
    })?;
    let id = id.expect("the mandatory member id was read");
    Ok(IdPat {
        country,
        id,
        doc_type,
    })
}

/// PDB_seq_id_.cpp:88-98: `mol` mandatory (the class `PDB-mol-id`, an alias of VisibleString);
/// `chain` with the default 32; `rel` a reference member without a set flag.
const PDB_SEQ_ID: Class = Class {
    name: "PDB-seq-id",
    members: &["mol", "chain", "rel", "chain-id"],
    mandatory: 0b0001,
    set_flag: 0b1011,
};

fn decode_pdb_seq_id(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<PdbSeqId, BerError> {
    let (mut mol, mut chain, mut rel, mut chain_id) = (Vec::new(), None, None, None);
    read_class(ber, &PDB_SEQ_ID, |ber, member| {
        match member {
            0 => mol = strings.read(ber)?,
            1 => chain = Some(ber.int()?),
            2 => rel = Some(decode_date(ber, strings)?),
            _ => chain_id = Some(strings.read(ber)?),
        }
        Ok(())
    })?;
    Ok(PdbSeqId {
        mol,
        chain,
        rel,
        chain_id,
    })
}

/// The variants of `Date` (Date_.cpp:134-141).
const DATE: [&str; 2] = ["str", "std"];

/// Date_std_.cpp:67-80: `year` mandatory; every member with a set flag.
const DATE_STD: Class = Class {
    name: "Date-std",
    members: &["year", "month", "day", "season", "hour", "minute", "second"],
    mandatory: 0b000_0001,
    set_flag: 0b111_1111,
};

fn decode_date(ber: &mut Reader<'_>, strings: &mut IdStrings) -> Result<Date, BerError> {
    read_choice(ber, &DATE, |ber, alternative| {
        Ok(match alternative {
            0 => Date::Str(strings.read(ber)?),
            _ => {
                let mut date = DateStd {
                    year: 0,
                    month: None,
                    day: None,
                    season: None,
                    hour: None,
                    minute: None,
                    second: None,
                };
                read_class(ber, &DATE_STD, |ber, member| {
                    match member {
                        0 => date.year = ber.int()?,
                        1 => date.month = Some(ber.int()?),
                        2 => date.day = Some(ber.int()?),
                        3 => date.season = Some(strings.read(ber)?),
                        4 => date.hour = Some(ber.int()?),
                        5 => date.minute = Some(ber.int()?),
                        _ => date.second = Some(ber.int()?),
                    }
                    Ok(())
                })?;
                Date::Std(date)
            }
        })
    })
}

#[cfg(test)]
#[path = "ber_tests.rs"]
mod tests;
