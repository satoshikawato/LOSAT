//! The steps of NCBI's binary ASN.1 object stream (`CObjectIStreamAsnBinary`) that the walk of
//! the type table (`read`) calls: BER as NCBI's reader takes it, over the bytes of a source
//! (`AsnSource`, the text reader's: a stream is kept from the last element on, so one object at a
//! time is in memory). The tags and lengths are the writer's (`report::serial::asn_ber`).
//!
//! NCBI's build reads with `VIRTUAL_MID_LEVEL_IO` (the binary stream's own `ReadClassSequential`,
//! `ReadChoiceSimple`, `ReadContainer` and skip functions), `USE_DEF_LEN` (a definite length is
//! taken where an indefinite one is expected, `ExpectIndefiniteLength`) and without
//! `CHECK_INSTREAM_STATE` and `CHECK_INSTREAM_LIMITS` (release build): the length of a primitive
//! value is not checked against the value that holds it. NCBI's types are automatically tagged
//! (`CTypeInfo::m_TagType`, typeinfo.cpp:108): each member and variant is an explicit context
//! tag `[n]` around its value, `n` its position (`Member::tag`).
//!
//! The end of the data is NCBI's `CEofException` (`ReadError::Eof`), which the walk turns into a
//! `CSerialException` when it leaves a frame with 2 or more frames left (`HandleEOF`): an
//! archive whose last byte is missing ends in the root class's `EndClass` and is no object at all
//! (rc 0, NCBI 2.17.0 probes P3.104 and P7.001-002); 2 to 16 bytes missing end inside the member
//! frame, an invalid input (P7.003-022).

use super::read::{read_delayed, Cx, Front, Post, ReadError, R};
use super::{EnumValue, Member, Storage, StringKind, Value};
use crate::ncbi_str::{string_to_double_posix, NCBI_NAN};
use crate::report::serial::asn_ber::{
    END_OF_CONTENTS, INDEFINITE_LENGTH, REAL_DECIMAL, REAL_ENCODING_MASK, REAL_MINUS_INFINITY,
    REAL_MINUS_ZERO, REAL_NOT_A_NUMBER, REAL_PLUS_INFINITY, TAG_APPLICATION_INTEGER,
    TAG_BIT_STRING, TAG_BOOLEAN, TAG_CLASS_MASK, TAG_CONSTRUCTED, TAG_CONTEXT, TAG_ENUMERATED,
    TAG_INTEGER, TAG_LONG, TAG_NULL, TAG_OCTET_STRING, TAG_REAL, TAG_SEQUENCE, TAG_STRING_STORE,
    TAG_UTF8_STRING, TAG_VALUE_MASK, TAG_VISIBLE_STRING,
};
use crate::report::serial::asn_text_reader::{AsnSource, SliceSource};
use crate::report::serial::encoder::StreamPosition;
use crate::report::serial::object_frame::FrameKind;

/// `kMaxDoubleLength`: the longest REAL NCBI's reader takes.
const MAX_DOUBLE_LENGTH: usize = 256;

/// The binary front end: NCBI's reader of one stream of BER objects.
pub(crate) struct Ber<S: AsnSource> {
    pub(crate) src: S,
    /// `GetStreamPosAsInt8()`: the position of the next byte.
    pub(crate) pos: usize,
    /// `m_CurrentTagLength`: the bytes of the identifier peeked (`PeekTag`), which `FlushTag`
    /// takes.
    tag_len: usize,
    /// `m_CurrentDataLimit`: where the value of a definite length being read ends (0: an
    /// indefinite length, ended by `00 00`). NCBI's `Int8`, its sum with a length wrapped as
    /// GCC converts a `size_t` to it.
    limit: i64,
    /// `m_DataLimits`: the limits of the values around it.
    limits: Vec<i64>,
    /// The warning post with subcode 10 (`ERR_POST_X_ONCE`) of `ExpectStringTag` was made.
    utf8_warned: bool,
    /// The start of a delay buffer being skipped: the source keeps its bytes.
    held: Option<usize>,
}

impl<S: AsnSource> Ber<S> {
    pub(crate) fn new(src: S) -> Self {
        Ber {
            src,
            pos: 0,
            tag_len: 0,
            limit: 0,
            limits: Vec::new(),
            utf8_warned: false,
            held: None,
        }
    }

    /// The source keeps no byte before the position (nor before a delay buffer being skipped).
    fn release(&mut self) {
        let keep = self.held.map_or(self.pos, |held| held.min(self.pos));
        self.src.release_before(keep);
    }

    /// `PeekTagByte` / `PeekChar` at an offset from the position: the end of the data is a
    /// `CEofException`.
    fn peek(&mut self, index: usize) -> R<u8> {
        self.src.byte_at(self.pos + index).ok_or(ReadError::Eof)
    }

    /// `GetChar` (`ReadByte`).
    fn get(&mut self) -> R<u8> {
        let byte = self.peek(0)?;
        self.pos += 1;
        Ok(byte)
    }

    /// `GetChars` / `ReadBytes` of `count` bytes: the next `count` bytes, taken; the end of the
    /// data before them is a `CEofException`.
    fn take(&mut self, count: usize) -> R<std::ops::Range<usize>> {
        let start = self.pos;
        let end = start.checked_add(count).ok_or(ReadError::Eof)?;
        if count > 0 {
            self.src.byte_at(end - 1).ok_or(ReadError::Eof)?;
        }
        self.pos = end;
        Ok(start..end)
    }

    /// `ExpectSysTagByte`: the identifier must be `byte`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:114-123
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:374-379
    fn expect_sys_tag(&mut self, byte: u8) -> R<()> {
        let got = self.peek(0)?;
        if got != byte {
            return Err(unexpected_sys_tag(got, byte));
        }
        self.tag_len = 1;
        Ok(())
    }

    /// `ExpectIntegerTag`: an INTEGER, universal or application (a big integer). NCBI's code also
    /// flags the application tag as the special case `eReadAsBigInt`, which changes nothing in
    /// the value.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:145-164
    fn expect_integer_tag(&mut self) -> R<()> {
        let got = self.peek(0)?;
        if got != TAG_INTEGER && got != TAG_APPLICATION_INTEGER {
            return Err(unexpected_sys_tag(got, TAG_INTEGER));
        }
        self.tag_len = 1;
        Ok(())
    }

    /// `PeekTag` with a class and constructed check (`class`: the class and constructed bits):
    /// the tag number of the identifier that starts with `first`, its length in `tag_len`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:166-200
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:174-180
    fn peek_tag(&mut self, first: u8, class: u8) -> R<i32> {
        if first & (TAG_CLASS_MASK | TAG_CONSTRUCTED) != class {
            return Err(ReadError::Serial(format!(
                "unexpected tag: {}, should be {}",
                tag_to_string(first),
                tag_to_string(class)
            )));
        }
        if first & TAG_VALUE_MASK != TAG_LONG {
            self.tag_len = 1;
            return Ok(i32::from(first & TAG_VALUE_MASK));
        }
        self.long_tag()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:151-172
    fn long_tag(&mut self) -> R<i32> {
        let mut tag = 0i32;
        let mut i = 1;
        loop {
            if tag >= 1 << 24 {
                return Err(ReadError::Serial(format!("tag number is too big: {tag}")));
            }
            let byte = self.peek(i)?;
            i += 1;
            tag = (tag << 7) | i32::from(byte & 0x7f);
            if byte & 0x80 == 0 {
                break;
            }
        }
        self.tag_len = i;
        Ok(tag)
    }

    /// `ExpectTag` of a universal class and the constructed form: the identifier of a SEQUENCE
    /// (16) or SET (17).
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:209-219
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:205-217
    fn expect_constructed(&mut self, number: i32) -> R<()> {
        let first = self.peek(0)?;
        let tag = self.peek_tag(first, TAG_CONSTRUCTED)?;
        if tag != number {
            return Err(ReadError::Serial(format!(
                "Unexpected tag: {tag}, expected: {number}"
            )));
        }
        Ok(())
    }

    /// `FlushTag`: the identifier taken, and the first byte of the length.
    fn flush_tag(&mut self) -> R<u8> {
        self.pos += self.tag_len;
        self.tag_len = 0;
        self.get()
    }

    /// `ReadShortLength`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:337-355
    // (`UnexpectedLongLength` is "ShortLength expected", objistrasnb.cpp:200-203; `StartTagData`
    // checks nothing without `CHECK_INSTREAM_LIMITS`.)
    fn short_length(&mut self) -> R<usize> {
        let byte = self.flush_tag()?;
        if byte >= 0x80 {
            return Err(ReadError::Serial("ShortLength expected".into()));
        }
        Ok(usize::from(byte))
    }

    /// `ReadLength`.
    fn length(&mut self) -> R<u64> {
        let byte = self.flush_tag()?;
        if byte >= 0x80 {
            return self.long_length(byte);
        }
        Ok(u64::from(byte))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:394-411
    fn long_length(&mut self, byte: u8) -> R<u64> {
        let mut count = usize::from(byte - 0x80);
        if count == 0 {
            return Err(ReadError::Serial("unexpected indefinite length".into()));
        }
        if count > 8 {
            return Err(ReadError::Serial("length overflow".into()));
        }
        let mut length = u64::from(self.get()?);
        if length == 0 {
            return Err(ReadError::Serial("illegal length start".into()));
        }
        while {
            count -= 1;
            count > 0
        } {
            length = (length << 8) | u64::from(self.get()?);
        }
        Ok(length)
    }

    /// `ReadLength` of a value whose bytes are then read: a length past the memory of the
    /// machine is past the end of the data (see `take`).
    fn data_length(&mut self) -> R<usize> {
        let length = self.length()?;
        usize::try_from(length).map_err(|_| ReadError::Eof)
    }

    /// `ExpectIndefiniteLength` (`USE_DEF_LEN`): the length of a constructed value, whose end
    /// becomes the limit (`HaveMoreElements`, `ExpectEndOfContent`).
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:276-312
    // The order of the two operands of the long form is not set by C++ (unsequenced); NCBI
    // 2.17.0's build takes the position before the length's bytes: impl D's probe (2026-11-20,
    // `F3/tmp/d/probe/probe_d.py`) gave a root `30 82 HH LL` whose length counts only its content
    // "Invalid input format for BLAST Archive." (rc 1, `def_len`) and formatted the one whose
    // length counts the two length bytes too (rc 0, `def_len_plus2`).
    fn expect_indefinite_length(&mut self) -> R<()> {
        let len = self.flush_tag()?;
        self.limits.push(self.limit);
        self.limit = if len == INDEFINITE_LENGTH {
            0
        } else if len > INDEFINITE_LENGTH {
            let at = self.pos as i64;
            let length = self.long_length(len)?;
            at.wrapping_add(length as i64)
        } else {
            self.pos as i64 + i64::from(len)
        };
        Ok(())
    }

    /// `HaveMoreElements`: before the limit, or before `00`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:242-251
    fn have_more_elements(&mut self) -> R<bool> {
        if self.limit != 0 {
            return Ok(self.limit > self.pos as i64);
        }
        Ok(self.peek(0)? != END_OF_CONTENTS)
    }

    /// `ExpectEndOfContent`: at the limit, or `00 00`; the limit of the value around it back.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:373-392
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:72-88
    // (`UnexpectedContinuation` is "EndOfContentsByte expected", objistrasnb.cpp:195-198.)
    fn expect_end_of_content(&mut self) -> R<()> {
        let continuation = || ReadError::Serial("EndOfContentsByte expected".into());
        if self.limit != 0 {
            if self.limit != self.pos as i64 {
                return Err(continuation());
            }
        } else {
            for _ in 0..2 {
                if self.peek(0)? != END_OF_CONTENTS {
                    return Err(continuation());
                }
                self.pos += 1;
            }
        }
        self.limit = self
            .limits
            .pop()
            .expect("each end of content closes a length expected before it");
        Ok(())
    }

    /// `SkipTagData`: the value's length, its bytes skipped (as `GetChars` takes them: the end of
    /// the data before them is a `CEofException`).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1785-1789
    fn skip_tag_data(&mut self) -> R<()> {
        let length = self.data_length()?;
        self.take(length)?;
        Ok(())
    }

    /// `ReadStdSigned` of a C++ integer of `size` bytes: the value's bytes, two's complement,
    /// more than `size` of them only as the sign's extension.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:475-512
    // (`ReadInt8` reads into a `Uint8`: the same bits.)
    fn std_signed(&mut self, size: usize) -> R<i64> {
        let overflow = || ReadError::Serial("overflow error".into());
        let mut length = self.short_length()?;
        if length == 0 {
            return Err(ReadError::Serial("zero length of number".into()));
        }
        let mut n;
        if length > size {
            length -= 1;
            let c = self.get()? as i8;
            if c != 0 && c != -1 {
                return Err(overflow());
            }
            while length > size {
                length -= 1;
                if self.get()? as i8 != c {
                    return Err(overflow());
                }
            }
            length -= 1;
            n = i64::from(self.get()? as i8);
            if (n ^ i64::from(c)) & 0x80 != 0 {
                return Err(overflow());
            }
        } else {
            length -= 1;
            n = i64::from(self.get()? as i8);
        }
        while length > 0 {
            length -= 1;
            n = n.wrapping_shl(8) | i64::from(self.get()?);
        }
        Ok(n)
    }

    /// `ReadStdUnsigned` of a `Uint4`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:514-549
    fn std_unsigned4(&mut self) -> R<u32> {
        let overflow = || ReadError::Serial("overflow error".into());
        let mut length = self.short_length()?;
        if length == 0 {
            return Err(ReadError::Serial("zero length of number".into()));
        }
        let mut n: u32;
        if length > 4 {
            while length > 4 {
                length -= 1;
                if self.get()? != 0 {
                    return Err(overflow());
                }
            }
            length -= 1;
            n = u32::from(self.get()?);
        } else if length == 4 {
            length -= 1;
            n = u32::from(self.get()?);
            if n & 0x80 != 0 {
                return Err(overflow());
            }
        } else {
            n = 0;
        }
        while length > 0 {
            length -= 1;
            n = n.wrapping_shl(8) | u32::from(self.get()?);
        }
        Ok(n)
    }

    /// `ExpectStringTag`: a VisibleString or a UTF8String takes the other's tag too
    /// (`SERIAL_READ_ANY_VISIBLESTRING_TAG` 1 and `SERIAL_READ_ANY_UTF8STRING_TAG` true, the
    /// defaults); a UTF8String's tag for a VisibleString is posted once.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:281-316
    // (`ERR_POST_X_ONCE` posts once in the process; LOSAT reads one archive or strategy stream
    // in a process, and once in it, its delay buffers included.)
    fn expect_string_tag(&mut self, utf8: bool, cx: &mut Cx<'_>) -> R<()> {
        let got = self.peek(0)?;
        if utf8 {
            if got == TAG_VISIBLE_STRING {
                return self.expect_sys_tag(got);
            }
            return self.expect_sys_tag(TAG_UTF8_STRING);
        }
        if got == TAG_UTF8_STRING {
            if !self.utf8_warned {
                self.utf8_warned = true;
                let trace = cx.trace();
                cx.post(Post::Utf8ForVisible(trace));
            }
            return self.expect_sys_tag(got);
        }
        self.expect_sys_tag(TAG_VISIBLE_STRING)
    }
}

/// `UnexpectedSysTagByte`'s message.
fn unexpected_sys_tag(got: u8, expected: u8) -> ReadError {
    ReadError::Serial(format!(
        "unexpected tag: {}, should be: {}",
        tag_to_string(got),
        tag_to_string(expected)
    ))
}

/// `UnexpectedMember`: the tag read, and every member or variant with its tag.
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1058-1067
fn unexpected_member(tag: i32, items: &[Member]) -> ReadError {
    let mut message = format!("unexpected member: [{tag}], should be one of: ");
    for item in items {
        message.push_str(&format!("{}[{}] ", item.name(), item.tag));
    }
    ReadError::Serial(message)
}

/// `TagToString`: the class, `constructed/`, the universal type's name or the number, and the
/// byte.
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:318-372
fn tag_to_string(byte: u8) -> String {
    let class = match byte & TAG_CLASS_MASK {
        0x40 => "application/",
        0x80 => "contextspecific/",
        0xc0 => "private/",
        _ => "",
    };
    let constructed = if byte & TAG_CONSTRUCTED != 0 {
        "constructed/"
    } else {
        ""
    };
    let value = byte & TAG_VALUE_MASK;
    let name = if byte & TAG_CLASS_MASK != 0 {
        value.to_string()
    } else {
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
    };
    format!("{class}{constructed}{name} ({byte})")
}

impl<S: AsnSource> Front for Ber<S> {
    // The lines are those of `END_OBJECT_FRAME()` in the functions that push the frames:
    // `CObjectIStream::Read` (objistr.cpp:959-969) and the binary stream's
    // `ReadClassSequential` (objistrasnb.cpp:1209-1238, the member frame ends in
    // `ReadClassSequentialContentsEnd` at 1231), `ReadChoiceSimple` (1380-1412),
    // `ReadContainer` (955-1000) and the skip functions `SkipClassSequential` (1269-1297),
    // `SkipChoiceSimple` (1414-1445) and `SkipContainer` (1002-1026).
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
        const BIN: &str = "src/serial/objistrasnb.cpp";
        match (kind, skip) {
            (FrameKind::Named | FrameKind::NamedType, _) => {
                ("src/serial/objistr.cpp", 968, "ncbi::CObjectIStream::Read")
            }
            (FrameKind::Class, false) => (
                BIN,
                1237,
                "ncbi::CObjectIStreamAsnBinary::ReadClassSequential",
            ),
            (FrameKind::ClassMember, false) => (
                BIN,
                1231,
                "ncbi::CObjectIStreamAsnBinary::ReadClassSequential",
            ),
            (FrameKind::Choice, false) => {
                (BIN, 1411, "ncbi::CObjectIStreamAsnBinary::ReadChoiceSimple")
            }
            (FrameKind::ChoiceVariant, false) => {
                (BIN, 1410, "ncbi::CObjectIStreamAsnBinary::ReadChoiceSimple")
            }
            (FrameKind::Array, false) => (BIN, 998, "ncbi::CObjectIStreamAsnBinary::ReadContainer"),
            (FrameKind::ArrayElement, false) => {
                (BIN, 990, "ncbi::CObjectIStreamAsnBinary::ReadContainer")
            }
            (FrameKind::Class, true) => (
                BIN,
                1296,
                "ncbi::CObjectIStreamAsnBinary::SkipClassSequential",
            ),
            (FrameKind::ClassMember, true) => (
                BIN,
                1290,
                "ncbi::CObjectIStreamAsnBinary::SkipClassSequential",
            ),
            (FrameKind::Choice, true) => {
                (BIN, 1444, "ncbi::CObjectIStreamAsnBinary::SkipChoiceSimple")
            }
            (FrameKind::ChoiceVariant, true) => {
                (BIN, 1443, "ncbi::CObjectIStreamAsnBinary::SkipChoiceSimple")
            }
            (FrameKind::Array, true) => (BIN, 1024, "ncbi::CObjectIStreamAsnBinary::SkipContainer"),
            (FrameKind::ArrayElement, true) => {
                (BIN, 1018, "ncbi::CObjectIStreamAsnBinary::SkipContainer")
            }
        }
    }

    fn start_object(&mut self) {
        self.held = None;
        self.release();
    }

    /// `CObjectIStream::ReadFileHeader`: a byte must follow; no type name.
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1121-1126
    fn file_header(&mut self, _cx: &mut Cx<'_>) -> R<Vec<u8>> {
        self.peek(0)?;
        Ok(Vec::new())
    }

    /// `BeginClass`: a SEQUENCE's identifier and length.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1029-1044
    fn begin_class(&mut self, _ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_constructed(i32::from(TAG_SEQUENCE & TAG_VALUE_MASK))?;
        self.expect_indefinite_length()
    }

    /// `BeginClassMember` of an automatically tagged class: the next member's context tag, one
    /// of the members from `pos`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1122-1174
    // (`Find` of a context-specific tag from the position `pos`, memberlist.cpp:462-473: the
    // member whose tag it is, when its position is `pos` or later; NCBI's readers skip no
    // unknown member by default.)
    fn next_member(
        &mut self,
        members: &'static [Member],
        pos: usize,
        _block: &mut bool,
        _cx: &mut Cx<'_>,
    ) -> R<Option<usize>> {
        self.release();
        if !self.have_more_elements()? {
            return Ok(None);
        }
        let first = self.peek(0)?;
        let tag = self.peek_tag(first, TAG_CONTEXT)?;
        self.expect_indefinite_length()?;
        members[pos..]
            .iter()
            .position(|member| i32::from(member.tag) == tag)
            .map(|at| Some(pos + at))
            .ok_or_else(|| unexpected_member(tag, members))
    }

    /// `EndClassMember`: the member's `00 00`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1176-1182
    fn end_member(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_end_of_content()
    }

    /// `EndClass`: the class's `00 00`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1046-1056
    fn end_class(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_end_of_content()
    }

    /// `ReadChoiceSimple` of an automatically tagged CHOICE: the variant's context tag.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1380-1412
    fn begin_variant(
        &mut self,
        _ty: u16,
        variants: &'static [Member],
        _cx: &mut Cx<'_>,
    ) -> R<usize> {
        let first = self.peek(0)?;
        let tag = self.peek_tag(first, TAG_CONTEXT)?;
        self.expect_indefinite_length()?;
        variants
            .iter()
            .position(|variant| i32::from(variant.tag) == tag)
            .ok_or_else(|| unexpected_member(tag, variants))
    }

    fn end_variant(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_end_of_content()
    }

    /// `ReadContainer`'s start: the SEQUENCE OF's (16) or SET OF's (17) identifier and length.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:955-1000
    fn begin_list(&mut self, set: bool, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_constructed(if set { 17 } else { 16 })?;
        self.expect_indefinite_length()
    }

    /// `HaveMoreElements`.
    fn next_element(&mut self, _element: u16, _block: &mut bool, _cx: &mut Cx<'_>) -> R<bool> {
        self.release();
        self.have_more_elements()
    }

    fn end_list(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_end_of_content()
    }

    /// `ReadInt4`, `ReadUint4`, `ReadInt8` (`ExpectIntegerTag` and `ReadStdSigned` /
    /// `ReadStdUnsigned`) and `ReadInt2` (an Int4 that must fit), or the skip (`SkipSNumber`,
    /// `SkipUNumber`: the tag and the bytes, no range).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:569-591
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2077-2087
    // (`CObjectIStream::ReadInt2`, objistr.cpp:1697-1704: "integer overflow".)
    fn integer(
        &mut self,
        _ty: u16,
        storage: Storage,
        skip: bool,
        _cx: &mut Cx<'_>,
    ) -> R<Option<i64>> {
        self.expect_integer_tag()?;
        if skip {
            self.skip_tag_data()?;
            return Ok(None);
        }
        let value = match storage {
            Storage::Int4 => i64::from(self.std_signed(4)? as i32),
            Storage::Uint4 => i64::from(self.std_unsigned4()?),
            Storage::Int8 => self.std_signed(8)?,
            Storage::Int2 => {
                let value = self.std_signed(4)? as i32;
                if i16::try_from(value).is_err() {
                    return Err(ReadError::Serial("integer overflow".into()));
                }
                i64::from(value)
            }
        };
        Ok(Some(value))
    }

    /// `ReadEnum` (the skip, `SkipEnum`, reads the same): an INTEGER with named values takes any
    /// number, an ENUMERATED one of its values.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2023-2039
    fn enumerated(
        &mut self,
        _ty: u16,
        values: &'static [EnumValue],
        integer: bool,
        _skip: bool,
        _cx: &mut Cx<'_>,
    ) -> R<i32> {
        if integer {
            self.expect_integer_tag()?;
            return Ok(self.std_signed(4)? as i32);
        }
        self.expect_sys_tag(TAG_ENUMERATED)?;
        let value = self.std_signed(4)? as i32;
        if !values.iter().any(|named| named.value == value) {
            return Err(ReadError::Serial("invalid value of enumerated type".into()));
        }
        Ok(value)
    }

    /// `ReadBool` / `SkipBool`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:551-558
    fn boolean(&mut self, _ty: u16, _skip: bool, _cx: &mut Cx<'_>) -> R<bool> {
        self.expect_sys_tag(TAG_BOOLEAN)?;
        self.expect_short_length(1)?;
        Ok(self.get()? != 0)
    }

    /// `ReadNull` / `SkipNull`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1505-1510
    fn null(&mut self, _ty: u16, _skip: bool, _cx: &mut Cx<'_>) -> R<()> {
        self.expect_sys_tag(TAG_NULL)?;
        self.expect_short_length(0)
    }

    /// `ReadDouble` / `SkipFNumber`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:603-646
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2089-2103
    fn real(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<Option<f64>> {
        let too_long =
            || ReadError::Serial(format!("too long REAL data: length > {MAX_DOUBLE_LENGTH}"));
        self.expect_sys_tag(TAG_REAL)?;
        let length = self.length()?;
        if skip {
            if length != 0 {
                if length > MAX_DOUBLE_LENGTH as u64 {
                    return Err(too_long());
                }
                self.get()?;
                self.take(length as usize - 1)?;
            }
            return Ok(None);
        }
        if length < 2 {
            if length == 0 {
                return Ok(Some(0.0));
            }
            return match self.get()? {
                REAL_PLUS_INFINITY => Ok(Some(f64::INFINITY)),
                REAL_MINUS_INFINITY => Ok(Some(f64::NEG_INFINITY)),
                REAL_NOT_A_NUMBER => Ok(Some(NCBI_NAN)),
                REAL_MINUS_ZERO => Ok(Some(-0.0)),
                _ => Err(ReadError::Serial("Unrecognized REAL data".into())),
            };
        }
        if length > MAX_DOUBLE_LENGTH as u64 {
            return Err(too_long());
        }
        let form = self.get()?;
        if form & REAL_ENCODING_MASK != REAL_DECIMAL {
            return Err(ReadError::Serial(format!(
                "Unsupported encoding of REAL data: encoding = {form}"
            )));
        }
        let digits = self.take(length as usize - 1)?;
        let digits = self.src.bytes(digits);
        // The C string ends at its first 0 byte.
        let text = digits.split(|&byte| byte == 0).next().unwrap_or(&[]);
        let read = string_to_double_posix(text, true);
        if read.end < text.len() {
            return Err(ReadError::Serial("bad REAL data string".into()));
        }
        Ok(Some(read.value))
    }

    /// `ReadString` (a VisibleString's bytes outside 0x20..0x7E fixed) / `ReadStringStore`
    /// (application tag 1, no fix), or the skips (`SkipString`, `SkipStringStore`: the tag and
    /// the bytes).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:785-827
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2105-2115
    // (A string the reader makes is empty, so its bytes go through the string form of
    // `FixVisibleChars`; the buffer form, for a string of the same length as the old one, fixes
    // the same bytes.)
    fn string(&mut self, _ty: u16, kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        match kind {
            StringKind::Store => self.expect_sys_tag(TAG_STRING_STORE)?,
            StringKind::Utf8 => self.expect_string_tag(true, cx)?,
            StringKind::Visible => self.expect_string_tag(false, cx)?,
        }
        if skip {
            self.skip_tag_data()?;
            return Ok(Vec::new());
        }
        let length = self.data_length()?;
        let range = self.take(length)?;
        let mut string = self.src.bytes(range).to_vec();
        if kind == StringKind::Visible {
            self.fix_visible_chars(&mut string, cx);
        }
        Ok(string)
    }

    /// `ReadBytes` of an OCTET STRING (or a BIT STRING, without its first byte), or the skip
    /// (`SkipByteBlock`: an OCTET STRING's tag and bytes).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1448-1473
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2124-2128
    // (A BIT STRING of length 0 has `SIZE_MAX` bytes: its bytes reach the end of the data.)
    fn octets(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<Option<Vec<u8>>> {
        if skip {
            self.expect_sys_tag(TAG_OCTET_STRING)?;
            self.skip_tag_data()?;
            return Ok(None);
        }
        let length = match self.peek(0)? {
            TAG_OCTET_STRING => {
                self.expect_sys_tag(TAG_OCTET_STRING)?;
                self.data_length()?
            }
            TAG_BIT_STRING => {
                self.expect_sys_tag(TAG_BIT_STRING)?;
                let length = self.data_length()?.wrapping_sub(1);
                self.get()?;
                length
            }
            _ => {
                return Err(ReadError::Serial(
                    "Unable to identify the type of byte block".into(),
                ))
            }
        };
        let range = self.take(length)?;
        Ok(Some(self.src.bytes(range).to_vec()))
    }

    /// `CObjectIStream::StartDelayBuffer`: the buffer starts at the member's identifier.
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    // (The member's value is skipped into the buffer.)
    fn delay_start(&mut self, _cx: &mut Cx<'_>) -> R<usize> {
        self.held = Some(self.pos);
        Ok(self.pos)
    }

    /// `EndDelayBuffer`, and the buffer parsed by a new binary stream over its bytes
    /// (`CDelayBuffer::Update`, delaybuf.cpp:63-84: `UpdateDelayedBuffer` reads the member's type
    /// under a frame of its own): LOSAT's value of a delayed member is read text
    /// (`Value::Delayed`), so a binary buffer is parsed at once rather than when it is first
    /// used.
    fn delay_end(&mut self, start: usize, ty: u16, cx: &mut Cx<'_>) -> R<Value> {
        self.held = None;
        let bytes = self.src.bytes(start..self.pos).to_vec();
        let mut buffer = Ber::new(SliceSource::new(&bytes));
        buffer.utf8_warned = self.utf8_warned;
        let value = read_delayed(&mut buffer, ty, cx);
        self.utf8_warned = buffer.utf8_warned;
        value
    }
}

impl<S: AsnSource> Ber<S> {
    /// `ExpectShortLength`.
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:364-371
    fn expect_short_length(&mut self, length: usize) -> R<()> {
        let got = self.short_length()?;
        if got != length {
            return Err(ReadError::Serial(format!(
                "unexpected length: {got}, should be: {length}"
            )));
        }
        Ok(())
    }

    /// `FixVisibleChars` with the method `eFNP_ReplaceAndWarn`: each byte outside 0x20..0x7E
    /// posted (the string as fixed so far, the position after its bytes) and replaced by `#`.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:703-749
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:815-819
    // (`CObjectIStreamAsnBinary` keeps `Location` of its base: `GetPosition`, objistr.cpp:807-813,
    // is "byte <position>".)
    fn fix_visible_chars(&mut self, string: &mut [u8], cx: &mut Cx<'_>) {
        for i in 0..string.len() {
            let byte = string[i];
            if !(b' '..=b'~').contains(&byte) {
                cx.post_bad_char(byte, string, StreamPosition::Byte(self.pos as u64));
                string[i] = b'#';
            }
        }
    }
}
