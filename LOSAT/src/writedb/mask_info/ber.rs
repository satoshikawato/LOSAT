//! ASN.1 binary (`MSerial_AsnBinary`, BER): the primitives of `CObjectIStreamAsnBinary` on a byte
//! slice that holds the rest of the file.
//!
//! The types of a mask file are read with "automatic tags": a SEQUENCE is `30 <len>`, each of
//! its members is `[n]` (`a0 | n`, constructed, context specific) around the value, a CHOICE is
//! its variant `[n]` around the value, and every constructed item ends with `00 00` when
//! its length is indefinite (`80`); a definite length is accepted as well.

use super::reader::{reserve_string, throw_error, Fail, Input, StringPost, R};
use super::schema::{ChoiceDef, ClassDef, EnumDef};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_text_reader::{PointerType, FIND_NAME};
use crate::report::serial::encoder::StreamPosition;
use crate::seqid::ber::StdSignedError;

/// A `ThrowError` of `CObjectIStreamAsnBinary` (objistrasnb.cpp) at `line` of `function` (`None`:
/// the template functions of the numbers, whose report names none).
const fn asnb_site(line: u32, kind: &'static str, function: Option<&'static str>) -> ThrowSite {
    ThrowSite {
        file: "src/serial/objistrasnb.cpp",
        line,
        kind,
        function,
    }
}

const FORMAT_ERROR: &str = "CSerialException::eFormatError";
const OVERFLOW: &str = "CSerialException::eOverflow";
const PEEK_LONG_TAG: ThrowSite = asnb_site(
    161,
    OVERFLOW,
    Some("ncbi::CObjectIStreamAsnBinary::PeekLongTag"),
);
const UNEXPECTED_TAG_CLASS_BYTE: ThrowSite = asnb_site(
    177,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedTagClassByte"),
);
const UNEXPECTED_SHORT_LENGTH: ThrowSite = asnb_site(
    185,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedShortLength"),
);
const UNEXPECTED_CONTINUATION: ThrowSite = asnb_site(
    197,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedContinuation"),
);
const UNEXPECTED_LONG_LENGTH: ThrowSite = asnb_site(
    202,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedLongLength"),
);
const UNEXPECTED_TAG_VALUE: ThrowSite = asnb_site(
    216,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedTagValue"),
);
const UNEXPECTED_SYS_TAG_BYTE: ThrowSite = asnb_site(
    377,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedSysTagByte"),
);
const READ_LENGTH_LONG: &str = "ncbi::CObjectIStreamAsnBinary::ReadLengthLong";
const PEEK_CLASS_TAG: ThrowSite = asnb_site(
    231,
    OVERFLOW,
    Some("ncbi::CObjectIStreamAsnBinary::PeekClassTag"),
);
const UNEXPECTED_MEMBER: ThrowSite = asnb_site(
    1066,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedMember"),
);

pub(super) struct BerIn<'a> {
    data: &'a [u8],
    pos: usize,
    /// `m_CurrentTagLength`: the bytes of the tag that was peeked, before its length.
    tag_len: usize,
    /// `m_CurrentDataLimit`: the end of the definite-length item being read, 0 for an
    /// indefinite one.
    limit: i64,
    /// `m_DataLimits`
    limits: Vec<i64>,
}

// NCBI reference (598d8ae6): c++/include/serial/impl/objstrasnb.hpp:60-130
const UNIVERSAL_CONSTRUCTED: u8 = 0x20;
const CONTEXT_CONSTRUCTED: u8 = 0xa0;
const CLASS_AND_CONSTRUCTED_MASK: u8 = 0xe0;
const TAG_BOOLEAN: u8 = 0x01;
const TAG_INTEGER: u8 = 0x02;
const TAG_APPLICATION_INTEGER: u8 = 0x42;
const TAG_NULL: u8 = 0x05;
const TAG_ENUMERATED: u8 = 0x0a;
const TAG_UTF8_STRING: u8 = 0x0c;
const TAG_VISIBLE_STRING: u8 = 0x1a;
/// `MakeTagByte` of `eApplication`, `ePrimitive` and `eObjectReference`: the object index of `@n`.
const TAG_OBJECT_REFERENCE: u8 = 0x5e;
/// `MakeTagByte` of `eApplication`, `eConstructed` and `eLongTag`: a class name follows, then its
/// value.
const TAG_OTHER_POINTER: u8 = 0x7f;
const TAG_SEQUENCE: i32 = 16;
const TAG_SET: i32 = 17;
const LONG_TAG: u8 = 0x1f;

impl<'a> BerIn<'a> {
    /// Bytes read so far (where the next object of the file starts).
    pub(super) fn offset(&self) -> usize {
        self.pos
    }

    pub(super) fn new(data: &'a [u8]) -> Self {
        BerIn {
            data,
            pos: 0,
            tag_len: 0,
            limit: 0,
            limits: Vec::new(),
        }
    }

    /// `ThrowError` with `fFormatError` or `fOverflow` and a message, at `site` at the current
    /// position (returned as the [`Fail`]: `return Err(self.throw(..))`, see [`throw_error`]).
    #[cold]
    #[inline(never)]
    fn throw(&self, site: ThrowSite, message: impl std::fmt::Display) -> Fail {
        throw_error(site, &self.position(), message.to_string())
    }

    /// `PeekChar` of `m_Input` at `i`: the end of the data throws `CEofException`.
    #[inline]
    fn peek(&self, i: usize) -> R<u8> {
        match self.data.get(self.pos + i) {
            Some(&byte) => Ok(byte),
            None => Err(Fail::Eof),
        }
    }

    /// `m_Input.GetChar()`
    fn byte(&mut self) -> R<u8> {
        let b = self.peek(0)?;
        self.pos += 1;
        Ok(b)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:132-171
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:158-170,175-188
    /// `PeekTag`: the tag number of the tag at the current position; `tag_len` is its length.
    fn peek_tag(&mut self) -> R<i32> {
        let first = self.peek(0)?;
        if first & LONG_TAG != LONG_TAG {
            self.tag_len = 1;
            return Ok(i32::from(first & LONG_TAG));
        }
        let mut tag: i32 = 0;
        let mut i = 1;
        loop {
            if tag >= (1 << 24) {
                return Err(self.throw(PEEK_LONG_TAG, format_args!("tag number is too big: {tag}")));
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

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:138-156,167-199
    fn peek_tag_of_class(&mut self, class_and_constructed: u8) -> R<i32> {
        let first = self.peek(0)?;
        if first & CLASS_AND_CONSTRUCTED_MASK != class_and_constructed {
            return Err(self.unexpected_tag_class_byte(first, class_and_constructed));
        }
        self.peek_tag()
    }

    /// `UnexpectedTagClassByte` of the first tag byte and the expected class byte.
    #[cold]
    #[inline(never)]
    fn unexpected_tag_class_byte(&self, first: u8, class_and_constructed: u8) -> Fail {
        self.throw(
            UNEXPECTED_TAG_CLASS_BYTE,
            format_args!(
                "unexpected tag: {}, should be {}",
                tag_to_string(first),
                tag_to_string(class_and_constructed)
            ),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:171-178
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:158-185,209-217
    fn expect_tag(&mut self, class_and_constructed: u8, value: i32) -> R<()> {
        let tag = self.peek_tag_of_class(class_and_constructed)?;
        if tag != value {
            return Err(self.unexpected_tag_value(tag, value));
        }
        Ok(())
    }

    /// `UnexpectedTagValue` of the tag class, the tag and the tag value. The class was checked
    /// (universal or context specific): no "Application " or "Private " before the number.
    #[cold]
    #[inline(never)]
    fn unexpected_tag_value(&self, tag: i32, value: i32) -> Fail {
        self.throw(
            UNEXPECTED_TAG_VALUE,
            format_args!("Unexpected tag: {tag}, expected: {value}"),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:394-410
    fn read_length_long(&mut self, byte: u8) -> R<u64> {
        let mut length_length = usize::from(byte - 0x80);
        if length_length == 0 {
            return Err(self.throw(
                asnb_site(398, FORMAT_ERROR, Some(READ_LENGTH_LONG)),
                "unexpected indefinite length",
            ));
        }
        if length_length > 8 {
            return Err(self.throw(
                asnb_site(401, OVERFLOW, Some(READ_LENGTH_LONG)),
                "length overflow",
            ));
        }
        let mut length = u64::from(self.byte()?);
        if length == 0 {
            return Err(self.throw(
                asnb_site(405, FORMAT_ERROR, Some(READ_LENGTH_LONG)),
                "illegal length start",
            ));
        }
        loop {
            length_length -= 1;
            if length_length == 0 {
                break;
            }
            length = (length << 8) | u64::from(self.byte()?);
        }
        Ok(length)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:262-311
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:236-246,253-263
    /// The tag that `peek_tag` found is followed by the length of a constructed item.
    fn expect_indefinite_length(&mut self) -> R<()> {
        self.pos += self.tag_len;
        let len = self.byte()?;
        self.limits.push(self.limit);
        self.limit = if len == 0x80 {
            0
        } else if len > 0x80 {
            // The stream position (`GetStreamPosAsInt8`) plus the length that `ReadLengthLong`
            // reads: C++ does not say which operand is evaluated first. The compiler of the 2.17.0
            // binaries evaluates the position first, so the limit is counted from the first length
            // byte: the item ends `n` bytes (the length of the length) before where its length
            // says. makeblastdb 2.17.0 fails with "EndOfContentsByte expected" on such an item
            // (checked on the `long_len_*` cases of the fixtures).
            //
            // The sum is an `Int8` added to a `size_t`: computed in 64 unsigned bits, it wraps, and
            // is stored in the `Int8` `m_CurrentDataLimit` (objistrasnb.hpp:267). A limit of 0 is
            // then an indefinite length (`HaveMoreElements`, `ExpectEndOfContent`), a negative one
            // ends the members at once (checked with makeblastdb 2.17.0, the b_length_long_wrap_*
            // fixtures).
            let base = self.pos as i64;
            let long_length = self.read_length_long(len)?;
            (base as u64).wrapping_add(long_length) as i64
        } else {
            self.pos as i64 + i64::from(len)
        };
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:174-182,242-250
    fn have_more_elements(&mut self) -> R<bool> {
        if self.limit != 0 {
            return Ok(self.limit > self.pos as i64);
        }
        Ok(self.peek(0)? != 0)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:358-405
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:70-90
    fn expect_end_of_content(&mut self) -> R<()> {
        if self.limit != 0 {
            if self.limit != self.pos as i64 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
        } else {
            if self.peek(0)? != 0 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
            self.pos += 1;
            if self.peek(0)? != 0 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
            self.pos += 1;
        }
        self.limit = self.limits.pop().unwrap_or(0);
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:309-329,338-352
    fn read_short_length(&mut self) -> R<usize> {
        self.pos += self.tag_len;
        self.tag_len = 0;
        let byte = self.byte()?;
        if byte >= 0x80 {
            return Err(self.throw(UNEXPECTED_LONG_LENGTH, "ShortLength expected"));
        }
        Ok(usize::from(byte))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:374-379
    /// `UnexpectedSysTagByte` of `byte`: the tag at the current position is not `byte`.
    #[cold]
    #[inline(never)]
    fn unexpected_sys_tag_byte(&self, got: u8, byte: u8) -> Fail {
        self.throw(
            UNEXPECTED_SYS_TAG_BYTE,
            format_args!(
                "unexpected tag: {}, should be: {}",
                tag_to_string(got),
                tag_to_string(byte)
            ),
        )
    }

    /// `ExpectSysTagByte` of `byte` (objistrasnb.inl:115-123): the exact tag byte.
    fn expect_sys_tag_byte(&mut self, byte: u8) -> R<()> {
        let got = self.peek(0)?;
        if got != byte {
            return Err(self.unexpected_sys_tag_byte(got, byte));
        }
        self.tag_len = 1;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:145-160
    fn expect_integer_tag(&mut self) -> R<()> {
        match self.peek(0)? {
            // The application tag only sets `eReadAsBigInt`, which has no effect on the number
            // read here (checked with the oracle on `b_small_application_integer.asnb`).
            TAG_INTEGER | TAG_APPLICATION_INTEGER => {
                self.tag_len = 1;
                Ok(())
            }
            got => Err(self.unexpected_sys_tag_byte(got, TAG_INTEGER)),
        }
    }

    /// `ReadStdSigned` (`crate::seqid::ber::read_std_signed`, the one port of it): `size` is
    /// the `sizeof` of `T`; the result is the value of `T` widened to `i64`.
    fn read_std_signed(&mut self, size: usize) -> R<i64> {
        let length = self.read_short_length()?;
        let result = crate::seqid::ber::read_std_signed(length, size, || self.byte());
        result.map_err(|error| match error {
            StdSignedError::ZeroLength => {
                self.throw(asnb_site(480, FORMAT_ERROR, None), "zero length of number")
            }
            StdSignedError::Overflow { line } => {
                self.throw(asnb_site(line, OVERFLOW, None), "overflow error")
            }
            StdSignedError::Read(fail) => fail,
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:514-549
    fn read_std_unsigned32(&mut self) -> R<u32> {
        let mut length = self.read_short_length()?;
        if length == 0 {
            return Err(self.throw(asnb_site(519, FORMAT_ERROR, None), "zero length of number"));
        }
        let mut n: u32;
        if length > 4 {
            while length > 4 {
                length -= 1;
                if self.byte()? != 0 {
                    return Err(self.throw(asnb_site(527, OVERFLOW, None), "overflow error"));
                }
            }
            length -= 1;
            n = u32::from(self.byte()?);
        } else if length == 4 {
            length -= 1;
            n = u32::from(self.byte()?);
            if n & 0x80 != 0 {
                return Err(self.throw(asnb_site(537, OVERFLOW, None), "overflow error"));
            }
        } else {
            n = 0;
        }
        while length > 0 {
            length -= 1;
            n = (n << 8) | u32::from(self.byte()?);
        }
        Ok(n)
    }
}

impl Input for BerIn<'_> {
    /// `GetPosition()` of the binary stream: the bytes read from the start of the object
    /// (`m_Input.GetStreamPosAsInt8()`; each object of the file is read by a stream of its own).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:118-121
    fn position(&self) -> String {
        format!("byte {}", self.pos)
    }

    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1029-1042
    // (A SEQUENCE: universal tag 16.)
    fn begin_class(&mut self) -> R<()> {
        self.expect_tag(UNIVERSAL_CONSTRUCTED, TAG_SEQUENCE)?;
        self.expect_indefinite_length()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1096-1179
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:462-473
    // (The classes of a mask file have the "automatic" tag type: member `n` has the tag `[n]`.)
    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        if !self.have_more_elements()? {
            return Ok(None);
        }
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        match usize::try_from(tag) {
            Ok(index) if index >= pos && index < def.members.len() => Ok(Some(index)),
            _ => Err(self.unexpected_member(tag, def.members.iter().map(|m| m.name))),
        }
    }

    fn end_member(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    fn end_class(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    fn begin_choice(&mut self) -> R<()> {
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1306-1332,1344-1353
    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize> {
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        match usize::try_from(tag) {
            Ok(index) if index < def.variants.len() => Ok(index),
            _ => Err(self.unexpected_member(tag, def.variants.iter().map(|v| v.name))),
        }
    }

    fn end_variant(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    fn end_choice(&mut self) -> R<()> {
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:921-934,956-1000
    // (`SEQUENCE OF`: universal tag 16, `SET OF`: 17.)
    fn begin_container(&mut self, set: bool) -> R<()> {
        self.expect_tag(
            UNIVERSAL_CONSTRUCTED,
            if set { TAG_SET } else { TAG_SEQUENCE },
        )?;
        self.expect_indefinite_length()
    }

    fn begin_element(&mut self) -> R<bool> {
        self.have_more_elements()
    }

    fn end_container(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    fn read_i32(&mut self) -> R<i32> {
        self.expect_integer_tag()?;
        Ok(self.read_std_signed(4)? as i32)
    }

    fn read_i64(&mut self) -> R<i64> {
        self.expect_integer_tag()?;
        self.read_std_signed(8)
    }

    fn read_u32(&mut self) -> R<u32> {
        self.expect_integer_tag()?;
        self.read_std_unsigned32()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:551-557
    fn read_bool(&mut self) -> R<bool> {
        self.expect_sys_tag_byte(TAG_BOOLEAN)?;
        let length = self.read_short_length()?;
        if length != 1 {
            return Err(self.throw(
                UNEXPECTED_SHORT_LENGTH,
                format_args!("unexpected length: {length}, should be: 1"),
            ));
        }
        Ok(self.byte()? != 0)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1502-1508
    fn read_null(&mut self) -> R<()> {
        self.expect_sys_tag_byte(TAG_NULL)?;
        let length = self.read_short_length()?;
        if length != 0 {
            return Err(self.throw(
                UNEXPECTED_SHORT_LENGTH,
                format_args!("unexpected length: {length}, should be: 0"),
            ));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:281-314
    // (`SERIAL_READ_ANY_VISIBLESTRING_TAG` is 1 by default: a UTF8String is accepted for a
    // VisibleString, with a warning.)
    fn read_string(&mut self, posts: &mut Vec<StringPost>) -> R<Vec<u8>> {
        match self.peek(0)? {
            TAG_VISIBLE_STRING => {}
            TAG_UTF8_STRING => posts.push(StringPost::Utf8Tag),
            got => return Err(self.unexpected_sys_tag_byte(got, TAG_VISIBLE_STRING)),
        }
        self.pos += 1;
        self.tag_len = 0;
        let byte = self.byte()?;
        let length = if byte >= 0x80 {
            self.read_length_long(byte)?
        } else {
            u64::from(byte)
        };
        // `ReadBytes` of `s` and `length`: `GetChars` of `str` and `count`, which reserves
        // `count` bytes first when they are not all in the input buffer, and reads to the end of
        // the file.
        let remaining = (self.data.len() - self.pos) as u64;
        if length > remaining {
            reserve_string(length)?;
            self.pos = self.data.len();
            return Err(Fail::Eof);
        }
        let length = length as usize;
        let mut s = self.data[self.pos..self.pos + length].to_vec();
        self.pos += length;
        fix_visible_chars(&mut s, self.pos, posts);
        Ok(s)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2023-2039
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        if def.integer {
            self.expect_integer_tag()?;
            return Ok(self.read_std_signed(4)? as i32);
        }
        self.expect_sys_tag_byte(TAG_ENUMERATED)?;
        let value = self.read_std_signed(4)? as i32;
        if !def.values.iter().any(|&(_, v)| v == value) {
            return Err(Fail::serial(NcbiException::new(
                FIND_NAME,
                format!("invalid value of enumerated type: {value}"),
            )));
        }
        Ok(value)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1861-1930
    // (The binary stream's own `ReadPointer`: the same cases as `CObjectIStream::ReadPointer` with
    // `ReadPointerType`, objistrasnb.cpp:1797-1823, inlined.)
    /// `ReadPointerType`: the tag byte, a NULL (`05 00`) taken as `ReadNull` takes it.
    fn read_pointer_type(&mut self) -> R<PointerType> {
        Ok(match self.peek(0)? {
            TAG_NULL => {
                self.read_null()?;
                PointerType::Null
            }
            TAG_OBJECT_REFERENCE => PointerType::Object,
            TAG_OTHER_POINTER => PointerType::Other,
            _ => PointerType::This,
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1831-1836
    // (`ReadPointer` does not take the tag: `m_CurrentTagLength` is still 0 from the end of the
    // item before (`EndOfTag`, `ExpectIndefiniteLength`, `ExpectEndOfContent`), so `ReadShortLength`
    // reads the tag byte itself as the length, 94, and `ReadStdSigned` of the 8-byte
    // `TObjectIndex` the 94 bytes after it. checked with makeblastdb 2.17.0, the b_pointer_*_at
    // fixtures.)
    fn read_object_pointer(&mut self) -> R<()> {
        self.read_std_signed(8).map(|_| ())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1838-1843
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:219-240
    /// `ReadOtherPointer`: the class name in the bytes of the long tag (7 bits each, the last one
    /// without the high bit), then the length of the value (`ExpectIndefiniteLength`).
    fn read_other_pointer(&mut self) -> R<Vec<u8>> {
        let mut name = Vec::new();
        let mut i = 1;
        loop {
            let c = self.peek(i)?;
            i += 1;
            name.push(c & 0x7f);
            if c & 0x80 == 0 {
                break;
            }
            if i > 1024 {
                return Err(self.throw(PEEK_CLASS_TAG, "tag number is too big (greater than 1024)"));
            }
        }
        self.tag_len = i;
        self.expect_indefinite_length()?;
        Ok(name)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1845-1848
    fn read_other_pointer_end(&mut self) -> R<()> {
        self.expect_end_of_content()
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1058-1067
impl BerIn<'_> {
    /// `UnexpectedMember` of the tag and the items, of a SEQUENCE or CHOICE whose members or
    /// variants are `names` (their automatic tags are 0, 1, ... in order).
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, tag: i32, names: impl Iterator<Item = &'static str>) -> Fail {
        let mut message = format!("unexpected member: [{tag}], should be one of: ");
        for (tag, name) in names.enumerate() {
            message.push_str(&format!("{name}[{tag}] "));
        }
        self.throw(UNEXPECTED_MEMBER, message)
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:316-372
/// `TagToString` of `byte`: the class, `constructed/`, the universal type's name or the number,
/// and the byte.
fn tag_to_string(byte: u8) -> String {
    let class = match byte & 0xc0 {
        0x40 => "application/",
        0x80 => "contextspecific/",
        0xc0 => "private/",
        _ => "",
    };
    let constructed = if byte & 0x20 != 0 { "constructed/" } else { "" };
    let value = byte & 0x1f;
    if byte & 0xc0 != 0 {
        return format!("{class}{constructed}{value} ({byte})");
    }
    let name = match value {
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
    };
    format!("{constructed}{name} ({byte})")
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:703-749
// (`ReadStringValue`, objistrasnb.cpp:803-826, fixes a new string after `ReadBytes`, before
// `EndOfTag`: the position is the byte after the string. `SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_
// READING` is 0.)
/// `FixVisibleChars` of the string `s`, read up to the byte `end` of the object: each byte
/// outside `' '..='~'` becomes `#`, posted with the string as fixed so far (the bytes before it
/// already `#`).
fn fix_visible_chars(s: &mut [u8], end: usize, posts: &mut Vec<StringPost>) {
    for i in 0..s.len() {
        let byte = s[i];
        if !(b' '..=b'~').contains(&byte) {
            posts.push(StringPost::BadChar {
                position: StreamPosition::Byte(end as u64),
                byte,
                string: s.to_vec(),
            });
            s[i] = b'#';
        }
    }
}
