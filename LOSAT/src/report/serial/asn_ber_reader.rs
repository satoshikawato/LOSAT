//! NCBI's ASN.1 binary object input stream (`CObjectIStreamAsnBinary`, BER): the primitives that
//! the object stacks of LOSAT's readers call through `asn_input::Input` (the table reader of
//! `asn_types` and the mask reader of `writedb::mask_info`) on any `AsnSource`: tags, lengths and
//! the data limits of definite lengths, the end-of-contents octets, numbers, strings, REAL, OCTET
//! STRING, BOOLEAN, NULL, NCBI's pointer forms, and what NCBI's stream throws, at its place in
//! NCBI's source. The header reader of a BLAST database (`seqid::ber`) shares its port of
//! `ReadStdSigned` (`read_std_signed`) and its tag bytes.
//!
//! NCBI's modules are read with "automatic tags": a SEQUENCE is `30 <len>`, each of its members is
//! `[n]` (`a0 | n`, constructed, context specific) around the value, a CHOICE is its variant `[n]`
//! around the value, and every constructed item ends with `00 00` when its length is indefinite
//! (`80`); a definite length is accepted as well (`USE_DEF_LEN`).
//!
//! Each object of a file is read by an object stream of its own in NCBI (`s_ReadObject` of
//! makeblastdb, the stream read as `MSerial_AsnBinary` into the object): the positions of its reports (`byte N`) and
//! its data limits count from the object's first byte (`BerReader::new_object`).

use super::asn_input::{
    invalid_enum, reserve, throw_error, Encoding, EnumItems, Fail, Input, Items, Names, Reserve,
    StringKind, StringPost, R,
};
use super::asn_text_reader::{AsnSource, PointerType, FIND_NAME};
use super::encoder::StreamPosition;
use crate::ncbi_exception::ThrowSite;
use crate::ncbi_str::string_to_double_posix;

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
const PEEK_CLASS_TAG: ThrowSite = asnb_site(
    231,
    OVERFLOW,
    Some("ncbi::CObjectIStreamAsnBinary::PeekClassTag"),
);
const UNEXPECTED_SYS_TAG_BYTE: ThrowSite = asnb_site(
    377,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedSysTagByte"),
);
const READ_LENGTH_LONG: &str = "ncbi::CObjectIStreamAsnBinary::ReadLengthLong";
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:603-646
const READ_DOUBLE: &str = "ncbi::CObjectIStreamAsnBinary::ReadDouble";
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2089-2103
const SKIP_F_NUMBER: ThrowSite = asnb_site(
    2095,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::SkipFNumber"),
);
const UNEXPECTED_MEMBER: ThrowSite = asnb_site(
    1066,
    FORMAT_ERROR,
    Some("ncbi::CObjectIStreamAsnBinary::UnexpectedMember"),
);
// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1448-1473
// (`ThrowError` with `fUnknownValue`: `ThrowError1` has no case for the flag, whose `default` is
// `eIoError`, objistr.cpp:828-852.)
const BEGIN_BYTES: ThrowSite = asnb_site(
    1470,
    "CSerialException::eIoError",
    Some("ncbi::CObjectIStreamAsnBinary::BeginBytes"),
);

// NCBI reference (598d8ae6): c++/include/serial/impl/objstrasnb.hpp:60-151
/// `MakeTagByte` of `eUniversal`, `eConstructed` and 0.
pub(crate) const UNIVERSAL_CONSTRUCTED: u8 = 0x20;
/// `MakeTagByte` of `eContextSpecific`, `eConstructed` and 0: a member or variant `[n]` is `0xa0 | n`.
pub(crate) const CONTEXT_CONSTRUCTED: u8 = 0xa0;
const CLASS_AND_CONSTRUCTED_MASK: u8 = 0xe0;
/// `eLongTag`: the tag number bits of a first tag byte, all set when the number follows.
pub(crate) const LONG_TAG: u8 = 0x1f;
pub(crate) const TAG_BOOLEAN: u8 = 0x01;
pub(crate) const TAG_INTEGER: u8 = 0x02;
const TAG_BIT_STRING: u8 = 0x03;
pub(crate) const TAG_OCTET_STRING: u8 = 0x04;
pub(crate) const TAG_NULL: u8 = 0x05;
const TAG_REAL: u8 = 0x09;
pub(crate) const TAG_ENUMERATED: u8 = 0x0a;
pub(crate) const TAG_UTF8_STRING: u8 = 0x0c;
pub(crate) const TAG_VISIBLE_STRING: u8 = 0x1a;
/// `MakeTagByte` of `eUniversal`, `eConstructed` and `eSequence`: the first byte of a SEQUENCE and of a
/// SEQUENCE OF.
pub(crate) const TAG_SEQUENCE: u8 = 0x30;
/// `MakeTagByte` of `eApplication`, `ePrimitive` and `eInteger`, which `ExpectIntegerTag` also takes.
pub(crate) const TAG_APPLICATION_INTEGER: u8 = 0x42;
/// `MakeTagByte` of `eApplication`, `ePrimitive` and `eStringStore`.
pub(crate) const TAG_STRING_STORE: u8 = 0x41;
/// `MakeTagByte` of `eApplication`, `ePrimitive` and `eObjectReference`: the object index of `@n`.
const TAG_OBJECT_REFERENCE: u8 = 0x5e;
/// `MakeTagByte` of `eApplication`, `eConstructed` and `eLongTag`: a class name follows, then its value.
const TAG_OTHER_POINTER: u8 = 0x7f;
/// The tag numbers of SEQUENCE (and SEQUENCE OF) and SET OF.
pub(crate) const SEQUENCE_NUMBER: i32 = 16;
pub(crate) const SET_NUMBER: i32 = 17;
/// The first content byte of a REAL: the encoding bits (`eDecimalEncoding`) and the four special
/// values of a REAL of one content byte (the binary writer's too, `asn_ber`).
const REAL_ENCODING_MASK: u8 = 0xc0;
pub(crate) const REAL_DECIMAL: u8 = 0x00;
pub(crate) const REAL_PLUS_INFINITY: u8 = 0x40;
pub(crate) const REAL_MINUS_INFINITY: u8 = 0x41;
pub(crate) const REAL_NOT_A_NUMBER: u8 = 0x42;
const REAL_MINUS_ZERO: u8 = 0x43;
/// `kMaxDoubleLength`.
const MAX_DOUBLE_LENGTH: u64 = 256;

/// The state of one `CObjectIStreamAsnBinary` over the bytes of `src`.
pub(crate) struct BerReader<S: AsnSource> {
    pub(crate) src: S,
    /// The position in `src` of the first byte of the object being read (`new_object`).
    base: usize,
    /// The position in `src` of the next byte; `at - base` is the position of NCBI's input
    /// buffer (the bytes of the object read so far, [`pos`](Self::pos)). The cursor is the
    /// source's position, so a byte is read without adding the object's start to it.
    at: usize,
    /// `m_CurrentTagLength`: the bytes of the tag that was peeked, before its length.
    tag_len: usize,
    /// `m_CurrentDataLimit`: the end of the definite-length item being read, 0 for an
    /// indefinite one.
    limit: i64,
    /// `m_DataLimits`
    limits: Vec<i64>,
    /// What LOSAT's refusal of a length calls the input (`asn_input::reserve`).
    names: Names,
}

impl<S: AsnSource> BerReader<S> {
    /// The binary input over `src`, which LOSAT's own messages call by `names`.
    pub(crate) fn new(src: S, names: Names) -> Self {
        BerReader {
            src,
            base: 0,
            at: 0,
            tag_len: 0,
            limit: 0,
            limits: Vec::new(),
            names,
        }
    }

    /// The position in the source of the next byte (where the next object of a file starts).
    pub(crate) fn offset(&self) -> usize {
        self.at
    }

    /// The bytes of the object read so far (`GetStreamPosAsInt8()`).
    #[inline(always)]
    fn pos(&self) -> usize {
        self.at - self.base
    }

    /// The reader will not ask for a byte before the current one again.
    #[inline]
    pub(crate) fn release(&mut self) {
        self.src.release_before(self.at);
    }

    /// The bytes that the source keeps from the cursor on.
    ///
    /// Three forms make most items of an object: a short tag and an indefinite length (the start
    /// of an item), the end-of-contents bytes (its end) and an INTEGER of a short length. When all
    /// the bytes of such a form are kept and NCBI reads them without an error, they are read from
    /// this slice in locals (`kept_indefinite_item`, `expect_end_of_content`, `kept_integer`). Any
    /// other bytes (a long tag or length, a form that goes past the slice: the end of a block or
    /// of the input, an error) are read one at a time (`peek`), which asks the source for more
    /// past its slice. Both readings leave the same state after the same bytes (the tests read
    /// the same files from memory and one byte per fill).
    #[inline(always)]
    fn kept(&self) -> &[u8] {
        self.src.rest(self.at)
    }

    /// A short tag whose first byte has `first & mask == want` and an indefinite length, both
    /// kept: `peek_tag` and `expect_indefinite_length` of them (the tag number). `None` reads
    /// nothing.
    #[inline(always)]
    fn kept_indefinite_item(&mut self, mask: u8, want: u8) -> Option<i32> {
        let first = match *self.kept() {
            [first, 0x80, ..] if first & mask == want && first & LONG_TAG != LONG_TAG => first,
            _ => return None,
        };
        self.at += 2;
        self.limits.push(self.limit);
        self.limit = 0;
        self.tag_len = 0;
        Some(i32::from(first & LONG_TAG))
    }

    /// An INTEGER tag (`expect_integer_tag`), a short length of 1 to `width` (at most 8) bytes
    /// (`read_short_length`) and that many content bytes, which `accept` takes (given the bytes
    /// as one big-endian number and their count), with at least 8 bytes kept after the length:
    /// that number and count, the bytes read. `None` reads nothing.
    #[inline(always)]
    fn kept_integer(
        &mut self,
        width: usize,
        accept: impl FnOnce(u64, usize) -> bool,
    ) -> Option<(u64, usize)> {
        let (n, length) = match *self.kept() {
            [TAG_INTEGER | TAG_APPLICATION_INTEGER, length, a, b, c, d, e, f, g, h, ..]
                if length != 0 && usize::from(length) <= width =>
            {
                // (The content bytes are the first `length` of the eight.)
                let eight = u64::from_be_bytes([a, b, c, d, e, f, g, h]);
                (eight >> (64 - 8 * u32::from(length)), usize::from(length))
            }
            _ => return None,
        };
        if !accept(n, length) {
            return None;
        }
        self.at += 2 + length;
        self.tag_len = 0;
        Some((n, length))
    }

    /// `ReadStdSigned` of a kept INTEGER of a signed type of `width` bytes whose content bytes
    /// are no more than the type's (`kept_integer`): the bytes, sign-extended.
    #[inline(always)]
    fn kept_signed(&mut self, width: usize) -> Option<i64> {
        let (n, length) = self.kept_integer(width, |_, _| true)?;
        let shift = 64 - 8 * length as u32;
        Some(((n << shift) as i64) >> shift)
    }

    /// `ThrowError(fFormatError or fOverflow, message)` at `site` at the current position
    /// (returned as the [`Fail`]: `return Err(self.throw(..))`, see [`throw_error`]).
    #[cold]
    #[inline(never)]
    fn throw(&self, site: ThrowSite, message: impl std::fmt::Display) -> Fail {
        throw_error(site, self.position(), message.to_string())
    }

    /// `PeekChar` of the input buffer, `i` bytes ahead: the end of the data throws
    /// `CEofException`. A byte the source keeps is read from its slice; the source reads more
    /// only past its end (`peek_after_fill`).
    #[inline(always)]
    fn peek(&mut self, i: usize) -> R<u8> {
        let at = self.at + i;
        match self.src.rest(at).first() {
            Some(&byte) => Ok(byte),
            None => self.peek_after_fill(at),
        }
    }

    /// [`peek`](Self::peek) of a byte that the source does not keep yet.
    #[cold]
    #[inline(never)]
    fn peek_after_fill(&mut self, at: usize) -> R<u8> {
        self.src.byte_after_fill(at).ok_or(Fail::Eof)
    }

    /// `GetChar` of the input buffer.
    #[inline(always)]
    fn byte(&mut self) -> R<u8> {
        let b = self.peek(0)?;
        self.at += 1;
        Ok(b)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:132-171
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:158-170,175-188
    /// `PeekTag`: the tag number of the tag at the current position; `tag_len` is its length.
    #[inline(always)]
    fn peek_tag(&mut self) -> R<i32> {
        let first = self.peek(0)?;
        if first & LONG_TAG != LONG_TAG {
            self.tag_len = 1;
            return Ok(i32::from(first & LONG_TAG));
        }
        self.peek_long_tag()
    }

    /// [`peek_tag`](Self::peek_tag) of a tag whose number follows its first byte.
    #[inline(never)]
    fn peek_long_tag(&mut self) -> R<i32> {
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
    /// `PeekTag` of the first tag byte, the class and the form: the tag number of a tag of that class and
    /// form.
    #[inline(always)]
    pub(crate) fn peek_tag_of_class(&mut self, class_and_constructed: u8) -> R<i32> {
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
    /// `ExpectTag` of the class, the form and the value.
    #[inline(always)]
    pub(crate) fn expect_tag(&mut self, class_and_constructed: u8, value: i32) -> R<()> {
        let tag = self.peek_tag_of_class(class_and_constructed)?;
        if tag != value {
            return Err(self.unexpected_tag_value(tag, value));
        }
        Ok(())
    }

    /// `UnexpectedTagValue` of the tag class, the tag and the tag value. The class was checked (universal or
    /// context specific): no "Application " or "Private " before the number.
    #[cold]
    #[inline(never)]
    fn unexpected_tag_value(&self, tag: i32, value: i32) -> Fail {
        self.throw(
            UNEXPECTED_TAG_VALUE,
            format_args!("Unexpected tag: {tag}, expected: {value}"),
        )
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:394-410
    #[inline(never)]
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
    /// `ExpectIndefiniteLength`: the tag that `peek_tag` found is followed by the length of a
    /// constructed item (a definite length sets the item's data limit).
    #[inline(always)]
    pub(crate) fn expect_indefinite_length(&mut self) -> R<()> {
        self.at += self.tag_len;
        let len = self.byte()?;
        self.limits.push(self.limit);
        self.limit = if len == 0x80 {
            0
        } else if len > 0x80 {
            // NCBI adds the stream position and the length that the length bytes give, and C++
            // does not say which operand is evaluated first. The compiler of the 2.17.0 binaries evaluates the
            // position first, so the limit is counted from the first length byte: the item ends
            // `n` bytes (the length of the length) before where its length says. makeblastdb 2.17.0
            // fails with "EndOfContentsByte expected" on such an item (checked on the
            // `long_len_*` cases of the mask fixtures).
            //
            // The sum is an `Int8` added to a `size_t`: computed in 64 unsigned bits, it wraps, and is stored in
            // the `Int8` `m_CurrentDataLimit` (objistrasnb.hpp:267). A limit of 0 is then an
            // indefinite length (`HaveMoreElements`, `ExpectEndOfContent`), a negative one ends the
            // members at once (checked with makeblastdb 2.17.0, the b_length_long_wrap_* fixtures).
            let base = self.pos() as i64;
            let long_length = self.read_length_long(len)?;
            (base as u64).wrapping_add(long_length) as i64
        } else {
            self.pos() as i64 + i64::from(len)
        };
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:174-182,242-250
    /// `HaveMoreElements`; the bytes before the next member or element are not read again
    /// (`AsnSource::release_before`).
    #[inline(always)]
    pub(crate) fn have_more_elements(&mut self) -> R<bool> {
        self.release();
        if self.limit != 0 {
            return Ok(self.limit > self.pos() as i64);
        }
        Ok(self.peek(0)? != 0)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:358-405
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:70-90
    /// `ExpectEndOfContent`: the end of a constructed item.
    #[inline(always)]
    pub(crate) fn expect_end_of_content(&mut self) -> R<()> {
        if self.limit == 0 && matches!(self.kept(), [0, 0, ..]) {
            self.at += 2;
            self.limit = self.limits.pop().unwrap_or(0);
            self.tag_len = 0;
            return Ok(());
        }
        self.expect_end_of_content_bytes()
    }

    /// [`expect_end_of_content`](Self::expect_end_of_content) one byte at a time.
    #[inline(never)]
    fn expect_end_of_content_bytes(&mut self) -> R<()> {
        if self.limit != 0 {
            if self.limit != self.pos() as i64 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
        } else {
            if self.peek(0)? != 0 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
            self.at += 1;
            if self.peek(0)? != 0 {
                return Err(self.throw(UNEXPECTED_CONTINUATION, "EndOfContentsByte expected"));
            }
            self.at += 1;
        }
        self.limit = self.limits.pop().unwrap_or(0);
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:309-329,338-352
    /// `ReadShortLength` after the tag.
    #[inline(always)]
    fn read_short_length(&mut self) -> R<usize> {
        self.at += self.tag_len;
        self.tag_len = 0;
        let byte = self.byte()?;
        if byte >= 0x80 {
            return Err(self.throw(UNEXPECTED_LONG_LENGTH, "ShortLength expected"));
        }
        Ok(usize::from(byte))
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:331-337
    /// `ReadLength` after the tag (`ReadLengthInlined`): a short or a long form.
    #[inline(always)]
    fn read_length(&mut self) -> R<u64> {
        self.at += self.tag_len;
        self.tag_len = 0;
        let byte = self.byte()?;
        if byte >= 0x80 {
            return self.read_length_long(byte);
        }
        Ok(u64::from(byte))
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
    #[inline(always)]
    pub(crate) fn expect_sys_tag_byte(&mut self, byte: u8) -> R<()> {
        let got = self.peek(0)?;
        if got != byte {
            return Err(self.unexpected_sys_tag_byte(got, byte));
        }
        self.tag_len = 1;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:145-160
    /// `ExpectIntegerTag`.
    #[inline(always)]
    pub(crate) fn expect_integer_tag(&mut self) -> R<()> {
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

    /// `ReadStdSigned` ([`read_std_signed`], the one port of it): `size` is the size of `T`; the result
    /// is the value of `T` widened to `i64`.
    pub(crate) fn read_std_signed(&mut self, size: usize) -> R<i64> {
        let length = self.read_short_length()?;
        let result = read_std_signed(length, size, || self.byte());
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
    /// `ReadStdUnsigned` of a `Uint4`.
    #[inline(always)]
    pub(crate) fn read_std_unsigned32(&mut self) -> R<u32> {
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

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1785-1789
    /// `SkipTagData`: the length after the tag, and that many bytes (`SkipBytes`).
    pub(crate) fn skip_tag_data(&mut self) -> R<()> {
        let length = self.read_length()?;
        self.skip_bytes(length)
    }

    /// `SkipBytes` of `count`: the end of the data inside them is a `CEofException`. (Values are
    /// skipped only inside a delayed member, whose bytes the source keeps anyway.)
    fn skip_bytes(&mut self, count: u64) -> R<()> {
        if count == 0 {
            return Ok(());
        }
        let last = usize::try_from(count)
            .ok()
            .and_then(|count| self.at.checked_add(count - 1));
        match last {
            Some(last) if self.src.byte_at(last).is_some() => {
                self.at += count as usize;
                Ok(())
            }
            _ => Err(Fail::Eof),
        }
    }

    /// `GetChars` of `count` bytes into `out`: the bytes of the source are appended as they come
    /// (`AsnSource::read_into`: a long value is not held twice); the end of the data before them is
    /// a `CEofException`.
    fn read_bytes(&mut self, count: u64, out: &mut Vec<u8>) -> R<()> {
        let read = self.src.read_into(self.at, count, out);
        self.at += read as usize;
        if read < count {
            return Err(Fail::Eof);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:281-314
    // (`SERIAL_READ_ANY_VISIBLESTRING_TAG` is 1 by default: a UTF8String is accepted for a
    // VisibleString, with a warning once per process; `SERIAL_READ_ANY_UTF8STRING_TAG` is true: a
    // VisibleString is accepted for a UTF8String.)
    /// `ExpectStringTag` of the string type, for a VisibleString (`visible`) or a UTF8String.
    fn expect_string_tag(&mut self, visible: bool, posts: &mut Vec<StringPost>) -> R<()> {
        let (tag, other) = if visible {
            (TAG_VISIBLE_STRING, TAG_UTF8_STRING)
        } else {
            (TAG_UTF8_STRING, TAG_VISIBLE_STRING)
        };
        let got = self.peek(0)?;
        if got == other {
            if visible {
                posts.push(StringPost::Utf8Tag);
            }
        } else if got != tag {
            return Err(self.unexpected_sys_tag_byte(got, tag));
        }
        self.tag_len = 1;
        Ok(())
    }

    /// The tag of a string of `kind` (`ExpectStringTag`, or the StringStore's `ExpectSysTag`).
    fn expect_string_kind(&mut self, kind: StringKind, posts: &mut Vec<StringPost>) -> R<()> {
        match kind {
            StringKind::Visible => self.expect_string_tag(true, posts),
            StringKind::Utf8 => self.expect_string_tag(false, posts),
            StringKind::Store => self.expect_sys_tag_byte(TAG_STRING_STORE),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1058-1067
    /// `UnexpectedMember` of the tag and the items, of a SEQUENCE or CHOICE whose members or variants are
    /// `items` (their automatic tags are 0, 1, ... in order).
    #[cold]
    #[inline(never)]
    fn unexpected_member(&self, tag: i32, items: &(impl Items + ?Sized)) -> Fail {
        let mut message = format!("unexpected member: [{tag}], should be one of: ");
        for i in 0..items.count() {
            message.push_str(&format!("{}[{i}] ", items.name(i)));
        }
        self.throw(UNEXPECTED_MEMBER, message)
    }

    /// The member or variant of tag number `tag` (from the `pos`-th item on), read after its
    /// length: an unknown one is an error at the byte after its length.
    #[inline(always)]
    fn member_index(&self, tag: i32, items: &(impl Items + ?Sized), pos: usize) -> R<usize> {
        match usize::try_from(tag) {
            Ok(index) if index >= pos && index < items.count() => Ok(index),
            _ => Err(self.unexpected_member(tag, items)),
        }
    }

    // The forms of the `Input` primitives below one byte at a time (see `kept`).

    /// `begin_class` and `begin_container`: a SEQUENCE (or SET) tag and its length.
    #[inline(never)]
    fn begin_constructed_bytes(&mut self, number: i32) -> R<()> {
        self.expect_tag(UNIVERSAL_CONSTRUCTED, number)?;
        self.expect_indefinite_length()
    }

    /// `begin_member`.
    #[inline(never)]
    fn begin_member_bytes(
        &mut self,
        items: &(impl Items + ?Sized),
        pos: usize,
    ) -> R<Option<usize>> {
        if !self.have_more_elements()? {
            return Ok(None);
        }
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        self.member_index(tag, items, pos).map(Some)
    }

    /// `begin_variant`.
    #[inline(never)]
    fn begin_variant_bytes(&mut self, items: &(impl Items + ?Sized)) -> R<usize> {
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        self.member_index(tag, items, 0)
    }

    /// `read_i32` and `read_i64`: a signed type of `width` bytes.
    #[inline(never)]
    fn read_signed_bytes(&mut self, width: usize) -> R<i64> {
        self.expect_integer_tag()?;
        self.read_std_signed(width)
    }

    /// `read_u32`.
    #[inline(never)]
    fn read_u32_bytes(&mut self) -> R<u32> {
        self.expect_integer_tag()?;
        self.read_std_unsigned32()
    }
}

impl<S: AsnSource> Input for BerReader<S> {
    const ENCODING: Encoding = Encoding::Binary;

    /// `GetPosition()` of the binary stream: the bytes read from the start of the object.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:118-121
    fn position(&self) -> StreamPosition {
        StreamPosition::Byte(self.pos() as u64)
    }

    /// The state of a new stream at the next byte: its positions count from there.
    fn new_object(&mut self) {
        self.base = self.at;
        self.tag_len = 0;
        self.limit = 0;
        self.limits.clear();
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:588-597
    #[cfg(test)]
    fn end_of_data(&mut self) -> bool {
        self.src.byte_at(self.at).is_none()
    }

    /// `ReadFileHeader` of the binary stream: `CObjectIStream`'s, which reads nothing (an empty
    /// input ends at the first tag: NCBI's `CEofException`).
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1029-1042
    // (A SEQUENCE: universal tag 16.)
    fn begin_class(&mut self) -> R<()> {
        if self.kept_indefinite_item(u8::MAX, TAG_SEQUENCE).is_some() {
            return Ok(());
        }
        self.begin_constructed_bytes(SEQUENCE_NUMBER)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1096-1179
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:462-473
    // (The length is read before the tag is looked up: the error of an unknown member is at the
    // byte after its length.)
    fn begin_member(&mut self, items: &(impl Items + ?Sized), pos: usize) -> R<Option<usize>> {
        // (In an item of an indefinite length, `have_more_elements` is the first byte of a
        // member's tag; a 0 ends the members.)
        if self.limit == 0 {
            self.release();
            if let Some(tag) =
                self.kept_indefinite_item(CLASS_AND_CONSTRUCTED_MASK, CONTEXT_CONSTRUCTED)
            {
                return self.member_index(tag, items, pos).map(Some);
            }
            if self.kept().first() == Some(&0) {
                return Ok(None);
            }
        }
        self.begin_member_bytes(items, pos)
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
    fn begin_variant(&mut self, items: &(impl Items + ?Sized)) -> R<usize> {
        match self.kept_indefinite_item(CLASS_AND_CONSTRUCTED_MASK, CONTEXT_CONSTRUCTED) {
            Some(tag) => self.member_index(tag, items, 0),
            None => self.begin_variant_bytes(items),
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
        let number = if set { SET_NUMBER } else { SEQUENCE_NUMBER };
        if self
            .kept_indefinite_item(u8::MAX, UNIVERSAL_CONSTRUCTED | number as u8)
            .is_some()
        {
            return Ok(());
        }
        self.begin_constructed_bytes(number)
    }

    fn begin_element(&mut self) -> R<bool> {
        self.have_more_elements()
    }

    fn end_container(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    #[inline]
    fn read_i32(&mut self) -> R<i32> {
        match self.kept_signed(4) {
            Some(n) => Ok(n as i32),
            None => self.read_signed_bytes(4).map(|n| n as i32),
        }
    }

    #[inline]
    fn read_i64(&mut self) -> R<i64> {
        match self.kept_signed(8) {
            Some(n) => Ok(n),
            None => self.read_signed_bytes(8),
        }
    }

    #[inline]
    fn read_u32(&mut self) -> R<u32> {
        // (`ReadStdUnsigned` takes a first byte of four with its high bit set as an overflow.)
        match self.kept_integer(4, |n, length| length < 4 || n >> 31 == 0) {
            Some((n, _)) => Ok(n as u32),
            None => self.read_u32_bytes(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2077-2087
    /// `SkipSNumber` and `SkipUNumber`: the tag and its bytes, whatever their value.
    fn skip_number(&mut self, _signed: bool) -> R<()> {
        self.expect_integer_tag()?;
        self.skip_tag_data()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:551-557
    /// `ReadBool`.
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
    /// `ReadNull`.
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

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:603-646
    /// `ReadDouble`: no content (0), one of the four special bytes, or a decimal number in text
    /// after its first byte (`NStr::StringToDoublePosix` with `fDecimalPosixFinite`, which stops at
    /// a 0 byte as at the end).
    fn read_real(&mut self) -> R<f64> {
        self.expect_sys_tag_byte(TAG_REAL)?;
        let length = self.read_length()?;
        if length < 2 {
            if length == 0 {
                return Ok(0.0);
            }
            return match self.byte()? {
                REAL_PLUS_INFINITY => Ok(f64::INFINITY),
                REAL_MINUS_INFINITY => Ok(f64::NEG_INFINITY),
                REAL_NOT_A_NUMBER => Ok(crate::ncbi_str::NCBI_NAN),
                REAL_MINUS_ZERO => Ok(-0.0),
                _ => Err(self.throw(
                    asnb_site(623, FORMAT_ERROR, Some(READ_DOUBLE)),
                    "Unrecognized REAL data",
                )),
            };
        }
        if length > MAX_DOUBLE_LENGTH {
            return Err(self.throw(
                asnb_site(626, FORMAT_ERROR, Some(READ_DOUBLE)),
                format_args!("too long REAL data: length > {MAX_DOUBLE_LENGTH}"),
            ));
        }
        let kind = self.byte()?;
        if kind & REAL_ENCODING_MASK != REAL_DECIMAL {
            return Err(self.throw(
                asnb_site(632, "CSerialException::eNotImplemented", Some(READ_DOUBLE)),
                format_args!("Unsupported encoding of REAL data: encoding = {kind}"),
            ));
        }
        let mut text = Vec::with_capacity(length as usize);
        self.read_bytes(length - 1, &mut text)?;
        // NCBI ends the text with a 0 byte after its length: the number ends at the first 0 byte.
        let end = text.iter().position(|&c| c == 0).unwrap_or(text.len());
        let read = string_to_double_posix(&text[..end], true);
        if read.end != end {
            return Err(self.throw(
                asnb_site(643, FORMAT_ERROR, Some(READ_DOUBLE)),
                "bad REAL data string",
            ));
        }
        Ok(read.value)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2079-2091
    /// `SkipFNumber`: the length checked, the bytes skipped.
    fn skip_real(&mut self) -> R<()> {
        self.expect_sys_tag_byte(TAG_REAL)?;
        let length = self.read_length()?;
        if length != 0 {
            if length > MAX_DOUBLE_LENGTH {
                return Err(self.throw(
                    SKIP_F_NUMBER,
                    format_args!("too long REAL data: length > {MAX_DOUBLE_LENGTH}"),
                ));
            }
            self.byte()?;
            self.skip_bytes(length - 1)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:785-826
    // (`ReadString` of a VisibleString fixes the string with `eFNP_ReplaceAndWarn`, the default
    // of `SERIAL_WRONG_CHARS_READ`; of a UTF8String, and `ReadStringStore`, with `eFNP_Allow`: a
    // StringStore keeps its bytes in binary, unlike the text stream.)
    /// A string (`kind`) appended to `out`, with what NCBI posts while it reads it
    /// in `posts` (also when it then fails): its tag, its length, its bytes (`ReadBytes`:
    /// `GetChars` of `str` and `count`, which reserves `count` bytes first when they are not all in the
    /// input buffer, and reads to the end of the data), and for a VisibleString each byte outside
    /// `' '..='~'` read as `#`.
    fn read_string(
        &mut self,
        kind: StringKind,
        out: &mut Vec<u8>,
        posts: &mut Vec<StringPost>,
    ) -> R<()> {
        self.expect_string_kind(kind, posts)?;
        let length = self.read_length()?;
        let start = out.len();
        if let Err(fail) = self.read_bytes(length, out) {
            // The end of the data inside the string: the reserve of its length came first.
            reserve(length, Reserve::String, self.names)?;
            return Err(fail);
        }
        if kind == StringKind::Visible {
            fix_visible_chars(&mut out[start..], self.position(), posts);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2093-2104
    /// `SkipString` of the string type / `SkipStringStore`: the tag, then the bytes skipped.
    fn skip_string(&mut self, kind: StringKind, posts: &mut Vec<StringPost>) -> R<()> {
        self.expect_string_kind(kind, posts)?;
        self.skip_tag_data()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1448-1490
    // NCBI reference (598d8ae6): c++/src/serial/stdtypes.cpp:1538-1548
    // (`BeginBytes` gives the block its length, so the `vector<char>` reader reserves it before
    // the bytes; the end of the data inside them is a `CEofException`.)
    /// An OCTET STRING (`BeginBytes`, the bytes, `EndBytes`); a BIT STRING is read as its bytes
    /// after the byte of its unused bits. The reserve of the length is the error when the data
    /// ends inside the bytes and NCBI's reserve throws (`asn_input::reserve`).
    fn read_octets(&mut self) -> R<Vec<u8>> {
        let length = match self.peek(0)? {
            TAG_OCTET_STRING => {
                self.tag_len = 1;
                self.read_length()?
            }
            TAG_BIT_STRING => {
                self.tag_len = 1;
                // The length less 1, in `size_t`: a length of 0 wraps to 2^64 - 1 bytes, which
                // pass the end of the data and every reserve.
                let length = self.read_length()?.wrapping_sub(1);
                self.byte()?;
                length
            }
            _ => {
                return Err(self.throw(BEGIN_BYTES, "Unable to identify the type of byte block"));
            }
        };
        let mut out = Vec::new();
        if let Err(fail) = self.read_bytes(length, &mut out) {
            reserve(length, Reserve::Bytes, self.names)?;
            return Err(fail);
        }
        Ok(out)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2124-2128
    /// `SkipByteBlock`: an OCTET STRING only.
    fn skip_octets(&mut self) -> R<()> {
        self.expect_sys_tag_byte(TAG_OCTET_STRING)?;
        self.skip_tag_data()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2023-2039
    /// `ReadEnum` of the enum's `values`: an INTEGER for a type with any number, else an ENUMERATED that must be
    /// one of the values (`FindName`).
    fn read_enum(&mut self, values: &(impl EnumItems + ?Sized)) -> R<i32> {
        if values.integer() {
            self.expect_integer_tag()?;
            return Ok(self.read_std_signed(4)? as i32);
        }
        self.expect_sys_tag_byte(TAG_ENUMERATED)?;
        let value = self.read_std_signed(4)? as i32;
        if !values.has_value(value) {
            return Err(invalid_enum(FIND_NAME, value.to_string()));
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

    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    /// `StartDelayBuffer` at the member's first byte: the source keeps the bytes from there.
    fn delay_start(&mut self) -> R<usize> {
        let start = self.offset();
        self.src.hold(Some(start));
        Ok(start)
    }

    fn delayed(&mut self, start: usize) -> Vec<u8> {
        let bytes = self.src.bytes(start..self.offset()).to_vec();
        self.src.hold(None);
        bytes
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:316-372
/// `TagToString` of `byte`: the class, `constructed/`, the universal type's name or the number,
/// and the byte.
pub(crate) fn tag_to_string(byte: u8) -> String {
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
/// `FixVisibleChars` of the string `s`, read up to `end`: each byte outside `' '..='~'`
/// becomes `#`, posted with the string as fixed so far (the bytes before it already `#`). The one
/// port of its posts: the binary object stream's and the header reader of a BLAST database's
/// (`seqid::ber::title_replaced_chars`).
pub(crate) fn fix_visible_chars(s: &mut [u8], end: StreamPosition, posts: &mut Vec<StringPost>) {
    for i in 0..s.len() {
        let byte = s[i];
        if !(b' '..=b'~').contains(&byte) {
            posts.push(StringPost::BadChar {
                position: end,
                byte,
                string: s.to_vec(),
            });
            s[i] = b'#';
        }
    }
}

/// Why `read_std_signed` stops: one of NCBI's two errors, or the byte source's own.
pub(crate) enum StdSignedError<E> {
    /// `fFormatError`, "zero length of number".
    ZeroLength,
    /// `fOverflow`, "overflow error", thrown at `line` of objistrasnb.cpp (488: the first byte
    /// beyond the width is not a sign; 493: a later one differs from it; 499: the value's sign
    /// differs from it).
    Overflow { line: u32 },
    /// The byte source failed (the data ends inside the value).
    Read(E),
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:475-512
/// NCBI's `ReadStdSigned` after its `ReadShortLength`: the INTEGER of `length` content bytes,
/// taken one at a time from `byte`, read into a signed type of `width` bytes (the size of `T`: 4 for
/// `Int4`, 8 for `Int8`, whose `Uint8` holds the same bits) and widened to `i64`. Content bytes
/// beyond the width are taken only as the sign extension of the value; anything else is an
/// overflow. The one port of the function: the binary object stream (`BerReader`) and the header
/// reader of a BLAST database (`seqid::ber::Reader::int`) read their INTEGERs with it.
pub(crate) fn read_std_signed<E>(
    mut length: usize,
    width: usize,
    mut byte: impl FnMut() -> Result<u8, E>,
) -> Result<i64, StdSignedError<E>> {
    if length == 0 {
        return Err(StdSignedError::ZeroLength);
    }
    let mut next = || byte().map_err(StdSignedError::Read);
    let mut n: i64;
    if length > width {
        length -= 1;
        let c = next()? as i8;
        if c != 0 && c != -1 {
            return Err(StdSignedError::Overflow { line: 488 });
        }
        while length > width {
            length -= 1;
            if next()? as i8 != c {
                return Err(StdSignedError::Overflow { line: 493 });
            }
        }
        length -= 1;
        n = i64::from(next()? as i8);
        if (n ^ i64::from(c)) & 0x80 != 0 {
            return Err(StdSignedError::Overflow { line: 499 });
        }
    } else {
        length -= 1;
        n = i64::from(next()? as i8);
    }
    while length > 0 {
        length -= 1;
        n = (n << 8) | i64::from(next()?);
    }
    Ok(n)
}
