//! NCBI's ASN.1 binary object stream (`CObjectOStreamAsnBinary`, outfmt 9) for the calls of
//! `Encoder`: BER with definite lengths for the primitive values, indefinite lengths for the
//! constructed ones, and an explicit context tag around every SEQUENCE member and CHOICE
//! variant.
//!
//! A document is the bytes of one top-level Seq-annot; documents follow each other with
//! nothing between them (`CObjectOStreamAsnBinary` has no separator, and `MSerial_AsnBinary`
//! writes no newline).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:382-386, 418-431

use super::encoder::{Encoder, Replaced, SerialError, Stream, StreamPosition};
use super::real::double_to_string_posix;

/// `DBL_DIG`: the precision of `WriteDouble` (`WriteDouble2` of `data` with `DBL_DIG`).
const DBL_DIG: u32 = 15;

/// The identifier octets of the primitive types and of the constructed types that NCBI
/// writes (`CObjectStreamAsnBinaryDefs::ETagValue`, include/serial/impl/objstrasnb.hpp, with
/// `eConstructed` = 0x20 and `eContextSpecific` = 0x80), which NCBI's binary reader expects
/// (`report::serial::asn_types::read_ber`).
pub(crate) const TAG_BOOLEAN: u8 = 0x01;
pub(crate) const TAG_INTEGER: u8 = 0x02;
pub(crate) const TAG_BIT_STRING: u8 = 0x03;
pub(crate) const TAG_REAL: u8 = 0x09;
pub(crate) const TAG_ENUMERATED: u8 = 0x0a;
pub(crate) const TAG_OCTET_STRING: u8 = 0x04;
pub(crate) const TAG_UTF8_STRING: u8 = 0x0c;
pub(crate) const TAG_VISIBLE_STRING: u8 = 0x1a;
pub(crate) const TAG_STRING_STORE: u8 = 0x41;
/// An INTEGER with the application class (`MakeTagByte` of `eApplication`, `ePrimitive` and
/// `eInteger`), which NCBI's reader takes as an INTEGER (`ExpectIntegerTag`).
pub(crate) const TAG_APPLICATION_INTEGER: u8 = 0x42;
pub(crate) const TAG_SEQUENCE: u8 = 0x30;
pub(crate) const TAG_SET: u8 = 0x31;
pub(crate) const TAG_NULL: u8 = 0x05;
/// A constructed context-specific tag (`[i]`, i < 31): `0xa0 | i`.
pub(crate) const TAG_CONTEXT: u8 = 0xa0;
/// `eTagClassMask`, `eTagConstructedMask` and `eTagValueMask` of an identifier octet.
pub(crate) const TAG_CLASS_MASK: u8 = 0xc0;
pub(crate) const TAG_CONSTRUCTED: u8 = 0x20;
pub(crate) const TAG_VALUE_MASK: u8 = 0x1f;
/// `eLongTag`: the tag number follows in base-128 bytes.
pub(crate) const TAG_LONG: u8 = 0x1f;
/// `eIndefiniteLengthByte`.
pub(crate) const INDEFINITE_LENGTH: u8 = 0x80;
/// `eEndOfContentsByte` (and `eZeroLengthByte`): the two bytes `00 00` that close an
/// indefinite length.
pub(crate) const END_OF_CONTENTS: u8 = 0x00;
/// The first content octet of a REAL written in decimal form (`eDecimal`, the NR form is
/// not selected, CXX-5248).
pub(crate) const REAL_DECIMAL: u8 = 0x00;
/// `eDecimalEncoding`: the bits of a REAL's first content octet that name its form.
pub(crate) const REAL_ENCODING_MASK: u8 = 0xc0;
/// The one-byte REALs (`ESpecialRealValues`): `ePositiveInfinity`, `eNegativeInfinity`,
/// `eNotANumber`, `eNegativeZero`.
pub(crate) const REAL_PLUS_INFINITY: u8 = 0x40;
pub(crate) const REAL_MINUS_INFINITY: u8 = 0x41;
pub(crate) const REAL_NOT_A_NUMBER: u8 = 0x42;
pub(crate) const REAL_MINUS_ZERO: u8 = 0x43;

/// The ASN.1 binary encoder, appending one document to `out`.
pub struct AsnBinaryEncoder<'b> {
    out: &'b mut Vec<u8>,
    error: Option<SerialError>,
    /// The bytes that `visible_str_replacing` replaced (`replace`).
    replaced: Vec<Replaced>,
}

impl<'b> AsnBinaryEncoder<'b> {
    /// An encoder that appends a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        Self {
            out,
            error: None,
            replaced: Vec::new(),
        }
    }

    /// `ReplaceVisibleChar` of the byte `byte`, at its offset in the buffer (`Replaced`;
    /// `encoder::Traced` makes it the offset in the document); `write_string` then writes `#`.
    #[cold]
    #[inline(never)]
    fn replace(&mut self, byte: u8) {
        self.replaced.push(Replaced {
            byte,
            position: StreamPosition::Byte(self.out.len() as u64),
            with_string: true,
            at: self.out.len(),
        });
    }

    /// The bytes of `WriteString` of a VisibleString: a byte outside 0x20..0x7E is written as
    /// `#`, the error `NonVisibleChar` without `REPLACE`, recorded at its offset with it
    /// (`replace`).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:841-865
    #[inline(always)]
    fn write_string<const REPLACE: bool>(&mut self, value: &[u8]) {
        self.out.push(TAG_VISIBLE_STRING);
        self.length(value.len());
        for &byte in value {
            if !Stream::AsnBinary.replaces(byte) {
                self.out.push(byte);
            } else {
                if REPLACE {
                    self.replace(byte);
                } else if self.error.is_none() {
                    self.error = Some(SerialError::NonVisibleChar(byte));
                }
                self.out.push(b'#');
            }
        }
    }

    /// `WriteLength`: one byte up to 127, else `0x80 + count` and the `count` bytes of the
    /// length, big-endian.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:393-416
    fn length(&mut self, length: usize) {
        if length <= 127 {
            self.out.push(length as u8);
            return;
        }
        let bytes = (length as u64).to_be_bytes();
        let count = if length <= 0xff {
            1
        } else if length <= 0xffff {
            2
        } else if length <= 0xff_ffff {
            3
        } else {
            // The `sizeof` of `length`: eight bytes on x86-64 Linux, not the shortest four.
            8
        };
        self.out.push(0x80 + count as u8);
        self.out.extend_from_slice(&bytes[8 - count..]);
    }

    /// `WriteNumberValue` of an `Int4`: the length (1 to 4) and the two's-complement bytes, the
    /// shortest of the ranges of one, two and three bytes, else four.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:533-554
    fn number(&mut self, value: i32) {
        let length = if (-0x80..=0x7f).contains(&value) {
            1
        } else if (-0x8000..=0x7fff).contains(&value) {
            2
        } else if (-0x80_0000..=0x7f_ffff).contains(&value) {
            3
        } else {
            4
        };
        self.out.push(length as u8);
        self.out
            .extend_from_slice(&value.to_be_bytes()[4 - length..]);
    }

    /// `BeginClassMember` / `BeginChoiceVariant`: the explicit context tag `[tag]` with an
    /// indefinite length.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1346-1375
    /// (`[i]`, an explicit tag; the tags of Seq-annot's modules are below 31, so the identifier is
    /// one byte.)
    fn context_open(&mut self, tag: u8) {
        debug_assert!(tag < 31, "a one-byte context tag");
        self.out.push(TAG_CONTEXT | tag);
        self.out.push(INDEFINITE_LENGTH);
    }

    /// `WriteEndOfContent`: the two zero bytes that close an indefinite-length value.
    fn end_of_content(&mut self) {
        self.out.extend_from_slice(&[0, 0]);
    }
}

impl Encoder for AsnBinaryEncoder<'_> {
    const STREAM: Option<super::encoder::Stream> = Some(super::encoder::Stream::AsnBinary);

    /// `WriteFileHeader` writes nothing in the binary stream (the type name `Seq-annot ::=`
    /// belongs to the text stream).
    fn begin_document(&mut self, _type_name: &str) {}

    /// No separator follows a binary document (`SetAutoSeparator` is not set).
    fn end_document(&mut self) {}

    /// `BeginClass`: `30 80`.
    fn begin_sequence(&mut self, _type_name: &str) {
        self.out
            .extend_from_slice(&[TAG_SEQUENCE, INDEFINITE_LENGTH]);
    }

    /// `EndClass`: `00 00`.
    fn end_sequence(&mut self) {
        self.end_of_content();
    }

    fn member(&mut self, _name: &str, tag: u8) {
        self.context_open(tag);
    }

    fn end_member(&mut self) {
        self.end_of_content();
    }

    /// `BeginContainer`: `30 80` for a SEQUENCE OF, `31 80` for a SET OF
    /// (`MakeContainerTagByte`).
    fn begin_list(&mut self, set: bool) {
        self.out
            .extend_from_slice(&[if set { TAG_SET } else { TAG_SEQUENCE }, INDEFINITE_LENGTH]);
    }

    /// `BeginContainerElement` writes nothing: the element writes its own tag.
    fn element(&mut self) {}

    fn end_list(&mut self) {
        self.end_of_content();
    }

    fn variant(&mut self, _type_name: &str, _name: &str, tag: u8) {
        self.context_open(tag);
    }

    fn end_variant(&mut self) {
        self.end_of_content();
    }

    /// `WriteInt4`: INTEGER.
    fn int(&mut self, value: i32) {
        self.out.push(TAG_INTEGER);
        self.number(value);
    }

    /// `WriteInt8`: INTEGER in the fewest big-endian two's-complement bytes of the value.
    fn big_int(&mut self, value: i64) {
        let bytes = value.to_be_bytes();
        let mut start = 0;
        while start < 7
            && ((bytes[start] == 0x00 && bytes[start + 1] & 0x80 == 0)
                || (bytes[start] == 0xff && bytes[start + 1] & 0x80 != 0))
        {
            start += 1;
        }
        self.out.push(TAG_INTEGER);
        self.out.push((8 - start) as u8);
        self.out.extend_from_slice(&bytes[start..]);
    }

    /// `WriteEnum`: ENUMERATED with the number of the value.
    fn enumerated(&mut self, _name: &str, value: i32) {
        self.out.push(TAG_ENUMERATED);
        self.number(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1006-1017
    fn named_int(&mut self, _name: &str, value: i32) {
        self.out.push(TAG_INTEGER);
        self.number(value);
    }

    /// `WriteBool`: `01 01 01` or `01 01 00`.
    fn boolean(&mut self, value: bool) {
        self.out
            .extend_from_slice(&[TAG_BOOLEAN, 1, u8::from(value)]);
    }

    /// `WriteString` of a VisibleString: the tag, the length and the bytes. A byte outside
    /// 0x20..0x7E is replaced by NCBI (`eFNP_ReplaceAndWarn`, as in the text stream): here an
    /// error (`write_string`).
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string::<false>(value);
    }

    /// `visible_str` that records each replaced byte (`replace`) instead of the error.
    fn visible_str_replacing(&mut self, value: &[u8]) {
        self.write_string::<true>(value);
    }

    /// An OCTET STRING: `04`, the length and the bytes.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1704-1708
    fn octets(&mut self, value: &[u8]) {
        self.out.push(TAG_OCTET_STRING);
        self.length(value.len());
        self.out.extend_from_slice(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:874-879
    fn string_store(&mut self, value: &[u8]) {
        self.out.push(TAG_STRING_STORE);
        self.length(value.len());
        self.out.extend_from_slice(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:433-437
    fn null(&mut self) {
        self.out.extend_from_slice(&[TAG_NULL, 0]);
    }

    /// `WriteString` of a UTF8String: the tag of its parameter `SERIAL_WRITE_UTF8STRING_TAG`,
    /// whose default (false, the one LOSAT writes) is a VisibleString's, the length and the
    /// bytes as they are (no byte is checked: only a VisibleString's are). NCBI BLAST+ 2.17.0's
    /// `psiblast -in_pssm d3a_c10_user -outfmt 11` with `ARCHIVE_FORMAT=bin` writes the
    /// `User-field` string `é\x01` as `1a 03 c3 a9 01` (probe XB of line RG-F3 PR 2, 2026-11-23).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:281-314
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:831-873
    fn utf8_str(&mut self, value: &[u8]) {
        self.out.push(TAG_VISIBLE_STRING);
        self.length(value.len());
        self.out.extend_from_slice(value);
    }

    /// NCBI parses a delay buffer read from text before it writes it in binary
    /// (`CMemberInfoFunctions::WriteLongMember`, member.cpp:1084-1119): no BLAST document that
    /// LOSAT writes in ASN.1 binary holds one.
    fn delayed(&mut self, _text: &[u8]) {
        if self.error.is_none() {
            self.error = Some(SerialError::Unwritten("a member read into a delay buffer"));
        }
    }

    fn buffered(&self) -> usize {
        self.out.len()
    }

    fn flush_to(&mut self, sink: &mut dyn std::io::Write, _all: bool) -> std::io::Result<()> {
        sink.write_all(self.out)?;
        self.out.clear();
        Ok(())
    }

    /// `WriteDouble2` of `data` with `DBL_DIG`: `09`, the length, the decimal form byte 0 and the
    /// digits of `DoubleToStringPosix`; a NaN or an infinity is its one special content octet
    /// (`09 01 42`, `40`, `41`; NCBI's branch for a negative zero is compiled out, so it takes the
    /// decimal form). NCBI BLAST+ 2.17.0's `psiblast -in_pssm` of a binary PSSM with
    /// `PLUS-INFINITY` in its frequency ratios, or `MINUS-INFINITY` or `NOT-A-NUMBER` as its
    /// lambda, `-outfmt 11` with `ARCHIVE_FORMAT=bin` (probe real of line RG-F3 PR 2,
    /// 2026-11-23).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:738-816
    fn real(&mut self, value: f64) {
        if !value.is_finite() {
            let special = if value.is_nan() {
                REAL_NOT_A_NUMBER
            } else if value > 0.0 {
                REAL_PLUS_INFINITY
            } else {
                REAL_MINUS_INFINITY
            };
            self.out.extend_from_slice(&[TAG_REAL, 1, special]);
            return;
        }
        let mut digits = Vec::with_capacity(24);
        double_to_string_posix(value, DBL_DIG, &mut digits);
        self.out.push(TAG_REAL);
        self.length(digits.len() + 1);
        self.out.push(REAL_DECIMAL);
        self.out.extend_from_slice(&digits);
    }

    /// The whole buffer (`flush_to` hands it all on, with what preceded the document).
    fn held(&self) -> &[u8] {
        self.out
    }

    fn hand_on(&mut self, sink: &mut dyn std::io::Write, count: usize) -> std::io::Result<()> {
        sink.write_all(&self.out[..count])?;
        self.out.drain(..count);
        Ok(())
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }

    fn take_replaced(&mut self) -> Vec<Replaced> {
        std::mem::take(&mut self.replaced)
    }
}
