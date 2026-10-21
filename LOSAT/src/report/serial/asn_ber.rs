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
//! ```c++
//! inline
//! void CObjectOStreamAsnBinary::WriteIndefiniteLength(void)
//! {
//!     WriteByte(eIndefiniteLengthByte);
//! }
//! ...
//! inline
//! void CObjectOStreamAsnBinary::WriteLength(size_t length)
//! {
//!     if ( length <= 127 )
//!         WriteShortLength(length);
//!     else
//!         WriteLongLength(length);
//! }
//!
//! inline
//! void CObjectOStreamAsnBinary::WriteEndOfContent(void)
//! {
//!     WriteSysTag(eNone);
//!     WriteShortLength(0);
//! }
//! ```

use super::encoder::{Encoder, SerialError};
use super::real::double_to_string_posix;

/// `DBL_DIG`: the precision of `WriteDouble` (`WriteDouble2(data, DBL_DIG)`).
const DBL_DIG: u32 = 15;

/// The identifier octets of the primitive types and of the constructed types that NCBI
/// writes (`CObjectStreamAsnBinaryDefs::ETagValue`, include/serial/impl/objstrasnb.hpp, with
/// `eConstructed` = 0x20 and `eContextSpecific` = 0x80).
const TAG_BOOLEAN: u8 = 0x01;
const TAG_INTEGER: u8 = 0x02;
const TAG_REAL: u8 = 0x09;
const TAG_ENUMERATED: u8 = 0x0a;
const TAG_VISIBLE_STRING: u8 = 0x1a;
const TAG_SEQUENCE: u8 = 0x30;
const TAG_SET: u8 = 0x31;
/// A constructed context-specific tag (`[i]`, i < 31): `0xa0 | i`.
const TAG_CONTEXT: u8 = 0xa0;
/// `eIndefiniteLengthByte`.
const INDEFINITE_LENGTH: u8 = 0x80;
/// The first content octet of a REAL written in decimal form (`eDecimal`, the NR form is
/// not selected, CXX-5248).
const REAL_DECIMAL: u8 = 0x00;

/// The ASN.1 binary encoder, appending one document to `out`.
pub struct AsnBinaryEncoder<'b> {
    out: &'b mut Vec<u8>,
    error: Option<SerialError>,
}

impl<'b> AsnBinaryEncoder<'b> {
    /// An encoder that appends a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        Self { out, error: None }
    }

    /// `WriteLength`: one byte up to 127, else `0x80 + count` and the `count` bytes of the
    /// length, big-endian.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:393-416
    /// ```c++
    /// void CObjectOStreamAsnBinary::WriteLongLength(size_t length)
    /// {
    ///     // long form
    ///     size_t count;
    ///     if ( length <= 0xffU )
    ///         count = 1;
    ///     else if ( length <= 0xffffU )
    ///         count = 2;
    ///     else if ( length <= 0xffffffU )
    ///         count = 3;
    ///     else {
    ///         count = sizeof(length);
    ///         ...
    ///     }
    ///     WriteByte(TByte(0x80 + count));
    ///     WriteBytesOf(length, count);
    /// }
    /// ```
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
            // `sizeof(length)`: eight bytes on x86-64 Linux, not the shortest four.
            8
        };
        self.out.push(0x80 + count as u8);
        self.out.extend_from_slice(&bytes[8 - count..]);
    }

    /// `WriteNumberValue(Int4)`: the length (1 to 4) and the two's-complement bytes, the
    /// shortest of the ranges of one, two and three bytes, else four.
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:533-554
    /// ```c++
    /// void CObjectOStreamAsnBinary::WriteNumberValue(Int4 data)
    /// {
    ///     size_t length;
    ///     if ( data >= Int4(-0x80) && data <= Int4(0x7f) ) {
    ///         // one byte
    ///         length = 1;
    ///     }
    ///     else if ( data >= Int4(-0x8000) && data <= Int4(0x7fff) ) {
    ///         // two bytes
    ///         length = 2;
    ///     }
    ///     else if ( data >= Int4(-0x800000) && data <= Int4(0x7fffff) ) {
    ///         // three bytes
    ///         length = 3;
    ///     }
    ///     else {
    ///         // full length signed
    ///         length = sizeof(data);
    ///     }
    ///     WriteShortLength(length);
    ///     WriteBytesOf(data, length);
    /// }
    /// ```
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
    /// ```c++
    /// void CObjectOStreamAsnBinary::BeginClassMember(const CMemberId& id)
    /// {
    ///     ...
    ///     WriteTag(eContextSpecific, eConstructed, id.GetTag());  // [i], explicit
    ///     WriteIndefiniteLength();
    /// }
    /// ```
    /// (the tags of Seq-annot's modules are below 31, so the identifier is one byte.)
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

    fn variant(&mut self, _name: &str, tag: u8) {
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

    /// `WriteBool`: `01 01 01` or `01 01 00`.
    fn boolean(&mut self, value: bool) {
        self.out
            .extend_from_slice(&[TAG_BOOLEAN, 1, u8::from(value)]);
    }

    /// `WriteString` of a VisibleString: the tag, the length and the bytes. A byte outside
    /// 0x20..0x7E is replaced by NCBI (`eFNP_ReplaceAndWarn`, as in the text stream); LOSAT
    /// refuses it (`NonVisibleChar`).
    fn visible_str(&mut self, value: &[u8]) {
        self.out.push(TAG_VISIBLE_STRING);
        self.length(value.len());
        for &byte in value {
            if (0x20..=0x7e).contains(&byte) {
                self.out.push(byte);
            } else {
                if self.error.is_none() {
                    self.error = Some(SerialError::NonVisibleChar(byte));
                }
                self.out.push(b'#');
            }
        }
    }

    /// `WriteDouble2(data, DBL_DIG)`: `09`, the length, the decimal form byte 0 and the
    /// digits of `DoubleToStringPosix`. (NCBI marks a NaN or an infinity by another form
    /// byte, which LOSAT does not write: BLAST reports have none.)
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:797-814
    /// ```c++
    ///     WriteLength(width + 1);
    ///     if (width) {
    ///         WriteByte(type);
    ///         WriteBytes(buffer, width);
    ///     } else {
    ///         WriteByte(type);
    ///     }
    /// ```
    fn real(&mut self, value: f64) {
        if !value.is_finite() {
            if self.error.is_none() {
                self.error = Some(SerialError::NonFiniteReal);
            }
            return;
        }
        let mut digits = Vec::with_capacity(24);
        double_to_string_posix(value, DBL_DIG, &mut digits);
        self.out.push(TAG_REAL);
        self.length(digits.len() + 1);
        self.out.push(REAL_DECIMAL);
        self.out.extend_from_slice(&digits);
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }
}
