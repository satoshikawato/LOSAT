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

use super::encoder::{Encoder, Replaced, SerialError};
use super::real::double_to_string_posix;
use super::stream_buffer::{StreamBuffer, ASN_BINARY_BACK_LIMIT};
use super::visible::{good_visible_char, Position};

/// `DBL_DIG`: the precision of `WriteDouble` (`WriteDouble2(data, DBL_DIG)`).
const DBL_DIG: u32 = 15;

/// The identifier octets of the primitive types and of the constructed types that NCBI
/// writes (`CObjectStreamAsnBinaryDefs::ETagValue`, include/serial/impl/objstrasnb.hpp, with
/// `eConstructed` = 0x20 and `eContextSpecific` = 0x80).
const TAG_BOOLEAN: u8 = 0x01;
const TAG_INTEGER: u8 = 0x02;
const TAG_REAL: u8 = 0x09;
const TAG_ENUMERATED: u8 = 0x0a;
const TAG_OCTET_STRING: u8 = 0x04;
const TAG_VISIBLE_STRING: u8 = 0x1a;
const TAG_STRING_STORE: u8 = 0x41;
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
pub struct AsnBinaryEncoder<'b, const MODEL: bool = false> {
    out: &'b mut Vec<u8>,
    /// Where the document starts in `out`.
    start: usize,
    error: Option<SerialError>,
    /// NCBI's form: the bytes written as `#` (`replace`).
    replaced: Vec<Replaced>,
    /// NCBI's form: NCBI's output buffer (back limit 0), and the bytes of the document handed to
    /// a sink (`flush_to`, `flush_through`).
    ncbi: StreamBuffer,
    drained: usize,
}

impl<'b> AsnBinaryEncoder<'b, false> {
    /// The plain encoder, appending a document to `out`.
    pub fn new(out: &'b mut Vec<u8>) -> Self {
        Self::with_out(out)
    }
}

impl<'b> AsnBinaryEncoder<'b, true> {
    /// NCBI's form of the encoder, appending a document to `out` (`encoder::Traced`): it writes
    /// `#` for a VisibleString byte outside 0x20..0x7E and records it at NCBI's stream position,
    /// with the bytes that NCBI's buffer has handed out then.
    pub fn ncbi(out: &'b mut Vec<u8>) -> Self {
        Self::with_out(out)
    }
}

impl<'b, const MODEL: bool> AsnBinaryEncoder<'b, MODEL> {
    fn with_out(out: &'b mut Vec<u8>) -> Self {
        let start = out.len();
        Self {
            out,
            start,
            error: None,
            replaced: Vec::new(),
            ncbi: StreamBuffer::new(ASN_BINARY_BACK_LIMIT),
            drained: 0,
        }
    }

    /// `WriteByte` (`m_Output.PutChar`, objostrasnb.cpp:161-225): one byte, one reservation in
    /// NCBI's form.
    #[inline(always)]
    fn byte(&mut self, byte: u8) {
        if MODEL {
            self.ncbi.skip(1);
        }
        self.out.push(byte);
    }

    /// Bytes that NCBI writes one `WriteByte` each (`WriteBytesOf`, a tag and its length).
    #[inline(always)]
    fn bytes_of(&mut self, bytes: &[u8]) {
        if MODEL {
            for _ in bytes {
                self.ncbi.skip(1);
            }
        }
        self.out.extend_from_slice(bytes);
    }

    /// `WriteBytes(bytes, size)` (`m_Output.PutString`, objostrasnb.cpp:230-247).
    #[inline(always)]
    fn put_bytes(&mut self, bytes: &[u8]) {
        if MODEL {
            self.ncbi.put_string(bytes.len());
        }
        self.out.extend_from_slice(bytes);
    }

    /// `GetStreamPos()`: the bytes of the document written so far.
    fn stream_pos(&self) -> usize {
        self.drained + (self.out.len() - self.start)
    }

    /// `ReplaceVisibleChar` of the byte `byte`, at its offset in the buffer (`Replaced`;
    /// `encoder::Traced` makes it the offset in the document); `write_string` then writes `#`.
    #[cold]
    #[inline(never)]
    fn replace(&mut self, index: usize) {
        self.replaced.push(Replaced {
            index,
            position: Position::Byte(self.stream_pos()),
            after: self.ncbi.flushed(),
            with_string: true,
        });
    }

    /// The bytes of `WriteString` of a VisibleString: a byte outside 0x20..0x7E is written as
    /// `#`, the error `NonVisibleChar` without `REPLACE`, recorded at its offset with it
    /// (`replace`).
    ///
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:841-865
    /// ```c++
    ///     WriteLength(length - fixed);
    ///     if ( type == eStringTypeVisible && x_FixCharsMethod() != eFNP_Allow ) {
    ///         size_t done = 0;
    ///         for ( size_t i = 0; i < length; ++i ) {
    ///             char c = str[i];
    ///             if ( !GoodVisibleChar(c) ) {
    ///                 ...
    ///                 if ( i > done ) {
    ///                     WriteBytes(str.data() + done, i - done);
    ///                 }
    ///                 c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, str, x_FixCharsSubst());
    ///                 if (c != 0) {
    ///                     WriteByte(c);
    ///                 }
    ///                 done = i + 1;
    ///             }
    ///         }
    /// ```
    #[inline(always)]
    fn write_string<const REPLACE: bool>(&mut self, value: &[u8]) {
        self.byte(TAG_VISIBLE_STRING);
        self.length(value.len());
        if MODEL {
            // NCBI's calls: `WriteBytes` of the good bytes before each bad one, the message, then
            // `WriteByte` of `#`; `WriteBytes` of the rest.
            let mut done = 0;
            for (index, &byte) in value.iter().enumerate() {
                if !good_visible_char(byte) {
                    if index > done {
                        self.put_bytes(&value[done..index]);
                    }
                    if REPLACE {
                        self.replace(index);
                    } else if self.error.is_none() {
                        self.error = Some(SerialError::NonVisibleChar(byte));
                    }
                    self.byte(b'#');
                    done = index + 1;
                }
            }
            if done < value.len() {
                self.put_bytes(&value[done..]);
            }
            return;
        }
        // The plain form (`REPLACE` is NCBI's form's only): each byte, `#` for a bad one, which is
        // the error `NonVisibleChar`.
        debug_assert!(!REPLACE);
        for &byte in value {
            if good_visible_char(byte) {
                self.out.push(byte);
            } else {
                if self.error.is_none() {
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
            self.byte(length as u8);
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
        self.byte(0x80 + count as u8);
        self.bytes_of(&bytes[8 - count..]);
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
        self.byte(length as u8);
        self.bytes_of(&value.to_be_bytes()[4 - length..]);
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
        self.byte(TAG_CONTEXT | tag);
        self.byte(INDEFINITE_LENGTH);
    }

    /// `WriteEndOfContent`: the two zero bytes that close an indefinite-length value.
    fn end_of_content(&mut self) {
        self.bytes_of(&[0, 0]);
    }
}

impl<const MODEL: bool> Encoder for AsnBinaryEncoder<'_, MODEL> {
    /// `WriteFileHeader` writes nothing in the binary stream (the type name `Seq-annot ::=`
    /// belongs to the text stream).
    fn begin_document(&mut self, _type_name: &'static str) {}

    /// No separator follows a binary document (`SetAutoSeparator` is not set).
    fn end_document(&mut self) {}

    /// `BeginClass`: `30 80`.
    fn begin_sequence(&mut self, _type_name: &str) {
        self.bytes_of(&[TAG_SEQUENCE, INDEFINITE_LENGTH]);
    }

    /// `EndClass`: `00 00`.
    fn end_sequence(&mut self) {
        self.end_of_content();
    }

    fn member(&mut self, _name: &'static str, tag: u8) {
        self.context_open(tag);
    }

    fn end_member(&mut self) {
        self.end_of_content();
    }

    /// `BeginContainer`: `30 80` for a SEQUENCE OF, `31 80` for a SET OF
    /// (`MakeContainerTagByte`).
    fn begin_list(&mut self, set: bool) {
        self.bytes_of(&[if set { TAG_SET } else { TAG_SEQUENCE }, INDEFINITE_LENGTH]);
    }

    /// `BeginContainerElement` writes nothing: the element writes its own tag.
    fn element(&mut self) {}

    fn end_list(&mut self) {
        self.end_of_content();
    }

    fn variant(&mut self, _type_name: &'static str, _name: &'static str, tag: u8) {
        self.context_open(tag);
    }

    fn end_variant(&mut self) {
        self.end_of_content();
    }

    /// `WriteInt4`: INTEGER.
    fn int(&mut self, value: i32) {
        self.byte(TAG_INTEGER);
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
        self.byte(TAG_INTEGER);
        self.byte((8 - start) as u8);
        self.bytes_of(&bytes[start..]);
    }

    /// `WriteEnum`: ENUMERATED with the number of the value.
    fn enumerated(&mut self, _name: &str, value: i32) {
        self.byte(TAG_ENUMERATED);
        self.number(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1006-1017
    /// ```c++
    /// void CObjectOStreamAsnBinary::WriteEnum(const CEnumeratedTypeValues& values,
    ///                                         TEnumValueType value)
    /// {
    ///     if ( values.IsInteger() ) {
    ///         WriteSysTag(eInteger);
    ///     }
    ///     ...
    ///     WriteNumberValue(value);
    /// }
    /// ```
    fn named_int(&mut self, _name: &str, value: i32) {
        self.byte(TAG_INTEGER);
        self.number(value);
    }

    /// `WriteBool`: `01 01 01` or `01 01 00`.
    fn boolean(&mut self, value: bool) {
        self.bytes_of(&[TAG_BOOLEAN, 1, u8::from(value)]);
    }

    /// `WriteString` of a VisibleString: the tag, the length and the bytes. A byte outside
    /// 0x20..0x7E is replaced by NCBI (`eFNP_ReplaceAndWarn`, as in the text stream): here an
    /// error (`write_string`).
    fn visible_str(&mut self, value: &[u8]) {
        self.write_string::<false>(value);
    }

    /// `visible_str` that records each replaced byte (`replace`) instead of the error.
    fn visible_str_replacing(&mut self, value: &[u8]) {
        // Only NCBI's form replaces (`encoder::Traced` writes through it).
        if MODEL {
            self.write_string::<true>(value);
        } else {
            self.write_string::<false>(value);
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
    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1704-1708
    /// ```c++
    /// void CObjectOStreamAsnBinary::BeginBytes(const ByteBlock& block)
    /// {
    ///     WriteSysTag(eOctetString);
    ///     WriteLength(block.GetLength());
    /// }
    /// ```
    fn octets(&mut self, value: &[u8]) {
        self.byte(TAG_OCTET_STRING);
        self.length(value.len());
        self.put_bytes(value);
    }

    /// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:874-879
    /// ```c++
    /// void CObjectOStreamAsnBinary::WriteStringStore(const string& str)
    /// {
    ///     WriteShortTag(eApplication, ePrimitive, eStringStore);
    ///     size_t length = str.size();
    ///     WriteLength(length);
    ///     WriteBytes(str.data(), length);
    /// }
    /// ```
    fn string_store(&mut self, value: &[u8]) {
        self.byte(TAG_STRING_STORE);
        self.length(value.len());
        self.put_bytes(value);
    }

    fn buffered(&self) -> usize {
        self.out.len() - self.start
    }

    /// The plain form hands out every byte (the binary stream looks back at nothing); NCBI's form
    /// what NCBI's buffer has handed out.
    fn flush_to(&mut self, sink: &mut dyn std::io::Write, all: bool) -> std::io::Result<()> {
        let through = if MODEL && !all {
            self.ncbi.flushed()
        } else {
            self.stream_pos()
        };
        self.flush_through(sink, through)
    }

    fn flush_keeping(&mut self, sink: &mut dyn std::io::Write, keep: usize) -> std::io::Result<()> {
        if MODEL {
            return self.flush_to(sink, false);
        }
        let through = self.stream_pos().saturating_sub(keep);
        self.flush_through(sink, through)
    }

    fn flush_through(
        &mut self,
        sink: &mut dyn std::io::Write,
        through: usize,
    ) -> std::io::Result<()> {
        if through <= self.drained {
            return Ok(());
        }
        let end = self.start + (through - self.drained);
        sink.write_all(&self.out[self.start..end])?;
        self.out.drain(self.start..end);
        self.drained = through;
        Ok(())
    }

    fn real(&mut self, value: f64) {
        if !value.is_finite() {
            if self.error.is_none() {
                self.error = Some(SerialError::NonFiniteReal);
            }
            return;
        }
        let mut digits = Vec::with_capacity(24);
        double_to_string_posix(value, DBL_DIG, &mut digits);
        self.byte(TAG_REAL);
        self.length(digits.len() + 1);
        self.byte(REAL_DECIMAL);
        self.put_bytes(&digits);
    }

    fn error(&self) -> Option<SerialError> {
        self.error
    }

    fn take_replaced(&mut self) -> Vec<Replaced> {
        std::mem::take(&mut self.replaced)
    }
}
