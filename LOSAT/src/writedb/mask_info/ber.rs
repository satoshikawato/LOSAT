//! ASN.1 binary (`MSerial_AsnBinary`, BER): the primitives of `CObjectIStreamAsnBinary` on a byte
//! slice that holds the rest of the file.
//!
//! The types of a mask file are read with "automatic tags": a SEQUENCE is `30 <len>`, each of
//! its members is `[n]` (`a0 | n`, constructed, context specific) around the value, a CHOICE is
//! its variant `[n]` around the value, and every constructed item ends with `00 00` when
//! its length is indefinite (`80`); a definite length is accepted as well.

use super::reader::{fix_visible_chars, format_err, Fail, Input, R};
use super::schema::{ChoiceDef, ClassDef, EnumDef};

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
// ```c++
// enum ETagClass {
//     eUniversal       = 0 << 6,
//     eApplication     = 1 << 6,
//     eContextSpecific = 2 << 6,
//     ePrivate         = 3 << 6
// };
// enum ETagConstructed {
//     ePrimitive       = 0 << 5,
//     eConstructed     = 1 << 5
// };
// enum ETagValue {
//     eNone = 0, eBoolean = 1, eInteger = 2, ... eNull = 5, ... eEnumerated = 10, ...
//     eSequence = 16, eSet = 17, ... eVisibleString = 26, ... eUTF8String = 12 ...
//     eLongTag = 31
// };
// ```
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
const TAG_SEQUENCE: i32 = 16;
const TAG_SET: i32 = 17;
const LONG_TAG: u8 = 0x1f;

impl<'a> BerIn<'a> {
    pub(super) fn new(data: &'a [u8]) -> Self {
        BerIn {
            data,
            pos: 0,
            tag_len: 0,
            limit: 0,
            limits: Vec::new(),
        }
    }

    /// `m_Input.PeekChar(i)`: the end of the data throws `CEofException`.
    fn peek(&self, i: usize) -> R<u8> {
        self.data.get(self.pos + i).copied().ok_or(Fail::Eof)
    }

    /// `m_Input.GetChar()`
    fn byte(&mut self) -> R<u8> {
        let b = self.peek(0)?;
        self.pos += 1;
        Ok(b)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:132-155
    // ```c++
    // CObjectIStreamAsnBinary::TLongTag
    // CObjectIStreamAsnBinary::PeekLongTag(void)
    // {
    //     TLongTag tag = 0;
    //     size_t i = 1;
    //     const size_t KBitsInByte = 8; // ?
    //     const size_t KTagBits = sizeof(tag) * KBitsInByte - 1;
    //     TByte byte;
    //     do {
    //         if ( tag >= (1 << (KTagBits - 7)) ) {
    //             ThrowError(fOverflow,
    //                        "tag number is too big: "+NStr::IntToString(tag));
    //         }
    //         byte = PeekTagByte(i++);
    //         tag = (tag << 7) | (byte & 0x7f);
    //     } while ( (byte & 0x80) != 0 );
    //     ...
    //     m_CurrentTagLength = i;
    //     return tag;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:158-170
    // ```c++
    // inline
    // CObjectIStreamAsnBinary::TLongTag
    // CObjectIStreamAsnBinary::PeekTag(TByte first_tag_byte)
    // {
    //     TByte byte = StartTag(first_tag_byte);
    //     ETagValue sysTag = GetTagValue(byte);
    //     if ( sysTag == eLongTag ) {
    //         return PeekLongTag();
    //     }
    //     ...
    //     m_CurrentTagLength = 1;
    //     return sysTag;
    // }
    // ```
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
                return format_err(format!("tag number is too big: {tag}"));
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

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:138-156
    // ```c++
    // inline
    // void CObjectIStreamAsnBinary::ExpectTagClassByte(TByte first_tag_byte,
    //                                                  TByte expected_class_byte)
    // {
    //     if ( GetTagClassAndConstructed(first_tag_byte) != expected_class_byte ) {
    //         UnexpectedTagClassByte(first_tag_byte, expected_class_byte);
    //     }
    // }
    // ...
    // inline
    // CObjectIStreamAsnBinary::TLongTag
    // CObjectIStreamAsnBinary::PeekTag(TByte first_tag_byte,
    //                                  ETagClass tag_class,
    //                                  ETagConstructed tag_constructed)
    // {
    //     ExpectTagClassByte(first_tag_byte,
    //                        MakeTagClassAndConstructed(tag_class, tag_constructed));
    //     return PeekTag(first_tag_byte);
    // }
    // ```
    fn peek_tag_of_class(&mut self, class_and_constructed: u8) -> R<i32> {
        let first = self.peek(0)?;
        if first & CLASS_AND_CONSTRUCTED_MASK != class_and_constructed {
            return format_err(format!(
                "unexpected tag: {first}, should be {class_and_constructed}"
            ));
        }
        self.peek_tag()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:171-178
    // ```c++
    // void CObjectIStreamAsnBinary::UnexpectedTagClassByte(TByte first_tag_byte,
    //                                                      TByte expected_class_byte)
    // ...
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:158-185
    // ```c++
    // inline
    // void CObjectIStreamAsnBinary::ExpectTag(TByte first_tag_byte,
    //     ETagClass tag_class, ETagConstructed tag_constructed, TLongTag tag_value)
    // {
    //     ExpectTagClassByte(first_tag_byte,
    //                        MakeTagClassAndConstructed(tag_class, tag_constructed));
    //     TLongTag tag = PeekTag(first_tag_byte);
    //     if (tag != tag_value) {
    //         UnexpectedTagValue( tag_class, tag, tag_value);
    //     }
    // }
    // ```
    fn expect_tag(&mut self, class_and_constructed: u8, value: i32) -> R<()> {
        let tag = self.peek_tag_of_class(class_and_constructed)?;
        if tag != value {
            return format_err(format!("Unexpected tag: {tag}, expected: {value}"));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:237-256
    // ```c++
    // size_t CObjectIStreamAsnBinary::ReadLengthLong(TByte byte)
    // {
    //     size_t lengthLength = byte - 0x80;
    //     if ( lengthLength == 0 ) {
    //         ThrowError(fFormatError, "unexpected indefinite length");
    //     }
    //     if ( lengthLength > sizeof(size_t) ) {
    //         ThrowError(fOverflow, "length overflow");
    //     }
    //     size_t length = m_Input.GetChar() & 0xff;
    //     if ( length == 0 ) {
    //         ThrowError(fFormatError, "illegal length start");
    //     }
    //     while ( --lengthLength > 0 ) {
    //         length = (length << 8) | (m_Input.GetChar() & 0xff);
    //     }
    //     return StartTagData(length);
    // }
    // ```
    fn read_length_long(&mut self, byte: u8) -> R<u64> {
        let mut length_length = usize::from(byte - 0x80);
        if length_length == 0 {
            return format_err("unexpected indefinite length");
        }
        if length_length > 8 {
            return format_err("length overflow");
        }
        let mut length = u64::from(self.byte()?);
        if length == 0 {
            return format_err("illegal length start");
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

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:262-308
    // ```c++
    // inline
    // void CObjectIStreamAsnBinary::ExpectIndefiniteLength(void)
    // {
    //     ...
    // #if USE_DEF_LEN
    //     TByte len = FlushTag();
    //     m_DataLimits.push_back(m_CurrentDataLimit);
    //     if (len == eIndefiniteLengthByte) {
    //         m_CurrentDataLimit = 0;
    //     } else  if (len > eIndefiniteLengthByte) {
    //         m_CurrentDataLimit = m_Input.GetStreamPosAsInt8() + ReadLengthLong(len);
    //     ...
    //     } else {
    //         m_CurrentDataLimit = m_Input.GetStreamPosAsInt8() + len;
    //     }
    // ...
    //     m_CurrentTagLength = 0;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:236-246
    // ```c++
    // inline
    // CObjectIStreamAsnBinary::TByte
    // CObjectIStreamAsnBinary::FlushTag(void)
    // {
    //     ...
    //     m_Input.SkipChars(m_CurrentTagLength);
    //     return TByte(m_Input.GetChar());
    // }
    // ```
    /// The tag that `peek_tag` found is followed by the length of a constructed item.
    fn expect_indefinite_length(&mut self) -> R<()> {
        self.pos += self.tag_len;
        let len = self.byte()?;
        self.limits.push(self.limit);
        self.limit = if len == 0x80 {
            0
        } else if len > 0x80 {
            // `m_Input.GetStreamPosAsInt8() + ReadLengthLong(len)`: C++ does not say which
            // operand is evaluated first. The compiler of the 2.17.0 binaries evaluates the
            // position first, so the limit is counted from the first length byte: the item ends
            // `n` bytes (the length of the length) before where its length says. makeblastdb 2.17.0
            // fails with "EndOfContentsByte expected" on such an item (checked on the
            // `long_len_*` cases of the fixtures).
            let base = self.pos as i64;
            let long_length = self.read_length_long(len)?;
            base.saturating_add(i64::try_from(long_length).unwrap_or(i64::MAX))
        } else {
            self.pos as i64 + i64::from(len)
        };
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:174-182
    // ```c++
    // inline
    // bool CObjectIStreamAsnBinary::HaveMoreElements(void)
    // {
    // #if USE_DEF_LEN
    //     if (m_CurrentDataLimit != 0) {
    //         return m_CurrentDataLimit > m_Input.GetStreamPosAsInt8();
    //     }
    // #endif
    //     return PeekTagByte() != eEndOfContentsByte;
    // }
    // ```
    fn have_more_elements(&mut self) -> R<bool> {
        if self.limit != 0 {
            return Ok(self.limit > self.pos as i64);
        }
        Ok(self.peek(0)? != 0)
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:358-405
    // ```c++
    // inline
    // void CObjectIStreamAsnBinary::ExpectEndOfContent(void)
    // {
    //     ...
    // #if USE_DEF_LEN
    //     if (m_CurrentDataLimit != 0) {
    //         if (m_CurrentDataLimit != m_Input.GetStreamPosAsInt8()) {
    //             UnexpectedContinuation();
    //         }
    //     } else {
    //         if ( !m_Input.SkipExpectedChars(char(eEndOfContentsByte),
    //                                         char(eZeroLengthByte)) ) {
    //             UnexpectedContinuation();
    //         }
    //     }
    //     m_CurrentDataLimit = m_DataLimits.back();
    //     m_DataLimits.pop_back();
    // ```
    // NCBI reference (598d8ae6): c++/include/util/strbuffer.inl:70-90
    // ```c++
    // inline
    // bool CIStreamBuffer::SkipExpectedChar(char c, size_t offset)
    // {
    //     const char* pos = m_CurrentPos+offset;
    //     if ( pos >= m_DataEndPos )
    //         pos = FillBuffer(pos);
    //     if ( *pos != c )
    //         return false;
    //     m_CurrentPos = pos+1;
    //     return true;
    // }
    //
    // inline
    // bool CIStreamBuffer::SkipExpectedChars(char c1, char c2, size_t offset)
    // {
    //     return SkipExpectedChar(c1, offset) && SkipExpectedChar(c2, 0);
    // }
    // ```
    fn expect_end_of_content(&mut self) -> R<()> {
        if self.limit != 0 {
            if self.limit != self.pos as i64 {
                return format_err("EndOfContentsByte expected");
            }
        } else {
            if self.peek(0)? != 0 {
                return format_err("EndOfContentsByte expected");
            }
            self.pos += 1;
            if self.peek(0)? != 0 {
                return format_err("EndOfContentsByte expected");
            }
            self.pos += 1;
        }
        self.limit = self.limits.pop().unwrap_or(0);
        self.tag_len = 0;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:309-329
    // ```c++
    // inline
    // size_t CObjectIStreamAsnBinary::ReadShortLength(void)
    // {
    //     TByte byte = FlushTag();
    //     if ( byte >= 0x80 ) {
    //         UnexpectedLongLength();
    //     }
    //     return StartTagData(byte);
    // }
    //
    // inline
    // size_t CObjectIStreamAsnBinary::ReadLengthInlined(void)
    // {
    //     TByte byte = FlushTag();
    //     if ( byte >= 0x80 ) {
    //         return ReadLengthLong(byte);
    //     }
    //     return StartTagData(byte);
    // }
    // ```
    fn read_short_length(&mut self) -> R<usize> {
        self.pos += self.tag_len;
        self.tag_len = 0;
        let byte = self.byte()?;
        if byte >= 0x80 {
            return format_err("ShortLength expected");
        }
        Ok(usize::from(byte))
    }

    /// `ExpectSysTagByte(byte)` (objistrasnb.inl:90-100): the exact tag byte.
    fn expect_sys_tag_byte(&mut self, byte: u8) -> R<()> {
        let got = self.peek(0)?;
        if got != byte {
            return format_err(format!("unexpected tag: {got}, should be: {byte}"));
        }
        self.tag_len = 1;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/include/serial/impl/objistrasnb.inl:112-132
    // ```c++
    // inline
    // void CObjectIStreamAsnBinary::ExpectIntegerTag(void)
    // {
    //     if (m_SkipNextTag) {
    //         m_SkipNextTag = false;
    //         return;
    //     }
    //     Uint1 tag = StartTag(PeekTagByte());
    //     Uint1 exp = MakeTagByte(  eUniversal, ePrimitive, eInteger);
    //     if (tag != exp) {
    //         if (tag != MakeTagByte(eApplication, ePrimitive, eInteger)) {
    //             UnexpectedSysTagByte(exp);
    //         }
    //         SetSpecialCaseUsed(CObjectIStream::eReadAsBigInt);
    //     }
    //     m_CurrentTagLength = 1;
    // }
    // ```
    fn expect_integer_tag(&mut self) -> R<()> {
        match self.peek(0)? {
            TAG_INTEGER => {
                self.tag_len = 1;
                Ok(())
            }
            // The tag sets `eReadAsBigInt`, which `ReadWithDefaultMemberX` takes for "not set"
            // (member.cpp:807-812); LOSAT does not model it.
            TAG_APPLICATION_INTEGER => Err(Fail::Unsupported(
                "an INTEGER with the application tag (big integer)".into(),
            )),
            got => format_err(format!("unexpected tag: {got}, should be: {TAG_INTEGER}")),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:396-428
    // ```c++
    // template<typename T>
    // void ReadStdSigned(CObjectIStreamAsnBinary& in, T& data)
    // {
    //     size_t length = in.ReadShortLength();
    //     if ( length == 0 ) {
    //         in.ThrowError(in.fFormatError, "zero length of number");
    //     }
    //     T n;
    //     if ( length > sizeof(data) ) {
    //         // skip
    //         --length;
    //         Int1 c = in.ReadSByte();
    //         if ( c != 0 && c != -1 ) {
    //             in.ThrowError(in.fOverflow, "overflow error");
    //         }
    //         while ( length > sizeof(data) ) {
    //             --length;
    //             if ( in.ReadSByte() != c ) {
    //                 in.ThrowError(in.fOverflow, "overflow error");
    //             }
    //         }
    //         --length;
    //         n = in.ReadSByte();
    //         if ( ((n ^ c) & 0x80) != 0 ) {
    //             in.ThrowError(in.fOverflow, "overflow error");
    //         }
    //     }
    //     else {
    //         --length;
    //         n = in.ReadSByte();
    //     }
    //     while ( length > 0 ) {
    //         --length;
    //         n = (n << 8) | in.ReadByte();
    //     }
    //     data = n;
    //     in.EndOfTag();
    // }
    // ```
    /// `size` is `sizeof(T)`; the result is the value of `T` widened to `i64`.
    fn read_std_signed(&mut self, size: usize) -> R<i64> {
        let mut length = self.read_short_length()?;
        if length == 0 {
            return format_err("zero length of number");
        }
        let mut n: i64;
        if length > size {
            length -= 1;
            let c = self.byte()? as i8;
            if c != 0 && c != -1 {
                return format_err("overflow error");
            }
            while length > size {
                length -= 1;
                if self.byte()? as i8 != c {
                    return format_err("overflow error");
                }
            }
            length -= 1;
            n = i64::from(self.byte()? as i8);
            if ((n ^ i64::from(c)) & 0x80) != 0 {
                return format_err("overflow error");
            }
        } else {
            length -= 1;
            n = i64::from(self.byte()? as i8);
        }
        while length > 0 {
            length -= 1;
            n = (n << 8) | i64::from(self.byte()?);
        }
        Ok(if size == 4 { i64::from(n as i32) } else { n })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:430-466
    // ```c++
    // template<typename T>
    // void ReadStdUnsigned(CObjectIStreamAsnBinary& in, T& data)
    // {
    //     size_t length = in.ReadShortLength();
    //     if ( length == 0 ) {
    //         in.ThrowError(in.fFormatError, "zero length of number");
    //     }
    //     T n;
    //     if ( length > sizeof(data) ) {
    //         // skip
    //         while ( length > sizeof(data) ) {
    //             --length;
    //             if ( in.ReadSByte() != 0 ) {
    //                 in.ThrowError(in.fOverflow, "overflow error");
    //             }
    //         }
    //         --length;
    //         n = in.ReadByte();
    //     }
    //     else if ( length == sizeof(data) ) {
    //         --length;
    //         n = in.ReadByte();
    //         if ( (n & 0x80) != 0 ) {
    //             in.ThrowError(in.fOverflow, "overflow error");
    //         }
    //     }
    //     else {
    //         n = 0;
    //     }
    //     while ( length > 0 ) {
    //         --length;
    //         n = (n << 8) | in.ReadByte();
    //     }
    //     data = n;
    //     in.EndOfTag();
    // }
    // ```
    fn read_std_unsigned32(&mut self) -> R<u32> {
        let mut length = self.read_short_length()?;
        if length == 0 {
            return format_err("zero length of number");
        }
        let mut n: u32;
        if length > 4 {
            while length > 4 {
                length -= 1;
                if self.byte()? != 0 {
                    return format_err("overflow error");
                }
            }
            length -= 1;
            n = u32::from(self.byte()?);
        } else if length == 4 {
            length -= 1;
            n = u32::from(self.byte()?);
            if n & 0x80 != 0 {
                return format_err("overflow error");
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
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        Ok(None)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1051-1065
    // ```c++
    // void CObjectIStreamAsnBinary::BeginClass(const CClassTypeInfo* classType)
    // {
    //     ...
    //     bool need_eoc = !m_SkipNextTag;
    //     if (!m_SkipNextTag) {
    //         ExpectTag(classType->GetTagClass(), eConstructed, classType->GetTag());
    //         ExpectIndefiniteLength();
    //     }
    //     m_SkipNextTag = classType->IsTagImplicit();
    //     TopFrame().SetNoEOC(!need_eoc);
    // ```
    // (A SEQUENCE: universal tag 16.)
    fn begin_class(&mut self) -> R<()> {
        self.expect_tag(UNIVERSAL_CONSTRUCTED, TAG_SEQUENCE)?;
        self.expect_indefinite_length()
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1096-1159
    // ```c++
    // TMemberIndex
    // CObjectIStreamAsnBinary::BeginClassMember(const CClassTypeInfo* classType,
    //                                           TMemberIndex pos)
    // {
    // #if USE_DEF_LEN
    //     if (!HaveMoreElements()) {
    //         return kInvalidMember;
    //     }
    //     TByte first_tag_byte = PeekTagByte();
    // ...
    //     TLongTag tag = PeekTag(first_tag_byte, eContextSpecific, eConstructed);
    //     ExpectIndefiniteLength();
    //     TMemberIndex index = classType->GetMembers().Find(tag, eContextSpecific, pos);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownMembers()) {
    //             ...
    //         }
    //         else {
    //             UnexpectedMember(tag, classType->GetItems());
    //         }
    //     }
    //     return index;
    // }
    //
    // void CObjectIStreamAsnBinary::EndClassMember(void)
    // {
    //     m_SkipNextTag = false;
    //     if (!TopFrame().GetNoEOC()) {
    //         ExpectEndOfContent();
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:462-473
    // ```c++
    //         TMemberIndex index = tag + zero_index;
    //         if ( index < pos || index > LastIndex() )
    //             return kInvalidMember;
    //         return index;
    // ```
    // (The classes of a mask file have the "automatic" tag type: member `n` has the tag `[n]`.)
    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        if !self.have_more_elements()? {
            return Ok(None);
        }
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        match usize::try_from(tag) {
            Ok(index) if index >= pos && index < def.members.len() => Ok(Some(index)),
            _ => format_err(format!("unexpected member: [{tag}]")),
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

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1306-1332
    // ```c++
    //     TLongTag tag = PeekTag(PeekTagByte(), eContextSpecific, eConstructed);
    //     ExpectIndefiniteLength();
    //     TMemberIndex index = choiceType->GetVariants().Find(tag, eContextSpecific);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownVariants()) {
    //             SetFailFlags(fUnknownValue);
    //         } else {
    //             UnexpectedMember(tag, choiceType->GetItems());
    //         }
    //         return index;
    //     }
    // ```
    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize> {
        let tag = self.peek_tag_of_class(CONTEXT_CONSTRUCTED)?;
        self.expect_indefinite_length()?;
        match usize::try_from(tag) {
            Ok(index) if index < def.variants.len() => Ok(index),
            _ => format_err(format!("unexpected member: [{tag}]")),
        }
    }

    fn end_variant(&mut self) -> R<()> {
        self.expect_end_of_content()
    }

    fn end_choice(&mut self) -> R<()> {
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:956-1000
    // ```c++
    // void CObjectIStreamAsnBinary::BeginContainer(const CContainerTypeInfo* cType)
    // {
    //     ...
    //     bool need_eoc = !m_SkipNextTag;
    //     if (!m_SkipNextTag) {
    //         ExpectTag(cType->GetTagClass(), eConstructed, cType->GetTag());
    //         ExpectIndefiniteLength();
    //     }
    //     m_SkipNextTag = cType->IsTagImplicit();
    //     TopFrame().SetNoEOC(!need_eoc);
    // ```
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

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:468-475
    // ```c++
    // bool CObjectIStreamAsnBinary::ReadBool(void)
    // {
    //     ExpectSysTag(eBoolean);
    //     ExpectShortLength(1);
    //     bool ret = ReadByte() != 0;
    //     EndOfTag();
    //     return ret;
    // }
    // ```
    fn read_bool(&mut self) -> R<bool> {
        self.expect_sys_tag_byte(TAG_BOOLEAN)?;
        let length = self.read_short_length()?;
        if length != 1 {
            return format_err(format!("unexpected length: {length}, should be: 1"));
        }
        Ok(self.byte()? != 0)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:1485-1491
    // ```c++
    // void CObjectIStreamAsnBinary::ReadNull(void)
    // {
    //     ExpectSysTag(eNull);
    //     ExpectShortLength(0);
    //     EndOfTag();
    // }
    // ```
    fn read_null(&mut self) -> R<()> {
        self.expect_sys_tag_byte(TAG_NULL)?;
        let length = self.read_short_length()?;
        if length != 0 {
            return format_err(format!("unexpected length: {length}, should be: 0"));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:218-251,820-835
    // ```c++
    // void CObjectIStreamAsnBinary::ExpectStringTag(EStringType type)
    // {
    //     ...
    //     else {
    //         static CSafeStatic<NCBI_PARAM_TYPE(SERIAL, READ_ANY_VISIBLESTRING_TAG)> sx_ReadAny;
    //         if ( sx_ReadAny->Get() ) {
    //             // may be eUTF8String
    //             TByte alt_tag = MakeTagByte(eUniversal,ePrimitive,eUTF8String);
    //             if ( PeekTagByte() == alt_tag ) {
    //                 // optionally issue a warning
    //                 ...
    //                 ExpectSysTagByte(alt_tag);
    //                 return;
    //             }
    //         }
    //         ExpectSysTag(eUniversal, ePrimitive, eVisibleString);
    //     }
    // }
    // ...
    // void CObjectIStreamAsnBinary::ReadString(string& s, EStringType type)
    // {
    //     ExpectStringTag(type);
    //     ReadStringValue(ReadLength(), s,
    //                     type == eStringTypeVisible? x_FixCharsMethod(): eFNP_Allow);
    // }
    // ```
    // (`SERIAL_READ_ANY_VISIBLESTRING_TAG` is 1 by default: a UTF8String is accepted for a
    // VisibleString, with a warning.)
    fn read_string(&mut self) -> R<Vec<u8>> {
        match self.peek(0)? {
            TAG_VISIBLE_STRING | TAG_UTF8_STRING => self.tag_len = 1,
            got => {
                return format_err(format!(
                    "unexpected tag: {got}, should be: {TAG_VISIBLE_STRING}"
                ))
            }
        }
        self.pos += self.tag_len;
        self.tag_len = 0;
        let byte = self.byte()?;
        let length = if byte >= 0x80 {
            self.read_length_long(byte)?
        } else {
            u64::from(byte)
        };
        // `m_Input.GetChars(str, count)` reads `count` bytes or throws at the end of the file;
        // the length is not trusted for an allocation.
        let remaining = (self.data.len() - self.pos) as u64;
        if length > remaining {
            return Err(Fail::Eof);
        }
        let length = length as usize;
        let mut s = self.data[self.pos..self.pos + length].to_vec();
        self.pos += length;
        fix_visible_chars(&mut s);
        Ok(s)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasnb.cpp:2023-2039
    // ```c++
    // TEnumValueType
    // CObjectIStreamAsnBinary::ReadEnum(const CEnumeratedTypeValues& values)
    // {
    //     TEnumValueType value;
    //     if ( values.IsInteger() ) {
    //         // allow any integer
    //         ExpectIntegerTag();
    //         ReadStdSigned(*this, value);
    //     }
    //     else {
    //         // enum element by value
    //         ExpectSysTag(eEnumerated);
    //         ReadStdSigned(*this, value);
    //         values.FindName(value, false); // check value
    //     }
    //     return value;
    // }
    // ```
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        if def.integer {
            self.expect_integer_tag()?;
            return Ok(self.read_std_signed(4)? as i32);
        }
        self.expect_sys_tag_byte(TAG_ENUMERATED)?;
        let value = self.read_std_signed(4)? as i32;
        if !def.values.iter().any(|&(_, v)| v == value) {
            return format_err(format!("invalid value of enumerated type: {value}"));
        }
        Ok(value)
    }

    fn offset(&self) -> usize {
        self.pos
    }
}
