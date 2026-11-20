//! The steps of NCBI's ASN.1 text object stream (`CObjectIStreamAsn`) that the walk of the type
//! table (`read`) calls, over the lexical layer of `report::serial::asn_text_reader` (`Reader`):
//! the file header `Type ::=`, the blocks `{ ... }` and their elements, the members by name or by
//! tag in brackets, INTEGER, ENUMERATED, REAL, strings, octets, and the skip functions.

use super::read::{Cx, FrameKind, Front, ReadError, R};
use super::{EnumValue, Member, Storage, StringKind, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_text_reader::{real_of_parts, AsnSource, MemberId, Reader};
use crate::report::serial::encoder::StreamPosition;

/// The text front end: the lexical reader of one text (or of a stream of texts).
pub(crate) struct Text<S: AsnSource> {
    pub(crate) r: Reader<S>,
}

impl<S: AsnSource> Text<S> {
    pub(crate) fn new(r: Reader<S>) -> Self {
        Text { r }
    }
}

// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:570-585
// ```c++
// void CIStreamBuffer::BadNumber(void)
//     THROWS1((CUtilException))
// {
//     m_Error = "bad number";
//     NCBI_THROW_FMT(CUtilException, eWrongData,
//                    "bad number in line " << GetLine());
// }
//
//
// void CIStreamBuffer::NumberOverflow(void)
//     THROWS1((CUtilException))
// {
//     m_Error = "number overflow";
//     NCBI_THROW_FMT(CUtilException, eWrongData,
//                    "number overflow in line " << GetLine());
// }
// ```
const BAD_NUMBER: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 574,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::BadNumber"),
};
const NUMBER_OVERFLOW: ThrowSite = ThrowSite {
    file: "src/util/strbuffer.cpp",
    line: 583,
    kind: "CUtilException::eWrongData",
    function: Some("ncbi::CIStreamBuffer::NumberOverflow"),
};

impl<S: AsnSource> Front for Text<S> {
    // The lines are those of `END_OBJECT_FRAME()` in the functions that push the frames:
    // `CObjectIStream::Read` (objistr.cpp:959-969), `CObjectIStreamAsn::ReadClassSequential`
    // (objistrasn.cpp:1380-1397, the member frame ends in `ReadClassSequentialContentsEnd` at
    // 1393), `CObjectIStream::ReadChoiceSimple` (objistr.cpp:1482-1501),
    // `CObjectIStreamAsn::ReadContainer` (objistrasn.cpp:1245-1274), and the skip functions
    // `SkipClassSequential` (objistrasn.cpp:1417-1433), `SkipChoiceSimple` (objistr.cpp:1503-1522)
    // and `SkipContainer` (objistrasn.cpp:1276-1292).
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
        const ASN: &str = "src/serial/objistrasn.cpp";
        const ISTR: &str = "src/serial/objistr.cpp";
        match (kind, skip) {
            (FrameKind::Named, _) => (ISTR, 968, "ncbi::CObjectIStream::Read"),
            (FrameKind::Class, false) => {
                (ASN, 1396, "ncbi::CObjectIStreamAsn::ReadClassSequential")
            }
            (FrameKind::ClassMember, false) => {
                (ASN, 1393, "ncbi::CObjectIStreamAsn::ReadClassSequential")
            }
            (FrameKind::Choice, false) => (ISTR, 1500, "ncbi::CObjectIStream::ReadChoiceSimple"),
            (FrameKind::ChoiceVariant, false) => {
                (ISTR, 1498, "ncbi::CObjectIStream::ReadChoiceSimple")
            }
            (FrameKind::Array, false) => (ASN, 1273, "ncbi::CObjectIStreamAsn::ReadContainer"),
            (FrameKind::ArrayElement, false) => {
                (ASN, 1270, "ncbi::CObjectIStreamAsn::ReadContainer")
            }
            (FrameKind::Class, true) => (ASN, 1432, "ncbi::CObjectIStreamAsn::SkipClassSequential"),
            (FrameKind::ClassMember, true) => {
                (ASN, 1429, "ncbi::CObjectIStreamAsn::SkipClassSequential")
            }
            (FrameKind::Choice, true) => (ISTR, 1521, "ncbi::CObjectIStream::SkipChoiceSimple"),
            (FrameKind::ChoiceVariant, true) => {
                (ISTR, 1519, "ncbi::CObjectIStream::SkipChoiceSimple")
            }
            (FrameKind::Array, true) => (ASN, 1291, "ncbi::CObjectIStreamAsn::SkipContainer"),
            (FrameKind::ArrayElement, true) => {
                (ASN, 1288, "ncbi::CObjectIStreamAsn::SkipContainer")
            }
        }
    }

    fn start_object(&mut self) {
        self.r.hold(None);
        self.r.release();
    }

    fn file_header(&mut self, _cx: &mut Cx<'_>) -> R<Vec<u8>> {
        let id = self.r.file_header()?;
        Ok(self.r.bytes(id).to_vec())
    }

    /// `StartBlock`.
    fn begin_class(&mut self, _ty: u16, _cx: &mut Cx<'_>) -> R<()> {
        self.r.start_block()?;
        Ok(())
    }

    /// `NextElement`, then the member named or tagged from `pos`.
    fn next_member(
        &mut self,
        members: &'static [Member],
        pos: usize,
        block: &mut bool,
        _cx: &mut Cx<'_>,
    ) -> R<Option<usize>> {
        if !self.r.next_element(block)? {
            return Ok(None);
        }
        self.item_index(members, pos, true).map(Some)
    }

    fn end_member(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }

    /// `EndBlock`.
    fn end_class(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(self.r.end_block()?)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1475 (`BeginChoiceVariant`,
    // quoted at `Reader::choice`: `NextElement` at the choice's start, then the variant's id)
    fn begin_variant(
        &mut self,
        _ty: u16,
        variants: &'static [Member],
        _cx: &mut Cx<'_>,
    ) -> R<usize> {
        let mut block_start = true;
        if !self.r.next_element(&mut block_start)? {
            return Err(ReadError::Serial("choice variant id expected".into()));
        }
        self.item_index(variants, 0, false)
    }

    fn end_variant(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(())
    }

    /// `StartBlock`.
    fn begin_list(&mut self, _set: bool, _cx: &mut Cx<'_>) -> R<()> {
        self.r.start_block()?;
        Ok(())
    }

    /// `NextElement`.
    fn next_element(&mut self, _element: u16, block: &mut bool, _cx: &mut Cx<'_>) -> R<bool> {
        Ok(self.r.next_element(block)?)
    }

    /// `EndBlock`.
    fn end_list(&mut self, _cx: &mut Cx<'_>) -> R<()> {
        Ok(self.r.end_block()?)
    }

    fn integer(
        &mut self,
        _ty: u16,
        storage: Storage,
        skip: bool,
        _cx: &mut Cx<'_>,
    ) -> R<Option<i64>> {
        self.integer_number(storage, skip)
    }

    fn enumerated(
        &mut self,
        _ty: u16,
        values: &'static [EnumValue],
        integer: bool,
        _skip: bool,
        _cx: &mut Cx<'_>,
    ) -> R<i32> {
        self.enum_value(values, integer)
    }

    fn boolean(&mut self, _ty: u16, _skip: bool, _cx: &mut Cx<'_>) -> R<bool> {
        Ok(self.r.boolean()?)
    }

    fn null(&mut self, _ty: u16, _skip: bool, _cx: &mut Cx<'_>) -> R<()> {
        Ok(self.r.null()?)
    }

    fn real(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<Option<f64>> {
        self.real_number(skip)
    }

    fn string(&mut self, _ty: u16, kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        self.string_value(kind, skip, cx)
    }

    fn octets(&mut self, _ty: u16, skip: bool, _cx: &mut Cx<'_>) -> R<Option<Vec<u8>>> {
        self.octet_string(skip)
    }

    /// The delayed member: `StartDelayBuffer` (white space skipped), then the value skipped
    /// (`Open::Delayed`), the text kept (`Value::Delayed`).
    // NCBI reference (598d8ae6): c++/src/serial/member.cpp:848-869
    // ```c++
    // void CMemberInfoFunctions::ReadLongMember(CObjectIStream& in,
    //                                           const CMemberInfo* memberInfo,
    //                                           TObjectPtr classPtr)
    // {
    //     if ( memberInfo->CanBeDelayed() ) {
    //         CDelayBuffer& buffer = memberInfo->GetDelayBuffer(classPtr);
    //         if ( !buffer ) {
    //             if (!in.ShouldParseDelayBuffer()) {
    //                 memberInfo->UpdateSetFlagYes(classPtr);
    //                 in.StartDelayBuffer();
    //                 memberInfo->GetTypeInfo()->SkipData(in);
    //                 in.EndDelayBuffer(buffer, memberInfo, classPtr);
    //                 return;
    //             }
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1648-1652
    // ```c++
    // void CObjectIStreamAsn::StartDelayBuffer(void)
    // {
    //     SkipWhiteSpace();
    //     CObjectIStream::StartDelayBuffer();
    // }
    // ```
    // (psiblast and tblastn set no hook on their stream: `ShouldParseDelayBuffer` is false,
    // objistr.cpp:745-760.)
    // (A stream keeps the text from its start until the value ends: `Reader::hold`.)
    fn delay_start(&mut self, _cx: &mut Cx<'_>) -> R<usize> {
        self.r.skip_white_space()?;
        self.r.hold(Some(self.r.pos));
        Ok(self.r.pos)
    }

    fn delay_end(&mut self, start: usize, _ty: u16, _cx: &mut Cx<'_>) -> R<Value> {
        let text = self.r.bytes(start..self.r.pos).to_vec();
        self.r.hold(None);
        Ok(Value::Delayed(text))
    }
}

impl<S: AsnSource> Text<S> {
    /// `ReadMemberId(SkipWhiteSpace())` and `GetMemberIndex(classType, id, pos)` (from `pos`,
    /// `in_class`) or `GetChoiceIndex(choiceType, id)`: the member or variant named or tagged.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:407-424
    // ```c++
    //     TMemberIndex idx;
    //     if (!id.empty()  &&  isdigit((unsigned char) id[0])) {
    //         idx = classType->GetMembers().Find
    //             (CMemberId::TTag(NStr::StringToInt(id)), CAsnBinaryDefs::eContextSpecific, pos);
    //     }
    //     else {
    //         idx = classType->GetMembers().Find(id, pos);
    //         if (idx == kInvalidMember) {
    //             idx = GetAltItemIndex(classType,id,pos);
    //         }
    //     }
    //     return idx;
    // ```
    // `GetAltItemIndex` (366-387) tries the name with its first letter in upper case and takes it
    // only for an item without a prefix, which no member of an ASN.1 module is (its names start
    // with a lower-case letter): an upper-case first letter is never a member here (c21 of
    // SURVEY_D3: `Biomol`).
    fn item_index(&mut self, items: &'static [Member], pos: usize, in_class: bool) -> R<usize> {
        let found = match self.r.member_or_tag()? {
            MemberId::Name(range) => {
                let id = self.r.bytes(range);
                if id.is_empty() && !in_class {
                    return Err(ReadError::Serial("choice variant id expected".into()));
                }
                items[pos..]
                    .iter()
                    .position(|item| item.name().as_bytes() == id)
                    .map(|at| at + pos)
            }
            MemberId::Tag(range) => {
                // NStr::StringToInt(id) with the default flags: its `CStringException` leaves
                // the reader (NCBI 2.17.0: `[3x]`, `[99999999999]`, probe e03, e07).
                let tag = crate::ncbi_str::string_to_int(self.r.bytes(range), 0)
                    .map_err(ReadError::Exception)?;
                items[pos..]
                    .iter()
                    .position(|item| i32::from(item.tag) == tag)
                    .map(|at| at + pos)
            }
        };
        found.ok_or_else(|| ReadError::Serial("unexpected member or variant".into()))
    }

    /// `CIStreamBuffer::SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
    fn get_after_ws(&mut self) -> R<u8> {
        loop {
            let c = self.r.peek_char(0)?;
            self.r.pos += 1;
            if !crate::ncbi_str::is_space(c) {
                return Ok(c);
            }
        }
    }

    fn util_error(&self, site: ThrowSite, what: &str) -> ReadError {
        ReadError::Exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.r.line),
        ))
    }

    /// `GetInt4`, `GetUint4` (`signed` false: no `-`) and `GetInt8` after `SkipWhiteSpace`, the
    /// digits checked against `max` as NCBI checks them.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-669
    // ```c++
    // Int4 CIStreamBuffer::GetInt4(void)
    //     THROWS1((CIOException,CUtilException))
    // {
    //     bool sign;
    //     char c = SkipWs();
    //     switch ( c ) {
    //     case '-':
    //         sign = true;
    //         c = GetChar();
    //         break;
    //     case '+':
    //         sign = false;
    //         c = GetChar();
    //         break;
    //     default:
    //         sign = false;
    //         break;
    //     }
    //
    //     Uint4 n = c - '0';
    //     if ( n > 9 )
    //         BadNumber();
    //
    //     // overflow limits
    //     const Uint4 kMaxBeforeMul = kMax_I4/10;
    //     const Uint1 kMaxLimitAdd = Uint1(kMax_I4%10 + sign);
    //
    //     for ( ;; ) {
    //         Uint1 d = (Uint1)(PeekCharNoEOF() - '0');
    //         if ( d > 9 )
    //             break;
    //         SkipChar();
    //
    //         // check multiplication overflow
    //         if ( n > kMaxBeforeMul || (n == kMaxBeforeMul && d > kMaxLimitAdd) )
    //             NumberOverflow();
    //
    //         n = n * 10 + d;
    //     }
    //     if ( sign )
    //         return -Int4(n);
    //     else
    //         return n;
    // }
    // ```
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        self.r.skip_white_space()?;
        let mut c = self.get_after_ws()?;
        let mut negative = false;
        if signed && c == b'-' {
            negative = true;
            c = self.r.peek_char(0)?;
            self.r.pos += 1;
        } else if c == b'+' {
            c = self.r.peek_char(0)?;
            self.r.pos += 1;
        }
        if !c.is_ascii_digit() {
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        loop {
            let d = self.r.peek_at(0);
            if !d.is_ascii_digit() {
                break;
            }
            self.r.pos += 1;
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
            }
            n = n * 10 + d;
        }
        // `-Int4(n)` of n = 2^31 is the least Int4 (and likewise in 64 bits).
        Ok(if negative {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    fn int4(&mut self) -> R<i32> {
        Ok(self.get_number(i32::MAX as u64, true)? as i32)
    }

    fn uint4(&mut self) -> R<u32> {
        Ok(self.get_number(u64::from(u32::MAX), false)? as u32)
    }

    /// An INTEGER of storage `storage` (`ReadStd`), or its skip (`SkipSNumber`, `SkipUNumber`):
    /// the number, `None` when skipped.
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
    // ```c++
    // Int2 CObjectIStream::ReadInt2(void)
    // {
    //     Int4 data = ReadInt4();
    //     Int2 ret = Int2(data);
    //     if ( ret != data )
    //         ThrowError(fOverflow, "integer overflow");
    //     return ret;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1780-1808 (the skips of every size:
    // `SkipSNumber` for the signed ones, `SkipUNumber` for the unsigned ones)
    fn integer_number(&mut self, storage: Storage, skip: bool) -> R<Option<i64>> {
        if skip {
            self.skip_number(storage != Storage::Uint4)?;
            return Ok(None);
        }
        let value = match storage {
            Storage::Int4 => i64::from(self.int4()?),
            Storage::Uint4 => i64::from(self.uint4()?),
            Storage::Int8 => self.get_number(i64::MAX as u64, true)?,
            Storage::Int2 => {
                let value = self.int4()?;
                if i16::try_from(value).is_err() {
                    return Err(ReadError::Serial("integer overflow".into()));
                }
                i64::from(value)
            }
        };
        Ok(Some(value))
    }

    /// `SkipSNumber` (`signed`) and `SkipUNumber`: an optional sign (`+` only for the unsigned)
    /// right before a digit, then digits, with no range.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1049
    // ```c++
    // void CObjectIStreamAsn::SkipSNumber(void)
    // {
    //     size_t i;
    //     char c = SkipWhiteSpace();
    //     switch ( c ) {
    //     case '-':
    //     case '+':
    //         c = m_Input.PeekChar(1);
    //         // next char
    //         i = 2;
    //         break;
    //     default:
    //         // next char
    //         i = 1;
    //         break;
    //     }
    //     if ( c < '0' || c > '9' ) {
    //         ThrowError(fFormatError, "bad signed integer in line "
    //             + NStr::SizetToString(m_Input.GetLine()));
    //     }
    //     while ( (c = m_Input.PeekChar(i)) >= '0' && c <= '9' ) {
    //         ++i;
    //     }
    //     m_Input.SkipChars(i);
    // }
    // ```
    // (`SkipUNumber`, 1050-1074, takes only `+` and peeks past the digits with `PeekCharNoEOF`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        let mut c = self.r.skip_white_space()?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = self.r.peek_char(1)?;
            i = 2;
        }
        if !c.is_ascii_digit() {
            return Err(ReadError::Serial(format!(
                "bad {} integer in line {}",
                if signed { "signed" } else { "unsigned" },
                self.r.line
            )));
        }
        loop {
            let c = if signed {
                self.r.peek_char(i)?
            } else {
                self.r.peek_at(i)
            };
            if !c.is_ascii_digit() {
                break;
            }
            i += 1;
        }
        self.r.pos += i;
        Ok(())
    }

    /// `ReadEnum` of the stream: a lower-case name of the list, or a number (one of the list for
    /// an ENUMERATED).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
    // ```c++
    // TEnumValueType CObjectIStreamAsn::ReadEnum(const CEnumeratedTypeValues& values)
    // {
    //     // not integer
    //     CTempString id = ReadLCaseId(SkipWhiteSpace());
    //     if ( !id.empty() ) {
    //         // enum element by name
    //         return values.FindValue(id);
    //     }
    //     // enum element by value
    //     TEnumValueType value = m_Input.GetInt4();
    //     if ( !values.IsInteger() ) // check value
    //         values.FindName(value, false);
    //
    //     return value;
    // }
    // ```
    fn enum_value(&mut self, values: &'static [EnumValue], integer: bool) -> R<i32> {
        let c = self.r.skip_white_space()?;
        let id = self.r.scan_end_of_id(c.is_ascii_lowercase())?;
        if !id.is_empty() {
            let id = self.r.bytes(id);
            return values
                .iter()
                .find(|value| value.name().as_bytes() == id)
                .map(|value| value.value)
                .ok_or_else(|| ReadError::Serial("invalid value of enumerated type".into()));
        }
        let value = self.int4()?;
        if !integer && !values.iter().any(|named| named.value == value) {
            return Err(ReadError::Serial("invalid value of enumerated type".into()));
        }
        Ok(value)
    }

    /// `ReadDouble` / `SkipFNumber`: a word, or `{ mantissa, base, exponent }` whose base is read
    /// with `ReadUint4` (a `CUtilException` for a `-`) and exponent with `ReadInt4`; the number,
    /// `None` when skipped.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (quoted at `Reader::real`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
    // ```c++
    // void CObjectIStreamAsn::SkipFNumber(void)
    // {
    //     if (PeekChar(true) != '{') {
    //         ScanEndOfId(true);
    //         return;
    //     }
    //     Expect('{', true);
    //     SkipSNumber();
    //     Expect(',', true);
    //     unsigned base = ReadUint4();
    //     Expect(',', true);
    //     SkipSNumber();
    //     Expect('}', true);
    //     if ( base != 2 && base != 10 )
    //         ThrowError(fFormatError, "illegal REAL base (must be 2 or 10)");
    // }
    // ```
    fn real_number(&mut self, skip: bool) -> R<Option<f64>> {
        if skip {
            if self.r.skip_white_space()? != b'{' {
                self.r.scan_end_of_id(true)?;
                return Ok(None);
            }
            self.r.pos += 1;
            self.skip_number(true)?;
            self.r.expect(b',')?;
            let base = self.uint4()?;
            self.r.expect(b',')?;
            self.skip_number(true)?;
            self.r.expect(b'}')?;
            if base != 2 && base != 10 {
                return Err(ReadError::Serial(
                    "illegal REAL base (must be 2 or 10)".into(),
                ));
            }
            return Ok(None);
        }
        if let Some(word) = self.r.real_word()? {
            return Ok(Some(word));
        }
        let (mantissa, is_negative) = self.r.real_mantissa()?;
        self.r.expect(b',')?;
        let base = self.uint4()?;
        self.r.expect(b',')?;
        let exp = self.int4()?;
        self.r.expect(b'}')?;
        real_of_parts(mantissa, is_negative, base, exp)
            .map(Some)
            .ok_or_else(|| ReadError::Serial("illegal REAL base (must be 2 or 10)".into()))
    }

    /// `ReadStringValue` with the fix of `kind` / `SkipString`: the string's lines joined, `""` as
    /// one quote, and for a VisibleString each byte outside 0x20..0x7E read as `#` with NCBI's
    /// message, which quotes the piece of the string being appended (`AppendStringData`: up to the
    /// end of the line, the quote, or 128 bytes; none when skipping).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:877-912 (`AppendStringData`)
    // ```c++
    //     const char* data = m_Input.GetCurrentPos();
    //     if ( fix_method != eFNP_Allow ) {
    //         size_t done = 0, valid = 0;
    //         for ( size_t i = 0; i < count; ++i ) {
    //             char c = data[i];
    //             if ( !GoodVisibleChar(c) ) {
    //                 ...
    //                 if (valid == 0) {
    //                     if ( i > done ) {
    //                         s.append(data + done, i - done);
    //                     }
    //                     c = ReplaceVisibleChar(c, fix_method, this, CTempString(data,count), x_FixCharsSubst());
    //                     if (c != 0) {
    //                         s += c;
    //                     }
    //                     done = i + 1;
    //                 }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-986 (`ReadStringValue`, quoted
    // at `Reader::string`: a piece is appended at an end of line, at the quote, and every 128
    // bytes; after `""` the next piece starts with the second quote)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1093-1155
    // ```c++
    // void CObjectIStreamAsn::SkipString(EStringType type)
    // {
    //     ...
    //             default:
    //                 if (type == eStringTypeVisible) {
    //                     if ( !GoodVisibleChar(c) ) {
    //                         ...
    //                         if (valid == 0) {
    //                             c = ReplaceVisibleChar(c, x_FixCharsMethod(), this, kEmptyStr, x_FixCharsSubst());
    //                         }
    // ```
    // (`SERIAL_ALLOW_UTF8_IN_VISIBLESTRING_ON_READING` is 0, objistr.hpp:1163: `valid` stays 0;
    // the fix is `eFNP_ReplaceAndWarn`, the default of `SERIAL_WRONG_CHARS_READ`, objistr.cpp:302.)
    fn string_value(&mut self, kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        let visible = kind != StringKind::Utf8;
        self.r.expect(b'"')?;
        let start_line = self.r.line;
        let unclosed = |line: usize| {
            ReadError::Serial(format!(
                "unclosed string starts at line {start_line} ({line})"
            ))
        };
        let mut out = Vec::new();
        let mut i = 0usize;
        loop {
            let Some(c) = self.r.src.byte_at(self.r.pos + i) else {
                return Err(unclosed(self.r.line));
            };
            match c {
                b'\r' | b'\n' => {
                    self.append(&mut out, i, visible, skip, cx);
                    self.r.pos += i + 1;
                    i = 0;
                    self.r.skip_end_of_line(c);
                }
                b'"' => {
                    self.append(&mut out, i, visible, skip, cx);
                    self.r.pos += i + 1;
                    if self.r.peek_at(0) != b'"' {
                        return Ok(out);
                    }
                    // `""`: the second quote starts the next piece (`SkipString` takes both).
                    if skip {
                        self.r.pos += 1;
                        i = 0;
                    } else {
                        i = 1;
                    }
                }
                _ => {
                    if skip && visible && !(b' '..=b'~').contains(&c) {
                        self.post_bad_char(c, &[], cx);
                    }
                    i += 1;
                    if i == 128 {
                        self.append(&mut out, i, visible, skip, cx);
                        self.r.pos += i;
                        i = 0;
                    }
                }
            }
        }
    }

    /// `AppendStringData` of the `count` bytes at the reader's position (not taken).
    fn append(
        &mut self,
        out: &mut Vec<u8>,
        count: usize,
        visible: bool,
        skip: bool,
        cx: &mut Cx<'_>,
    ) {
        if skip {
            return;
        }
        let data = self.r.bytes(self.r.pos..self.r.pos + count).to_vec();
        if !visible {
            out.extend_from_slice(&data);
            return;
        }
        let mut done = 0;
        for (i, &c) in data.iter().enumerate() {
            if !(b' '..=b'~').contains(&c) {
                out.extend_from_slice(&data[done..i]);
                self.post_bad_char(c, &data, cx);
                out.push(b'#');
                done = i + 1;
            }
        }
        out.extend_from_slice(&data[done..]);
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (quoted at
    // `report::serial::encoder::ReplacedChar`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:807-813
    // ```c++
    // string CObjectIStream::GetPosition(void) const
    // {
    //     string loc_type;
    //     size_t loc_pos;
    //     Location(loc_type, loc_pos);
    //     return loc_type+" "+NStr::Int8ToString(NcbiStreamposToInt8(loc_pos));
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:111-115
    // ```c++
    // void CObjectIStreamAsn::Location(string& loc_type, size_t& loc) const
    // {
    //     loc_type = "line";
    //     loc = m_Input.GetLine();
    // }
    // ```
    fn post_bad_char(&mut self, byte: u8, string: &[u8], cx: &mut Cx<'_>) {
        cx.post_bad_char(byte, string, StreamPosition::Line(self.r.line as u64));
    }

    /// `ReadBytes` / `SkipByteBlock` (the skip takes `H` or `B` after white space): the bytes,
    /// `None` when skipped.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
    // ```c++
    // void CObjectIStreamAsn::SkipByteBlock(void)
    // {
    //     Expect('\'', true);
    //     for ( ;; ) {
    //         char c = GetChar();
    //         if ( c >= '0' && c <= '9' ) {
    //             continue;
    //         }
    //         ...
    //         else if ( c == '\'' ) {
    //             break;
    //         }
    //         else if ( c == '\r' || c == '\n' ) {
    //             SkipEndOfLine(c);
    //         }
    //         else {
    //             m_Input.UngetChar(c);
    //             ThrowError(fFormatError, "bad char in octet string: #"
    //                 + NStr::IntToString(c));
    //         }
    //     }
    //     Expect('H', 'B', true);
    // }
    // ```
    fn octet_string(&mut self, skip: bool) -> R<Option<Vec<u8>>> {
        if !skip {
            return Ok(Some(self.r.octets()?));
        }
        self.r.expect(b'\'')?;
        loop {
            let c = self.r.peek_char(0)?;
            self.r.pos += 1;
            match c {
                b'0'..=b'9' | b'A'..=b'F' | b'a'..=b'f' => {}
                b'\'' => break,
                b'\r' | b'\n' => self.r.skip_end_of_line(c),
                _ => {
                    self.r.pos -= 1;
                    return Err(ReadError::Serial(format!(
                        "bad char in octet string: #{}",
                        c as i8
                    )));
                }
            }
        }
        match self.r.skip_white_space()? {
            b'H' | b'B' => {
                self.r.pos += 1;
                Ok(None)
            }
            _ => Err(ReadError::Serial("'H' or 'B' expected".into())),
        }
    }
}
