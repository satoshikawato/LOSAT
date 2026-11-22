//! The steps of NCBI's ASN.1 text object stream (`CObjectIStreamAsn`) that the walk of the type
//! table (`read`) calls, over the lexical layer of `report::serial::asn_text_reader` (`Reader`):
//! the file header `Type ::=`, the blocks `{ ... }` and their elements, the members by name or by
//! tag in brackets, INTEGER, ENUMERATED, REAL, strings, octets, and the skip functions.

use super::read::{Cx, Front, ReadError, R};
use super::{EnumValue, Member, Storage, StringKind, Value};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_text_reader::{
    real_of_parts, AsnSource, MemberId, Reader, BAD_NUMBER, NUMBER_OVERFLOW,
};
use crate::report::serial::encoder::StreamPosition;
use crate::report::serial::object_frame::{text_frame_site, FrameKind};

/// The text front end: the lexical reader of one text (or of a stream of texts).
pub(crate) struct Text<S: AsnSource> {
    pub(crate) r: Reader<S>,
}

impl<S: AsnSource> Text<S> {
    pub(crate) fn new(r: Reader<S>) -> Self {
        Text { r }
    }
}

impl<S: AsnSource> Front for Text<S> {
    /// The places of the text stream's frames (`object_frame::text_frame_site`, their one owner).
    fn frame_site(kind: FrameKind, skip: bool) -> (&'static str, u32, &'static str) {
        text_frame_site(kind, skip)
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

    /// `BeginChoiceVariant`: `NextElement` at the choice's start, then the variant's id.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1475 (`Reader::choice`
    // cites it)
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
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1648-1652
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
    /// `ReadMemberId` after `SkipWhiteSpace` and `GetMemberIndex` of `classType`, `id` and `pos`
    /// (from `pos`, `in_class`) or `GetChoiceIndex` of `choiceType` and `id`: the member or variant
    /// named or tagged.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:407-424
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
                // `NStr::StringToInt` with the default flags: its `CStringException` leaves
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

    fn util_error(&self, site: ThrowSite, what: &str) -> ReadError {
        ReadError::Exception(NcbiException::new(
            site,
            format!("{what} in line {}", self.r.line),
        ))
    }

    /// `GetInt4`, `GetUint4` (`signed` false: no `-`) and `GetInt8` after `SkipWhiteSpace`, the
    /// digits checked against `max` as NCBI checks them.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:626-668
    // (`GetUint4`, 672-700, takes only `+`, and `GetInt8`, 703-745, is `GetInt4` in 64 bits.)
    fn get_number(&mut self, max: u64, signed: bool) -> R<i64> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let r = &mut self.r;
        let (pos, c) = r.white_space_at(r.pos)?;
        // `SkipWs` then `GetChar`: the first byte that is not C's `isspace`, taken.
        let (pos, mut c) = r.ws_at(pos, c)?;
        let mut pos = pos + 1;
        let mut negative = false;
        if signed && c == b'-' {
            negative = true;
            c = r.char_at(pos)?;
            pos += 1;
        } else if c == b'+' {
            c = r.char_at(pos)?;
            pos += 1;
        }
        if !c.is_ascii_digit() {
            r.pos = pos;
            return Err(self.util_error(BAD_NUMBER, "bad number"));
        }
        let mut n = u64::from(c - b'0');
        let before_mul = max / 10;
        let limit_add = max % 10 + u64::from(negative);
        let mut overflow = false;
        let (end, _) = r.src.scan(pos, |d| {
            if !d.is_ascii_digit() {
                return false;
            }
            let d = u64::from(d - b'0');
            if n > before_mul || (n == before_mul && d > limit_add) {
                overflow = true;
                return false;
            }
            n = n * 10 + d;
            true
        });
        if overflow {
            // `GetInt4` takes the digit before it checks the overflow.
            r.pos = end + 1;
            return Err(self.util_error(NUMBER_OVERFLOW, "number overflow"));
        }
        r.pos = end;
        // The negation of n = 2^31 as an Int4 is the least Int4 (and likewise in 64 bits).
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
    /// the number, `None` when skipped. An Int2 is read as an Int4 and must fit in 16 bits
    /// (`ReadInt2`'s "integer overflow").
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1697-1704
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
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1024-1048
    // (`SkipUNumber`, 1050-1074, takes only `+` and peeks past the digits with `PeekCharNoEOF`.)
    fn skip_number(&mut self, signed: bool) -> R<()> {
        // The position is kept in a local and written once (`Reader::white_space_at`).
        let r = &mut self.r;
        let (start, mut c) = r.white_space_at(r.pos)?;
        let mut i = 1;
        if c == b'+' || (signed && c == b'-') {
            c = match r.src.byte_at(start + 1) {
                Some(c) => c,
                None => return Err(Self::end_of_number(r, start)),
            };
            i = 2;
        }
        if !c.is_ascii_digit() {
            r.pos = start;
            return Err(ReadError::Serial(format!(
                "bad {} integer in line {}",
                if signed { "signed" } else { "unsigned" },
                r.line
            )));
        }
        let (end, c) = r.src.scan(start + i, |c| c.is_ascii_digit());
        // `SkipSNumber` peeks past the digits with `PeekChar`, which fails at the end of the
        // text; `SkipUNumber` with `PeekCharNoEOF`.
        if signed && c.is_none() {
            return Err(Self::end_of_number(r, start));
        }
        r.pos = end;
        Ok(())
    }

    /// The end of the text inside a number that starts at `start`: `PeekChar`'s, at `start`.
    #[cold]
    fn end_of_number(r: &mut Reader<S>, start: usize) -> ReadError {
        r.pos = start;
        r.end_of_text().into()
    }

    /// `ReadEnum` of the stream: a lower-case name of the list, or a number (one of the list for
    /// an ENUMERATED).
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:723-737
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
    /// `None` when skipped. The skip checks the base as the read does ("illegal REAL base").
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:801-857 (cited at `Reader::real`)
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1076-1091
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

    /// `ReadStringValue` with the fix of `kind` / `SkipString` (`Reader::visible_string`): for a
    /// VisibleString each byte outside 0x20..0x7E read as `#`, with NCBI's message, which quotes
    /// the piece of the string being appended (none when skipping).
    fn string_value(&mut self, kind: StringKind, skip: bool, cx: &mut Cx<'_>) -> R<Vec<u8>> {
        let mut out = Vec::new();
        let mut replaced = Vec::new();
        let result = self.r.visible_string(
            &mut out,
            kind != StringKind::Utf8,
            skip,
            &mut |line, byte, piece| replaced.push((line, byte, piece.to_vec())),
        );
        for (line, byte, piece) in replaced {
            Self::post_bad_char(line, byte, &piece, cx);
        }
        result?;
        Ok(out)
    }

    /// `ReplaceVisibleChar`'s post of `byte` at the line `line`: the text stream's position is
    /// the line of its reader (`GetPosition` of the base stream over `Location`).
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (cited at
    // `report::serial::encoder::ReplacedChar`)
    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:807-813
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:111-115
    fn post_bad_char(line: usize, byte: u8, string: &[u8], cx: &mut Cx<'_>) {
        cx.post_bad_char(byte, string, StreamPosition::Line(line as u64));
    }

    /// `ReadBytes` / `SkipByteBlock`: the bytes, `None` when skipped. The skip takes hex digits
    /// up to the closing quote (an end of line skipped, any other byte a "bad char in octet
    /// string" error, put back), then `H` or `B` after white space.
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1170-1197
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
