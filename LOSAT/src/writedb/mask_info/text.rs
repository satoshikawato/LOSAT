//! ASN.1 text (`MSerial_AsnText`): the primitives of `CObjectIStreamAsn` and of the
//! `CIStreamBuffer` it reads from, on a byte slice that holds the rest of the file.

use super::reader::{format_err, Fail, Input, R};
use super::schema::{ChoiceDef, ClassDef, EnumDef};

pub(super) struct TextIn<'a> {
    data: &'a [u8],
    pos: usize,
    /// `m_BlockStart`
    block_start: bool,
}

// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:101-107
// ```c++
// inline
// bool CObjectIStreamAsn::FirstIdChar(char c)
// {
//     return isalpha((unsigned char) c) || c == '_';
// }
//
// inline
// bool CObjectIStreamAsn::IdChar(char c)
// {
//     return isalnum((unsigned char) c) || c == '_' || c == '.';
// }
// ```
// (The "C" locale: ASCII letters and digits only.)
fn first_id_char(c: u8) -> bool {
    c.is_ascii_alphabetic() || c == b'_'
}

fn id_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// C `isspace` in the "C" locale.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

impl<'a> TextIn<'a> {
    pub(super) fn new(data: &'a [u8]) -> Self {
        TextIn {
            data,
            pos: 0,
            block_start: false,
        }
    }

    /// `m_Input.PeekChar(i)`: the end of the data throws `CEofException`.
    fn peek(&self, i: usize) -> R<u8> {
        self.data.get(self.pos + i).copied().ok_or(Fail::Eof)
    }

    /// `m_Input.PeekCharNoEOF(i)`: NUL at the end of the data.
    fn peek_no_eof(&self, i: usize) -> u8 {
        self.data.get(self.pos + i).copied().unwrap_or(0)
    }

    /// `m_Input.GetChar()`
    fn get_char(&mut self) -> R<u8> {
        let c = self.peek(0)?;
        self.pos += 1;
        Ok(c)
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:526-536
    // ```c++
    // void CIStreamBuffer::SkipEndOfLine(char lastChar)
    //     THROWS1((CIOException))
    // {
    //     ...
    //     m_Line++;
    //     char nextChar = PeekCharNoEOF();
    //     // lastChar either '\r' or \n'
    //     // if nextChar is compliment, skip it
    //     if ( (lastChar + nextChar) == ('\r' + '\n') )
    //         SkipChar();
    // }
    // ```
    fn skip_end_of_line(&mut self, last: u8) {
        let next = self.peek_no_eof(0);
        if u32::from(last) + u32::from(next) == u32::from(b'\r') + u32::from(b'\n') {
            self.pos += 1;
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:207-243
    // ```c++
    // char CObjectIStreamAsn::SkipWhiteSpace(void)
    // {
    //     try { // catch CEofException
    //         for ( ;; ) {
    //             char c = m_Input.SkipSpaces();
    //             switch ( c ) {
    //             case '\t':
    //                 m_Input.SkipChar();
    //                 continue;
    //             case '\r':
    //             case '\n':
    //                 m_Input.SkipChar();
    //                 SkipEndOfLine(c);
    //                 continue;
    //             case '-':
    //                 // check for comments
    //                 if ( m_Input.PeekChar(1) != '-' ) {
    //                     return '-';
    //                 }
    //                 m_Input.SkipChars(2);
    //                 // skip comments
    //                 SkipComments();
    //                 continue;
    //             default:
    //                 return c;
    //             }
    //         }
    //     } catch (CEofException& e) {
    //         if (GetStackDepth() <= 2) {
    //             throw;
    //         } else {
    //             // There should be no eof here, report as error
    //             ThrowError(fEOF, e.what());
    //         }
    //     }
    //     return '\0';
    // }
    // ```
    // (`Rd::frame` makes the depth decision.)
    fn skip_white_space(&mut self) -> R<u8> {
        loop {
            // CIStreamBuffer::SkipSpaces: only ' ' is skipped.
            while self.peek(0)? == b' ' {
                self.pos += 1;
            }
            let c = self.peek(0)?;
            match c {
                b'\t' => self.pos += 1,
                b'\r' | b'\n' => {
                    self.pos += 1;
                    self.skip_end_of_line(c);
                }
                b'-' => {
                    if self.peek(1)? != b'-' {
                        return Ok(b'-');
                    }
                    self.pos += 2;
                    self.skip_comments();
                }
                _ => return Ok(c),
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:245-274
    // ```c++
    // void CObjectIStreamAsn::SkipComments(void)
    // {
    //     try {
    //         for ( ;; ) {
    //             char c = GetChar();
    //             switch ( c ) {
    //             case '\r':
    //             case '\n':
    //                 SkipEndOfLine(c);
    //                 return;
    //             case '-':
    //                 c = GetChar();
    //                 switch ( c ) {
    //                 case '\r':
    //                 case '\n':
    //                     SkipEndOfLine(c);
    //                     return;
    //                 case '-':
    //                     return;
    //                 }
    //                 continue;
    //             default:
    //                 continue;
    //             }
    //         }
    //     }
    //     catch ( CEofException& /* ignored */ ) {
    //         return;
    //     }
    // }
    // ```
    fn skip_comments(&mut self) {
        loop {
            let Ok(c) = self.get_char() else { return };
            match c {
                b'\r' | b'\n' => {
                    self.skip_end_of_line(c);
                    return;
                }
                b'-' => {
                    let Ok(c) = self.get_char() else { return };
                    match c {
                        b'\r' | b'\n' => {
                            self.skip_end_of_line(c);
                            return;
                        }
                        b'-' => return,
                        _ => {}
                    }
                }
                _ => {}
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:276-290
    // ```c++
    // CTempString CObjectIStreamAsn::ScanEndOfId(bool isId)
    // {
    //     if ( isId ) {
    //         for ( size_t i = 1; ; ++i ) {
    //             char c = m_Input.PeekCharNoEOF(i);
    //             if ( !IdChar(c) &&
    //                  (c != '-' || !IdChar(m_Input.PeekChar(i + 1))) ) {
    //                 const char* ptr = m_Input.GetCurrentPos();
    //                 m_Input.SkipChars(i);
    //                 return CTempString(ptr, i);
    //             }
    //         }
    //     }
    //     return CTempString();
    // }
    // ```
    fn scan_end_of_id(&mut self, is_id: bool) -> R<&'a [u8]> {
        if !is_id {
            return Ok(&[]);
        }
        let mut i = 1;
        loop {
            let c = self.peek_no_eof(i);
            if !id_char(c) && (c != b'-' || !id_char(self.peek(i + 1)?)) {
                let start = self.pos;
                self.pos += i;
                return Ok(&self.data[start..start + i]);
            }
            i += 1;
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:343-365
    // ```c++
    // inline
    // CTempString CObjectIStreamAsn::ReadMemberId(char c)
    // {
    //     if ( c == '[' ) {
    //         for ( size_t i = 1; ; ++i ) {
    //             switch ( m_Input.PeekChar(i) ) {
    //             case '\r':
    //             case '\n':
    //                 ThrowError(fFormatError, "end of line: expected ']'");
    //                 break;
    //             case ']':
    //                 {
    //                     const char* ptr = m_Input.GetCurrentPos();
    //                     m_Input.SkipChars(++i);
    //                     return CTempString(ptr + 1, i - 2);
    //                 }
    //             }
    //         }
    //     }
    //     else {
    //         return ScanEndOfId(islower((unsigned char) c) != 0);
    //     }
    // }
    // ```
    fn read_member_id(&mut self, c: u8) -> R<&'a [u8]> {
        if c == b'[' {
            let mut i = 1;
            loop {
                match self.peek(i)? {
                    b'\r' | b'\n' => return format_err("end of line: expected ']'"),
                    b']' => {
                        let start = self.pos;
                        self.pos += i + 1;
                        return Ok(&self.data[start + 1..start + i]);
                    }
                    _ => {}
                }
                i += 1;
            }
        }
        self.scan_end_of_id(c.is_ascii_lowercase())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:128-205
    // ```c++
    // inline
    // bool CObjectIStreamAsn::GetChar(char expect, bool skipWhiteSpace)
    // {
    //     if ( PeekChar(skipWhiteSpace) != expect ) {
    //         return false;
    //     }
    //     m_Input.SkipChar();
    //     return true;
    // }
    //
    // void CObjectIStreamAsn::Expect(char expect, bool skipWhiteSpace)
    // {
    //     if ( !GetChar(expect, skipWhiteSpace) ) {
    //         string msg("\'");
    //         msg += expect;
    //         msg += "' expected";
    //         ThrowError(fFormatError, msg);
    //     }
    // }
    // ```
    fn expect(&mut self, expect: u8, skip_white_space: bool) -> R<()> {
        let c = if skip_white_space {
            self.skip_white_space()?
        } else {
            self.peek(0)?
        };
        if c != expect {
            return format_err(format!("'{}' expected", expect as char));
        }
        self.pos += 1;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1199-1228
    // ```c++
    // void CObjectIStreamAsn::StartBlock(void)
    // {
    //     Expect('{', true);
    //     m_BlockStart = true;
    // }
    //
    // bool CObjectIStreamAsn::NextElement(void)
    // {
    //     char c = SkipWhiteSpace();
    //     if ( m_BlockStart ) {
    //         // first element
    //         m_BlockStart = false;
    //         return c != '}';
    //     }
    //     else {
    //         // next element
    //         if ( c == ',' ) {
    //             m_Input.SkipChar();
    //             return true;
    //         }
    //         else if ( c != '}' )
    //             ThrowError(fFormatError, "',' or '}' expected");
    //         return false;
    //     }
    // }
    //
    // void CObjectIStreamAsn::EndBlock(void)
    // {
    //     Expect('}');
    // }
    // ```
    fn start_block(&mut self) -> R<()> {
        self.expect(b'{', true)?;
        self.block_start = true;
        Ok(())
    }

    fn next_element(&mut self) -> R<bool> {
        let c = self.skip_white_space()?;
        if self.block_start {
            self.block_start = false;
            Ok(c != b'}')
        } else if c == b',' {
            self.pos += 1;
            Ok(true)
        } else if c != b'}' {
            format_err("',' or '}' expected")
        } else {
            Ok(false)
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:389-424
    // ```c++
    // TMemberIndex CObjectIStreamAsn::GetMemberIndex
    //     (const CClassTypeInfo* classType,
    //      const CTempString& id,
    //      const TMemberIndex pos)
    // {
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
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/memberlist.cpp:432-439,462-473
    // ```c++
    // TMemberIndex CItemsInfo::Find(const CTempString& name, TMemberIndex pos) const
    // {
    //     for ( CIterator i(*this, pos); i.Valid(); ++i ) {
    //         if ( name == GetItemInfo(i)->GetId().GetName() )
    //             return *i;
    //     }
    //     return kInvalidMember;
    // }
    // ...
    //         TMemberIndex index = tag + zero_index;
    //         if ( index < pos || index > LastIndex() )
    //             return kInvalidMember;
    //         return index;
    // ```
    // (The names of these types start with a lower case letter, so `GetAltItemIndex`, which looks up
    // the name with an upper case first letter, finds nothing.)
    fn find_member(
        names: impl Iterator<Item = &'static str>,
        id: &[u8],
        pos: usize,
    ) -> R<Option<usize>> {
        if id.first().is_some_and(u8::is_ascii_digit) {
            let tag = std::str::from_utf8(id)
                .ok()
                .and_then(|s| s.parse::<u32>().ok())
                .ok_or_else(|| Fail::Format("bad member number".into()))?
                as usize;
            return Ok(names
                .enumerate()
                .find(|&(i, _)| i == tag && i >= pos)
                .map(|(i, _)| i));
        }
        Ok(names
            .enumerate()
            .skip(pos)
            .find(|(_, name)| name.as_bytes() == id)
            .map(|(i, _)| i))
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1337-1362
    // ```c++
    // TMemberIndex
    // CObjectIStreamAsn::BeginClassMember(const CClassTypeInfo* classType,
    //                                     TMemberIndex pos)
    // {
    //     if ( !NextElement() )
    //         return kInvalidMember;
    //
    //     CTempString id = ReadMemberId(SkipWhiteSpace());
    //     TMemberIndex index = GetMemberIndex(classType, id, pos);
    //     if ( index == kInvalidMember ) {
    //         if (m_TypeAlias && classType->GetMembers().GetItemInfo(pos)->GetId().HasNotag()) {
    //             ...
    //         } else if (CanSkipUnknownMembers()) {
    //             ...
    //         } else {
    //             UnexpectedMember(id, classType->GetMembers());
    //         }
    //     }
    //     return index;
    // }
    // ```
    // (`CanSkipUnknownMembers()` is false by default, `m_TypeAlias` is null.)
    fn begin_member_impl(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        if !self.next_element()? {
            return Ok(None);
        }
        let c = self.skip_white_space()?;
        let id = self.read_member_id(c)?;
        match Self::find_member(def.members.iter().map(|m| m.name), id, pos)? {
            Some(index) => Ok(Some(index)),
            None => format_err(format!(
                "\"{}\": unexpected member",
                String::from_utf8_lossy(id)
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1453-1478
    // ```c++
    // TMemberIndex CObjectIStreamAsn::BeginChoiceVariant(const CChoiceTypeInfo* choiceType)
    // {
    //     bool skipname = !m_BlockStart;
    //     if ( !NextElement() )
    //         return kInvalidMember;
    //     CTempString id = ReadMemberId(SkipWhiteSpace());
    //     if (skipname) {
    //         id = ReadMemberId(SkipWhiteSpace());
    //     }
    //     if ( id.empty() )
    //         ThrowError(fFormatError, "choice variant id expected");
    //
    //     TMemberIndex index = GetChoiceIndex(choiceType, id);
    //     if ( index == kInvalidMember ) {
    //         if (CanSkipUnknownVariants()) {
    //             SetFailFlags(fUnknownValue);
    //         } else {
    //             UnexpectedMember(id, choiceType->GetVariants());
    //         }
    //     }
    //
    //     return index;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:1436-1451
    // ```c++
    // void CObjectIStreamAsn::BeginChoice(const CChoiceTypeInfo* choiceType)
    // {
    //     if (choiceType->GetVariantInfo(kFirstMemberIndex)->GetId().IsAttlist()) {
    //         TopFrame().SetNotag();
    //         StartBlock();
    //     }
    //     m_BlockStart = true;
    // }
    // void CObjectIStreamAsn::EndChoice(void)
    // {
    //     ...
    //     m_BlockStart = false;
    // }
    // ```
    fn begin_variant_impl(&mut self, def: &'static ChoiceDef) -> R<usize> {
        // BeginChoice has just set m_BlockStart, so `skipname` is false.
        if !self.next_element()? {
            return format_err("choice variant id expected");
        }
        let c = self.skip_white_space()?;
        let id = self.read_member_id(c)?;
        if id.is_empty() {
            return format_err("choice variant id expected");
        }
        let found = if id[0].is_ascii_digit() {
            Self::find_member(def.variants.iter().map(|v| v.name), id, 0)?
        } else {
            Self::find_member(def.variants.iter().map(|v| v.name), id, 0)?
        };
        match found {
            Some(index) => Ok(index),
            None => format_err(format!(
                "\"{}\": unexpected member",
                String::from_utf8_lossy(id)
            )),
        }
    }

    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:617-670
    // ```c++
    // char CIStreamBuffer::SkipWs(void)
    // {
    //     char c;
    //     while (isspace((unsigned char)(c = GetChar())))
    //     ;
    //     return c;
    // }
    //
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
    // `GetInt8` is the same with `Int8`/`Uint8`; `GetUint4` has no sign but accepts a leading '+'
    // (strbuffer.cpp:672-700, 703-736).
    fn skip_ws(&mut self) -> R<u8> {
        loop {
            let c = self.get_char()?;
            if !is_c_space(c) {
                return Ok(c);
            }
        }
    }

    /// The digits of a number whose first digit is `first`; `max_before_mul` and `max_limit_add`
    /// are the overflow limits.
    fn digits(&mut self, first: u8, max_before_mul: u64, max_limit_add: u8) -> R<u64> {
        let d0 = first.wrapping_sub(b'0');
        if d0 > 9 {
            return format_err("bad number");
        }
        let mut n = u64::from(d0);
        loop {
            let d = self.peek_no_eof(0).wrapping_sub(b'0');
            if d > 9 {
                return Ok(n);
            }
            self.pos += 1;
            if n > max_before_mul || (n == max_before_mul && d > max_limit_add) {
                return format_err("number overflow");
            }
            n = n * 10 + u64::from(d);
        }
    }

    fn get_int<const MAX: u64>(&mut self) -> R<(bool, u64)> {
        let mut c = self.skip_ws()?;
        let sign = match c {
            b'-' => {
                c = self.get_char()?;
                true
            }
            b'+' => {
                c = self.get_char()?;
                false
            }
            _ => false,
        };
        let n = self.digits(c, MAX / 10, (MAX % 10) as u8 + u8::from(sign))?;
        Ok((sign, n))
    }

    fn get_uint<const MAX: u64>(&mut self) -> R<u64> {
        let mut c = self.skip_ws()?;
        if c == b'+' {
            c = self.get_char()?;
        }
        self.digits(c, MAX / 10, (MAX % 10) as u8)
    }
}

impl Input for TextIn<'_> {
    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:292-309,709-721
    // ```c++
    // CTempString CObjectIStreamAsn::ReadTypeId(char c)
    // {
    //     if ( c == '[' ) {
    //         for ( size_t i = 1; ; ++i ) {
    //             switch ( m_Input.PeekChar(i) ) {
    //             case '\r':
    //             case '\n':
    //                 ThrowError(fFormatError, "end of line: expected ']'");
    //                 break;
    //             case ']':
    //                 {
    //                     const char* ptr = m_Input.GetCurrentPos();
    //                     m_Input.SkipChars(i);
    //                     return CTempString(ptr + 1, i - 2);
    //                 }
    //             }
    //         }
    //     }
    //     else {
    //         return ScanEndOfId(FirstIdChar(c));
    //     }
    // }
    // ...
    // string CObjectIStreamAsn::ReadFileHeader()
    // {
    //     CTempString id = ReadTypeId(SkipWhiteSpace());
    //     string s(id);
    //     if ( SkipWhiteSpace() == ':' &&
    //          m_Input.PeekCharNoEOF(1) == ':' &&
    //          m_Input.PeekCharNoEOF(2) == '=' ) {
    //         m_Input.SkipChars(3);
    //     }
    //     else
    //         ThrowError(fFormatError, "'::=' expected");
    //     return s;
    // }
    // ```
    fn read_file_header(&mut self) -> R<Option<Vec<u8>>> {
        let c = self.skip_white_space()?;
        let name: Vec<u8> = if c == b'[' {
            let mut i = 1;
            loop {
                match self.peek(i)? {
                    b'\r' | b'\n' => return format_err("end of line: expected ']'"),
                    b']' => {
                        // SkipChars(i) leaves the ']' in the input, so the "::=" test below fails.
                        self.pos += i;
                        break Vec::new();
                    }
                    _ => {}
                }
                i += 1;
            }
        } else {
            self.scan_end_of_id(first_id_char(c))?.to_vec()
        };
        if self.skip_white_space()? == b':'
            && self.peek_no_eof(1) == b':'
            && self.peek_no_eof(2) == b'='
        {
            self.pos += 3;
        } else {
            return format_err("'::=' expected");
        }
        Ok(Some(name))
    }

    fn begin_class(&mut self) -> R<()> {
        self.start_block()
    }

    fn begin_member(&mut self, def: &'static ClassDef, pos: usize) -> R<Option<usize>> {
        self.begin_member_impl(def, pos)
    }

    fn end_member(&mut self) -> R<()> {
        Ok(())
    }

    fn end_class(&mut self) -> R<()> {
        self.expect(b'}', false)
    }

    fn begin_choice(&mut self) -> R<()> {
        self.block_start = true;
        Ok(())
    }

    fn begin_variant(&mut self, def: &'static ChoiceDef) -> R<usize> {
        self.begin_variant_impl(def)
    }

    fn end_variant(&mut self) -> R<()> {
        Ok(())
    }

    fn end_choice(&mut self) -> R<()> {
        self.block_start = false;
        Ok(())
    }

    fn begin_container(&mut self, _set: bool) -> R<()> {
        self.start_block()
    }

    fn begin_element(&mut self) -> R<bool> {
        self.next_element()
    }

    fn end_container(&mut self) -> R<()> {
        self.expect(b'}', false)
    }

    fn read_i32(&mut self) -> R<i32> {
        self.skip_white_space()?;
        let (sign, n) = self.get_int::<{ i32::MAX as u64 }>()?;
        // `-Int4(n)` wraps for n = 2^31.
        Ok(if sign {
            (n as i32).wrapping_neg()
        } else {
            n as i32
        })
    }

    fn read_i64(&mut self) -> R<i64> {
        self.skip_white_space()?;
        let (sign, n) = self.get_int::<{ i64::MAX as u64 }>()?;
        Ok(if sign {
            (n as i64).wrapping_neg()
        } else {
            n as i64
        })
    }

    fn read_u32(&mut self) -> R<u32> {
        self.skip_white_space()?;
        Ok(self.get_uint::<{ u32::MAX as u64 }>()? as u32)
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:739-767
    // ```c++
    // bool CObjectIStreamAsn::ReadBool(void)
    // {
    //     switch ( SkipWhiteSpace() ) {
    //     case 'T':
    //         if ( m_Input.PeekCharNoEOF(1) == 'R' &&
    //              m_Input.PeekCharNoEOF(2) == 'U' &&
    //              m_Input.PeekCharNoEOF(3) == 'E' &&
    //              !IdChar(m_Input.PeekCharNoEOF(4)) ) {
    //             m_Input.SkipChars(4);
    //             return true;
    //         }
    //         break;
    //     case 'F':
    //         if ( m_Input.PeekCharNoEOF(1) == 'A' &&
    //              m_Input.PeekCharNoEOF(2) == 'L' &&
    //              m_Input.PeekCharNoEOF(3) == 'S' &&
    //              m_Input.PeekCharNoEOF(4) == 'E' &&
    //              !IdChar(m_Input.PeekCharNoEOF(5)) ) {
    //             m_Input.SkipChars(5);
    //             return false;
    //         }
    //         break;
    //     }
    //     ThrowError(fFormatError, "TRUE or FALSE expected");
    //     return false;
    // }
    // ```
    fn read_bool(&mut self) -> R<bool> {
        let c = self.skip_white_space()?;
        for (word, value) in [(&b"TRUE"[..], true), (&b"FALSE"[..], false)] {
            if c == word[0]
                && word[1..]
                    .iter()
                    .enumerate()
                    .all(|(i, &w)| self.peek_no_eof(i + 1) == w)
                && !id_char(self.peek_no_eof(word.len()))
            {
                self.pos += word.len();
                return Ok(value);
            }
        }
        format_err("TRUE or FALSE expected")
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:444-455
    // ```c++
    // void CObjectIStreamAsn::ReadNull(void)
    // {
    //     if ( SkipWhiteSpace() == 'N' &&
    //          m_Input.PeekCharNoEOF(1) == 'U' &&
    //          m_Input.PeekCharNoEOF(2) == 'L' &&
    //          m_Input.PeekCharNoEOF(3) == 'L' &&
    //          !IdChar(m_Input.PeekCharNoEOF(4)) ) {
    //         m_Input.SkipChars(4);
    //     }
    //     else
    //         ThrowError(fFormatError, "'NULL' expected");
    // }
    // ```
    fn read_null(&mut self) -> R<()> {
        if self.skip_white_space()? == b'N'
            && self.peek_no_eof(1) == b'U'
            && self.peek_no_eof(2) == b'L'
            && self.peek_no_eof(3) == b'L'
            && !id_char(self.peek_no_eof(4))
        {
            self.pos += 4;
            Ok(())
        } else {
            format_err("'NULL' expected")
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:938-994
    // ```c++
    // void CObjectIStreamAsn::ReadStringValue(string& s, EFixNonPrint fix_method)
    // {
    //     Expect('\"', true);
    //     size_t startLine = m_Input.GetLine();
    //     size_t i = 0;
    //     s.erase();
    //     try {
    //         for (;;) {
    //             char c = m_Input.PeekChar(i);
    //             switch ( c ) {
    //             case '\r':
    //             case '\n':
    //                 // flush string
    //                 AppendLongStringData(s, i, fix_method, startLine);
    //                 m_Input.SkipChar(); // '\r' or '\n'
    //                 i = 0;
    //                 // skip end of line
    //                 SkipEndOfLine(c);
    //                 break;
    //             case '\"':
    //                 s.reserve(s.size() + i);
    //                 AppendStringData(s, i, fix_method, startLine);
    //                 m_Input.SkipChar(); // quote
    //                 if ( m_Input.PeekCharNoEOF() != '\"' ) {
    //                     // end of string
    //                     return;
    //                 }
    //                 else {
    //                     // double quote -> one quote
    //                     i = 1;
    //                 }
    //                 break;
    //             default:
    //                 // ok: append char
    //                 if ( ++i == 128 ) {
    //                     // too long string -> flush it
    //                     AppendLongStringData(s, i, fix_method, startLine);
    //                     i = 0;
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    //     catch ( CEofException& ) {
    //         SetFailFlags(fEOF);
    //         UnendedString(startLine);
    //         throw;
    //     }
    // }
    // ```
    // (A line end inside the string is dropped; the quote that follows the closing quote stays in
    // the string; the end of the file inside the string is `UnendedString`, an error.)
    fn read_string(&mut self) -> R<(Vec<u8>, bool)> {
        self.expect(b'"', true)?;
        let mut s = Vec::new();
        loop {
            let Some(&c) = self.data.get(self.pos) else {
                return format_err("unclosed string");
            };
            self.pos += 1;
            match c {
                b'\r' | b'\n' => self.skip_end_of_line(c),
                b'"' => {
                    if self.peek_no_eof(0) == b'"' {
                        s.push(b'"');
                        self.pos += 1;
                    } else {
                        break;
                    }
                }
                _ => s.push(c),
            }
        }
        Ok((s, false))
    }

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
    // NCBI reference (598d8ae6): c++/src/serial/enumerated.cpp:123-137,139-157
    // ```c++
    // TEnumValueType CEnumeratedTypeValues::FindValue(const CTempString& name) const
    // {
    //     const TNameToValue& m = NameToValue();
    //     TNameToValue::const_iterator i = m.find(name);
    //     if ( i == m.end() ) {
    //         string name_alt = string(name);
    //         name_alt[0] = (char)toupper((unsigned char)name_alt[0]);
    //         i = m.find(name_alt);
    //         if ( i == m.end() ) {
    //             NCBI_THROW(CSerialException,eInvalidData,
    //                        "invalid value of enumerated type: " + string(name));
    //         }
    //     }
    //     return i->second;
    // }
    // ```
    // (No name of these enumerations starts with an upper case letter.)
    fn read_enum(&mut self, def: &'static EnumDef) -> R<i32> {
        let c = self.skip_white_space()?;
        let id = self.scan_end_of_id(c.is_ascii_lowercase())?;
        if !id.is_empty() {
            return match def.values.iter().find(|(name, _)| name.as_bytes() == id) {
                Some(&(_, v)) => Ok(v),
                None => format_err(format!(
                    "invalid value of enumerated type: {}",
                    String::from_utf8_lossy(id)
                )),
            };
        }
        let (sign, n) = self.get_int::<{ i32::MAX as u64 }>()?;
        let value = if sign {
            (n as i32).wrapping_neg()
        } else {
            n as i32
        };
        if !def.integer && !def.values.iter().any(|&(_, v)| v == value) {
            return format_err(format!("invalid value of enumerated type: {value}"));
        }
        Ok(value)
    }

    fn offset(&self) -> usize {
        self.pos
    }
}
