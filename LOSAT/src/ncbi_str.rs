//! NCBI's number conversions of `NStr` (`c++/src/corelib/ncbistr.cpp`) as LOSAT's ports call
//! them: base 10, the flags each caller passes, and the `CStringException` each throws, with the
//! place it is thrown from ([`ThrowSite`]). `NStr::StringToUInt8_DataSize` is in
//! `writedb::data_size`.
//!
//! Not ported: `fConvErr_NoThrow` and the errno bookkeeping of `CS2N_Guard` (no caller asks for
//! them), the bases other than 10, and the flags no caller passes (`fMandatorySign`,
//! `fAllowCommas`, `fAllowTrailingSymbols`, the leading symbols other than spaces).

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::seqid::parse::is_space;
use crate::writedb::data_size::printable_string;

/// A throw of `S2N_CONVERT_ERROR` in `NStr::<function>` at `line` of ncbistr.cpp.
const fn convert_site(line: u32, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/corelib/ncbistr.cpp",
        line,
        kind: "CStringException::eConvert",
        function: Some(function),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:544-558
// ```c++
// string CS2N_Guard::Message(const CTempString str, const char* to_type, const CTempString msg)
// {
//     string s;
//     s.reserve(str.length() + msg.length() + 50);
//     s += "Cannot convert string '";
//     s += NStr::PrintableString(str);
//     s += "' to ";
//     s += to_type;
//     if ( !msg.empty() ) {
//         s += ", ";
//         s += msg;
//     }
//     return s;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:569-575
// ```c++
// #define S2N_CONVERT_ERROR(to_type, msg, errcode, pos)                 \
//     do {                                                              \
//         err_guard.Set(errcode);                                       \
//         if ( !(flags & NStr::fConvErr_NoThrow) ) {                    \
//             err_guard.Throw();                                        \
//             NCBI_THROW2(CStringException, eConvert,                   \
//                         err_guard.Message(str, #to_type, msg), pos);  \
// ```
/// `S2N_CONVERT_ERROR(to_type, msg, errcode, pos)` at `site` for the string `str`.
fn convert_error(
    str: &[u8],
    to_type: &str,
    msg: &str,
    site: ThrowSite,
    pos: usize,
) -> NcbiException {
    let mut s = format!(
        "Cannot convert string '{}' to {to_type}",
        printable_string(str)
    );
    if !msg.is_empty() {
        s.push_str(", ");
        s.push_str(msg);
    }
    NcbiException::with_pos(site, s, pos)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-873
// ```c++
// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD(flags);
//
//     // Current position in the string
//     SIZE_TYPE pos = 0;
//
//     // Skip allowed leading symbols
//     if (flags & fAllowLeadingSymbols) {
//         bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
//         s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
//     }
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     ...
//     }
//     SIZE_TYPE pos0 = pos;
//     ...
//     Int8 n = 0;
//     Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
//     Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);
//     ...
//     while (char ch = str[pos]) {
//         int  delta;   // corresponding numeric value of 'ch'
//         ...
//         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
//             break;
//         }
//         // Overflow check
//         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
//             S2N_CONVERT_ERROR_OVERFLOW(Int8);
//         }
//         n *= base;
//         n += delta;
//         pos++;
//     }
//
//     // Last checks
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     ...
//     // Assign sign before the end pointer check
//     n = sign ? -n : n;
//     CHECK_ENDPTR(Int8);
//
//     return n;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:740-759
// ```c++
// void s_SkipAllowedSymbols(const CTempString       str,
//                           SIZE_TYPE&              pos,
//                           ESkipMode               skip_mode,
//                           NStr::TStringToNumFlags flags)
// {
//     ...
//     for ( SIZE_TYPE len = str.length();  pos < len;  ++pos ) {
//         unsigned char ch = str[pos];
//         if ( isdigit(ch)  ||  ch == '+' ||  ch == '-'  ||  s_IsDecimalPoint(ch,flags) ) {
//             break;
//         }
//         if ( (skip_mode == eSkipSpacesOnly)  &&  !isspace(ch) ) {
//             break;
//         }
//     }
// }
// ```
// The string is a `std::string`'s or a C string's: the loop ends at its terminating NUL, and at
// an embedded one (whose end check then passes). The errors are thrown at lines 852
// (overflow), 861 (no digit) and 870 (`CHECK_ENDPTR`).
const STRING_TO_INT8_OVERFLOW: ThrowSite = convert_site(852, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_INVAL: ThrowSite = convert_site(861, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_ENDPTR: ThrowSite = convert_site(870, "ncbi::NStr::StringToInt8");

/// `NStr::StringToInt8(str, flags)` in base 10; `leading_spaces` is `fAllowLeadingSpaces`,
/// the only flag a caller passes.
pub(crate) fn string_to_int8(s: &[u8], leading_spaces: bool) -> Result<i64, NcbiException> {
    let at = |pos: usize| s.get(pos).copied().unwrap_or(0);
    let error = |site: ThrowSite, msg: &str, pos: usize| convert_error(s, "Int8", msg, site, pos);
    let mut pos = 0;
    if leading_spaces {
        while pos < s.len() {
            let ch = s[pos];
            if ch.is_ascii_digit() || ch == b'+' || ch == b'-' || !is_space(ch) {
                break;
            }
            pos += 1;
        }
    }
    let sign = match at(pos) {
        b'-' => {
            pos += 1;
            true
        }
        b'+' => {
            pos += 1;
            false
        }
        _ => false,
    };
    let pos0 = pos;
    let limdiv = i64::MAX / 10;
    let limoff = i64::MAX % 10 + i64::from(sign);
    let mut n: i64 = 0;
    loop {
        let ch = at(pos);
        if ch == 0 || !ch.is_ascii_digit() {
            break;
        }
        let delta = i64::from(ch - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(STRING_TO_INT8_OVERFLOW, "overflow", pos));
        }
        // `i64::MIN` passes the check as its magnitude, one past `i64::MAX`.
        n = n.wrapping_mul(10).wrapping_add(delta);
        pos += 1;
    }
    if pos == pos0 {
        return Err(error(STRING_TO_INT8_INVAL, "", pos));
    }
    let n = if sign { n.wrapping_neg() } else { n };
    if at(pos) != 0 {
        return Err(error(STRING_TO_INT8_ENDPTR, "", pos));
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
// ```c++
// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Int8 value = StringToInt8(str, flags, base);
//     if ( value < kMin_Int  ||  value > kMax_Int ) {
//         S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
//     }
//     return (int) value;
// }
// ```
const STRING_TO_INT_OVERFLOW: ThrowSite = convert_site(640, "ncbi::NStr::StringToInt");

/// `NStr::StringToInt(str, flags)` in base 10; `leading_spaces` as for [`string_to_int8`].
pub(crate) fn string_to_int(s: &[u8], leading_spaces: bool) -> Result<i32, NcbiException> {
    let value = string_to_int8(s, leading_spaces)?;
    i32::try_from(value).map_err(|_| convert_error(s, "int", "overflow", STRING_TO_INT_OVERFLOW, 0))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:876-917
// ```c++
// Uint8 NStr::StringToUInt8(const CTempString str,
//                           TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD(flags);
//
//     const TStringToNumFlags slow_flags =
//         fMandatorySign|fAllowCommas|fAllowLeadingSymbols|fAllowTrailingSymbols;
//
//     if ( base == 10  &&  (flags & slow_flags) == 0 ) {
//         // fast conversion
//
//         // Current position in the string
//         CTempString::const_iterator ptr = str.begin(), end = str.end();
//
//         // Determine sign
//         if ( ptr != end && *ptr == '+' ) {
//             ++ptr;
//         }
//         if ( ptr == end ) {
//             S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
//         }
//
//         // Begin conversion
//         Uint8 n = 0;
//
//         const Uint8 limdiv = kMax_UI8/10;
//         const int   limoff = int(kMax_UI8 % 10);
//
//         do {
//             char ch = *ptr;
//             int  delta = ch - '0';
//             if ( unsigned(delta) >= 10 ) {
//                 S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
//             }
//             // Overflow check
//             if ( n >= limdiv && (n > limdiv || delta > limoff) ) {
//                 S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, ptr-str.begin());
//             }
//             n = n*10+delta;
//         } while ( ++ptr != end );
//
//         return n;
//     }
// ```
// With the default flags only the fast conversion runs. Its overflow error has no message
// after the type (`kEmptyStr`, unlike `S2N_CONVERT_ERROR_OVERFLOW` of the other conversions):
// makeblastdb 2.17.0 with `BLASTDB_LMDB_MAP_SIZE=99999999999999999999` prints "Cannot convert
// string '99999999999999999999' to Uint8 (m_Pos = 19)". A byte of 0x80 or more is a negative
// `char`: not a digit either.
const STRING_TO_UINT8_EMPTY: ThrowSite = convert_site(895, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_DIGIT: ThrowSite = convert_site(908, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_OVERFLOW: ThrowSite = convert_site(912, "ncbi::NStr::StringToUInt8");

/// `NStr::StringToUInt8(str)` with the default flags and base (the fast conversion).
pub(crate) fn string_to_uint8(s: &[u8]) -> Result<u64, NcbiException> {
    let error = |site: ThrowSite, pos: usize| convert_error(s, "Uint8", "", site, pos);
    let mut ptr = 0;
    if ptr != s.len() && s[ptr] == b'+' {
        ptr += 1;
    }
    if ptr == s.len() {
        return Err(error(STRING_TO_UINT8_EMPTY, ptr));
    }
    let limdiv = u64::MAX / 10;
    let limoff = u64::MAX % 10;
    let mut n: u64 = 0;
    while ptr != s.len() {
        let ch = s[ptr];
        if !ch.is_ascii_digit() {
            return Err(error(STRING_TO_UINT8_DIGIT, ptr));
        }
        let delta = u64::from(ch - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(STRING_TO_UINT8_OVERFLOW, ptr));
        }
        n = n * 10 + delta;
        ptr += 1;
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:646-655
// ```c++
// unsigned int
// NStr::StringToUInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Uint8 value = StringToUInt8(str, flags, base);
//     if ( value > kMax_UInt ) {
//         S2N_CONVERT_ERROR(unsigned int, "overflow", ERANGE, 0);
//     }
//     return (unsigned int) value;
// }
// ```
const STRING_TO_UINT_OVERFLOW: ThrowSite = convert_site(652, "ncbi::NStr::StringToUInt");

/// `NStr::StringToUInt(str)` with the default flags and base.
pub(crate) fn string_to_uint(s: &[u8]) -> Result<u32, NcbiException> {
    let value = string_to_uint8(s)?;
    u32::try_from(value)
        .map_err(|_| convert_error(s, "unsigned int", "overflow", STRING_TO_UINT_OVERFLOW, 0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The message, the line and the position of an error.
    fn failure(e: NcbiException) -> (String, u32, String) {
        let what = e.what();
        let pos = what
            .rsplit_once("(m_Pos = ")
            .map(|(_, rest)| rest.trim_end_matches(")\n").to_string())
            .unwrap_or_default();
        (e.get_msg().to_string(), e.site().line, pos)
    }

    // makeblastdb 2.17.0 (the official Linux x64 release), probe of 2026-11-08 (cases
    // mapsz_*_2seq, txn_*_2seq; taxid_map_*): the line and position of each conversion error.
    #[test]
    fn errors_are_thrown_where_ncbi_throws_them() {
        let cases: [(Result<(), NcbiException>, &str, u32, &str); 9] = [
            (
                string_to_uint8(b"abc").map(drop),
                "Cannot convert string 'abc' to Uint8",
                908,
                "0",
            ),
            (
                string_to_uint8(b"1x").map(drop),
                "Cannot convert string '1x' to Uint8",
                908,
                "1",
            ),
            (
                string_to_uint8(b"+").map(drop),
                "Cannot convert string '+' to Uint8",
                895,
                "1",
            ),
            (
                string_to_uint8(b"99999999999999999999").map(drop),
                "Cannot convert string '99999999999999999999' to Uint8",
                912,
                "19",
            ),
            (
                string_to_int(b"abc", false).map(drop),
                "Cannot convert string 'abc' to Int8",
                861,
                "0",
            ),
            (
                string_to_int(b"-1x", false).map(drop),
                "Cannot convert string '-1x' to Int8",
                870,
                "2",
            ),
            (
                string_to_int(b"99999999999", false).map(drop),
                "Cannot convert string '99999999999' to int, overflow",
                640,
                "0",
            ),
            (
                string_to_int(b"99999999999999999999", false).map(drop),
                "Cannot convert string '99999999999999999999' to Int8, overflow",
                852,
                "18",
            ),
            (
                string_to_uint(b"4294967296").map(drop),
                "Cannot convert string '4294967296' to unsigned int, overflow",
                652,
                "0",
            ),
        ];
        for (result, message, line, pos) in cases {
            assert_eq!(
                failure(result.unwrap_err()),
                (message.to_string(), line, pos.to_string()),
                "{message}"
            );
        }
        let e = string_to_uint8(b"abc").unwrap_err();
        assert_eq!(
            e.site(),
            ThrowSite {
                file: "src/corelib/ncbistr.cpp",
                line: 908,
                kind: "CStringException::eConvert",
                function: Some("ncbi::NStr::StringToUInt8"),
            }
        );
    }

    // NStr::StringToInt / StringToUInt / StringToUInt8 with the default flags (makeblastdb
    // 2.17.0: MAX_LMDB_TXN_ENTRY="+3" is 3, " 3" and "abc" stop it, "3000000000" overflows
    // int, BLASTDB_LMDB_MAP_SIZE="-5" stops it).
    #[test]
    fn values_convert_as_ncbi_converts_them() {
        assert_eq!(string_to_int(b"+3", false), Ok(3));
        assert_eq!(string_to_int(b"-1", false), Ok(-1));
        assert_eq!(string_to_int(b"-2147483648", false), Ok(i32::MIN));
        assert!(string_to_int(b" 3", false).is_err());
        assert_eq!(string_to_int(b" 3", true), Ok(3));
        assert!(string_to_int(b"3 ", true).is_err());
        assert!(string_to_int(b"", false).is_err());
        assert_eq!(string_to_int8(b"-9223372036854775808", false), Ok(i64::MIN));
        assert_eq!(string_to_int8(b"12\x0034", false), Ok(12));
        assert_eq!(string_to_uint(b"4294967295"), Ok(u32::MAX));
        assert_eq!(string_to_uint(b"+0"), Ok(0));
        assert!(string_to_uint8(b"-5").is_err());
        assert!(string_to_uint8(b"").is_err());
        assert!(string_to_uint8(b" 1").is_err());
        assert!(string_to_uint8(b"\xb11").is_err());
        assert_eq!(string_to_uint8(b"+40960"), Ok(40960));
        assert_eq!(string_to_uint8(b"18446744073709551615"), Ok(u64::MAX));
        assert_eq!(
            failure(string_to_uint8(b"18446744073709551616").unwrap_err()),
            (
                "Cannot convert string '18446744073709551616' to Uint8".to_string(),
                912,
                "19".to_string()
            )
        );
    }
}
