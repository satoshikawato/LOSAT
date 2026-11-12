//! The functions of NCBI's `NStr` (`c++/src/corelib/ncbistr.cpp`) that makeblastdb and
//! rpsbproc's data loader (`rpsbproc::data`) use to convert numbers, as LOSAT's ports call them:
//! base 10, the flags each caller passes, and the `CStringException` each throws, with the place
//! it is thrown from ([`ThrowSite`]); with them `NStr::PrintableString`, which their messages
//! quote the string with, `NStr::StringToUInt8_DataSize` (`-max_file_sz`, module `data_size`)
//! and `NStr::StringToDouble` with glibc's `strtod` (rpsbproc's bit score thresholds, module
//! `double`). Also the string functions that more than one port calls: `NStr::ParseEscapes`
//! (makeblastdb's titles and the registry reader) and `NStr::TruncateSpaces`.
//!
//! Not ported: `fConvErr_NoThrow` and the errno bookkeeping of `CS2N_Guard` (no caller asks for
//! them), the bases other than 10 (but 16 of `StringToInt8`, which the argument reader tries,
//! `ncbi_app::args`), and the flags no caller passes (`fMandatorySign`, `fAllowCommas`, the
//! leading and trailing symbols other than spaces). Other modules keep copies
//! of some of these conversions (`blastinput::query_batch`, `blastinput::value_parsers`,
//! `seqdb::isam`, `seqdb::mask`, `blastdbcmd`), not yet folded here.

mod data_size;
mod double;

pub use data_size::string_to_uint8_data_size;
pub(crate) use double::string_to_double;

use crate::ncbi_exception::{NcbiException, ThrowSite};

/// `NStr::TStringToNumFlags`, of the flags LOSAT's callers pass.
pub(crate) type Flags = u32;

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:294-298
// ```c++
//         fAllowLeadingSpaces      = (1 << 19),     ///< Ignore leading whitespace characters in converted string.
//         fAllowLeadingSymbols     = (1 << 20) | fAllowLeadingSpaces,
//         ...
//         fAllowTrailingSpaces     = (1 << 21),     ///< Ignore trailing whitespace characters.
//         fAllowTrailingSymbols    = (1 << 22) | fAllowTrailingSpaces,
// ```
/// `NStr::fAllowLeadingSpaces`.
pub(crate) const ALLOW_LEADING_SPACES: Flags = 1 << 19;
/// `NStr::fAllowTrailingSpaces`.
pub(crate) const ALLOW_TRAILING_SPACES: Flags = 1 << 21;

/// `isspace` in the C locale.
pub(crate) fn is_space(c: u8) -> bool {
    c == b' ' || (b'\t'..=b'\r').contains(&c)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:106-115
// ```c++
// bool NStr::IsBlank(const CTempString str, SIZE_TYPE pos)
// {
//     SIZE_TYPE len = str.length();
//     for (SIZE_TYPE idx = pos; idx < len; ++idx) {
//         if (!isspace((unsigned char) str[idx])) {
//             return false;
//         }
//     }
//     return true;
// }
// ```
/// `NStr::IsBlank(s)`: every byte is a C `isspace` space (an empty string included).
pub(crate) fn is_blank(s: &[u8]) -> bool {
    s.iter().all(|&c| is_space(c))
}

/// A throw of `S2N_CONVERT_ERROR` in `NStr::<function>` at `line` of ncbistr.cpp.
const fn convert_site(line: u32, function: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/corelib/ncbistr.cpp",
        line,
        kind: "CStringException::eConvert",
        function: Some(function),
    }
}

/// `NStr::PrintableString(str)` with its default mode `fNewLine_Quote | fNonAscii_Passthru` (both 0) and language C,
/// which is what `CS2N_Guard::Message` uses on the input string.
/// Ported branches: `s_PrintableString` with `mode == 0`, `lang == eLanguage_C`
/// (so `fNewLine_Passthru`, `fNonAscii_Quote`, `fPrintable_Full` and `eLanguage_Javascript` are not ported).
pub(crate) fn printable_string(str_: &[u8]) -> String {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3736-3742
    // ```c++
    //     unique_ptr<CNcbiOstrstream> out;
    //     SIZE_TYPE i, j = 0;
    //
    //     for (i = 0;  i < str.size();  ++i) {
    //         bool octal = false;
    //         char c = str[i];
    //         switch (c) {
    // ```
    let mut out: Option<Vec<u8>> = None;
    let mut j: usize = 0;
    for i in 0..str_.len() {
        let mut octal = false;
        let mut c = str_[i];
        match c {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3743-3771
            // ```c++
            //         case '\a':
            //             if (lang == eLanguage_C)
            //                 c = 'a';
            //             else
            //                 octal = true;
            //             break;
            //         case '\b':
            //             c = 'b';
            //             break;
            //         case '\f':
            //             c = 'f';
            //             break;
            //         case '\r':
            //             c = 'r';
            //             break;
            //         case '\t':
            //             c = 't';
            //             break;
            //         case '\v':
            //             c = 'v';
            //             break;
            //         case '\n':
            //             if (!(mode & NStr::fNewLine_Passthru))
            //                 c = 'n';
            //             /*FALLTHRU*/
            //         case '\\':
            //         case '\'':
            //         case '"':
            //             break;
            // ```
            // '\a' (lang is C)
            0x07 => c = b'a',
            0x08 => c = b'b',
            0x0c => c = b'f',
            b'\r' => c = b'r',
            b'\t' => c = b't',
            0x0b => c = b'v',
            // '\n': fNewLine_Passthru is not set, so c becomes 'n'
            b'\n' => c = b'n',
            b'\\' | b'\'' | b'"' => {}
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3772-3775
            // ```c++
            //         case '&':
            //             if (lang == eLanguage_Javascript)
            //                 break;
            //             continue;
            // ```
            b'&' => continue,
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3776-3783
            // ```c++
            //         case '?':
            //             if (lang == eLanguage_C) {
            //                 if (i  &&  str[i - 1] == '?')
            //                     break;
            //                 if (i < str.size() - 1  &&  str[i + 1] == '?')
            //                     break;
            //             }
            //             continue;
            // ```
            b'?' => {
                if !((i != 0 && str_[i - 1] == b'?') || (i < str_.len() - 1 && str_[i + 1] == b'?'))
                {
                    continue;
                }
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3784-3795
            // ```c++
            //         default:
            //             if (!isascii((unsigned char) c)) {
            //                 if (mode & NStr::fNonAscii_Quote) {
            //                     octal = true;
            //                     break;
            //                 }
            //             }
            //             if (!isprint((unsigned char) c)) {
            //                 octal = true;
            //                 break;
            //             }
            //             continue;
            // ```
            _ => {
                // fNonAscii_Quote is not set; isprint of the "C" locale is 0x20..=0x7e
                if !(0x20..=0x7e).contains(&c) {
                    octal = true;
                } else {
                    continue;
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3797-3803
        // ```c++
        //         if (!out.get()) {
        //             out.reset(new CNcbiOstrstream);
        //         }
        //         if (i > j) {
        //             out->write(str.data() + j, i - j);
        //         }
        //         out->put('\\');
        // ```
        let out_buf = out.get_or_insert_with(Vec::new);
        if i > j {
            out_buf.extend_from_slice(&str_[j..i]);
        }
        out_buf.push(b'\\');
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3804-3831
        // ```c++
        //         if (c == '\n') {
        //             out->write("n\\\n", 3);
        //         } else if (octal) {
        //             bool reduce;
        //             if (!(mode & NStr::fPrintable_Full)) {
        //                 reduce = (i == str.size() - 1  ||
        //                           str[i + 1] < '0' || '7' < str[i + 1] ? true : false);
        //             } else {
        //                 reduce = false;
        //             }
        //             unsigned char v;
        //             char val[3];
        //             int k = 0;
        //             v =  (unsigned char) c >> 6;
        //             if (v  ||  !reduce) {
        //                 val[k++] = char('0' + v);
        //                 reduce = false;
        //             }
        //             v = ((unsigned char) c >> 3) & 7;
        //             if (v  ||  !reduce) {
        //                 val[k++] = char('0' + v);
        //             }
        //             v =  (unsigned char) c       & 7;
        //             val[k++] = char('0' + v);
        //             out->write(val, k);
        //         } else {
        //             out->put(c);
        //         }
        // ```
        // `c == '\n'` cannot be true here: c was changed to 'n' above (no fNewLine_Passthru).
        if octal {
            // fPrintable_Full is not set
            let mut reduce = i == str_.len() - 1 || str_[i + 1] < b'0' || b'7' < str_[i + 1];
            let mut val: Vec<u8> = Vec::with_capacity(3);
            let mut v = c >> 6;
            if v != 0 || !reduce {
                val.push(b'0' + v);
                reduce = false;
            }
            v = (c >> 3) & 7;
            if v != 0 || !reduce {
                val.push(b'0' + v);
            }
            v = c & 7;
            val.push(b'0' + v);
            out_buf.extend_from_slice(&val);
        } else {
            out_buf.push(c);
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3832-3832
        // ```c++
        //         j = i + 1;
        // ```
        j = i + 1;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3834-3844
    // ```c++
    //     if (j  &&  i > j) {
    //         _ASSERT(out.get());
    //         out->write(str.data() + j, i - j);
    //     }
    //     if (out.get()) {
    //         // Return encoded string
    //         return CNcbiOstrstreamToString(*out);
    //     }
    //
    //     // All characters are good - return (a copy of) the original string
    //     return str;
    // ```
    let i = str_.len();
    if j != 0 && i > j {
        out.as_mut()
            .expect("out is set whenever j != 0")
            .extend_from_slice(&str_[j..i]);
    }
    match out {
        Some(buf) => String::from_utf8_lossy(&buf).into_owned(),
        None => String::from_utf8_lossy(str_).into_owned(),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:545-558
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
//     // Skip allowed trailing symbols
//     if (flags & fAllowTrailingSymbols) {
//         bool spaces = ((flags & fAllowTrailingSymbols) == fAllowTrailingSpaces);
//         s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
//     }
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
// `CTempString::operator[]` gives NUL past the end (tempstr.hpp:341-347): the loop ends there and
// at an embedded NUL, whose end check (`if ( str[pos] )`, lines 603-606) then passes, so
// "12\0x" is 12. A decimal point stops the skipping of trailing spaces, and the end check then
// fails on it, as on any other byte that is not a space. The errors are thrown at lines 852
// (overflow), 861 (no digit) and 870 (`CHECK_ENDPTR`).
const STRING_TO_INT8_OVERFLOW: ThrowSite = convert_site(852, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_INVAL: ThrowSite = convert_site(861, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_ENDPTR: ThrowSite = convert_site(870, "ncbi::NStr::StringToInt8");

/// `s_SkipAllowedSymbols(str, pos, eSkipSpacesOnly, flags)`: the position after the white space
/// from `pos`.
fn skip_spaces(s: &[u8], mut pos: usize) -> usize {
    while pos < s.len() {
        let ch = s[pos];
        if ch.is_ascii_digit() || ch == b'+' || ch == b'-' || !is_space(ch) {
            break;
        }
        pos += 1;
    }
    pos
}

/// `NStr::StringToInt8(str, flags)` in base 10, `flags` [`ALLOW_LEADING_SPACES`] and
/// [`ALLOW_TRAILING_SPACES`] or 0.
pub(crate) fn string_to_int8(s: &[u8], flags: Flags) -> Result<i64, NcbiException> {
    string_to_int8_base(s, flags, 10)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:765-795
// ```c++
// static inline
// bool s_CheckRadix(const CTempString str, SIZE_TYPE& pos, int& base)
// {
//     if ( base == 10  ||  base == 8 ) {
//         // shortcut for most frequent case
//         return true;
//     }
//     ...
//     unsigned char ch   = str[pos];
//     unsigned char next = str[pos+1];
//     ...
//     // Remove leading '0x' for hex numbers
//     if ( base == 16 ) {
//         if (ch == '0'  &&  (next == 'x' || next == 'X')) {
//             pos += 2;
//         }
//     }
//     return true;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:684-712
// ```c++
// bool s_IsGoodCharForRadix(char ch, int base, int* value = 0)
// {
//     if ( base <= 10 ) {
//         // shortcut for most frequent case
//         int delta = ch-'0';
//         if ( unsigned(delta) < unsigned(base) ) {
//             ...
//             return true;
//         }
//         return false;
//     }
//     if (!isalnum((unsigned char) ch)) {
//         return false;
//     }
//     // Corresponding numeric value of *endptr
//     int delta;
//     if (isdigit((unsigned char) ch)) {
//         delta = ch - '0';
//     } else {
//         ch = (char) tolower((unsigned char) ch);
//         delta = ch - 'a' + 10;
//     }
//     ...
//     return delta < base;
//  }
// ```
// `pos0` is set before the radix check: a `0x` with no digit after it is 0 (`pos != pos0`).
/// `NStr::StringToInt8(str, flags, base)` in base 10 or 16 (`s_StringToInt8` of the argument
/// reader tries 16 after 10, `ncbi_app::args`), `flags` as for [`string_to_int8`].
pub(crate) fn string_to_int8_base(s: &[u8], flags: Flags, base: u8) -> Result<i64, NcbiException> {
    debug_assert!(base == 10 || base == 16, "only bases 10 and 16 are ported");
    let at = |pos: usize| s.get(pos).copied().unwrap_or(0);
    let error = |site: ThrowSite, msg: &str, pos: usize| convert_error(s, "Int8", msg, site, pos);
    let mut pos = 0;
    if flags & ALLOW_LEADING_SPACES != 0 {
        pos = skip_spaces(s, pos);
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
    if base == 16 && at(pos) == b'0' && matches!(at(pos + 1), b'x' | b'X') {
        pos += 2;
    }
    let radix = i64::from(base);
    let limdiv = i64::MAX / radix;
    let limoff = i64::MAX % radix + i64::from(sign);
    let mut n: i64 = 0;
    loop {
        let ch = at(pos);
        let digit = match ch {
            0 => None,
            b'0'..=b'9' => Some(ch - b'0'),
            _ if base > 10 && ch.is_ascii_alphabetic() => Some(ch.to_ascii_lowercase() - b'a' + 10),
            _ => None,
        };
        let Some(delta) = digit.filter(|&delta| delta < base) else {
            break;
        };
        let delta = i64::from(delta);
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(STRING_TO_INT8_OVERFLOW, "overflow", pos));
        }
        // `i64::MIN` passes the check as its magnitude, one past `i64::MAX`.
        n = n.wrapping_mul(radix).wrapping_add(delta);
        pos += 1;
    }
    if pos == pos0 {
        return Err(error(STRING_TO_INT8_INVAL, "", pos));
    }
    if flags & ALLOW_TRAILING_SPACES != 0 {
        pos = skip_spaces(s, pos);
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

/// `NStr::StringToInt(str, flags)` in base 10, `flags` as for [`string_to_int8`].
pub(crate) fn string_to_int(s: &[u8], flags: Flags) -> Result<i32, NcbiException> {
    let value = string_to_int8(s, flags)?;
    i32::try_from(value).map_err(|_| convert_error(s, "int", "overflow", STRING_TO_INT_OVERFLOW, 0))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:876-918
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

/// `NStr::ParseEscapes(str)` with its default mode (`eEscSeqRange_Standard`), as makeblastdb's
/// title (`writedb::fasta_source::title`) and NCBI's registry reader
/// (`blastinput::ncbi_registry`) call it. An error is the message of NCBI's `CStringException`
/// (an escape that ends the string, `\x` without a hexadecimal digit, or one whose value
/// `StringToUInt` cannot hold).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892
/// ```c++
/// string NStr::ParseEscapes(const CTempString str, EEscSeqRange mode, char user_char)
/// {
///     string out;
///     out.reserve(str.size());  // result string can only be smaller
///     SIZE_TYPE pos = 0;
///     bool is_error = false;
///
///     while (pos < str.size()  ||  !is_error) {
///         SIZE_TYPE pos2 = str.find('\\', pos);
///         if (pos2 == NPOS) {
///             //~ out += str.substr(pos);
///             CTempString sub(str, pos);
///             out += sub;
///             break;
///         }
///         //~ out += str.substr(pos, pos2 - pos);
///         CTempString sub(str, pos, pos2-pos);
///         out += sub;
///         if (++pos2 == str.size()) {
///             NCBI_THROW2(CStringException, eFormat,
///                         "Unterminated escape sequence", pos2);
///         }
///         switch (str[pos2]) {
///         case 'a':  out += '\a';  break;
///         case 'b':  out += '\b';  break;
///         case 'f':  out += '\f';  break;
///         case 'n':  out += '\n';  break;
///         case 'r':  out += '\r';  break;
///         case 't':  out += '\t';  break;
///         case 'v':  out += '\v';  break;
///         case 'x':
///             {{
///                 pos = ++pos2;
///                 while (pos < str.size()
///                        &&  isxdigit((unsigned char) str[pos])) {
///                     pos++;
///                 }
///                 if (pos > pos2) {
///                     SIZE_TYPE len = pos-pos2;
///                     if ((mode == eEscSeqRange_FirstByte) && (len > 2)) {
///                         // Take only 2 first hex-digits
///                         len = 2;
///                         pos = pos2 + 2;
///                     }
///                     unsigned int value =
///                         StringToUInt(CTempString(str, pos2, len), 0, 16);
///                     if ((mode != eEscSeqRange_Standard)  &&  (value > 255)) {
///                         // eEscSeqRange_Standard -- by default
///                         switch (mode) {
///                         case eEscSeqRange_FirstByte:
///                             // Already have right value
///                             break;
///                         case eEscSeqRange_Throw:
///                             NCBI_THROW2(CStringException, eFormat,
///                                 "Escape sequence '" + NStr::PrintableString(CTempString(str, pos2, len)) +
///                                 "' is out of range [0-255]", pos2);
///                             break;
///                         case eEscSeqRange_Errno:
///                             CNcbiError::SetErrno(errno = ERANGE, str);
///                             is_error = true;
///                             continue;
///                         case eEscSeqRange_User:
///                             value = (unsigned)user_char;
///                             break;
///                         default:
///                             NCBI_THROW2(CStringException, eFormat, "Wrong set of flags", pos2);
///                         }
///                     }
///                     out += static_cast<char>(value);
///                 } else {
///                     NCBI_THROW2(CStringException, eFormat,
///                                 "\\x followed by no hexadecimal digits", pos);
///                 }
///             }}
///             continue;
///         case '0':  case '1':  case '2':  case '3':
///         case '4':  case '5':  case '6':  case '7':
///             {{
///                 pos = pos2;
///                 unsigned char c = (unsigned char)(str[pos++] - '0');
///                 while (pos < pos2 + 3  &&  pos < str.size()
///                        &&  str[pos] >= '0'  &&  str[pos] <= '7') {
///                     c = (unsigned char)((c << 3) | (str[pos++] - '0'));
///                 }
///                 out += c;
///             }}
///             continue;
///         case '\n':
///             // quoted EOL means no EOL
///             break;
///         default:
///             out += str[pos2];
///             break;
///         }
///         pos = pos2 + 1;
///     }
///     if (mode == eEscSeqRange_Errno) {
///         if (is_error) {
///             return kEmptyStr;
///         }
///         errno = 0;
///     }
///     return out;
/// }
/// ```
pub(crate) fn parse_escapes(str: &[u8]) -> Result<Vec<u8>, &'static str> {
    let mut out = Vec::with_capacity(str.len());
    let mut pos = 0;
    loop {
        let Some(found) = str[pos..].iter().position(|&byte| byte == b'\\') else {
            out.extend_from_slice(&str[pos..]);
            break;
        };
        let mut pos2 = pos + found;
        out.extend_from_slice(&str[pos..pos2]);
        pos2 += 1;
        if pos2 == str.len() {
            return Err("Unterminated escape sequence");
        }
        match str[pos2] {
            b'a' => out.push(0x07),
            b'b' => out.push(0x08),
            b'f' => out.push(0x0c),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b't' => out.push(b'\t'),
            b'v' => out.push(0x0b),
            b'x' => {
                pos2 += 1;
                pos = pos2;
                while pos < str.len() && str[pos].is_ascii_hexdigit() {
                    pos += 1;
                }
                if pos == pos2 {
                    return Err("\\x followed by no hexadecimal digits");
                }
                // `StringToUInt(..., 0, 16)` throws for a value above 4,294,967,295.
                let mut value: u32 = 0;
                for &digit in &str[pos2..pos] {
                    let digit = (digit as char).to_digit(16).unwrap_or(0);
                    value = value
                        .checked_mul(16)
                        .and_then(|value| value.checked_add(digit))
                        .ok_or("Cannot convert string to unsigned int")?;
                }
                out.push(value as u8);
                continue;
            }
            b'0'..=b'7' => {
                pos = pos2;
                let mut c = str[pos] - b'0';
                pos += 1;
                while pos < pos2 + 3 && pos < str.len() && (b'0'..=b'7').contains(&str[pos]) {
                    c = (c << 3) | (str[pos] - b'0');
                    pos += 1;
                }
                out.push(c);
                continue;
            }
            b'\n' => {}
            other => out.push(other),
        }
        pos = pos2 + 1;
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3145-3180
// ```c++
// template <class TStr>
// TStr s_TruncateSpaces(const TStr& str, NStr::ETrunc where,
//                       const TStr& empty_str)
// {
//     SIZE_TYPE length = str.length();
//     if (length == 0) {
//         return empty_str;
//     }
//     SIZE_TYPE beg = 0;
//     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
//         _ASSERT(beg < length);
//         while ( isspace((unsigned char) str[beg]) ) {
//             if (++beg == length) {
//                 return empty_str;
//             }
//         }
//     }
//     SIZE_TYPE end = length;
//     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
//         _ASSERT(beg < end);
//         while (isspace((unsigned char) str[--end])) {
//             if (beg == end) {
//                 return empty_str;
//             }
//         }
//         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
//         ++end;
//     }
// ```
/// `NStr::TruncateSpaces(str)` (`eTrunc_Both`, as `TruncateSpacesInPlace`): `str` without the C
/// `isspace` spaces at its ends.
pub(crate) fn truncate_spaces(str: &[u8]) -> &[u8] {
    let Some(beg) = str.iter().position(|&c| !is_space(c)) else {
        return &[];
    };
    let end = str
        .iter()
        .rposition(|&c| !is_space(c))
        .map_or(beg, |end| end + 1);
    &str[beg..end]
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
                string_to_int(b"abc", 0).map(drop),
                "Cannot convert string 'abc' to Int8",
                861,
                "0",
            ),
            (
                string_to_int(b"-1x", 0).map(drop),
                "Cannot convert string '-1x' to Int8",
                870,
                "2",
            ),
            (
                string_to_int(b"99999999999", 0).map(drop),
                "Cannot convert string '99999999999' to int, overflow",
                640,
                "0",
            ),
            (
                string_to_int(b"99999999999999999999", 0).map(drop),
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
        assert_eq!(string_to_int(b"+3", 0), Ok(3));
        assert_eq!(string_to_int(b"-1", 0), Ok(-1));
        assert_eq!(string_to_int(b"-2147483648", 0), Ok(i32::MIN));
        assert!(string_to_int(b" 3", 0).is_err());
        assert_eq!(string_to_int(b" 3", ALLOW_LEADING_SPACES), Ok(3));
        assert!(string_to_int(b"3 ", ALLOW_LEADING_SPACES).is_err());
        assert!(string_to_int(b"", 0).is_err());
        assert_eq!(string_to_int8(b"-9223372036854775808", 0), Ok(i64::MIN));
        assert_eq!(string_to_int8(b"12\x0034", 0), Ok(12));
        // Base 16 (ncbistr.cpp:684-712, 765-795): `0x` or `0X` is skipped, a letter of either case
        // is a digit, and a `0x` with no digit after it is 0.
        assert_eq!(string_to_int8_base(b"0x10", 0, 16), Ok(16));
        assert_eq!(string_to_int8_base(b"0XfF", 0, 16), Ok(255));
        assert_eq!(string_to_int8_base(b"ff", 0, 16), Ok(255));
        assert_eq!(string_to_int8_base(b"0x", 0, 16), Ok(0));
        assert_eq!(string_to_int8_base(b"-0x10", 0, 16), Ok(-16));
        assert_eq!(
            string_to_int8_base(b"0x7fffffffffffffff", 0, 16),
            Ok(i64::MAX)
        );
        assert_eq!(
            failure(string_to_int8_base(b"0x8000000000000000", 0, 16).unwrap_err()),
            (
                "Cannot convert string '0x8000000000000000' to Int8, overflow".to_string(),
                852,
                "17".to_string()
            )
        );
        assert_eq!(
            failure(string_to_int8_base(b"0x1g", 0, 16).unwrap_err()),
            (
                "Cannot convert string '0x1g' to Int8".to_string(),
                870,
                "3".to_string()
            )
        );
        // fAllowLeadingSpaces | fAllowTrailingSpaces (a taxid list's entry); C `isspace`
        // includes the vertical tab.
        let both = ALLOW_LEADING_SPACES | ALLOW_TRAILING_SPACES;
        assert_eq!(string_to_int(b" \x0b-5\t\r ", both), Ok(-5));
        assert_eq!(string_to_int(b"7 \x00x", both), Ok(7));
        for bad in [&b"1 2"[..], b"9606.0", b"+", b" ", b"0x10", b"3 -"] {
            assert!(string_to_int(bad, both).is_err(), "{bad:?}");
        }
        assert!(string_to_int(b"3 ", ALLOW_TRAILING_SPACES).is_ok());
        assert!(string_to_int(b" 3", ALLOW_TRAILING_SPACES).is_err());
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

    // ncbistr.cpp:3736-3844 (`s_PrintableString`, quoted at `printable_string`): C escapes,
    // `?` escaped only next to another `?`, other bytes outside 0x20..0x7e in octal, as short as
    // the next byte allows (rpsbproc's data loader quotes its lines with it).
    #[test]
    fn printable_string_escapes_as_ncbi() {
        assert_eq!(printable_string(b"a\tb"), "a\\tb");
        assert_eq!(printable_string(b"'\"\\"), "\\'\\\"\\\\");
        assert_eq!(printable_string(b"??x?"), "\\?\\?x?");
        assert_eq!(printable_string(b"\x01x"), "\\1x");
        assert_eq!(printable_string(b"\x017"), "\\0017");
        assert_eq!(printable_string(&[0xc3, 0xa9]), "\\303\\251");
    }

    /// `NStr::IsBlank` is C `isspace` over the whole string: the six spaces of the C locale
    /// (vertical tab and form feed included), and an empty string, are blank; NUL is not.
    #[test]
    fn is_blank_takes_the_spaces_of_c_isspace() {
        assert!(is_blank(b""));
        assert!(is_blank(b" \t\n\x0b\x0c\r"));
        for text in [&b"\0"[..], b" 9 ", b"\x85", b"\xa0"] {
            assert!(!is_blank(text), "{text:?}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892 (`NStr::ParseEscapes`,
    // quoted at `parse_escapes`): the C escapes, `\x` with every hexadecimal digit after it (the
    // value cut to a byte), up to three octal digits (cut to a byte), an escaped line feed that
    // is no byte, any other byte for itself; and the three errors.
    #[test]
    fn parse_escapes_reads_ncbis_escapes() {
        assert_eq!(parse_escapes(b"a\\tb\\n").unwrap(), b"a\tb\n");
        assert_eq!(parse_escapes(b"a\\x41\\x4a4bz").unwrap(), b"aA\x4bz");
        assert_eq!(parse_escapes(b"\\x000000041").unwrap(), b"A");
        assert_eq!(parse_escapes(b"\\101\\1012\\8\\q").unwrap(), b"AA28q");
        assert_eq!(parse_escapes(b"\\777\\0\\08").unwrap(), b"\xff\0\08");
        assert_eq!(
            parse_escapes(b"\\a\\b\\f\\r\\v\\\\\\\"").unwrap(),
            b"\x07\x08\x0c\r\x0b\\\""
        );
        assert_eq!(parse_escapes(b"a\\\nb").unwrap(), b"ab");
        assert_eq!(parse_escapes(b"").unwrap(), b"");
        assert_eq!(parse_escapes(b"a\\"), Err("Unterminated escape sequence"));
        for bad in [&b"a\\xg"[..], b"\\xZZ", b"\\x"] {
            assert_eq!(
                parse_escapes(bad),
                Err("\\x followed by no hexadecimal digits"),
                "{bad:?}"
            );
        }
        assert_eq!(parse_escapes(b"\\xffffffff").unwrap(), b"\xff");
        assert!(parse_escapes(b"\\x100000000").is_err());
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3145-3180 (quoted at
    // `truncate_spaces`).
    #[test]
    fn truncate_spaces_cuts_the_c_spaces_at_both_ends() {
        assert_eq!(truncate_spaces(b""), b"");
        assert_eq!(truncate_spaces(b" \t\n\x0b\x0c\r"), b"");
        assert_eq!(truncate_spaces(b"\x0b a b\r\n"), b"a b");
        assert_eq!(truncate_spaces(b"a"), b"a");
        assert_eq!(truncate_spaces(b"\0 \xa0"), b"\0 \xa0");
    }
}
