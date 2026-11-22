//! The functions of NCBI's `NStr` (`c++/src/corelib/ncbistr.cpp`) that LOSAT's ports call to
//! convert numbers (the argument readers of the BLAST programs, makeblastdb and its FASTA
//! reader, the database readers, blastdbcmd, the Seq-id parser and its accession guide, the
//! taxid lists and rpsbproc's data loader `rpsbproc::data`): base 10, the flags each caller
//! passes, and the `CStringException` each throws, with the place it is thrown from
//! ([`ThrowSite`]); with them `NStr::PrintableString`, which their messages quote the string
//! with, `NStr::StringToUInt8_DataSize` (`-max_file_sz`, module `data_size`) and
//! `NStr::StringToDouble` with glibc's `strtod` and NCBI's own `NStr::StringToDoublePosix`
//! (module `double`: the real arguments, the SEG locut and hicut, rpsbproc's bit score
//! thresholds, the constraints of the arguments and the REAL values of the ASN.1 text and XML
//! readers). Also the string functions that more than one port calls: `NStr::ParseEscapes`
//! (makeblastdb's titles and the registry reader), `NStr::TruncateSpaces` (each `ETrunc`: the
//! FASTA readers, the argument readers, the format guess), `NStr::IsBlank`, `NStr::StringToBool`
//! and C's `isspace`.
//!
//! Not ported: `fConvErr_NoThrow` and the errno bookkeeping of `CS2N_Guard` (a caller that passes
//! `fConvErr_NoThrow` takes NCBI's 0 for an error itself), the bases other than 10 (but 16 of
//! `StringToInt8`, which the argument readers try, `ncbi_app::args` and
//! `blastinput::value_parsers`), and the flags no caller passes (`fMandatorySign`,
//! `fAllowCommas`, `fAllowLeadingSymbols` beyond its spaces).
//! Also here: `NStr::StringToNonNegativeInt`.
//!
//! The parsers of the frozen web ABI v1 (`blastinput::value_parsers::parse_seg_filtering` and
//! `tblastx_v1_outfmt`, `web_api`) keep conversions of their own.

mod data_size;
mod double;
pub(crate) mod utf8;

pub use data_size::string_to_uint8_data_size;
pub(crate) use double::{
    string_to_double, string_to_double_posix, string_to_double_posix_or_local, NCBI_NAN,
};

use std::borrow::Cow;

use crate::ncbi_exception::{NcbiException, ThrowSite};

/// `NStr::TStringToNumFlags`, of the flags LOSAT's callers pass.
pub(crate) type Flags = u32;

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:294-298
/// `NStr::fAllowLeadingSpaces`.
pub(crate) const ALLOW_LEADING_SPACES: Flags = 1 << 19;
/// `NStr::fAllowTrailingSpaces`.
pub(crate) const ALLOW_TRAILING_SPACES: Flags = 1 << 21;
/// `NStr::fAllowTrailingSymbols`: any bytes after the number (white space included).
pub(crate) const ALLOW_TRAILING_SYMBOLS: Flags = (1 << 22) | ALLOW_TRAILING_SPACES;

/// `isspace` in the C locale: the space, `\t`, `\n`, `\v`, `\f` and `\r` (a `const fn` for the
/// character tables built at compile time).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3156 (`s_TruncateSpaces`)
pub(crate) const fn is_space(c: u8) -> bool {
    matches!(c, b'\t'..=b'\r' | b' ')
}

/// `atoi` of glibc, which NCBI's Linux builds call: `strtol` of `s` in base 10 cast to `int`, in
/// the C locale (white space, a sign, the digits up to the first other byte; a value beyond `long`
/// is `LONG_MAX` or `LONG_MIN`), truncated to `int`.
pub(crate) fn c_atoi(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && is_space(s[i]) {
        i += 1;
    }
    let mut negative = false;
    if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        negative = s[i] == b'-';
        i += 1;
    }
    // The magnitude, up to that of LONG_MIN (2^63).
    let mut value: u64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        value = value
            .saturating_mul(10)
            .saturating_add(u64::from(s[i] - b'0'));
        i += 1;
    }
    let long = if negative {
        0i64.wrapping_sub_unsigned(value.min(1 << 63))
    } else {
        value.min(i64::MAX as u64) as i64
    };
    long as i32
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:106-115
/// `NStr::IsBlank`: every byte of `s` is a C `isspace` space (an empty string included).
pub(crate) fn is_blank(s: &[u8]) -> bool {
    s.iter().all(|&c| is_space(c))
}

/// `NStr::StringToNonNegativeInt`: the value, or -1 (the `errno` it sets is not kept; no
/// caller reads it).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:457-491
/// It reads `str.size()` bytes: a NUL is not a digit here (unlike `StringToInt8`).
pub(crate) fn string_to_non_negative_int(s: &[u8]) -> i32 {
    if s.is_empty() {
        return -1;
    }
    let digits = if s[0] == b'+' && s.len() > 1 {
        &s[1..]
    } else {
        s
    };
    let mut value: u32 = 0;
    for &c in digits {
        let d = u32::from(c.wrapping_sub(b'0'));
        if d > 9 {
            return -1;
        }
        let next = value.wrapping_mul(10).wrapping_add(d);
        const OVERFLOW_LIMIT: u32 = (i32::MAX as u32 - 9) / 10 + 1;
        if value >= OVERFLOW_LIMIT && (value > OVERFLOW_LIMIT || next > i32::MAX as u32) {
            return -1;
        }
        value = next;
    }
    value as i32
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

/// `NStr::PrintableString` with its default mode (the flags `fNewLine_Quote` and
/// `fNonAscii_Passthru`, both 0) and language C: what `CS2N_Guard::Message` quotes the input string
/// with, and the form in which `CFastaReader` keeps a definition line (fasta.cpp:600, read back by
/// makeblastdb, `writedb::fasta_source::title`, and written by the archive, `report::archive`).
/// Ported branches: `s_PrintableString` with a `mode` of 0 and a `lang` of `eLanguage_C` (so
/// `fNewLine_Passthru`, `fNonAscii_Quote`, `fPrintable_Full` and `eLanguage_Javascript` are not
/// ported). `isprint` is the C locale's (0x20 to 0x7E), so every byte of the result is printable
/// ASCII; it borrows `str_` when nothing is escaped (NCBI returns a copy of it).
pub(crate) fn printable_string(str_: &[u8]) -> Cow<'_, [u8]> {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3736-3742
    let mut out: Option<Vec<u8>> = None;
    let mut j = 0;
    for (i, &byte) in str_.iter().enumerate() {
        let mut octal = false;
        let mut c = byte;
        match byte {
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3743-3771
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
            b'&' => continue,
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3776-3783
            b'?' => {
                let after_one = i > 0 && str_[i - 1] == b'?';
                let before_one = i + 1 < str_.len() && str_[i + 1] == b'?';
                if !after_one && !before_one {
                    continue;
                }
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3784-3795
            // fNonAscii_Quote is not set; isprint of the "C" locale is 0x20..=0x7e
            0x20..=0x7e => continue,
            _ => octal = true,
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3797-3803
        let out = out.get_or_insert_with(|| Vec::with_capacity(str_.len() + 8));
        out.extend_from_slice(&str_[j..i]);
        out.push(b'\\');
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3804-3831
        // `c == '\n'` cannot be true here: c was changed to 'n' above (no fNewLine_Passthru).
        if octal {
            // fPrintable_Full is not set
            let mut reduce = i == str_.len() - 1 || !(b'0'..=b'7').contains(&str_[i + 1]);
            let v = c >> 6;
            if v != 0 || !reduce {
                out.push(b'0' + v);
                reduce = false;
            }
            let v = (c >> 3) & 7;
            if v != 0 || !reduce {
                out.push(b'0' + v);
            }
            out.push(b'0' + (c & 7));
        } else {
            out.push(c);
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3832-3832
        j = i + 1;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3834-3844
    // `out` is set exactly when `j != 0`; the rest is empty when `i == j`.
    match out {
        Some(mut out) => {
            out.extend_from_slice(&str_[j..]);
            Cow::Owned(out)
        }
        None => Cow::Borrowed(str_),
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:545-558
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:569-575
/// `S2N_CONVERT_ERROR` with `to_type`, `msg`, `errcode` and `pos`, at `site` for the string `str`.
fn convert_error(
    str: &[u8],
    to_type: &str,
    msg: &str,
    site: ThrowSite,
    pos: usize,
) -> NcbiException {
    let mut s = format!(
        "Cannot convert string '{}' to {to_type}",
        String::from_utf8_lossy(&printable_string(str))
    );
    if !msg.is_empty() {
        s.push_str(", ");
        s.push_str(msg);
    }
    NcbiException::with_pos(site, s, pos)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-873
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:740-759
// `CTempString::operator[]` gives NUL past the end (tempstr.hpp:341-347): the loop ends there and
// at an embedded NUL, whose end check (a test of `str[pos]`, lines 603-606) then passes, so
// "12\0x" is 12. A decimal point stops the skipping of trailing spaces, and the end check then
// fails on it, as on any other byte that is not a space. The errors are thrown at lines 852
// (overflow), 861 (no digit) and 870 (`CHECK_ENDPTR`).
const STRING_TO_INT8_OVERFLOW: ThrowSite = convert_site(852, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_INVAL: ThrowSite = convert_site(861, "ncbi::NStr::StringToInt8");
const STRING_TO_INT8_ENDPTR: ThrowSite = convert_site(870, "ncbi::NStr::StringToInt8");

/// `s_SkipAllowedSymbols` in the mode `eSkipSpacesOnly`: the position after the white space
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

/// `NStr::StringToInt8` with `flags` in base 10, `flags` [`ALLOW_LEADING_SPACES`] and
/// [`ALLOW_TRAILING_SPACES`] or [`ALLOW_TRAILING_SYMBOLS`], or 0.
pub(crate) fn string_to_int8(s: &[u8], flags: Flags) -> Result<i64, NcbiException> {
    string_to_int8_base(s, flags, 10)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:765-795
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:684-712
// `pos0` is set before the radix check: a `0x` with no digit after it is 0 (`pos != pos0`).
/// `NStr::StringToInt8` with `flags` and a `base` of 10 or 16 (`s_StringToInt8` of the argument
/// readers tries 16 after 10: `ncbi_app::args`, `blastinput::value_parsers::ncbi_int8`), `flags`
/// as for [`string_to_int8`].
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
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:740-748
    // `fAllowTrailingSymbols` skips to the end (`eSkipAll`), past any byte and any NUL.
    if flags & ALLOW_TRAILING_SYMBOLS != 0 {
        if flags & ALLOW_TRAILING_SYMBOLS == ALLOW_TRAILING_SPACES {
            pos = skip_spaces(s, pos);
        } else {
            pos = s.len();
        }
    }
    let n = if sign { n.wrapping_neg() } else { n };
    if at(pos) != 0 {
        return Err(error(STRING_TO_INT8_ENDPTR, "", pos));
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
const STRING_TO_INT_OVERFLOW: ThrowSite = convert_site(640, "ncbi::NStr::StringToInt");

/// `NStr::StringToInt` with `flags` in base 10, `flags` as for [`string_to_int8`].
pub(crate) fn string_to_int(s: &[u8], flags: Flags) -> Result<i32, NcbiException> {
    let value = string_to_int8(s, flags)?;
    i32::try_from(value).map_err(|_| convert_error(s, "int", "overflow", STRING_TO_INT_OVERFLOW, 0))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:876-918
// With the default flags only the fast conversion runs. Its overflow error has no message
// after the type (`kEmptyStr`, unlike `S2N_CONVERT_ERROR_OVERFLOW` of the other conversions):
// makeblastdb 2.17.0 with `BLASTDB_LMDB_MAP_SIZE=99999999999999999999` prints "Cannot convert
// string '99999999999999999999' to Uint8 (m_Pos = 19)". A byte of 0x80 or more is a negative
// `char`: not a digit either.
const STRING_TO_UINT8_EMPTY: ThrowSite = convert_site(895, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_DIGIT: ThrowSite = convert_site(908, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_OVERFLOW: ThrowSite = convert_site(912, "ncbi::NStr::StringToUInt8");

/// `NStr::StringToUInt8` with the default flags and base (the fast conversion).
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
const STRING_TO_UINT_OVERFLOW: ThrowSite = convert_site(652, "ncbi::NStr::StringToUInt");

/// `NStr::StringToUInt` with the default flags and base.
pub(crate) fn string_to_uint(s: &[u8]) -> Result<u32, NcbiException> {
    let value = string_to_uint8(s)?;
    u32::try_from(value)
        .map_err(|_| convert_error(s, "unsigned int", "overflow", STRING_TO_UINT_OVERFLOW, 0))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:920-982
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:766-795
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:684-713
// `pos0` is taken before the "0x" is skipped, so "0x" alone is 0. The loop ends at the end of
// the string and at an embedded NUL (`CTempString::operator[]`); the size check then fails on
// the NUL, unlike `StringToInt8`'s check of `str[pos]`. NCBI BLAST+ 2.17.0 (line NP M3d, the
// `>` lines of ascii unit counts) throws at lines 964 (overflow), 973 (no digit) and 980
// (`CHECK_ENDPTR_SIZE`).
const STRING_TO_UINT8_BASE_OVERFLOW: ThrowSite = convert_site(964, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_BASE_INVAL: ThrowSite = convert_site(973, "ncbi::NStr::StringToUInt8");
const STRING_TO_UINT8_BASE_ENDPTR: ThrowSite = convert_site(980, "ncbi::NStr::StringToUInt8");

/// `NStr::StringToUInt8` with flags 0 and `base`: `base` 0 (from the string: `0x` hexadecimal, a
/// leading `0` octal, else decimal) or 2 to 36; base 10 is the fast conversion
/// ([`string_to_uint8`]).
pub(crate) fn string_to_uint8_base(s: &[u8], base: u32) -> Result<u64, NcbiException> {
    debug_assert!(base == 0 || (2..=36).contains(&base));
    if base == 10 {
        return string_to_uint8(s);
    }
    let at = |pos: usize| s.get(pos).copied().unwrap_or(0);
    let error = |site: ThrowSite, msg: &str, pos: usize| convert_error(s, "Uint8", msg, site, pos);
    let mut pos = 0;
    if at(pos) == b'+' {
        pos += 1;
    }
    let pos0 = pos;
    let (ch, next) = (at(pos), at(pos + 1));
    let base = match base {
        0 if ch != b'0' => 10,
        0 if next == b'x' || next == b'X' => 16,
        0 => 8,
        base => base,
    };
    if base == 16 && ch == b'0' && (next == b'x' || next == b'X') {
        pos += 2;
    }
    let base = u64::from(base);
    let limdiv = u64::MAX / base;
    let limoff = u64::MAX % base;
    let mut n: u64 = 0;
    loop {
        let ch = at(pos);
        if ch == 0 {
            break;
        }
        let delta = match ch {
            b'0'..=b'9' => u64::from(ch - b'0'),
            b'a'..=b'z' if base > 10 => u64::from(ch - b'a') + 10,
            b'A'..=b'Z' if base > 10 => u64::from(ch - b'A') + 10,
            _ => break,
        };
        if delta >= base {
            break;
        }
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(STRING_TO_UINT8_BASE_OVERFLOW, "overflow", pos));
        }
        n = n * base + delta;
        pos += 1;
    }
    if pos == pos0 {
        return Err(error(STRING_TO_UINT8_BASE_INVAL, "", pos));
    }
    if pos < s.len() {
        return Err(error(STRING_TO_UINT8_BASE_ENDPTR, "", pos));
    }
    Ok(n)
}

/// `NStr::StringToUInt` with flags 0 and `base` (`base` as for [`string_to_uint8_base`]).
pub(crate) fn string_to_uint_base(s: &[u8], base: u32) -> Result<u32, NcbiException> {
    let value = string_to_uint8_base(s, base)?;
    u32::try_from(value)
        .map_err(|_| convert_error(s, "unsigned int", "overflow", STRING_TO_UINT_OVERFLOW, 0))
}

/// `NStr::ParseEscapes` with its default mode (`eEscSeqRange_Standard`), as makeblastdb's
/// title (`writedb::fasta_source::title`) and NCBI's registry reader
/// (`blastinput::ncbi_registry`) call it. An error is the message of NCBI's `CStringException`
/// (an escape that ends the string, `\x` without a hexadecimal digit, or one whose value
/// `StringToUInt` cannot hold).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892
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
                // `StringToUInt` (flags 0, base 16) throws for a value above 4,294,967,295.
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
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3187-3195
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3182-3185
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3197-3235 (`TruncateSpacesInPlace`,
// which takes a string and an `ETrunc` `where`: the loops of `s_TruncateSpaces` over `str.data()`,
// then) `NStr::TruncateSpaces`, `TruncateSpaces_Unsafe` and the two `TruncateSpacesInPlace` cut the
// same bytes; the four functions below are its `where`, as the slice of `str` that is left (empty
// when every byte is a space), `truncate_spaces_str` for a `string`.
/// `NStr::TruncateSpaces` (`eTrunc_Both`, as `TruncateSpacesInPlace`): `str` without the C
/// `isspace` spaces at its ends.
pub(crate) fn truncate_spaces(str: &[u8]) -> &[u8] {
    truncate_spaces_end(truncate_spaces_begin(str))
}

/// `NStr::TruncateSpaces` with `NStr::eTrunc_Begin`: `str` without the C `isspace` spaces at its
/// start.
pub(crate) fn truncate_spaces_begin(str: &[u8]) -> &[u8] {
    let beg = str.iter().position(|&c| !is_space(c)).unwrap_or(str.len());
    &str[beg..]
}

/// `NStr::TruncateSpaces` with `NStr::eTrunc_End`: `str` without the C `isspace` spaces at its
/// end.
pub(crate) fn truncate_spaces_end(str: &[u8]) -> &[u8] {
    let end = str
        .iter()
        .rposition(|&c| !is_space(c))
        .map_or(0, |end| end + 1);
    &str[..end]
}

/// [`truncate_spaces`] of a `str` (`NStr::TruncateSpaces` of a `string`): the spaces are ASCII
/// bytes, so what is left starts and ends at characters.
pub(crate) fn truncate_spaces_str(str: &str) -> &str {
    let beg = str.len() - truncate_spaces_begin(str.as_bytes()).len();
    let end = truncate_spaces_end(str.as_bytes()).len().max(beg);
    &str[beg..end]
}

/// `NStr::StringToBool`: `None` where NCBI throws.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2817-2839
pub(crate) fn string_to_bool(value: &str) -> Option<bool> {
    let is = |word: &str| value.eq_ignore_ascii_case(word);
    if value == "1" || is("true") || is("t") || is("yes") || is("y") || is("on") {
        Some(true)
    } else if value == "0" || is("false") || is("f") || is("no") || is("n") || is("off") {
        Some(false)
    } else {
        None
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1952-1993
// Each chunk after the first takes 11 characters, while its 9 digits need 3 commas once the
// count of the chunk before has run on (the comment's "2 commas" holds for the first chunk
// only): from the second chunk on, a chunk loses its top digit or digits, so a value of
// 4,294,967,296,000,000,000 (2^32 * 10^9) or more prints with digits lost (NCBI BLAST+
// 2.17.0: LENGTH 4294967296000000000 prints `494,967,296,000,000,000`; line RG D1a-3, probe
// p5_pr_L4295e15). The loop of the last characters is reproduced as it runs: it ends where the
// chunk's 11 characters end.
/// NCBI's `s_PrintBase10Uint8` with `NStr::fWithCommas`: `value` in decimal with commas, as
/// the text of `NStr::UInt8ToString` (and of `Int8ToString` after the sign).
fn print_base10_uint8_with_commas(mut value: u64) -> Vec<u8> {
    const PRINT_INT8_CHUNK: u64 = 1_000_000_000;
    const PRINT_INT8_CHUNK_SIZE: usize = 9;
    // NCBI's buffer is 64 characters (`CHAR_BIT` times the size of `value`); the text fills it from
    // the end.
    let mut buffer = [0u8; 64];
    let mut pos = buffer.len();
    let mut cnt: i32 = -1;
    let mut put = |pos: &mut usize, digit: u8, cnt: &mut i32| {
        *cnt += 1;
        if *cnt == 3 {
            *pos -= 1;
            buffer[*pos] = b',';
            *cnt = 0;
        }
        *pos -= 1;
        buffer[*pos] = digit;
    };
    while value & !u64::from(u32::MAX) != 0 {
        let mut chunk = value as u32;
        value /= PRINT_INT8_CHUNK;
        chunk = chunk.wrapping_sub((PRINT_INT8_CHUNK as u32).wrapping_mul(value as u32));
        let end = pos - PRINT_INT8_CHUNK_SIZE - 2;
        loop {
            put(&mut pos, b'0' + (chunk % 10) as u8, &mut cnt);
            chunk /= 10;
            // A step writes one or two characters; the count of each chunk's start makes the
            // 11 characters end on a digit, so `pos` meets `end` (an overshoot cannot occur).
            if pos <= end {
                break;
            }
        }
    }
    let mut chunk = value as u32;
    loop {
        put(&mut pos, b'0' + (chunk % 10) as u8, &mut cnt);
        chunk /= 10;
        if chunk == 0 {
            break;
        }
    }
    buffer[pos..].to_vec()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2055-2073
/// `NStr::UInt8ToString` with `NStr::fWithCommas`: `n` in decimal with a comma between groups of
/// three digits, with NCBI's loss of digits from 4,294,967,296,000,000,000
/// (`print_base10_uint8_with_commas`).
pub(crate) fn with_commas(n: u64) -> String {
    String::from_utf8(print_base10_uint8_with_commas(n)).expect("digits and commas")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2031-2048
// `-svalue` of INT64_MIN overflows; GCC's wrapped result, cast to `Uint8`, is 2^63, its magnitude
// (`unsigned_abs`), so the text is the value's under either reading (probe p5_neg_2p63 of line
// RG D1a-3: `-923,372,036,854,775,808`).
/// `NStr::Int8ToString` with `NStr::fWithCommas`: the sign, then the magnitude as
/// [`with_commas`] prints it.
pub(crate) fn int8_with_commas(n: i64) -> String {
    let digits = print_base10_uint8_with_commas(n.unsigned_abs());
    let mut text = String::with_capacity(digits.len() + 1);
    if n < 0 {
        text.push('-');
    }
    text.push_str(std::str::from_utf8(&digits).expect("digits and commas"));
    text
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1825-1862, 1867-1878
// The `int` is widened to `long` before it is negated: INT_MIN prints whole
// (`-2,147,483,648`, probe p5_pr_N2p31 of line RG D1a-3), and every magnitude of an `int` is
// below 2^32, where the digits are whole.
/// `NStr::IntToString` with `NStr::fWithCommas`: the sign, then the magnitude with commas.
pub(crate) fn int_with_commas(n: i32) -> String {
    let digits = print_base10_uint8_with_commas(u64::from(n.unsigned_abs()));
    let mut text = String::with_capacity(digits.len() + 1);
    if n < 0 {
        text.push('-');
    }
    text.push_str(std::str::from_utf8(&digits).expect("digits and commas"));
    text
}

#[cfg(test)]
#[path = "ncbi_str_tests.rs"]
mod tests;
