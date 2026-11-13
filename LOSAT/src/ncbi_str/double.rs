//! `NStr::StringToDouble(str)` with its default flags (corelib/ncbistr.cpp), which rpsbproc's
//! data loader calls (`NStr::StringToNumeric<double>`, ncbistr.hpp:3632-3637) for the bit score
//! thresholds of `bitscore_specific.txt`, and glibc's `strtod` that it reads the number with.
//!
//! Also `NStr::StringToDouble(str, NStr::fDecimalPosixOrLocal)`, which NCBI's argument reader
//! calls for an `eDouble` key (`CArg_Double`: segmasker's `-locut` and `-hicut`,
//! `ncbi_app::args`): [`string_to_double_posix_or_local`].
//!
//! Only these two flag sets are ported: the leading and trailing symbols, `fMandatorySign`,
//! `fDecimalPosix`, `fDecimalPosixFinite` and `fConvErr_NoThrow` are not (no caller passes them).
//! The decimal point of `strtod` is that of the C locale, `.` (`s_IsDecimalPoint` with no flag
//! asks `localeconv`; NCBI's applications keep the C locale).

use super::{convert_error, convert_site};
use crate::ncbi_exception::{NcbiException, ThrowSite};

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1312-1319
// ```c++
//     // For consistency make additional check on incorrect leading symbols.
//     // Because strtod() may just skip such symbols.
//     if (!(flags & NStr::fAllowLeadingSymbols)) {
//         char c = str[pos];
//         if ( !isdigit((unsigned char)c)  &&  !s_IsDecimalPoint(c,flags)  &&  c != '-'  &&  c != '+') {
//             S2N_CONVERT_ERROR_INVAL(double);
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1342-1349
// ```c++
//     if ( !endptr  ||  endptr == begptr ) {
//         S2N_CONVERT_ERROR(double, kEmptyStr, EINVAL, s_DiffPtr(endptr, begptr) + pos);
//     }
//     // some libs set ERANGE, others do not
//     // here, we do not consider ERANGE as error
//     if ( errno_ref && errno_ref != ERANGE ) {
//         S2N_CONVERT_ERROR(double, kEmptyStr, errno_ref, s_DiffPtr(endptr, begptr) + pos);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1367-1375
// ```c++
//     pos += s_DiffPtr(endptr, begptr);
//
//     // Skip allowed trailing symbols
//     if (flags & NStr::fAllowTrailingSymbols) {
//         ...
//     }
//     CHECK_ENDPTR(double);
//     return n;
// ```
// The errors are thrown in `s_StringToDouble` (`double ncbi::s_StringToDouble(const char*,
// size_t, ...)`, whose class part is `ncbi` for `CDiagCompileInfo::ParseCurrFunctName`,
// ncbidiag.cpp:826-888) at lines 1317 (the first byte), 1343 (no number) and 1374
// (`CHECK_ENDPTR`); glibc's `strtod` sets no errno other than `ERANGE`, so line 1348 is not
// reached.
const STRING_TO_DOUBLE_INVAL: ThrowSite = convert_site(1317, "ncbi::s_StringToDouble");
const STRING_TO_DOUBLE_NO_NUMBER: ThrowSite = convert_site(1343, "ncbi::s_StringToDouble");
const STRING_TO_DOUBLE_ENDPTR: ThrowSite = convert_site(1374, "ncbi::s_StringToDouble");

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1386-1403
// ```c++
// double NStr::StringToDouble(const CTempStringEx str, TStringToNumFlags flags)
// {
//     size_t size = str.size();
//     if ( str.HasZeroAtEnd() ) {
//         // string has zero at the end already
//         return s_StringToDouble(str.data(), size, flags);
//     }
//     char buf[256]; // small temporary buffer on stack for appending zero char
//     if ( size < sizeof(buf) ) {
//         memcpy(buf, str.data(), size);
//         buf[size] = '\0';
//         return s_StringToDouble(buf, size, flags);
//     }
//     else {
//         // use std::string() to allocate memory for appending zero char
//         return s_StringToDouble(string(str).c_str(), size, flags);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1283-1285
// ```c++
// static double s_StringToDouble(const char* str, size_t size,
//                                NStr::TStringToNumFlags flags)
// {
// ```
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:422-433
// ```c++
// CTempString::CTempString(const char* str)
// ...
//     m_String = str;
//     m_Length = strlen(str);
// ```
// `str` of `s_StringToDouble` is the C string, so `strtod` and the end check (`str[pos]`) stop
// at its first 0 byte ("12\0x" is 12), and the message of an error, `Message(str, ...)`, quotes
// the text before that byte (the other conversions quote the whole `CTempString`).
/// `NStr::StringToDouble(str)` with the default flags: the first byte must be a digit, the
/// decimal point or a sign; `strtod` reads the number, which must end the C string.
pub(crate) fn string_to_double(s: &[u8]) -> Result<f64, NcbiException> {
    let c_string = &s[..s.iter().position(|&b| b == 0).unwrap_or(s.len())];
    let error = |site: ThrowSite, pos: usize| convert_error(c_string, "double", "", site, pos);
    let at = |pos: usize| c_string.get(pos).copied().unwrap_or(0);
    let c = at(0);
    if !c.is_ascii_digit() && c != b'.' && c != b'-' && c != b'+' {
        return Err(error(STRING_TO_DOUBLE_INVAL, 0));
    }
    let (value, used) = strtod_prefix(c_string);
    if used == 0 {
        return Err(error(STRING_TO_DOUBLE_NO_NUMBER, 0));
    }
    if at(used) != 0 {
        return Err(error(STRING_TO_DOUBLE_ENDPTR, used));
    }
    Ok(value)
}

/// Why [`string_to_double_posix_or_local`] gave no number.
#[derive(Debug, PartialEq)]
pub(crate) enum PosixOrLocalError {
    /// NCBI's `CStringException` (`eConvert`).
    Ncbi(NcbiException),
    /// A number that NCBI computes in the x87 extended precision of `long double` in a way that
    /// LOSAT does not reproduce (see `posix_read_past_strtod`).
    ExtendedPrecision,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:724-737
// ```c++
// bool s_IsDecimalPoint(unsigned char ch, NStr::TStringToNumFlags  flags)
// {
//     if ( ch != '.' && ch != ',') {
//         return false;
//     }
//     if (flags & NStr::fDecimalPosix) {
//         return ch == '.';
//     }
//     else if (flags & NStr::fDecimalPosixOrLocal) {
//         return ch == '.' || ch == ',';
//     }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1331-1341
// ```c++
//     if (flags & NStr::fDecimalPosix) {
//         n = NStr::StringToDoublePosix(begptr, &endptr, flags);
//     } else {
//         n = strtod(begptr, &endptr);
//     }
//     if (flags & NStr::fDecimalPosixOrLocal) {
//         char* endptr2 = 0;
//         double n2 = NStr::StringToDoublePosix(begptr, &endptr2, flags);
//         if (!endptr || (endptr2 && endptr2 > endptr)) {
//             n = n2;
//             endptr = endptr2;
//         }
//     }
// ```
/// `NStr::StringToDouble(str, NStr::fDecimalPosixOrLocal)`: as [`string_to_double`], but a comma
/// passes the check of the first byte (no number then starts there, and the error is that of
/// no number), and NCBI's own reader `StringToDoublePosix` reads the text too; its number is
/// taken where it reads further than `strtod`.
pub(crate) fn string_to_double_posix_or_local(s: &[u8]) -> Result<f64, PosixOrLocalError> {
    let c_string = &s[..s.iter().position(|&b| b == 0).unwrap_or(s.len())];
    let error = |site: ThrowSite, pos: usize| {
        PosixOrLocalError::Ncbi(convert_error(c_string, "double", "", site, pos))
    };
    let at = |pos: usize| c_string.get(pos).copied().unwrap_or(0);
    let c = at(0);
    if !c.is_ascii_digit() && !matches!(c, b'.' | b',' | b'-' | b'+') {
        return Err(error(STRING_TO_DOUBLE_INVAL, 0));
    }
    let (mut value, mut used) = strtod_prefix(c_string);
    if let Some((posix, posix_used)) = posix_read_past_strtod(c_string) {
        if posix_used > used {
            value = posix?;
            used = posix_used;
        }
    }
    if used == 0 {
        return Err(error(STRING_TO_DOUBLE_NO_NUMBER, 0));
    }
    if at(used) != 0 {
        return Err(error(STRING_TO_DOUBLE_ENDPTR, used));
    }
    Ok(value)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1027-1069 (`StringToDoublePosix`, the
// mantissa)
// ```c++
//     for ( ; ; c = *ptr++ ) {
//         if (c >= '0' && c <= '9') {
//             // digits: accumulate
//             c = (char)(c - '0');
//             anydigits = true;
//             ++digits;
//             if (first == 0) {
//                 first = c;
//                 if ( first == 0 ) {
//                     // omit leading zeros
//                     --digits;
//                     if (dot) {
//                         --dot_position;
//                     }
//                 }
//             } else if (digits <= 9) {
//                 // first 9 digits come to 'first'
//                 first = first*10 + c;
//             } else if (digits <= 18) {
//                 // next 9 digits come to 'second'
//                 first_mul *= 10;
//                 second = second*10 + c;
//             } else {
//                 // other digits come to 'third'
//                 second_mul *= NCBI_CONST_LONGDOUBLE(10.);
//                 third = third * NCBI_CONST_LONGDOUBLE(10.) + c;
//             }
//         }
//         else if (c == '.') {
//             // dot
//             // if second dot, stop
//             if (dot) {
//                 --ptr;
//                 break;
//             }
//             dot_position = digits;
//             dot = true;
//         }
//         else if (c == 'e' || c == 'E') {
//             // if exponent, stop
//             if (!anydigits) {
//                 --ptr;
//                 break;
//             }
//             expn = true;
//             break;
//         }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1113-1116, 1142-1153, 1158-1169
// ```c++
//     int exponent = dot ? dot_position - digits : 0;
//
//     // read exponent
//     if (expn && *ptr) {
//     ...
//         // if no digits, rollback
//         if (!expdigits) {
//             // rollback sign
//             if (expsign) {
//                 --ptr;
//             }
//             // rollback exponent
//             if (expn) {
//                 --ptr;
//             }
//         }
//     ...
//     long double ret;
//     if ( first_mul > 1 ) {
//         _ASSERT(first);
//         ret = ((long double)first * first_mul + second)* second_mul + third;
//     }
//     else {
//         ...
//         ret = first;
//     }
//     // calculate exponent
//     if ( first && exponent ) {
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1218-1237
// ```c++
//                 exponent = -exponent;
//                 static const long double mul1[16] = {
//                     NCBI_CONST_LONGDOUBLE(1e-0),
//                     NCBI_CONST_LONGDOUBLE(1e-1),
//                     ...
//                 };
//                 ret *= mul1[exponent&15];
// ```
/// Where `StringToDoublePosix` reads further than glibc's `strtod`: a decimal mantissa (at least
/// one digit, at most one point, a sign before it) followed by `e` or `E` that ends the text.
/// `StringToDoublePosix` takes the mark into its number when nothing follows it (`*ptr` is 0),
/// and steps back to it otherwise (an exponent with digits `strtod` reads as well); its other
/// readings (`inf`, `nan`, a hexadecimal `0x`) end where `strtod`'s do or before. The number and
/// the end of `StringToDoublePosix` for those texts, `None` for the others.
///
/// The number is `long double` arithmetic on x87 (the 64-bit mantissa of GCC's x86-64
/// `long double`, rounded to nearest, then to the double): the digits as an integer, times the
/// constant `1e-k` for `k` digits after the point. LOSAT computes it exactly for at most 18
/// digits (the integer is exact) and at most 15 places (one multiplication); a longer mantissa
/// (`third`, more roundings) or more places (`mul2`) is `ExtendedPrecision`.
fn posix_read_past_strtod(text: &[u8]) -> Option<(Result<f64, PosixOrLocalError>, usize)> {
    let at = |pos: usize| text.get(pos).copied().unwrap_or(0);
    let negative = at(0) == b'-';
    let mut pos = usize::from(matches!(at(0), b'+' | b'-'));
    let (mut dot, mut anydigits) = (false, false);
    let (mut digits, mut dot_position) = (0i32, 0i32);
    // `first * first_mul + second` while it holds every digit (at most 18 after the leading
    // zeros); `None` once a digit goes to `third`.
    let mut mantissa: Option<u64> = Some(0);
    loop {
        let c = at(pos);
        if c.is_ascii_digit() {
            let digit = u64::from(c - b'0');
            anydigits = true;
            digits += 1;
            if mantissa == Some(0) {
                mantissa = Some(digit);
                if digit == 0 {
                    digits -= 1;
                    if dot {
                        dot_position -= 1;
                    }
                }
            } else if digits <= 18 {
                mantissa = mantissa.map(|m| m * 10 + digit);
            } else {
                mantissa = None;
            }
        } else if c == b'.' {
            if dot {
                return None;
            }
            dot_position = digits;
            dot = true;
        } else if matches!(c, b'e' | b'E') && anydigits && at(pos + 1) == 0 {
            break;
        } else {
            return None;
        }
        pos += 1;
    }
    let exponent = if dot { dot_position - digits } else { 0 };
    let value = match mantissa {
        Some(0) => Ok(0.0),
        // The integer is exact in the 64-bit mantissa, then rounded once to the double.
        Some(m) if exponent == 0 => Ok(m as f64),
        Some(m) if (-15..0).contains(&exponent) => Ok(x87_scaled_down(m, exponent.unsigned_abs())),
        _ => Err(PosixOrLocalError::ExtendedPrecision),
    };
    let value = value.map(|v| if negative { -v } else { v });
    Some((value, pos + 1))
}

/// `(double)((long double)m * 1e-kL)` on x87, for `m < 10^18` and `1 <= k <= 15`: the constant
/// `1e-kL` rounded to 64 bits (GCC rounds its constants correctly), the product rounded to 64
/// bits, then to 53; to nearest, ties to even, each time.
fn x87_scaled_down(m: u64, k: u32) -> f64 {
    let power = 10u128.pow(k);
    let bits = 128 - power.leading_zeros();
    // 1e-k = power_mantissa * 2^-(63 + bits), power_mantissa in [2^63, 2^64]
    let scaled = 1u128 << (63 + bits);
    let (quotient, remainder) = (scaled / power, scaled % power);
    // No tie: 5^k divides no power of two.
    let power_mantissa = quotient + u128::from(2 * remainder > power);
    let product = u128::from(m) * power_mantissa;
    let product_bits = 128 - product.leading_zeros();
    let mut exponent = -i64::from(63 + bits);
    let mut kept = product;
    if product_bits > 64 {
        let shift = product_bits - 64;
        let rest = product & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        kept = product >> shift;
        if rest > half || (rest == half && kept & 1 == 1) {
            kept += 1;
        }
        exponent += i64::from(shift);
    }
    round_binary(kept, exponent, false)
}

/// glibc `strtod` of the start of `text` (no leading white space here) in the C locale and the
/// default rounding: the value and the bytes read (0 when no number starts there).
fn strtod_prefix(text: &[u8]) -> (f64, usize) {
    let at = |pos: usize| text.get(pos).copied().unwrap_or(0);
    let mut pos = usize::from(matches!(at(0), b'+' | b'-'));
    let negative = at(0) == b'-';
    let lower = |pos: usize| at(pos).to_ascii_lowercase();
    if lower(pos) == b'0' && lower(pos + 1) == b'x' {
        if let Some((value, end)) = hex_prefix(text, pos + 2) {
            return (if negative { -value } else { value }, end);
        }
        // `0x` without a hexadecimal digit reads as the number 0, which ends before the `x`.
        return (if negative { -0.0 } else { 0.0 }, pos + 1);
    }
    if lower(pos) == b'i' && lower(pos + 1) == b'n' && lower(pos + 2) == b'f' {
        let mut end = pos + 3;
        if (0..5).all(|k| lower(end + k) == b"inity"[k]) {
            end += 5;
        }
        return (
            if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            end,
        );
    }
    if lower(pos) == b'n' && lower(pos + 1) == b'a' && lower(pos + 2) == b'n' {
        let mut end = pos + 3;
        let mut payload = 0;
        if at(end) == b'(' {
            let mut k = end + 1;
            while at(k).is_ascii_alphanumeric() || at(k) == b'_' {
                k += 1;
            }
            if at(k) == b')' {
                payload = nan_payload(&text[end + 1..k]);
                end = k + 1;
            }
        }
        // glibc keeps the sign of a NaN (`-nan` is printed back as `-nan`).
        let nan = f64::from_bits(f64::NAN.to_bits() | payload);
        return (if negative { -nan } else { nan }, end);
    }
    let mut digits = 0;
    while at(pos).is_ascii_digit() {
        pos += 1;
        digits += 1;
    }
    if at(pos) == b'.' {
        pos += 1;
        while at(pos).is_ascii_digit() {
            pos += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return (0.0, 0);
    }
    let mantissa_end = pos;
    if matches!(at(pos), b'e' | b'E') {
        let mut k = pos + 1;
        if matches!(at(k), b'+' | b'-') {
            k += 1;
        }
        if at(k).is_ascii_digit() {
            while at(k).is_ascii_digit() {
                k += 1;
            }
            pos = k;
        }
    }
    let number = std::str::from_utf8(&text[..pos]).expect("an ASCII number");
    let value: f64 = number
        .parse()
        .ok()
        .or_else(|| number[..mantissa_end].parse().ok())
        .expect("a decimal number parses");
    (value, pos)
}

/// The payload of glibc's `nan(n-char-sequence)` (`__strtod_nan`): the sequence read by
/// `strtoull` in base 0 (`0x` hexadecimal, a leading `0` octal, else decimal; the largest value
/// on overflow), its low 51 bits put below the quiet bit, when the whole sequence is that number;
/// 0 otherwise.
fn nan_payload(sequence: &[u8]) -> u64 {
    let (digits, radix) = match sequence {
        [b'0', b'x' | b'X', rest @ ..] if !rest.is_empty() => (rest, 16),
        [b'0', rest @ ..] if !rest.is_empty() => (rest, 8),
        _ => (sequence, 10),
    };
    if digits.is_empty() {
        return 0;
    }
    let mut value: u64 = 0;
    for &byte in digits {
        let Some(digit) = char::from(byte).to_digit(radix) else {
            return 0;
        };
        value = value
            .checked_mul(u64::from(radix))
            .and_then(|value| value.checked_add(u64::from(digit)))
            .unwrap_or(u64::MAX);
    }
    value & ((1 << 51) - 1)
}

/// The hexadecimal number of glibc's `strtod` after its `0x`, from `start`: hexadecimal digits
/// with an optional point, at least one digit, then an optional binary exponent `p` with an
/// optional sign and decimal digits (a `p` without digits is not read); the value rounded once
/// to the nearest double (ties to even), the infinity above the largest. `None` when no digit
/// follows (`strtod` then reads the `0` alone). Its absolute value.
// (The hexadecimal form of C's `strtod`, ISO C11 7.22.1.3, as glibc reads it in the C locale;
// the unit test `hex_numbers_read_as_glibcs_strtod` holds glibc 2.39's values and ends for these
// forms.)
fn hex_prefix(text: &[u8], start: usize) -> Option<(f64, usize)> {
    let at = |pos: usize| text.get(pos).copied().unwrap_or(0);
    let digit = |pos: usize| char::from(at(pos)).to_digit(16);
    // The digits as an integer of at most 116 bits, the binary exponent of its last bit, and
    // whether a digit that did not fit was not 0.
    let (mut mantissa, mut exponent, mut sticky) = (0u128, 0i64, false);
    let mut push = |value: u32, after_point: bool| {
        if mantissa == 0 && value == 0 {
            exponent -= if after_point { 4 } else { 0 };
        } else if mantissa < 1 << 112 {
            mantissa = mantissa * 16 + u128::from(value);
            exponent -= if after_point { 4 } else { 0 };
        } else {
            sticky |= value != 0;
            exponent += if after_point { 0 } else { 4 };
        }
    };
    let mut pos = start;
    let mut any = false;
    while let Some(value) = digit(pos) {
        push(value, false);
        any = true;
        pos += 1;
    }
    if at(pos) == b'.' && (any || digit(pos + 1).is_some()) {
        pos += 1;
        while let Some(value) = digit(pos) {
            push(value, true);
            any = true;
            pos += 1;
        }
    }
    if !any {
        return None;
    }
    if matches!(at(pos), b'p' | b'P') {
        let mut k = pos + 1;
        let negative = at(k) == b'-';
        if matches!(at(k), b'+' | b'-') {
            k += 1;
        }
        if at(k).is_ascii_digit() {
            let mut value = 0i64;
            while at(k).is_ascii_digit() {
                value = (value * 10 + i64::from(at(k) - b'0')).min(1 << 40);
                k += 1;
            }
            exponent += if negative { -value } else { value };
            pos = k;
        }
    }
    Some((round_binary(mantissa, exponent, sticky), pos))
}

/// `mantissa * 2^exponent` (plus less than its last bit when `sticky`) rounded to the nearest
/// double, ties to even: subnormal below 2^-1022, the infinity from 2^1024 on.
fn round_binary(mantissa: u128, exponent: i64, sticky: bool) -> f64 {
    if mantissa == 0 {
        return 0.0;
    }
    let bits = i64::from(128 - mantissa.leading_zeros());
    // The exponent of the leading bit, and the bits the double keeps from it.
    let lead = exponent + bits - 1;
    if lead > 1023 {
        return f64::INFINITY;
    }
    let keep = if lead >= -1022 {
        53
    } else {
        53 - (-1022 - lead)
    };
    if keep < 0 {
        return 0.0;
    }
    let shift = bits - keep;
    let mut kept = if shift <= 0 {
        mantissa << (-shift)
    } else {
        let kept = mantissa >> shift;
        let rest = mantissa & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        let up = rest > half || (rest == half && (sticky || kept & 1 == 1));
        kept + u128::from(up)
    };
    let mut lead = lead;
    if keep == 53 && kept == 1 << 53 {
        kept >>= 1;
        lead += 1;
        if lead > 1023 {
            return f64::INFINITY;
        }
    }
    if keep == 53 {
        f64::from_bits(((lead + 1023) as u64) << 52 | (kept as u64 & ((1 << 52) - 1)))
    } else {
        // A subnormal: its bits are the kept ones (2^52 is the smallest normal number).
        f64::from_bits(kept as u64)
    }
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

    // NCBI's rpsbproc printed `Cannot convert string 'abc' to double (m_Pos = 0)` for a bit score
    // of `bitscore_specific.txt` (PS16a probe 4), read `0x10` as 16 and `0X1p3` as 8 and refused
    // `0x` at position 1 (PS16a review oracle A4-1 cases `bs12`, `bs14`, `bs13`).
    #[test]
    fn values_and_errors_are_ncbis() {
        assert_eq!(string_to_double(b"55.5"), Ok(55.5));
        assert_eq!(string_to_double(b"-1e3"), Ok(-1000.0));
        assert_eq!(string_to_double(b"0x10"), Ok(16.0));
        assert_eq!(string_to_double(b"0X1p3"), Ok(8.0));
        assert_eq!(string_to_double(b"12\x00x"), Ok(12.0));
        let cases: [(&[u8], &str, u32, &str); 7] = [
            (b"abc", "Cannot convert string 'abc' to double", 1317, "0"),
            (b"", "Cannot convert string '' to double", 1317, "0"),
            (b"-x", "Cannot convert string '-x' to double", 1343, "0"),
            (
                b"12.5\r",
                "Cannot convert string '12.5\\r' to double",
                1374,
                "4",
            ),
            (b"1e", "Cannot convert string '1e' to double", 1374, "1"),
            (b"0x", "Cannot convert string '0x' to double", 1374, "1"),
            // the message quotes the C string, up to its first 0 byte
            (
                b"1x\x00y",
                "Cannot convert string '1x' to double",
                1374,
                "1",
            ),
        ];
        for (text, message, line, pos) in cases {
            assert_eq!(
                failure(string_to_double(text).unwrap_err()),
                (message.to_string(), line, pos.to_string()),
                "{text:?}"
            );
        }
        assert_eq!(
            string_to_double(b"abc").unwrap_err().site(),
            ThrowSite {
                file: "src/corelib/ncbistr.cpp",
                line: 1317,
                kind: "CStringException::eConvert",
                function: Some("ncbi::s_StringToDouble"),
            }
        );
    }

    // `LOSAT/tests/unit/string_to_double_posix_or_local_vectors.tsv`: NCBI's own
    // `NStr::StringToDouble(text, fDecimalPosixOrLocal)` (its corelib linked to a probe), for
    // the texts that end in an exponent mark (`StringToDoublePosix` reads them past `strtod`,
    // in x87 extended precision), the edges of glibc's `strtod` that `CArg_Double` reads, and
    // the errors. LOSAT refuses only the texts whose mark follows more than 18 digits or 15
    // places (`ExtendedPrecision`).
    #[test]
    fn posix_or_local_values_are_ncbis() {
        let table = include_str!("../../tests/unit/string_to_double_posix_or_local_vectors.tsv");
        let (mut rows, mut refused) = (0, Vec::new());
        for row in table.lines().filter(|row| !row.starts_with('#')) {
            let (text, bits) = row.rsplit_once('\t').unwrap();
            match (string_to_double_posix_or_local(text.as_bytes()), bits) {
                (Err(PosixOrLocalError::Ncbi(_)), "ERR") => {}
                (Ok(value), bits) if bits != "ERR" => {
                    assert_eq!(format!("{:016x}", value.to_bits()), bits, "{text:?}")
                }
                (Err(PosixOrLocalError::ExtendedPrecision), bits) if bits != "ERR" => {
                    refused.push(text)
                }
                (other, bits) => panic!("{text:?}: {other:?}, NCBI {bits}"),
            }
            rows += 1;
        }
        assert_eq!(rows, 422);
        assert_eq!(
            refused,
            [
                "0.0000000000000001e",
                "1234567890123456789e",
                "1.234567890123456789e",
                "9.99999999999999999e"
            ]
        );
    }

    // `LOSAT/tests/unit/strtod_hex_vectors.tsv`: glibc 2.39's `strtod` of hexadecimal texts
    // (edge cases of rounding, subnormals, overflow and the ends of the form, then random ones),
    // made with `ctypes` by `psiblast-20261025/PR/tmp/impl16a/strtodhex.py`.
    #[test]
    fn hex_numbers_read_as_glibcs_strtod() {
        let table = include_str!("../../tests/unit/strtod_hex_vectors.tsv");
        let mut rows = 0;
        for row in table.lines().filter(|row| !row.starts_with('#')) {
            let mut fields = row.split('\t');
            let (text, bits, used) = (
                fields.next().unwrap(),
                fields.next().unwrap(),
                fields.next().unwrap(),
            );
            let (value, end) = strtod_prefix(text.as_bytes());
            assert_eq!(
                (format!("{:016x}", value.to_bits()), end.to_string()),
                (bits.to_string(), used.to_string()),
                "{text:?}"
            );
            rows += 1;
        }
        assert_eq!(rows, 360);
    }
}
