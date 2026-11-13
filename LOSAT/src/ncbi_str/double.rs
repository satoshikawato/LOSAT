//! `NStr::StringToDouble(str)` with its default flags (corelib/ncbistr.cpp), which rpsbproc's
//! data loader calls (`NStr::StringToNumeric<double>`, ncbistr.hpp:3632-3637) for the bit score
//! thresholds of `bitscore_specific.txt`, the BLAST programs' argument readers for the SEG locut
//! and hicut (`blastinput::app::read_seg_option`) and the constraints of the arguments
//! (`blastinput::value_parsers`), and glibc's `strtod` that it reads the number with.
//!
//! Also `NStr::StringToDouble(str, NStr::fDecimalPosixOrLocal)`, which NCBI's argument reader
//! calls for an `eDouble` key (`CArg_Double`: every real argument of the BLAST programs and
//! rpsbproc, `blastinput::value_parsers::ncbi_double`, and segmasker's `-locut` and `-hicut`,
//! `ncbi_app::args`): [`string_to_double_posix_or_local`]; and NCBI's own decimal reader
//! `NStr::StringToDoublePosix` that it calls, which the ASN.1 text reader calls too for a REAL
//! (`report::serial::asn_text_reader`): [`string_to_double_posix`].
//!
//! Only these flag sets are ported: the leading and trailing symbols, `fMandatorySign`,
//! `fDecimalPosix` and `fConvErr_NoThrow` of `StringToDouble` are not (no caller passes them).
//! The decimal point of `strtod` is that of the C locale, `.` (`s_IsDecimalPoint` with no flag
//! asks `localeconv`; NCBI's applications keep the C locale).

use super::{convert_error, convert_site, is_space};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::utils::x87::X87;

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
    /// A number whose `int` arithmetic in `StringToDoublePosix` overflows ([`IntOverflow`]): a
    /// text of 2^31 bytes or more.
    IntOverflow,
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
// `StringToDoublePosix` reads further than glibc's `strtod` only for a decimal mantissa followed
// by an `e` or `E` that ends the text (`1e`, `2.5E`: it keeps the mark when nothing follows it and
// steps back to it before a sign or another byte; an exponent with digits `strtod` reads as
// well; a comma is never its decimal point; its `inf` and `nan` end where `strtod`'s do or
// before). Its `errno` (`EINVAL` without digits, `ERANGE`) never turns a number of `strtod` into
// an error: it is `EINVAL` only where `strtod` reads nothing either.
/// `NStr::StringToDouble(str, NStr::fDecimalPosixOrLocal)`: as [`string_to_double`], but a comma
/// passes the check of the first byte (no number then starts there, and the error is that of
/// no number), and NCBI's own reader `StringToDoublePosix` ([`string_to_double_posix`]) reads the
/// text too; its number is taken where it reads further than `strtod`.
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
    let posix = string_to_double_posix(c_string, false);
    if posix.end > used {
        value = posix
            .value
            .map_err(|IntOverflow| PosixOrLocalError::IntOverflow)?;
        used = posix.end;
    }
    if used == 0 {
        return Err(error(STRING_TO_DOUBLE_NO_NUMBER, 0));
    }
    if at(used) != 0 {
        return Err(error(STRING_TO_DOUBLE_ENDPTR, used));
    }
    Ok(value)
}

/// What [`string_to_double_posix`] read: NCBI's number, and the bytes of the text that NCBI's
/// reader took (`*endptr - ptr`), which do not depend on the number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PosixDouble {
    pub(crate) value: Result<f64, IntOverflow>,
    pub(crate) end: usize,
}

/// An `int` of `StringToDoublePosix`'s arithmetic that overflows, whose result is undefined
/// behaviour in NCBI's C++ (LOSAT refuses the number): `digits` or `dot_position` (a mantissa of
/// 2^31 digits or leading zeros), `dot_position - digits` or the exponent `exponent - expvalue`
/// (`1.55e-2147483647`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct IntOverflow;

/// NCBI's `HUGE_VAL/HUGE_VAL` (`NCBI_FAKE_WARNING`), which GCC leaves to the run time (the
/// division raises the invalid-operation flag): `divsd` of x86-64 gives its default NaN, the
/// "QNaN floating-point indefinite" with the sign bit set (NCBI's corelib, `nfix/probe/x87`:
/// `nan` reads as 0xfff8000000000000).
const HUGE_VAL_DIV_HUGE_VAL: f64 = f64::from_bits(0xfff8_0000_0000_0000);

/// The x87 `long double` constant `1e<n>L` of `StringToDoublePosix` (`NCBI_CONST_LONGDOUBLE`):
/// `1e-15` to `1e15` and the multiples of 16 up to 256 of either sign.
fn long_double_pow10(n: i32) -> X87 {
    crate::report::serial::real::long_double_pow10(n)
        .expect("the long double literals of StringToDoublePosix")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:985-1278
// ```c++
// double NStr::StringToDoublePosix(const char* ptr, char** endptr, TStringToNumFlags flags)
// {
//     S2N_CONVERT_GUARD(NStr::fConvErr_NoThrow);
//
//     const char* start = ptr;
//     char c = *ptr++;
//
//     // skip leading whitespace
//     while ( isspace((unsigned char)c) ) {
//         c = *ptr++;
//     }
//
//     int sign = 0;
//     if ( c == '-' ) {
//         sign = -1;
//         c = *ptr++;
//     }
//     else if ( c == '+' ) {
//         sign = +1;
//         c = *ptr++;
//     }
//
//     if (c == 0) {
//         if (endptr) {
//             *endptr = (char*)start;
//         }
//         err_guard.Set(EINVAL);
//         return 0.;
//     }
//
//     // short-cut - single digit
//     if ( !*ptr && c >= '0' && c <= '9' ) {
//         if (endptr) {
//             *endptr = (char*)ptr;
//         }
//         double result = c-'0';
//         // some compilers fail to negate zero
//         return sign < 0 ? (c == '0' ? -0. : -result) : result;
//     }
//
//     bool         dot = false, expn = false, anydigits = false;
//     int          digits = 0, dot_position = 0;
//     unsigned int first=0, second=0, first_mul=1;
//     long double  second_mul = NCBI_CONST_LONGDOUBLE(1.),
//                  third      = NCBI_CONST_LONGDOUBLE(0.);
//
//     // up to exponent
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
//         else {
//             --ptr;
//             if (!anydigits) {
//                 if ( !dot && (c == 'n' || c == 'N') &&
//                      NStr::strncasecmp(ptr,"nan",3)==0) {
//                     if (endptr) {
//                         *endptr = (char*)(ptr+3);
//                     }
//                     return HUGE_VAL/HUGE_VAL; /* NCBI_FAKE_WARNING */
//                 }
//                 if ( (c == 'i' || c == 'I') ) {
//                     if ( NStr::strncasecmp(ptr,"inf",3)==0) {
//                         ptr += 3;
//                         if ( NStr::strncasecmp(ptr,"inity",5)==0) {
//                             ptr += 5;
//                         }
//                         if (endptr) {
//                             *endptr = (char*)ptr;
//                         }
//                         return sign < 0 ? -HUGE_VAL : HUGE_VAL;
//                     }
//                 }
//             }
//             break;
//         }
//     }
//     // if no digits, stop now - error
//     if (!anydigits) {
//         if (endptr) {
//             *endptr = (char*)start;
//         }
//         err_guard.Set(EINVAL);
//         return 0.;
//     }
//     int exponent = dot ? dot_position - digits : 0;
//
//     // read exponent
//     if (expn && *ptr) {
//         int expvalue = 0;
//         bool expsign = false, expnegate= false;
//         int expdigits= 0;
//         for( ; ; ++ptr) {
//             c = *ptr;
//             // sign: should be no digits at this point
//             if (c == '-' || c == '+') {
//                 // if there was sign or digits, stop
//                 if (expsign || expdigits) {
//                     break;
//                 }
//                 expsign = true;
//                 expnegate = c == '-';
//             }
//             // digits: accumulate
//             else if (c >= '0' && c <= '9') {
//                 ++expdigits;
//                 int newexpvalue = expvalue*10 + (c-'0');
//                 if (newexpvalue > expvalue) {
//                     expvalue = newexpvalue;
//                 }
//             }
//             else {
//                 break;
//             }
//         }
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
//         else {
//             exponent = expnegate ? exponent - expvalue : exponent + expvalue;
//         }
//     }
//     long double ret;
//     if ( first_mul > 1 ) {
//         _ASSERT(first);
//         ret = ((long double)first * first_mul + second)* second_mul + third;
//     }
//     else {
//         _ASSERT(first_mul == 1);
//         _ASSERT(second == 0);
//         _ASSERT(second_mul == 1);
//         _ASSERT(third == 0);
//         ret = first;
//     }
//     // calculate exponent
//     if ( first && exponent ) {
//         // multiply by power of 10 only non-zero mantissa
//         if (exponent > 2*DBL_MAX_10_EXP) {
//             ret = (flags & fDecimalPosixFinite) ? DBL_MAX :  HUGE_VAL;
//             err_guard.Set(ERANGE);
//         } else if (exponent < 2*DBL_MIN_10_EXP) {
//             ret = (flags & fDecimalPosixFinite) ? DBL_MIN : 0.;
//             err_guard.Set(ERANGE);
//         } else {
//             if ( exponent > 0 ) {
//                 static const double mul1[16] = {
//                     1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7,
//                     1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15
//                 };
//                 ret *= mul1[exponent&15];
//                 if ( exponent >>= 4 ) {
//                     static const long double mul2[16] = {
//                         NCBI_CONST_LONGDOUBLE(1e0),
//                         NCBI_CONST_LONGDOUBLE(1e16),
//                         ...
//                         NCBI_CONST_LONGDOUBLE(1e240)
//                     };
//                     ret *= mul2[exponent&15];
//                     for ( exponent >>= 4; exponent; --exponent ) {
//                         ret *= NCBI_CONST_LONGDOUBLE(1e256);
//                     }
//                 }
//                 if (!finite(double(ret))) {
//                     if (flags & fDecimalPosixFinite) {
//                         ret = DBL_MAX;
//                     }
//                     err_guard.Set(ERANGE);
//                 }
//             }
//             else {
//                 exponent = -exponent;
//                 static const long double mul1[16] = {
//                     NCBI_CONST_LONGDOUBLE(1e-0),
//                     ...
//                     NCBI_CONST_LONGDOUBLE(1e-15)
//                 };
//                 ret *= mul1[exponent&15];
//                 if ( exponent >>= 4 ) {
//                     static const long double mul2[16] = {
//                         NCBI_CONST_LONGDOUBLE(1e-0),
//                         NCBI_CONST_LONGDOUBLE(1e-16),
//                         ...
//                         NCBI_CONST_LONGDOUBLE(1e-240)
//                     };
//                     ret *= mul2[exponent&15];
//                     for ( exponent >>= 4; exponent; --exponent ) {
//                         ret *= NCBI_CONST_LONGDOUBLE(1e-256);
//                     }
//                 }
//                 if ( ret < DBL_MIN ) {
//                     if (flags & fDecimalPosixFinite) {
//                         ret = DBL_MIN;
//                     }
//                     err_guard.Set(ERANGE);
//                 }
//             }
//         }
//     }
//     if ( sign < 0 ) {
//         ret = -ret;
//     }
//     // done
//     if (endptr) {
//         *endptr = (char*)ptr;
//     }
//     return (double)ret;
// }
// ```
/// `NStr::StringToDoublePosix(text, &endptr, flags)`, `finite` for `fDecimalPosixFinite`: NCBI's
/// own reader of a decimal number, for the C string `text` (its first 0 byte, or its end, ends
/// it). The `long double` arithmetic of NCBI's x86-64 build is the x87 unit's (`X87`, rounded to
/// 64 bits each time); a mantissa of more than 18 digits goes through `third` and `second_mul`,
/// whose x87 values pass the largest 80-bit value from some 4,900 digits on and become the
/// infinity (`None` here). The `errno` it sets is not kept: NCBI's callers do not read it
/// (`string_to_double_posix_or_local`, `CObjectIStreamAsn::ReadDouble`).
pub(crate) fn string_to_double_posix(text: &[u8], finite: bool) -> PosixDouble {
    // `*p` of the NUL-terminated string.
    let at = |p: usize| text.get(p).copied().unwrap_or(0);
    let starts_with_ignore_case = |p: usize, name: &[u8]| {
        text.get(p..p + name.len())
            .is_some_and(|bytes| bytes.eq_ignore_ascii_case(name))
    };
    let mut ptr = 0usize;
    let mut c = at(ptr);
    ptr += 1;
    while is_space(c) {
        c = at(ptr);
        ptr += 1;
    }
    let mut sign = 0;
    if c == b'-' {
        sign = -1;
        c = at(ptr);
        ptr += 1;
    } else if c == b'+' {
        sign = 1;
        c = at(ptr);
        ptr += 1;
    }
    if c == 0 {
        return PosixDouble {
            value: Ok(0.0),
            end: 0,
        };
    }
    if at(ptr) == 0 && c.is_ascii_digit() {
        let result = f64::from(c - b'0');
        let value = if sign < 0 {
            if c == b'0' {
                -0.0
            } else {
                -result
            }
        } else {
            result
        };
        return PosixDouble {
            value: Ok(value),
            end: ptr,
        };
    }
    let (mut dot, mut expn, mut anydigits) = (false, false, false);
    let (mut digits, mut dot_position): (i32, i32) = (0, 0);
    let (mut first, mut second, mut first_mul): (u32, u32, u32) = (0, 0, 1);
    // `second_mul` and `third`, nonnegative; `None` is the x87 infinity.
    let ten = X87::from_u32(10);
    let (mut second_mul, mut third) = (Some(X87::from_u32(1)), Some(X87::ZERO));
    // Whether an `int` overflowed (`IntOverflow`); the text is still read to NCBI's end.
    let mut overflow = false;
    loop {
        if c.is_ascii_digit() {
            let d = c - b'0';
            anydigits = true;
            let (next, wrapped) = digits.overflowing_add(1);
            digits = next;
            overflow |= wrapped;
            if overflow {
                // The number is refused; its digits are only read.
            } else if first == 0 {
                first = u32::from(d);
                if first == 0 {
                    digits -= 1;
                    if dot {
                        let (next, wrapped) = dot_position.overflowing_sub(1);
                        dot_position = next;
                        overflow |= wrapped;
                    }
                }
            } else if digits <= 9 {
                first = first * 10 + u32::from(d);
            } else if digits <= 18 {
                first_mul *= 10;
                second = second * 10 + u32::from(d);
            } else {
                second_mul = second_mul.and_then(|mul| mul.checked_mul(ten));
                third = third
                    .and_then(|third| third.checked_mul(ten))
                    .and_then(|third| third.checked_add(X87::from_u32(u32::from(d))));
            }
        } else if c == b'.' {
            if dot {
                ptr -= 1;
                break;
            }
            dot_position = digits;
            dot = true;
        } else if c == b'e' || c == b'E' {
            if !anydigits {
                ptr -= 1;
                break;
            }
            expn = true;
            break;
        } else {
            ptr -= 1;
            if !anydigits {
                if !dot && (c == b'n' || c == b'N') && starts_with_ignore_case(ptr, b"nan") {
                    return PosixDouble {
                        value: Ok(HUGE_VAL_DIV_HUGE_VAL),
                        end: ptr + 3,
                    };
                }
                if (c == b'i' || c == b'I') && starts_with_ignore_case(ptr, b"inf") {
                    let mut end = ptr + 3;
                    if starts_with_ignore_case(end, b"inity") {
                        end += 5;
                    }
                    let value = if sign < 0 {
                        f64::NEG_INFINITY
                    } else {
                        f64::INFINITY
                    };
                    return PosixDouble {
                        value: Ok(value),
                        end,
                    };
                }
            }
            break;
        }
        c = at(ptr);
        ptr += 1;
    }
    if !anydigits {
        return PosixDouble {
            value: Ok(0.0),
            end: 0,
        };
    }
    let mut exponent: i32 = 0;
    if dot {
        let (next, wrapped) = dot_position.overflowing_sub(digits);
        exponent = next;
        overflow |= wrapped;
    }
    if expn && at(ptr) != 0 {
        let (mut expvalue, mut expsign, mut expnegate, mut expdigits) = (0i32, false, false, 0);
        loop {
            let c = at(ptr);
            if c == b'-' || c == b'+' {
                if expsign || expdigits != 0 {
                    break;
                }
                expsign = true;
                expnegate = c == b'-';
            } else if c.is_ascii_digit() {
                expdigits += 1;
                // NCBI's guard against the overflow of `expvalue * 10`: the wrapped value is
                // not taken (NCBI's corelib keeps 214748364 for `1e2147483648`).
                let newexpvalue = expvalue.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
                if newexpvalue > expvalue {
                    expvalue = newexpvalue;
                }
            } else {
                break;
            }
            ptr += 1;
        }
        if expdigits == 0 {
            if expsign {
                ptr -= 1;
            }
            if expn {
                ptr -= 1;
            }
        } else {
            let (next, wrapped) = if expnegate {
                exponent.overflowing_sub(expvalue)
            } else {
                exponent.overflowing_add(expvalue)
            };
            exponent = next;
            overflow |= wrapped;
        }
    }
    if overflow {
        return PosixDouble {
            value: Err(IntOverflow),
            end: ptr,
        };
    }
    // `(long double)first * first_mul + second` is exact (below 10^18 < 2^64).
    let mut ret = if first_mul > 1 {
        let integer = X87::from_u64(u64::from(first) * u64::from(first_mul) + u64::from(second));
        second_mul
            .and_then(|mul| integer.checked_mul(mul))
            .and_then(|ret| third.and_then(|third| ret.checked_add(third)))
    } else {
        Some(X87::from_u32(first))
    };
    // DBL_MAX_10_EXP 308, DBL_MIN_10_EXP -307.
    if first != 0 && exponent != 0 {
        if exponent > 2 * 308 {
            ret = finite.then(|| X87::from_f64(f64::MAX));
        } else if exponent < 2 * -307 {
            ret = Some(if finite {
                X87::from_f64(f64::MIN_POSITIVE)
            } else {
                X87::ZERO
            });
        } else if exponent > 0 {
            const MUL1: [f64; 16] = [
                1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14,
                1e15,
            ];
            let mut e = exponent;
            ret = ret.and_then(|ret| ret.checked_mul(X87::from_f64(MUL1[(e & 15) as usize])));
            e >>= 4;
            if e != 0 {
                ret = ret.and_then(|ret| ret.checked_mul(long_double_pow10(16 * (e & 15))));
                e >>= 4;
                while e != 0 {
                    ret = ret.and_then(|ret| ret.checked_mul(long_double_pow10(256)));
                    e -= 1;
                }
            }
            // `!finite(double(ret))`
            if finite && !ret.is_some_and(|ret| ret.to_f64().is_finite()) {
                ret = Some(X87::from_f64(f64::MAX));
            }
        } else {
            let mut e = -exponent;
            ret = ret.and_then(|ret| ret.checked_mul(long_double_pow10(-(e & 15))));
            e >>= 4;
            if e != 0 {
                ret = ret.and_then(|ret| ret.checked_mul(long_double_pow10(-16 * (e & 15))));
                e >>= 4;
                while e != 0 {
                    ret = ret.and_then(|ret| ret.checked_mul(long_double_pow10(-256)));
                    e -= 1;
                }
            }
            // `ret < DBL_MIN`: the infinity is not.
            if finite && ret.is_some_and(|ret| ret.lt(X87::from_f64(f64::MIN_POSITIVE))) {
                ret = Some(X87::from_f64(f64::MIN_POSITIVE));
            }
        }
    }
    // `ret = -ret` (`fchs`) and `(double)ret` (the store to a double); the infinity keeps its sign.
    let value = match ret {
        Some(ret) if sign < 0 => ret.neg().to_f64(),
        Some(ret) => ret.to_f64(),
        None if sign < 0 => f64::NEG_INFINITY,
        None => f64::INFINITY,
    };
    PosixDouble {
        value: Ok(value),
        end: ptr,
    }
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
    // the errors.
    #[test]
    fn posix_or_local_values_are_ncbis() {
        let table = include_str!("../../tests/unit/string_to_double_posix_or_local_vectors.tsv");
        let mut rows = 0;
        for row in table.lines().filter(|row| !row.starts_with('#')) {
            let (text, bits) = row.rsplit_once('\t').unwrap();
            match (string_to_double_posix_or_local(text.as_bytes()), bits) {
                (Err(PosixOrLocalError::Ncbi(_)), "ERR") => {}
                (Ok(value), bits) if bits != "ERR" => {
                    assert_eq!(format!("{:016x}", value.to_bits()), bits, "{text:?}")
                }
                (other, bits) => panic!("{text:?}: {other:?}, NCBI {bits}"),
            }
            rows += 1;
        }
        assert_eq!(rows, 422);
    }

    /// The four conversions of a text as NCBI's probe prints them: `StringToDouble(text,
    /// fDecimalPosixOrLocal)`, `StringToDouble(text)`, `StringToDoublePosix(text, &end, 0)` and
    /// with `fDecimalPosixFinite`, each the bits of the double or `ERR`, the last two with
    /// `:end`.
    fn four_conversions(text: &[u8]) -> [String; 4] {
        let bits = |value: f64| format!("{:016x}", value.to_bits());
        let posix = |finite| {
            let read = string_to_double_posix(text, finite);
            format!(
                "{}:{}",
                read.value.map_or("OVERFLOW".into(), bits),
                read.end
            )
        };
        [
            string_to_double_posix_or_local(text).map_or("ERR".into(), bits),
            string_to_double(text).map_or("ERR".into(), bits),
            posix(false),
            posix(true),
        ]
    }

    // `LOSAT/tests/unit/string_to_double_posix_vectors.tsv`: NCBI's own `StringToDouble` with
    // both flag sets and `StringToDoublePosix` with both (its corelib linked to a probe,
    // `swebok-improve-20261105/sw8/nfix/probe/x87`): mantissas of more than 18 digits and more
    // than 15 places (`third`, `mul2`), exponents at the ends of NCBI's ranges, NaN (NCBI's
    // `HUGE_VAL/HUGE_VAL` has the sign bit), and random texts. The probe compared 900,184 texts
    // (these, 900,000 random ones and mantissas of up to 5,000 digits): all equal but the one
    // whose `int` arithmetic overflows (below).
    #[test]
    fn posix_values_and_ends_are_ncbis() {
        let table = include_str!("../../tests/unit/string_to_double_posix_vectors.tsv");
        let mut rows = 0;
        for row in table.lines().filter(|row| !row.starts_with('#')) {
            let mut fields = row.split('\t');
            let text = fields.next().unwrap();
            let ncbi: Vec<&str> = fields.collect();
            assert_eq!(
                four_conversions(text.as_bytes()),
                ncbi.as_slice(),
                "{text:?}"
            );
            rows += 1;
        }
        assert_eq!(rows, 564);
        // Texts of their own: white space first (`StringToDoublePosix` skips it, `StringToDouble`
        // refuses it), the empty text, and long mantissas, whose x87 value passes the largest
        // 80-bit value from 4,933 digits on (the infinity; scaled by at least 10^-614 it stays
        // above the largest double, so NCBI's result is the infinity or, below 10^-614, zero).
        let repeat = |byte: &str, count: usize| byte.repeat(count);
        let cases: [(String, [&str; 4]); 12] = [
            (
                "\x0b1.5".into(),
                ["ERR", "ERR", "3ff8000000000000:4", "3ff8000000000000:4"],
            ),
            (
                "\x0c-2.5e".into(),
                ["ERR", "ERR", "c004000000000000:6", "c004000000000000:6"],
            ),
            (
                "".into(),
                ["ERR", "ERR", "0000000000000000:0", "0000000000000000:0"],
            ),
            (
                format!("1{}", repeat("0", 4931)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "7ff0000000000000:4932",
                    "7ff0000000000000:4932",
                ],
            ),
            (
                format!("1{}", repeat("0", 4933)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "7ff0000000000000:4934",
                    "7ff0000000000000:4934",
                ],
            ),
            (
                format!("1{}e", repeat("0", 5000)),
                [
                    "7ff0000000000000",
                    "ERR",
                    "7ff0000000000000:5002",
                    "7ff0000000000000:5002",
                ],
            ),
            (
                format!("{}e1", repeat("9", 5000)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "7ff0000000000000:5002",
                    "7fefffffffffffff:5002",
                ],
            ),
            (
                format!("{}e-1", repeat("9", 5000)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "7ff0000000000000:5003",
                    "7ff0000000000000:5003",
                ],
            ),
            (
                format!("{}.{}e", repeat("9", 5000), repeat("9", 700)),
                [
                    "0000000000000000",
                    "ERR",
                    "0000000000000000:5702",
                    "0010000000000000:5702",
                ],
            ),
            (
                format!("{}.{}e", repeat("9", 5000), repeat("9", 600)),
                [
                    "7ff0000000000000",
                    "ERR",
                    "7ff0000000000000:5602",
                    "7ff0000000000000:5602",
                ],
            ),
            (
                format!("{}e-600", repeat("5", 4950)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "7ff0000000000000:4955",
                    "7ff0000000000000:4955",
                ],
            ),
            (
                format!("{}e-616", repeat("5", 4950)),
                [
                    "7ff0000000000000",
                    "7ff0000000000000",
                    "0000000000000000:4955",
                    "0010000000000000:4955",
                ],
            ),
        ];
        for (text, ncbi) in cases {
            assert_eq!(
                four_conversions(text.as_bytes()),
                ncbi,
                "{}",
                &text[..text.len().min(20)]
            );
        }
        // NCBI's `exponent - expvalue` overflows `int` (its build gives the infinity, undefined
        // behaviour): refused, at NCBI's end; `strtod` reads the whole text, so `CArg_Double`
        // takes its number (0).
        assert_eq!(
            four_conversions(b"1.55e-2147483647"),
            [
                "0000000000000000",
                "0000000000000000",
                "OVERFLOW:16",
                "OVERFLOW:16"
            ]
        );
        assert_eq!(
            string_to_double_posix(b"1.55e-2147483647", false).value,
            Err(IntOverflow)
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

    // glibc's `strtod` reads these to their end and stops before the end of these (the SEG
    // locut and hicut of the argument reader, `blastinput::app::read_seg_option`, test
    // `string_to_double` with it, which also needs a digit, a point or a sign first).
    #[test]
    fn strtod_reads_the_whole_of_these_only() {
        let whole = |text: &str| {
            let (_, used) = strtod_prefix(text.as_bytes());
            used != 0 && used == text.len()
        };
        for text in [
            "1",
            "+1",
            "-.5",
            "5.",
            "1e5",
            "1E-5",
            "0x10",
            "0X1p3",
            "0x.8p-1",
            "0x1.",
            "inf",
            "+INFINITY",
            "-nan",
            "nan(x_1)",
            "1e400",
        ] {
            assert!(whole(text), "{text}");
        }
        for text in [
            "1e", "1e+", "2.5x", "2,2", ".", "+.", ".e1", "0x", "0x.", "0x1p", "1_0", "--5",
            "2.5.5", "nan(", "nan(-)", "infx", "",
        ] {
            assert!(!whole(text), "{text}");
        }
    }
}
