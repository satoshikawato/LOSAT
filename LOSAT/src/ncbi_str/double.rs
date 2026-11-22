//! `NStr::StringToDouble` with its default flags (corelib/ncbistr.cpp), which rpsbproc's
//! data loader calls (`NStr::StringToNumeric<double>`, ncbistr.hpp:3632-3637) for the bit score
//! thresholds of `bitscore_specific.txt`, the BLAST programs' argument readers for the SEG locut
//! and hicut (`blastinput::app::read_seg_option`) and the constraints of the arguments
//! (`blastinput::value_parsers`), and glibc's `strtod` that it reads the number with.
//!
//! Also `NStr::StringToDouble` with `NStr::fDecimalPosixOrLocal`, which NCBI's argument reader
//! calls for an `eDouble` key (`CArg_Double`: every real argument of the BLAST programs and
//! rpsbproc, `blastinput::value_parsers::ncbi_double`, and segmasker's `-locut` and `-hicut`,
//! `ncbi_app::args`): [`string_to_double_posix_or_local`]; and NCBI's own decimal reader
//! `NStr::StringToDoublePosix` that it calls, which the ASN.1 text and XML readers call too for a
//! REAL (`report::serial::asn_text_reader`, `report::serial::xml_reader`):
//! [`string_to_double_posix`].
//!
//! Only these flag sets are ported: the leading and trailing symbols, `fMandatorySign`,
//! `fDecimalPosix` and `fConvErr_NoThrow` of `StringToDouble` are not (no caller passes them).
//! The decimal point of `strtod` is that of the C locale, `.` (`s_IsDecimalPoint` with no flag
//! asks `localeconv`; NCBI's applications keep the C locale).

use super::{convert_error, convert_site, is_space};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::utils::x87::X87;

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1312-1319
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1342-1349
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1367-1375
// The errors are thrown in `s_StringToDouble` (a function of the `ncbi` namespace that takes a
// `const char*`, a `size_t` and more, whose class part is `ncbi` for
// `CDiagCompileInfo::ParseCurrFunctName`, ncbidiag.cpp:826-888) at lines 1317 (the first byte),
// 1343 (no number) and 1374 (`CHECK_ENDPTR`); glibc's `strtod` sets no errno other than `ERANGE`,
// so line 1348 is not reached.
const STRING_TO_DOUBLE_INVAL: ThrowSite = convert_site(1317, "ncbi::s_StringToDouble");
const STRING_TO_DOUBLE_NO_NUMBER: ThrowSite = convert_site(1343, "ncbi::s_StringToDouble");
const STRING_TO_DOUBLE_ENDPTR: ThrowSite = convert_site(1374, "ncbi::s_StringToDouble");

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1386-1403
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1283-1285
// NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:422-433
// `str` of `s_StringToDouble` is the C string, so `strtod` and the end check (`str[pos]`) stop
// at its first 0 byte ("12\0x" is 12), and the message of an error, `Message`, quotes
// the text before that byte (the other conversions quote the whole `CTempString`).
/// `NStr::StringToDouble` with the default flags: the first byte must be a digit, the
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:724-737
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1329-1341
// `StringToDoublePosix` reads further than glibc's `strtod` for a decimal mantissa followed by an
// `e` or `E` that ends the text (`1e`, `2.5E`: it keeps the mark when nothing follows it and steps
// back to it before a sign or another byte; an exponent with digits `strtod` reads as well), and
// for an `inf` after a dot (`.inf`, `-.Infinity`: its `inf` is read after a dot, its `nan` is not;
// `strtod` reads neither). A comma is never its decimal point; its other `inf` and its `nan` end
// where `strtod`'s do or before. Its `errno` (`EINVAL` without digits, `ERANGE`) never turns a
// number of `strtod` into an error: it is `EINVAL` only where `strtod` reads nothing either.
/// `NStr::StringToDouble` with `NStr::fDecimalPosixOrLocal`: as [`string_to_double`], but a comma
/// passes the check of the first byte (no number then starts there, and the error is that of
/// no number), and NCBI's own reader `StringToDoublePosix` ([`string_to_double_posix`]) reads the
/// text too; its number is taken where it reads further than `strtod`.
pub(crate) fn string_to_double_posix_or_local(s: &[u8]) -> Result<f64, NcbiException> {
    let c_string = &s[..s.iter().position(|&b| b == 0).unwrap_or(s.len())];
    let error = |site: ThrowSite, pos: usize| convert_error(c_string, "double", "", site, pos);
    let at = |pos: usize| c_string.get(pos).copied().unwrap_or(0);
    let c = at(0);
    if !c.is_ascii_digit() && !matches!(c, b'.' | b',' | b'-' | b'+') {
        return Err(error(STRING_TO_DOUBLE_INVAL, 0));
    }
    let (mut value, mut used) = strtod_prefix(c_string);
    let posix = string_to_double_posix(c_string, false);
    if posix.end > used {
        value = posix.value;
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
/// reader took (from `ptr` to `*endptr`), which do not depend on the number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PosixDouble {
    pub(crate) value: f64,
    pub(crate) end: usize,
}

/// NCBI's `HUGE_VAL/HUGE_VAL` (`NCBI_FAKE_WARNING`), the NaN of NCBI's readers of a REAL (`nan`
/// of `StringToDoublePosix`, `NOT-A-NUMBER` of ASN.1 text): GCC leaves the division to the run
/// time (it raises the invalid-operation flag), and `divsd` of x86-64 gives its default NaN, the
/// "QNaN floating-point indefinite", whose sign bit is set, so that an `ostream` prints it `-nan`
/// (NCBI's corelib, `nfix/probe/x87`: `nan` reads as 0xfff8000000000000; NCBI rpsbproc 0.5.1,
/// both builds, PS16d review oracles O1 `x_nan` and `x_mnan`, `Parameters_expect` `nan` and
/// `-nan`; O2 `a_nan`, an archive's E-value threshold `NOT-A-NUMBER`; O1 `s_doc2nan`, a kept
/// `Hsp_evalue` `nan`: `-nan` on the SESSION and DOMAINS lines). `StringToDoublePosix` gives it
/// for `nan` whatever the sign before it.
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1082-1087
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:810-811
pub(crate) const NCBI_NAN: f64 = f64::from_bits(0xfff8_0000_0000_0000);

/// The x87 long double constant `1e<n>L` of `StringToDoublePosix` (`NCBI_CONST_LONGDOUBLE`):
/// `1e-15` to `1e15` and the multiples of 16 up to 256 of either sign.
fn long_double_pow10(n: i64) -> X87 {
    i32::try_from(n)
        .ok()
        .and_then(crate::report::serial::real::long_double_pow10)
        .expect("the long double literals of StringToDoublePosix")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:985-1278
/// `NStr::StringToDoublePosix` with `text`, an end pointer and `flags`, `finite` for
/// `fDecimalPosixFinite`: NCBI's own reader of a decimal number, for the C string `text` (its first
/// 0 byte, or its end, ends it). The long double arithmetic of NCBI's x86-64 build is the x87
/// unit's (`X87`, rounded to 64 bits each time); a mantissa of more than 18 digits goes through
/// `third` and `second_mul`, whose x87 values pass the largest 80-bit value from some 4,900 digits
/// on and become the infinity (`None` here). The `errno` it sets is not kept: no result depends on
/// it (`string_to_double_posix_or_local` reads it, where it is `EINVAL` only if `strtod` read
/// nothing too; `CObjectIStreamAsn::ReadDouble` and `CObjectIStreamXml::ReadDouble` do not).
///
/// NCBI's `int` counters and exponent (`digits`, `dot_position`, `expvalue`, `exponent`) are
/// `i64` here, and `expvalue * 10 + digit` and `exponent ± expvalue` saturate at its ends: where
/// NCBI's `int` overflows (a signed overflow, undefined in C++; NCBI's corelib wraps
/// `1.55e-2147483647`'s exponent to 2147483647 and gives the infinity, `DBL_MAX` with
/// `fDecimalPosixFinite`; PD-LOSAT-NCBI-DEFECTS exception 11), LOSAT computes the
/// exponent without the wrap, as the Owner's rule of 2026-10-10 has it for NCBI's signed
/// overflows (OWNER_DECISIONS 2026-10-10, NCBI parity DECISIONS 17 of 2026-10-22
/// 「広い型になおしましょう」): `1.55e-2147483647` is 0 (`DBL_MIN` with `fDecimalPosixFinite`). Every
/// exponent that NCBI's `int` holds is the same, and so is every result: a saturated exponent
/// (beyond 2^63) is far outside -614..616, as the exact one is, and NCBI's guard
/// (`newexpvalue > expvalue`) keeps every exponent that does not wrap.
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
        return PosixDouble { value: 0.0, end: 0 };
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
        return PosixDouble { value, end: ptr };
    }
    let (mut dot, mut expn, mut anydigits) = (false, false, false);
    // NCBI's `digits` and `dot_position`, both `int`, widened (see above).
    let (mut digits, mut dot_position): (i64, i64) = (0, 0);
    let (mut first, mut second, mut first_mul): (u32, u32, u32) = (0, 0, 1);
    // `second_mul` and `third`, nonnegative; `None` is the x87 infinity.
    let ten = X87::from_u32(10);
    let (mut second_mul, mut third) = (Some(X87::from_u32(1)), Some(X87::ZERO));
    loop {
        if c.is_ascii_digit() {
            let d = c - b'0';
            anydigits = true;
            digits += 1;
            if first == 0 {
                first = u32::from(d);
                if first == 0 {
                    digits -= 1;
                    if dot {
                        dot_position -= 1;
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
                        value: NCBI_NAN,
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
                    return PosixDouble { value, end };
                }
            }
            break;
        }
        c = at(ptr);
        ptr += 1;
    }
    if !anydigits {
        return PosixDouble { value: 0.0, end: 0 };
    }
    // NCBI's `exponent`, an `int`, widened: both counters are at most the length of the text.
    let mut exponent: i64 = if dot { dot_position - digits } else { 0 };
    if expn && at(ptr) != 0 {
        let (mut expvalue, mut expsign, mut expnegate, mut expdigits) = (0i64, false, false, 0u64);
        loop {
            let c = at(ptr);
            if c == b'-' || c == b'+' {
                if expsign || expdigits != 0 {
                    break;
                }
                expsign = true;
                expnegate = c == b'-';
            } else if c.is_ascii_digit() {
                expdigits = expdigits.saturating_add(1);
                // `expvalue*10 + (c-'0')`, widened and saturated: NCBI's guard keeps every value
                // that does not wrap (NCBI's corelib keeps 214748364 for `1e2147483648`, which
                // is the infinity either way).
                let newexpvalue = expvalue
                    .saturating_mul(10)
                    .saturating_add(i64::from(c - b'0'));
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
            exponent = if expnegate {
                exponent.saturating_sub(expvalue)
            } else {
                exponent.saturating_add(expvalue)
            };
        }
    }
    // `(long double)first * first_mul + second` is exact (below 10^18 < 2^64). Without digits
    // after the 18th (`digits <= 18`), `second_mul` is 1 and `third` 0, so NCBI's `* second_mul +
    // third` is exact and is skipped.
    let mut ret = if first_mul > 1 {
        let integer = X87::from_u64(u64::from(first) * u64::from(first_mul) + u64::from(second));
        if digits <= 18 {
            Some(integer)
        } else {
            second_mul
                .and_then(|mul| integer.checked_mul(mul))
                .and_then(|ret| third.and_then(|third| ret.checked_add(third)))
        }
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
            // `ret` below `DBL_MIN`: the infinity is not.
            if finite && ret.is_some_and(|ret| ret.lt(X87::from_f64(f64::MIN_POSITIVE))) {
                ret = Some(X87::from_f64(f64::MIN_POSITIVE));
            }
        }
    }
    // The negation of `ret` (`fchs`) and `ret` cast to `double` (the store to a double); the
    // infinity keeps its sign.
    let value = match ret {
        Some(ret) if sign < 0 => ret.neg().to_f64(),
        Some(ret) => ret.to_f64(),
        None if sign < 0 => f64::NEG_INFINITY,
        None => f64::INFINITY,
    };
    PosixDouble { value, end: ptr }
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

/// The payload of `nan(n-char-sequence)`: the sequence read as by `strtoull` with base 0 (C
/// standard: `0x` hexadecimal, a leading `0` octal, else decimal; the largest value on
/// overflow), its low 51 bits kept below the quiet bit, when the whole sequence is that number;
/// 0 otherwise. These are the bits of NCBI 2.17.0's `NStr::StringToDouble` as measured
/// (`LOSAT/tests/unit/string_to_double_posix_vectors.tsv`, which has every case above).
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
#[path = "double_tests.rs"]
mod tests;
