//! NCBI's real-number digits for the serialized formats: `NStr::DoubleToString_Ecvt`,
//! `NStr::DoubleToStringPosix` and the `WriteDouble` forms of the ASN.1 text, JSON and XML
//! object streams (outfmt 5, 8, 9, 11, 12, 13-16).
//!
//! NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:78-81
//! ```c++
//! NCBI_PARAM_DECL(bool, SERIAL, FastWriteDouble);
//! NCBI_PARAM_DEF(bool, SERIAL, FastWriteDouble, true);
//! typedef NCBI_PARAM_TYPE(SERIAL, FastWriteDouble) TFastWriteDouble;
//! static CSafeStatic<TFastWriteDouble> s_FastWriteDouble;
//! ```
//! NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:305-308
//! ```c++
//! void CObjectOStreamAsn::WriteDouble(double data)
//! {
//!     WriteDouble2(data, DBL_DIG);
//! }
//! ```
//!
//! Only the default `FastWriteDouble` branches exist here (a registry or environment value that
//! turns it off is rejected by the callers, outfmt-20261019 DECISIONS 13). The digits come from
//! 80-bit x87 `long double` arithmetic on NCBI's x86-64 build; [`crate::utils::x87`] reproduces it, so
//! no `f64` formatting of the Rust standard library is involved and every platform prints the
//! same bytes.

use crate::utils::x87::X87;

/// `DBL_DIG`: the digits every `WriteDouble` of the ASN.1 text, ASN.1 binary, JSON and XML
/// object streams passes to the digit generator.
pub const DBL_DIG: u32 = 15;

/// The `long double` literals of `DoubleToString_Ecvt` (`__NLG(1.e<n>)`): `(n, significand,
/// unbiased exponent)` of the 80-bit value gcc gives the literal `1.e<n>L` (correct rounding to
/// 64 significand bits, ties to even). Each entry is checked by exact big-integer arithmetic in
/// `real_tests.rs`.
const LONG_DOUBLE_POW10: [(i32, u64, i32); 71] = [
    (-256, 0xC031_4325_637A_193A, -851),
    (-240, 0xD560_5FCD_CF32_E1D7, -798),
    (-224, 0xECE5_3CEC_4A31_4EBE, -745),
    (-208, 0x8380_DEA9_3DA4_BC60, -691),
    (-192, 0x91FF_8377_5423_CC06, -638),
    (-176, 0xA217_27DB_38CB_0030, -585),
    (-160, 0xB3F4_E093_DB73_A093, -532),
    (-144, 0xC7CA_BA6E_7C53_82C9, -479),
    (-128, 0xDDD0_467C_64BC_E4A1, -426),
    (-112, 0xF643_35BC_F065_D37D, -373),
    (-96, 0x88B4_02F7_FD75_539B, -319),
    (-80, 0x97C5_60BA_6B09_19A6, -266),
    (-64, 0xA87F_EA27_A539_E9A5, -213),
    (-48, 0xBB12_7C53_B17E_C159, -160),
    (-32, 0xCFB1_1EAD_4539_94BA, -107),
    (-16, 0xE695_94BE_C44D_E15B, -54),
    (-15, 0x901D_7CF7_3AB0_ACD9, -50),
    (-14, 0xB424_DC35_095C_D80F, -47),
    (-13, 0xE12E_1342_4BB4_0E13, -44),
    (-12, 0x8CBC_CC09_6F50_88CC, -40),
    (-11, 0xAFEB_FF0B_CB24_AAFF, -37),
    (-10, 0xDBE6_FECE_BDED_D5BF, -34),
    (-9, 0x8970_5F41_36B4_A597, -30),
    (-8, 0xABCC_7711_8461_CEFD, -27),
    (-7, 0xD6BF_94D5_E57A_42BC, -24),
    (-6, 0x8637_BD05_AF6C_69B6, -20),
    (-5, 0xA7C5_AC47_1B47_8423, -17),
    (-4, 0xD1B7_1758_E219_652C, -14),
    (-3, 0x8312_6E97_8D4F_DF3B, -10),
    (-2, 0xA3D7_0A3D_70A3_D70A, -7),
    (-1, 0xCCCC_CCCC_CCCC_CCCD, -4),
    (0, 0x8000_0000_0000_0000, 0),
    (1, 0xA000_0000_0000_0000, 3),
    (2, 0xC800_0000_0000_0000, 6),
    (3, 0xFA00_0000_0000_0000, 9),
    (4, 0x9C40_0000_0000_0000, 13),
    (5, 0xC350_0000_0000_0000, 16),
    (6, 0xF424_0000_0000_0000, 19),
    (7, 0x9896_8000_0000_0000, 23),
    (8, 0xBEBC_2000_0000_0000, 26),
    (9, 0xEE6B_2800_0000_0000, 29),
    (10, 0x9502_F900_0000_0000, 33),
    (11, 0xBA43_B740_0000_0000, 36),
    (12, 0xE8D4_A510_0000_0000, 39),
    (13, 0x9184_E72A_0000_0000, 43),
    (14, 0xB5E6_20F4_8000_0000, 46),
    (15, 0xE35F_A931_A000_0000, 49),
    (16, 0x8E1B_C9BF_0400_0000, 53),
    (17, 0xB1A2_BC2E_C500_0000, 56),
    (18, 0xDE0B_6B3A_7640_0000, 59),
    (19, 0x8AC7_2304_89E8_0000, 63),
    (20, 0xAD78_EBC5_AC62_0000, 66),
    (21, 0xD8D7_26B7_177A_8000, 69),
    (22, 0x8786_7832_6EAC_9000, 73),
    (23, 0xA968_163F_0A57_B400, 76),
    (24, 0xD3C2_1BCE_CCED_A100, 79),
    (32, 0x9DC5_ADA8_2B70_B59E, 106),
    (48, 0xAF29_8D05_0E43_95D7, 159),
    (64, 0xC278_1F49_FFCF_A6D5, 212),
    (80, 0xD7E7_7A8F_87DA_F7FC, 265),
    (96, 0xEFB3_AB16_C59B_14A3, 318),
    (112, 0x850F_ADC0_9923_329E, 372),
    (128, 0x93BA_47C9_80E9_8CE0, 425),
    (144, 0xA402_B9C5_A8D3_A6E7, 478),
    (160, 0xB616_A12B_7FE6_17AA, 531),
    (176, 0xCA28_A291_859B_BF93, 584),
    (192, 0xE070_F78D_3927_556B, 637),
    (208, 0xF92E_0C35_3782_6146, 690),
    (224, 0x8A52_96FF_E33C_C930, 744),
    (240, 0x9991_A6F3_D6BF_1766, 797),
    (256, 0xAA7E_EBFB_9DF9_DE8E, 850),
];

/// `__NLG(1.e<n>)`. Used only in constant items, so a missing entry is a compile error.
const fn ld(n: i32) -> X87 {
    let mut i = 0;
    while i < LONG_DOUBLE_POW10.len() {
        let (k, significand, exponent) = LONG_DOUBLE_POW10[i];
        if k == n {
            return X87::from_parts(significand, exponent);
        }
        i += 1;
    }
    panic!("no long double literal 1e<n> in LONG_DOUBLE_POW10");
}

#[cfg(test)]
pub(super) fn long_double_pow10_table() -> &'static [(i32, u64, i32)] {
    &LONG_DOUBLE_POW10
}

const LD_1: X87 = ld(0);
const LD_1E8: X87 = ld(8);
const LD_1E9: X87 = ld(9);
const LD_1E16: X87 = ld(16);
const LD_1E256: X87 = ld(256);
const LD_1EM16: X87 = ld(-16);
const LD_1EM256: X87 = ld(-256);

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2486-2502
/// ```c++
///         if (value >= __NLG(1.e16)) {
///             if      (value>=__NLG(1.e240)) {value*=__NLG(1.e-240); exp+=240;}
///             else if (value>=__NLG(1.e224)) {value*=__NLG(1.e-224); exp+=224;}
///             ...
///             else if (value>=__NLG(1.e16))  {value*=__NLG(1.e-16);  exp+=16;}
///         }
/// ```
/// The `else if` lines in source order: `(threshold, factor, exp step)`.
const POSITIVE_SCALE: [(X87, X87, u32); 15] = [
    (ld(240), ld(-240), 240),
    (ld(224), ld(-224), 224),
    (ld(208), ld(-208), 208),
    (ld(192), ld(-192), 192),
    (ld(176), ld(-176), 176),
    (ld(160), ld(-160), 160),
    (ld(144), ld(-144), 144),
    (ld(128), ld(-128), 128),
    (ld(112), ld(-112), 112),
    (ld(96), ld(-96), 96),
    (ld(80), ld(-80), 80),
    (ld(64), ld(-64), 64),
    (ld(48), ld(-48), 48),
    (ld(32), ld(-32), 32),
    (ld(16), ld(-16), 16),
];

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2503-2520
/// ```c++
///         if      (value<   __NLG(1.)) {mult=__NLG(1.e+9); exp-= 1;}
///         else if (value<  __NLG(10.)) {mult=__NLG(1.e+8);         }
///         else if (value< __NLG(1.e2)) {mult=__NLG(1.e+7); exp+= 1;}
///         ...
///         else if (value<__NLG(1.e16)) {mult=__NLG(1.e-7); exp+=15;}
///         else                         {mult=__NLG(1.e-8); exp+=16;}
/// ```
/// The `if` / `else if` lines in source order: `(bound, mult, exp step)`; the first line
/// (`exp -= 1`) and the final `else` are written out in the port.
const POSITIVE_MULT: [(X87, X87, u32); 16] = [
    (ld(1), ld(8), 0),
    (ld(2), ld(7), 1),
    (ld(3), ld(6), 2),
    (ld(4), ld(5), 3),
    (ld(5), ld(4), 4),
    (ld(6), ld(3), 5),
    (ld(7), ld(2), 6),
    (ld(8), ld(1), 7),
    (ld(9), ld(0), 8),
    (ld(10), ld(-1), 9),
    (ld(11), ld(-2), 10),
    (ld(12), ld(-3), 11),
    (ld(13), ld(-4), 12),
    (ld(14), ld(-5), 13),
    (ld(15), ld(-6), 14),
    (ld(16), ld(-7), 15),
];

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2524-2540
/// ```c++
///         if (value <= __NLG(1.e-16)) {
///             if      (value<=__NLG(1.e-240)) {value*=__NLG(1.e240); exp+=240;}
///             ...
///             else if (value<=__NLG(1.e-16))  {value*=__NLG(1.e16);  exp+=16;}
///         }
/// ```
const NEGATIVE_SCALE: [(X87, X87, u32); 15] = [
    (ld(-240), ld(240), 240),
    (ld(-224), ld(224), 224),
    (ld(-208), ld(208), 208),
    (ld(-192), ld(192), 192),
    (ld(-176), ld(176), 176),
    (ld(-160), ld(160), 160),
    (ld(-144), ld(144), 144),
    (ld(-128), ld(128), 128),
    (ld(-112), ld(112), 112),
    (ld(-96), ld(96), 96),
    (ld(-80), ld(80), 80),
    (ld(-64), ld(64), 64),
    (ld(-48), ld(48), 48),
    (ld(-32), ld(32), 32),
    (ld(-16), ld(16), 16),
];

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2541-2557
/// ```c++
///         if      (value<__NLG(1.e-15)) {mult=__NLG(1.e24); exp+=16;}
///         else if (value<__NLG(1.e-14)) {mult=__NLG(1.e23); exp+=15;}
///         ...
///         else if (value<__NLG(1.))     {mult=__NLG(1.e9);  exp+=1;}
///         else                          {mult=__NLG(1.e8);         }
/// ```
const NEGATIVE_MULT: [(X87, X87, u32); 16] = [
    (ld(-15), ld(24), 16),
    (ld(-14), ld(23), 15),
    (ld(-13), ld(22), 14),
    (ld(-12), ld(21), 13),
    (ld(-11), ld(20), 12),
    (ld(-10), ld(19), 11),
    (ld(-9), ld(18), 10),
    (ld(-8), ld(17), 9),
    (ld(-7), ld(16), 8),
    (ld(-6), ld(15), 7),
    (ld(-5), ld(14), 6),
    (ld(-4), ld(13), 5),
    (ld(-3), ld(12), 4),
    (ld(-2), ld(11), 3),
    (ld(-1), ld(10), 2),
    (ld(0), ld(9), 1),
];

/// The result of `NStr::DoubleToString_Ecvt`: the digit buffer, `*dec` and `*sign`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcvtDigits {
    buffer: [u8; 32],
    len: usize,
    dec: i32,
    sign: i32,
}

impl EcvtDigits {
    /// The returned digits (`buffer[0..return value]`): significant digits without trailing
    /// zeros; `"0"` for +0 and `"-0"` for -0; empty for precision 0.
    pub fn digits(&self) -> &[u8] {
        &self.buffer[..self.len]
    }

    /// `*dec`: the decimal exponent of the first digit.
    pub fn dec(&self) -> i32 {
        self.dec
    }

    /// `*sign`: -1 or 1; 0 for +0 and for precision 0.
    pub fn sign(&self) -> i32 {
        self.sign
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2422-2442
/// ```c++
/// static char* s_ncbi_append_int2str(char* buffer, unsigned int value, size_t digits, bool zeros)
/// {
///     char* buffer_start = buffer;
///     char* buffer_end = (buffer += digits-1);
///     if (zeros) {
///         do {
///             *buffer-- = (char)('0' + (value % 10));
///             value /= 10;
///         } while (--digits);
///     } else {
///         do {
///             *buffer-- = (char)('0' + (value % 10));
///         } while (value /= 10);
///
///         if (++buffer != buffer_start) {
///             memmove(buffer_start, buffer, buffer_end-buffer+1);
///             buffer_end -= buffer - buffer_start;
///         }
///     }
///     return ++buffer_end;
/// }
/// ```
/// Writes at `buffer[start..]` and returns the index after the last character. Without
/// `zeros`, the callers' values never have more than `digits` digits (NCBI would write before
/// `buffer_start`).
fn s_ncbi_append_int2str(
    buffer: &mut [u8],
    start: usize,
    mut value: u32,
    digits: usize,
    zeros: bool,
) -> usize {
    let mut buffer_end = start + digits - 1;
    if zeros {
        for i in 0..digits {
            buffer[buffer_end - i] = b'0' + (value % 10) as u8;
            value /= 10;
        }
    } else {
        let mut written = 0;
        loop {
            debug_assert!(written < digits);
            buffer[buffer_end - written] = b'0' + (value % 10) as u8;
            written += 1;
            value /= 10;
            if value == 0 {
                break;
            }
        }
        let first = buffer_end + 1 - written;
        if first != start {
            buffer.copy_within(first..=buffer_end, start);
            buffer_end -= first - start;
        }
    }
    buffer_end + 1
}

/// `NStr::DoubleToString_Ecvt(val, precision, buffer, bufsize, &dec, &sign)` for a finite
/// `val` (NCBI's scaling loop never ends for an infinity; its writers test NaN and infinity
/// first). The buffer is the external one (`bufsize > 20`, as every serial caller passes).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2447-2654 (quoted piece by piece
/// below; `__NLG` = `NCBI_CONST_LONGDOUBLE`, a `long double` literal)
pub fn double_to_string_ecvt(val: f64, precision: u32) -> EcvtDigits {
    assert!(
        val.is_finite(),
        "DoubleToString_Ecvt of a NaN or an infinity"
    );
    let mut out = EcvtDigits {
        buffer: [0; 32],
        len: 0,
        dec: 0,
        sign: 0,
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2452-2474
    // ```c++
    //     *dec = *sign = 0;
    //     if (precision==0) {
    //         return 0;
    //     }
    //     if (precision > DBL_DIG) {
    //         precision = DBL_DIG;
    //     }
    //     if (val == 0.) {
    //         double zero = 0.;
    //         if (memcmp(&val, &zero, sizeof(double)) == 0) {
    //             *buffer='0';
    //             return 1;
    //         }
    //         *buffer++='-';
    //         *buffer='0';
    //         *sign = -1;
    //         return 2;
    //     }
    //     *sign = val < 0. ? -1 : 1;
    //     if (*sign < 0) {
    //         val = -val;
    //     }
    //     bool high_precision = precision > 9;
    // ```
    if precision == 0 {
        return out;
    }
    let precision = precision.min(DBL_DIG) as usize;
    if val == 0.0 {
        if val.to_bits() == 0 {
            out.buffer[0] = b'0';
            out.len = 1;
            return out;
        }
        out.buffer[0] = b'-';
        out.buffer[1] = b'0';
        out.sign = -1;
        out.len = 2;
        return out;
    }
    out.sign = if val < 0.0 { -1 } else { 1 };
    let val = if out.sign < 0 { -val } else { val };
    let high_precision = precision > 9;

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2476-2485
    // ```c++
    // // calculate exponent
    //     unsigned int exp=0;
    //     bool exp_positive = val >= 1.;
    //     unsigned int first, second=0;
    //     long double mult = __NLG(1.);
    //     long double value = val;
    //
    //     if (exp_positive) {
    //         while (value>=__NLG(1.e256))
    //             {value*=__NLG(1.e-256); exp+=256;}
    // ```
    let mut exp: u32 = 0;
    let mut exp_positive = val >= 1.0;
    let first: u32;
    let mut second: u32 = 0;
    let mut mult = LD_1;
    let mut value = X87::from_f64(val);

    if exp_positive {
        while value.ge(LD_1E256) {
            value = value.mul(LD_1EM256);
            exp += 256;
        }
        // ncbistr.cpp:2486-2502 (POSITIVE_SCALE)
        if value.ge(LD_1E16) {
            for &(threshold, factor, step) in &POSITIVE_SCALE {
                if value.ge(threshold) {
                    value = value.mul(factor);
                    exp += step;
                    break;
                }
            }
        }
        // ncbistr.cpp:2503-2520 (POSITIVE_MULT). `value < 1` only after a scaling step
        // (val >= 1), so `exp -= 1` never wraps.
        if value.lt(LD_1) {
            mult = LD_1E9;
            exp -= 1;
        } else {
            let mut matched = false;
            for &(bound, m, step) in &POSITIVE_MULT {
                if value.lt(bound) {
                    mult = m;
                    exp += step;
                    matched = true;
                    break;
                }
            }
            if !matched {
                mult = ld(-8);
                exp += 16;
            }
        }
    } else {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2521-2523
        // ```c++
        //     } else {
        //         while (value<=__NLG(1.e-256))
        //             {value*=__NLG(1.e256); exp+=256;}
        // ```
        while value.le(LD_1EM256) {
            value = value.mul(LD_1E256);
            exp += 256;
        }
        // ncbistr.cpp:2524-2540 (NEGATIVE_SCALE)
        if value.le(LD_1EM16) {
            for &(threshold, factor, step) in &NEGATIVE_SCALE {
                if value.le(threshold) {
                    value = value.mul(factor);
                    exp += step;
                    break;
                }
            }
        }
        // ncbistr.cpp:2541-2557 (NEGATIVE_MULT)
        let mut matched = false;
        for &(bound, m, step) in &NEGATIVE_MULT {
            if value.lt(bound) {
                mult = m;
                exp += step;
                matched = true;
                break;
            }
        }
        if !matched {
            mult = LD_1E8;
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2560-2577
    // ```c++
    // // get all digits
    //     long double t1 = value * mult;
    //     if (t1 >= __NLG(1.e9)) {
    //         first = 999999999;
    //     } else if (t1 < __NLG(1.e8)) {
    //         first = 100000000;
    //         t1 = first;
    //     } else {
    //         first = (unsigned int)t1;
    //     }
    //     if (high_precision) {
    //         long double t2 = (t1-first) * __NLG(1.e8);
    //         if (t2 >= __NLG(1.e8)) {
    //             second = 99999999;
    //         } else {
    //             second = (unsigned int)t2;
    //         }
    //     }
    // ```
    let mut t1 = value.mul(mult);
    if t1.ge(LD_1E9) {
        first = 999_999_999;
    } else if t1.lt(LD_1E8) {
        first = 100_000_000;
        t1 = X87::from_u32(first);
    } else {
        first = t1.to_u32_trunc();
    }
    if high_precision {
        let t2 = t1.sub_u32_exact(first).mul(LD_1E8);
        if t2.ge(LD_1E8) {
            second = 99_999_999;
        } else {
            second = t2.to_u32_trunc();
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2579-2589
    // ```c++
    // // convert them into string
    //     bool use_ext_buffer = bufsize > 20;
    //     char tmp[32];
    //     char *digits = use_ext_buffer ? buffer : tmp;
    //     char *digits_end = s_ncbi_append_int2str(digits,first,9,false);
    //     if (high_precision) {
    //         digits_end = s_ncbi_append_int2str(digits_end,second,8,true);
    //     }
    //     size_t digits_len = digits_end - digits;
    //     size_t digits_got = digits_len;
    //     size_t digits_expected = high_precision ? 17 : 9;
    // ```
    let digits = &mut out.buffer;
    let mut digits_end = s_ncbi_append_int2str(digits, 0, first, 9, false);
    if high_precision {
        digits_end = s_ncbi_append_int2str(digits, digits_end, second, 8, true);
    }
    let mut digits_len = digits_end;
    let digits_got = digits_len;
    let digits_expected = if high_precision { 17 } else { 9 };

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2591-2614
    // ```c++
    // // get significant digits according to requested precision
    //     size_t pos = precision;
    //     if (digits_len > precision) {
    //         digits_len = precision;
    //
    //         // this is questionable, but in fact,
    //         // improves the result (on average)
    // #if 1
    //         if (high_precision) {
    //             if (digits[pos] == '4') {
    //                 size_t pt = pos-1;
    //                 while (pt != 0 && digits[--pt] == '9')
    //                     ;
    //                 if (pt != 0 && (pos-pt) > precision/2)
    //                     digits[pos]='5';
    //             } else if (digits[pos] == '5') {
    //                 size_t pt = pos;
    //                 while (pt != 0 && digits[--pt] == '0')
    //                     ;
    //                 if (pt != 0 && (pos-pt) > precision/2)
    //                     digits[pos]='4';
    //             }
    //         }
    // #endif
    // ```
    let mut pos = precision;
    if digits_len > precision {
        digits_len = precision;

        if high_precision {
            if digits[pos] == b'4' {
                let mut pt = pos - 1;
                while pt != 0 {
                    pt -= 1;
                    if digits[pt] != b'9' {
                        break;
                    }
                }
                if pt != 0 && (pos - pt) > precision / 2 {
                    digits[pos] = b'5';
                }
            } else if digits[pos] == b'5' {
                let mut pt = pos;
                while pt != 0 {
                    pt -= 1;
                    if digits[pt] != b'0' {
                        break;
                    }
                }
                if pt != 0 && (pos - pt) > precision / 2 {
                    digits[pos] = b'4';
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2616-2637
        // ```c++
        //         if (digits[pos] >= '5') {
        //             do {
        //                 if (digits[--pos] < '9') {
        //                     ++digits[pos++];
        //                     break;
        //                 }
        //                 digits[pos]='0';
        //             } while (pos > 0);
        //             if (pos == 0) {
        //                 if (digits_expected <= digits_got) {
        //                     if (exp_positive) {
        //                        ++exp;
        //                     } else {
        // // exp cannot be 0, by design
        //                         exp_positive = --exp == 0;
        //                     }
        //                 }
        //                 *digits = '1';
        //                 digits_len = 1;
        //             }
        //         }
        //     }
        // ```
        if digits[pos] >= b'5' {
            loop {
                pos -= 1;
                if digits[pos] < b'9' {
                    digits[pos] += 1;
                    pos += 1;
                    break;
                }
                digits[pos] = b'0';
                if pos == 0 {
                    break;
                }
            }
            if pos == 0 {
                if digits_expected <= digits_got {
                    if exp_positive {
                        exp += 1;
                    } else {
                        exp -= 1;
                        exp_positive = exp == 0;
                    }
                }
                digits[0] = b'1';
                digits_len = 1;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2639-2643
    // ```c++
    // // truncate trailing zeros
    //     for (pos = digits_len; pos-- > 0 && digits[pos] == '0';)
    //         --digits_len;
    //
    //     *dec = exp_positive ? int(exp) : -int(exp);
    // ```
    pos = digits_len;
    while pos > 0 {
        pos -= 1;
        if digits[pos] != b'0' {
            break;
        }
        digits_len -= 1;
    }
    out.dec = if exp_positive {
        exp as i32
    } else {
        -(exp as i32)
    };
    out.len = digits_len;
    out
}

/// `NStr::DoubleToStringPosix(val, precision, buffer, bufsize)` appended to `out`, for a
/// finite `val`. (NCBI throws when `bufsize < precision + 8`; `out` grows instead, and the
/// serial callers pass 64- and 512-byte buffers with precision 15.)
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2658-2681
/// ```c++
/// SIZE_TYPE NStr::DoubleToStringPosix(double val, unsigned int precision,
///                                     char* buffer, SIZE_TYPE bufsize)
/// {
///     ...
///     int dec=0, sign=0;
///     char digits[32];
///     size_t digits_len = DoubleToString_Ecvt(
///         val, precision, digits, sizeof(digits), &dec, &sign);
///     if (digits_len == 0) {
///         errno = 0;
///         return 0;
///     }
///     if (val == 0.) {
///         strncpy(buffer,digits, digits_len);
///         return digits_len;
///     }
///     if (digits_len == 1 && dec == 0 && sign >=0) {
///         *buffer = digits[0];
///         errno = 0;
///         return 1;
///     }
/// ```
pub fn double_to_string_posix(val: f64, precision: u32, out: &mut Vec<u8>) {
    let ecvt = double_to_string_ecvt(val, precision);
    let digits = ecvt.digits();
    let mut digits_len = digits.len();
    if digits_len == 0 {
        return;
    }
    if val == 0.0 {
        out.extend_from_slice(digits);
        return;
    }
    if digits_len == 1 && ecvt.dec == 0 && ecvt.sign >= 0 {
        out.push(digits[0]);
        return;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2682-2692
    // ```c++
    //     bool exp_positive = dec >= 0;
    //     unsigned int exp= (unsigned int)(exp_positive ? dec : (-dec));
    //
    //     // assemble the result
    //     char *buffer_pos = buffer;
    // //    char *buffer_end = buffer + bufsize;
    //     char *digits_pos = digits;
    //
    //     if (sign < 0) {
    //         *buffer_pos++ = '-';
    //     }
    // ```
    let exp_positive = ecvt.dec >= 0;
    let mut exp: u32 = ecvt.dec.unsigned_abs();
    if ecvt.sign < 0 {
        out.push(b'-');
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2693-2715
    // ```c++
    //     // The 'e' format is used when the exponent of the value is less than -4
    //     //  or greater than or equal to the precision argument
    //     if ((exp_positive && exp >= precision) || (!exp_positive && exp > 4)) {
    //         *buffer_pos++ = *digits_pos++;
    //         --digits_len;
    //         if (digits_len != 0) {
    //             *buffer_pos++ = '.';
    //             strncpy(buffer_pos,digits_pos,digits_len);
    //             buffer_pos += digits_len;
    //         }
    //         *buffer_pos++ = 'e';
    //         *buffer_pos++ = exp_positive ? '+' : '-';
    //
    // //#if defined(NCBI_OS_MSWIN)
    // #if NCBI_COMPILER_MSVC && _MSC_VER < 1900
    //         ...
    // #else
    //         bool need_zeros = exp < 10 ? true : false;
    //         size_t need_digits = exp < 100 ? 2 : 3;
    // #endif
    //         // assuming exp < 1000
    //         buffer_pos = s_ncbi_append_int2str(buffer_pos, exp, need_digits,need_zeros);
    // ```
    // `precision` here is the caller's value, not the one DoubleToString_Ecvt clamps to
    // DBL_DIG. The oracle is gcc's build, so the non-MSVC exponent width applies.
    if (exp_positive && exp >= precision) || (!exp_positive && exp > 4) {
        out.push(digits[0]);
        digits_len -= 1;
        if digits_len != 0 {
            out.push(b'.');
            out.extend_from_slice(&digits[1..1 + digits_len]);
        }
        out.push(b'e');
        out.push(if exp_positive { b'+' } else { b'-' });
        let need_zeros = exp < 10;
        let need_digits = if exp < 100 { 2 } else { 3 };
        let mut tmp = [0u8; 3];
        let end = s_ncbi_append_int2str(&mut tmp, 0, exp, need_digits, need_zeros);
        out.extend_from_slice(&tmp[..end]);
    } else if exp_positive {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2716-2732
        // ```c++
        //     } else if (exp_positive) {
        //         *buffer_pos++ = *digits_pos++;
        //         --digits_len;
        //         if (digits_len > exp) {
        //             strncpy(buffer_pos,digits_pos,exp);
        //             buffer_pos += exp;
        //             *buffer_pos++ = '.';
        //             strncpy(buffer_pos,digits_pos+exp,digits_len-exp);
        //             buffer_pos += digits_len-exp;
        //         } else {
        //             strncpy(buffer_pos,digits_pos,digits_len);
        //             buffer_pos += digits_len;
        //             exp -= (unsigned int)digits_len;
        //             while (exp--) {
        //                 *buffer_pos++ = '0';
        //             }
        //         }
        // ```
        out.push(digits[0]);
        digits_len -= 1;
        let rest = &digits[1..1 + digits_len];
        let exp_digits = exp as usize;
        if digits_len > exp_digits {
            out.extend_from_slice(&rest[..exp_digits]);
            out.push(b'.');
            out.extend_from_slice(&rest[exp_digits..]);
        } else {
            out.extend_from_slice(rest);
            exp -= digits_len as u32;
            for _ in 0..exp {
                out.push(b'0');
            }
        }
    } else {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2733-2741
        // ```c++
        //     } else {
        //         *buffer_pos++ = '0';
        //         *buffer_pos++ = '.';
        //         for (--exp; exp--;) {
        //             *buffer_pos++ = '0';
        //         }
        //         strncpy(buffer_pos,digits_pos, digits_len);
        //         buffer_pos += digits_len;
        //     }
        // ```
        out.push(b'0');
        out.push(b'.');
        for _ in 1..exp {
            out.push(b'0');
        }
        out.extend_from_slice(digits);
    }
}

/// `PutInt4`: plain decimal with a leading `-`.
///
/// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:951-971
/// ```c++
///     Int4 n = v;
///     if ( n < 0 ) {
///         n = -n;
///     }
///     char* pos = b + BSIZE;
///     do {
///         Int4 a = '0'+n;
///         *--pos = char(a-10*(n/=Uint4(10)));
///     } while ( n );
///     if ( v < 0 ) {
///         *--pos = '-';
///     }
/// ```
fn put_int4(out: &mut Vec<u8>, v: i32) {
    let mut b = [0u8; 12];
    let mut pos = b.len();
    let mut n = v.unsigned_abs();
    loop {
        pos -= 1;
        b[pos] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    if v < 0 {
        pos -= 1;
        b[pos] = b'-';
    }
    out.extend_from_slice(&b[pos..]);
}

/// `CObjectOStreamAsn::WriteDouble(data)` (ASN.1 text, outfmt 8 and 11): `{ m, 10, e }` with
/// `DBL_DIG` digits, appended to `out`.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:175-221
/// ```c++
/// void CObjectOStreamAsn::WriteDouble2(double data, unsigned digits)
/// {
///     ...
///     if (isnan(data)) {
///         m_Output.PutString("NOT-A-NUMBER");
///         return;
///     }
///     if (!finite(data)) {
///         if (data > 0) {
///             m_Output.PutString("PLUS-INFINITY");
///         } else {
///             m_Output.PutString("MINUS-INFINITY");
///         }
///         return;
///     }
/// #endif
///     if ( data == 0.0 ) {
///         double zero = 0.;
///         if (memcmp(&data, &zero, sizeof(double)) == 0) {
///             m_Output.PutString("{ 0, 10, 0 }");
///         } else {
///             m_Output.PutString("{ -0, 10, 0 }");
///         }
///         return;
///     }
///
///     char buffer[128];
///     if (m_FastWriteDouble) {
///         int dec, sign;
///         size_t len = NStr::DoubleToString_Ecvt(
///             data, digits, buffer, sizeof(buffer), &dec, &sign);
///         _ASSERT(len > 0);
///         m_Output.PutString("{ ");
///         if (sign < 0) {
///             m_Output.PutString("-");
///         }
///         m_Output.PutString(buffer,len);
///         m_Output.PutString(", 10, ");
///         m_Output.PutInt4(dec - (int)(len-1));
///         m_Output.PutString(" }");
/// ```
pub fn write_asn_text_double(data: f64, out: &mut Vec<u8>) {
    if data.is_nan() {
        out.extend_from_slice(b"NOT-A-NUMBER");
        return;
    }
    if !data.is_finite() {
        if data > 0.0 {
            out.extend_from_slice(b"PLUS-INFINITY");
        } else {
            out.extend_from_slice(b"MINUS-INFINITY");
        }
        return;
    }
    if data == 0.0 {
        if data.to_bits() == 0 {
            out.extend_from_slice(b"{ 0, 10, 0 }");
        } else {
            out.extend_from_slice(b"{ -0, 10, 0 }");
        }
        return;
    }
    let ecvt = double_to_string_ecvt(data, DBL_DIG);
    let digits = ecvt.digits();
    out.extend_from_slice(b"{ ");
    if ecvt.sign < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(digits);
    out.extend_from_slice(b", 10, ");
    put_int4(out, ecvt.dec - (digits.len() as i32 - 1));
    out.extend_from_slice(b" }");
}

/// `CObjectOStreamJson::WriteDouble2` refuses NaN and infinities (`ThrowError(fInvalidData,
/// ...)`, a `CSerialException`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum JsonDoubleError {
    /// "invalid double: not a number"
    NotANumber,
    /// "invalid double: infinite"
    Infinite,
}

impl JsonDoubleError {
    /// NCBI's message text.
    pub fn message(self) -> &'static str {
        match self {
            JsonDoubleError::NotANumber => "invalid double: not a number",
            JsonDoubleError::Infinite => "invalid double: infinite",
        }
    }
}

/// `CObjectOStreamJson::WriteDouble(data)` (outfmt 12, 13, 15): `DoubleToStringPosix(data,
/// DBL_DIG)` appended to `out`.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrjson.cpp:205-225
/// ```c++
/// void CObjectOStreamJson::WriteDouble(double data)
/// {
///     WriteDouble2(data,DBL_DIG);
/// }
///
/// void CObjectOStreamJson::WriteDouble2(double data, unsigned digits)
/// {
///     if (isnan(data)) {
///         ThrowError(fInvalidData, "invalid double: not a number");
///     }
///     if (!finite(data)) {
///         ThrowError(fInvalidData, "invalid double: infinite");
///     }
///     if (m_FastWriteDouble) {
///         char buffer[64];
///         WriteKeywordValue( string(buffer,
///             NStr::DoubleToStringPosix(data, digits, buffer, sizeof(buffer))));
///     } else {
/// ```
pub fn write_json_double(data: f64, out: &mut Vec<u8>) -> Result<(), JsonDoubleError> {
    if data.is_nan() {
        return Err(JsonDoubleError::NotANumber);
    }
    if !data.is_finite() {
        return Err(JsonDoubleError::Infinite);
    }
    double_to_string_posix(data, DBL_DIG, out);
    Ok(())
}

/// `CObjectOStreamXml::WriteDouble(data)` with the default scientific real format (outfmt 5,
/// 14, 16): `DoubleToStringPosix(data, DBL_DIG)` appended to `out`.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrxml.cpp:661-713
/// ```c++
/// void CObjectOStreamXml::WriteDouble2(double data, unsigned digits)
/// {
///     if (m_SpecialCaseWrite && x_SpecialCaseWrite()) {
///         return;
///     }
///     if (isnan(data)) {
///         m_Output.PutString("NaN", 3);
///         return;
///     }
///     if (!finite(data)) {
///         if (data < 0) {
///             m_Output.PutChar('-');
///         }
///         m_Output.PutString("INF", 3);
///         return;
///     }
///     char buffer[512];
///     SIZE_TYPE width;
///     if (m_RealFmt == eRealFixedFormat) {
///         ...
///     } else {
///         if (m_FastWriteDouble) {
///             width = NStr::DoubleToStringPosix(data, digits, buffer, sizeof(buffer));
///         } else {
///             ...
///         }
///     }
///     m_Output.PutString(buffer, width);
/// }
///
/// void CObjectOStreamXml::WriteDouble(double data)
/// {
///     WriteDouble2(data, DBL_DIG);
/// }
/// ```
/// (`m_RealFmt` starts as `eRealScientificFormat`, objostrxml.cpp:80; the special-case write
/// belongs to the element writer that calls this.)
pub fn write_xml_double(data: f64, out: &mut Vec<u8>) {
    if data.is_nan() {
        out.extend_from_slice(b"NaN");
        return;
    }
    if !data.is_finite() {
        if data < 0.0 {
            out.push(b'-');
        }
        out.extend_from_slice(b"INF");
        return;
    }
    double_to_string_posix(data, DBL_DIG, out);
}
