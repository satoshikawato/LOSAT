//! `NStr::StringToUInt8_DataSize(str, flags = 0)` (corelib/ncbistr.cpp), which makeblastdb calls for
//! the value of `-max_file_sz`.
//!
//! Only the default-flags behaviour of `NStr::StringToUInt8_DataSize(str)` is ported (the two-argument overload, `flags = 0`).
//! The deprecated three-argument overload (`flags`, `base`, ncbistr.cpp:1460) is not what makeblastdb calls and is not ported.
//! With `flags == 0` these flag-dependent branches of the C++ never run and are therefore not ported:
//! `fConvErr_NoThrow` (errors always throw), `fMandatorySign`, `fAllowCommas`, `fAllowLeadingSymbols`,
//! `fAllowTrailingSymbols`, `fDS_ForceBinary`, `fDS_ProhibitFractions` (so a '.' is always a decimal point) and
//! `fDS_ProhibitSpaceBeforeSuffix` (so one space before the suffix is always accepted).
//! Also not ported: the "Wrong set of flags" check and the errno bookkeeping of `CS2N_Guard`.
//!
//! C `isdigit`, `isprint` and `toupper` are the ones of the "C" locale (ASCII only); the C++ works on bytes, so this port works on the bytes of the `&str`.

use super::printable_string;
use crate::ncbi_exception::{NcbiException, ThrowSite};

/// The `CStringException` (`eConvert`) NCBI throws from `StringToUInt8_DataSize`.
///
/// `site` is where it is thrown, `message` the text NCBI builds (`CS2N_Guard::Message`, the
/// exception's `GetMsg()`), `pos` the position given to `NCBI_THROW2`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataSizeError {
    pub site: ThrowSite,
    pub message: String,
    pub pos: usize,
}

impl DataSizeError {
    /// The exception: makeblastdb reports it through its general handler (`Error: ` and
    /// `what()`, exit status 255).
    pub fn exception(&self) -> NcbiException {
        NcbiException::with_pos(self.site, self.message.as_str(), self.pos)
    }
}

/// A throw of `S2N_CONVERT_ERROR(Uint8, ...)` at `line` of `NStr::StringToUInt8_DataSize`.
const fn data_size_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/corelib/ncbistr.cpp",
        line,
        kind: "CStringException::eConvert",
        function: Some("ncbi::NStr::StringToUInt8_DataSize"),
    }
}

/// `err_guard.Message(str, "Uint8", kEmptyStr)` as thrown by `S2N_CONVERT_ERROR(Uint8, kEmptyStr, errcode, pos)`
/// at `site`: a `CStringException` with this message and position `pos` (the `msg` argument is always empty in the
/// ported paths).
fn convert_error(str_: &[u8], site: ThrowSite, pos: usize) -> DataSizeError {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:488-558
    // ```c++
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
    let mut s = String::new();
    s += "Cannot convert string '";
    s += &String::from_utf8_lossy(&printable_string(str_));
    s += "' to ";
    s += "Uint8";
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
    DataSizeError {
        site,
        message: s,
        pos,
    }
}

/// `NStr::StringToUInt8_DataSize(text)` with default flags (0).
pub fn string_to_uint8_data_size(text: &str) -> Result<u64, DataSizeError> {
    let str_ = text.as_bytes();
    let str_end = str_.len();
    let err = |site: ThrowSite, pos: usize| convert_error(str_, site, pos);

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1552-1553
    // ```c++
    //     const char* str_ptr = str.data();
    //     const char* str_end = str_ptr + str.size();
    // ```
    let mut str_ptr: usize = 0;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1572-1579
    // ```c++
    //     if (str_ptr < str_end  &&  *str_ptr == '+') {
    //         ++str_ptr;
    //     }
    //     else if ((str_ptr < str_end  &&  *str_ptr == '-')
    //              ||  (flags & fMandatorySign))
    //     {
    //         S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    //     }
    // ```
    const SIGN: ThrowSite = data_size_site(1578);
    if str_ptr < str_end && str_[str_ptr] == b'+' {
        str_ptr += 1;
    } else if str_ptr < str_end && str_[str_ptr] == b'-' {
        // (flags & fMandatorySign) is 0
        return Err(err(SIGN, str_ptr));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1581-1585
    // ```c++
    //     const char* num_start = str_ptr;
    //     bool have_dot = false;
    //     bool allow_commas = (flags & fAllowCommas) != 0;
    //     bool allow_dot = (flags & fDS_ProhibitFractions) == 0;
    //     Uint4 digs_pre_dot = 0, digs_post_dot = 0;
    // ```
    let num_start = str_ptr;
    let mut have_dot = false;
    let mut digs_pre_dot: u32 = 0;
    let mut digs_post_dot: u32 = 0;

    // The `c == ','` branch (lines 1604-1611) needs fAllowCommas and is not ported.
    // allow_dot is true because fDS_ProhibitFractions is 0.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1587-1618
    // ```c++
    //     for (; str_ptr < str_end; ++str_ptr) {
    //         char c = *str_ptr;
    //         if (isdigit(c)) {
    //             if (have_dot)
    //                 ++digs_post_dot;
    //             else
    //                 ++digs_pre_dot;
    //         }
    //         else if (c == '.'  &&  allow_dot) {
    //             if (have_dot  ||  str_ptr == num_start)
    //                 break;
    //             if (*(str_ptr - 1) == ',') {
    //                 --str_ptr;
    //                 break;
    //             }
    //             have_dot = true;
    //         }
    // ...
    //         else
    //             break;
    //     }
    //     if (have_dot  &&  digs_post_dot == 0)
    //         --str_ptr;
    //     else if (str_ptr > num_start  &&  *(str_ptr - 1) == ',')
    //         --str_ptr;
    // ```
    while str_ptr < str_end {
        let c = str_[str_ptr];
        if c.is_ascii_digit() {
            if have_dot {
                digs_post_dot += 1;
            } else {
                digs_pre_dot += 1;
            }
        } else if c == b'.' {
            if have_dot || str_ptr == num_start {
                break;
            }
            if str_[str_ptr - 1] == b',' {
                str_ptr -= 1;
                break;
            }
            have_dot = true;
        } else {
            break;
        }
        str_ptr += 1;
    }
    if have_dot && digs_post_dot == 0 {
        str_ptr -= 1;
    }
    // (the `else if (... *(str_ptr - 1) == ',')` arm cannot fire without fAllowCommas)

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1620-1631
    // ```c++
    //     const char* num_end = str_ptr;
    //     if (num_start == num_end) {
    //         S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    //     }
    //     if (str_ptr < str_end  &&  *str_ptr == ' '
    //         &&  !(flags & fDS_ProhibitSpaceBeforeSuffix))
    //     {
    //         ++str_ptr;
    //     }
    //     char suff_c = 0;
    //     if (str_ptr < str_end)
    //         suff_c = (char)toupper(*str_ptr);
    // ```
    const NO_DIGIT: ThrowSite = data_size_site(1622);
    let num_end = str_ptr;
    if num_start == num_end {
        return Err(err(NO_DIGIT, str_ptr));
    }
    // fDS_ProhibitSpaceBeforeSuffix is 0
    if str_ptr < str_end && str_[str_ptr] == b' ' {
        str_ptr += 1;
    }
    let mut suff_c: u8 = 0;
    if str_ptr < str_end {
        suff_c = str_[str_ptr].to_ascii_uppercase();
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1633-1645
    // ```c++
    //     static const char s_Suffixes[] = {'K', 'M', 'G', 'T', 'P', 'E'};
    //     static const char* const s_BinCoefs[] = {"1024", "1048576", "1073741824",
    //                                              "1099511627776",
    //                                              "1125899906842624",
    //                                              "1152921504606846976"};
    //     static const Uint4 s_NumSuffixes = (Uint4)(sizeof(s_Suffixes) / sizeof(s_Suffixes[0]));
    //
    //     bool binary_suff = (flags & fDS_ForceBinary) != 0;
    //     Uint4 suff_idx = 0;
    //     for (; suff_idx < s_NumSuffixes; ++suff_idx) {
    //         if (suff_c == s_Suffixes[suff_idx])
    //             break;
    //     }
    // ```
    const S_SUFFIXES: [u8; 6] = [b'K', b'M', b'G', b'T', b'P', b'E'];
    const S_BIN_COEFS: [&str; 6] = [
        "1024",
        "1048576",
        "1073741824",
        "1099511627776",
        "1125899906842624",
        "1152921504606846976",
    ];
    const S_NUM_SUFFIXES: u32 = S_SUFFIXES.len() as u32;

    // fDS_ForceBinary is 0
    let mut binary_suff = false;
    let mut suff_idx: u32 = 0;
    while suff_idx < S_NUM_SUFFIXES {
        if suff_c == S_SUFFIXES[suff_idx as usize] {
            break;
        }
        suff_idx += 1;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1646-1661
    // ```c++
    //     if (suff_idx < s_NumSuffixes) {
    //         ++str_ptr;
    //         if (str_ptr + 1 < str_end  &&  toupper(*str_ptr) == 'I'
    //             &&  toupper(*(str_ptr + 1)) == 'B')
    //         {
    //             str_ptr += 2;
    //             binary_suff = true;
    //         }
    //         else if (str_ptr < str_end  &&  toupper(*str_ptr) == 'B')
    //             ++str_ptr;
    //     }
    //     else if (suff_c == 'B') {
    //         ++str_ptr;
    //     }
    //     else if (*(str_ptr - 1) == ' ')
    //         --str_ptr;
    // ```
    if suff_idx < S_NUM_SUFFIXES {
        str_ptr += 1;
        if str_ptr + 1 < str_end
            && str_[str_ptr].eq_ignore_ascii_case(&b'I')
            && str_[str_ptr + 1].eq_ignore_ascii_case(&b'B')
        {
            str_ptr += 2;
            binary_suff = true;
        } else if str_ptr < str_end && str_[str_ptr].eq_ignore_ascii_case(&b'B') {
            str_ptr += 1;
        }
    } else if suff_c == b'B' {
        str_ptr += 1;
    } else if str_[str_ptr - 1] == b' ' {
        str_ptr -= 1;
    }

    // fAllowTrailingSymbols is 0 (lines 1663-1672 are not ported)
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1673-1675
    // ```c++
    //     if (str_ptr != str_end) {
    //         S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    //     }
    // ```
    const TRAILING: ThrowSite = data_size_site(1674);
    if str_ptr != str_end {
        return Err(err(TRAILING, str_ptr));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1677-1688
    // ```c++
    //     Uint4 orig_digs = digs_pre_dot + digs_post_dot;
    //     AutoArray<Uint1> orig_num(orig_digs);
    //     str_ptr = num_start;
    //     for (Uint4 i = 0; str_ptr < num_end; ++str_ptr) {
    //         if (*str_ptr == ','  ||  *str_ptr == '.')
    //             continue;
    //         orig_num[i++] = Uint1(*str_ptr - '0');
    //     }
    //
    //     Uint1* num_to_conv = orig_num.get();
    //     Uint4 digs_to_conv = digs_pre_dot;
    //     AutoArray<Uint1> mul_num;
    // ```
    let orig_digs: u32 = digs_pre_dot + digs_post_dot;
    let mut orig_num: Vec<u8> = vec![0; orig_digs as usize];
    str_ptr = num_start;
    let mut k: usize = 0;
    while str_ptr < num_end {
        if str_[str_ptr] == b',' || str_[str_ptr] == b'.' {
            str_ptr += 1;
            continue;
        }
        orig_num[k] = str_[str_ptr] - b'0';
        k += 1;
        str_ptr += 1;
    }

    // num_to_conv is `mul_num[num_off..]` when `use_mul` is set, else `orig_num[num_off..]`
    // (the C++ points `num_to_conv` into one of the two arrays).
    let mut use_mul = false;
    let mut num_off: usize = 0;
    let mut digs_to_conv: u32 = digs_pre_dot;
    let mut mul_num: Vec<u8> = Vec::new();
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1689-1725
    // ```c++
    //     if (binary_suff  &&  suff_idx < s_NumSuffixes) {
    //         const char* coef = s_BinCoefs[suff_idx];
    //         Uint4 coef_size = Uint4(strlen(coef));
    //         mul_num = new Uint1[orig_digs + coef_size];
    //         memset(mul_num.get(), 0, orig_digs + coef_size);
    //         for (Uint4 coef_i = 0; coef_i < coef_size; ++coef_i) {
    //             Uint1 coef_d = Uint1(coef[coef_i] - '0');
    //             Uint1 carry = 0;
    //             Uint4 res_idx = orig_digs + coef_i;
    //             for (int orig_i = orig_digs - 1; orig_i >= 0; --orig_i, --res_idx) {
    //                 Uint1 orig_d = orig_num[orig_i];
    //                 Uint1 res_d = Uint1(coef_d * orig_d + carry + mul_num[res_idx]);
    //                 carry = 0;
    //                 while (res_d >= 10) {
    //                     res_d = (Uint1)(res_d - 10); // res_d -= 10;
    //                     ++carry;
    //                 }
    //                 mul_num[res_idx] = res_d;
    //             }
    //             _ASSERT(carry <= 9);
    //             for (; carry != 0; --res_idx) {
    //                 Uint1 res_d = Uint1(mul_num[res_idx] + carry);
    //                 carry = 0;
    //                 while (res_d >= 10) {
    //                     res_d = (Uint1)(res_d - 10); // res_d -= 10;
    //                     ++carry;
    //                 }
    //                 mul_num[res_idx] = res_d;
    //             }
    //         }
    //         digs_to_conv = orig_digs + coef_size - digs_post_dot;
    //         num_to_conv = mul_num.get();
    //         while (digs_to_conv > 1  &&  *num_to_conv == 0) {
    //             --digs_to_conv;
    //             ++num_to_conv;
    //         }
    //     }
    // ```
    if binary_suff && suff_idx < S_NUM_SUFFIXES {
        let coef = S_BIN_COEFS[suff_idx as usize].as_bytes();
        let coef_size = coef.len() as u32;
        mul_num = vec![0; (orig_digs + coef_size) as usize];
        for coef_i in 0..coef_size {
            let coef_d: u8 = coef[coef_i as usize] - b'0';
            let mut carry: u8 = 0;
            let mut res_idx: u32 = orig_digs + coef_i;
            let mut orig_i: i32 = orig_digs as i32 - 1;
            while orig_i >= 0 {
                let orig_d: u8 = orig_num[orig_i as usize];
                let mut res_d: u8 = coef_d * orig_d + carry + mul_num[res_idx as usize];
                carry = 0;
                while res_d >= 10 {
                    res_d -= 10;
                    carry += 1;
                }
                mul_num[res_idx as usize] = res_d;
                orig_i -= 1;
                res_idx = res_idx.wrapping_sub(1);
            }
            // _ASSERT(carry <= 9);
            while carry != 0 {
                let mut res_d: u8 = mul_num[res_idx as usize] + carry;
                carry = 0;
                while res_d >= 10 {
                    res_d -= 10;
                    carry += 1;
                }
                mul_num[res_idx as usize] = res_d;
                res_idx = res_idx.wrapping_sub(1);
            }
        }
        digs_to_conv = orig_digs + coef_size - digs_post_dot;
        use_mul = true;
        num_off = 0;
        while digs_to_conv > 1 && mul_num[num_off] == 0 {
            digs_to_conv -= 1;
            num_off += 1;
        }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1726-1742
    // ```c++
    //     else if (suff_idx < s_NumSuffixes) {
    //         Uint4 coef_size = (suff_idx + 1) * 3;
    //         if (coef_size <= digs_post_dot) {
    //             digs_to_conv += coef_size;
    //             digs_post_dot -= coef_size;
    //         }
    //         else {
    //             digs_to_conv += digs_post_dot;
    //             coef_size -= digs_post_dot;
    //             digs_post_dot = 0;
    //             mul_num = new Uint1[digs_to_conv + coef_size];
    //             memmove(mul_num.get(), num_to_conv, digs_to_conv);
    //             memset(mul_num.get() + digs_to_conv, 0, coef_size);
    //             num_to_conv = mul_num.get();
    //             digs_to_conv += coef_size;
    //         }
    //     }
    // ```
    } else if suff_idx < S_NUM_SUFFIXES {
        let mut coef_size: u32 = (suff_idx + 1) * 3;
        if coef_size <= digs_post_dot {
            digs_to_conv += coef_size;
            digs_post_dot -= coef_size;
        } else {
            digs_to_conv += digs_post_dot;
            coef_size -= digs_post_dot;
            digs_post_dot = 0;
            mul_num = vec![0; (digs_to_conv + coef_size) as usize];
            mul_num[..digs_to_conv as usize].copy_from_slice(&orig_num[..digs_to_conv as usize]);
            use_mul = true;
            num_off = 0;
            digs_to_conv += coef_size;
        }
    }
    let num_to_conv: &[u8] = if use_mul { &mul_num } else { &orig_num };
    let num_to_conv = &num_to_conv[num_off..];

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1744-1754
    // ```c++
    //     const Uint8 limdiv = kMax_UI8/10;
    //     const int   limoff = int(kMax_UI8 % 10);
    //     Uint8 n = 0;
    //     for (Uint4 i = 0; i < digs_to_conv; ++i) {
    //         Uint1 d = num_to_conv[i];
    //         if (n >= limdiv  &&  (n > limdiv  ||  d > limoff)) {
    //             S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, i);
    //         }
    //         n *= 10;
    //         n += d;
    //     }
    // ```
    const LIMDIV: u64 = u64::MAX / 10;
    const LIMOFF: i32 = (u64::MAX % 10) as i32;
    const OVERFLOW: ThrowSite = data_size_site(1750);
    let mut n: u64 = 0;
    for i in 0..digs_to_conv {
        let d: u8 = num_to_conv[i as usize];
        if n >= LIMDIV && (n > LIMDIV || d as i32 > LIMOFF) {
            return Err(err(OVERFLOW, i as usize));
        }
        n *= 10;
        n += d as u64;
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1755-1762
    // ```c++
    //     if (digs_post_dot != 0  &&  num_to_conv[digs_to_conv] >= 5) {
    //         if (n == kMax_UI8) {
    //             S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, digs_to_conv);
    //         }
    //         ++n;
    //     }
    //     return n;
    // }
    // ```
    const ROUNDING_OVERFLOW: ThrowSite = data_size_site(1757);
    if digs_post_dot != 0 && num_to_conv[digs_to_conv as usize] >= 5 {
        if n == u64::MAX {
            return Err(err(ROUNDING_OVERFLOW, digs_to_conv as usize));
        }
        n += 1;
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values the source makes certain (decimal digits scaled by the suffix, 'K' = 1000, "KiB" = 1024).
    #[test]
    fn parsed_values() {
        let cases: &[(&str, u64)] = &[
            ("0", 0),
            ("300", 300),
            ("+5", 5),
            ("00010", 10),
            ("300B", 300),
            ("1KB", 1000),
            ("1K", 1000),
            ("1k", 1000),
            ("1 KB", 1000),
            ("1kb", 1000),
            ("1KiB", 1024),
            ("1kib", 1024),
            ("1 KiB", 1024),
            ("1MiB", 1048576),
            ("1MB", 1_000_000),
            ("3GB", 3_000_000_000),
            ("4GB", 4_000_000_000),
            ("4GiB", 4 << 30),
            ("1TB", 1_000_000_000_000),
            ("1TiB", 1 << 40),
            ("1PB", 1_000_000_000_000_000),
            ("1PiB", 1 << 50),
            ("1EB", 1_000_000_000_000_000_000),
            ("1EiB", 1 << 60),
            ("15EiB", 15 << 60),
            ("4294967295", 4294967295),
            ("18446744073709551615", u64::MAX),
            ("18446744073709551.615K", u64::MAX),
            // a fraction is rounded half up on its first dropped digit
            ("1.4", 1),
            ("1.5", 2),
            ("0.5", 1),
            ("0.4", 0),
            ("1.5KB", 1500),
            ("1.5KiB", 1536),
            ("1.5MiB", 1572864),
            ("2.5GiB", 2684354560),
            ("0.0005K", 1),
            ("0.0004999K", 0),
            ("1.2345K", 1235),
            ("1.2344K", 1234),
            ("10239.5", 10240),
            ("10239.4", 10239),
            ("18446744073709551615.4", u64::MAX),
        ];
        for (input, want) in cases {
            assert_eq!(
                string_to_uint8_data_size(input),
                Ok(*want),
                "input {:?}",
                input
            );
        }
    }

    #[test]
    fn overflow_errors() {
        // ERANGE: m_Pos is the index of the digit that overflows
        let e = string_to_uint8_data_size("18446744073709551616").unwrap_err();
        assert_eq!(
            e.message,
            "Cannot convert string '18446744073709551616' to Uint8"
        );
        assert_eq!(e.pos, 19);
        assert_eq!(
            string_to_uint8_data_size("18446744073709551615.5")
                .unwrap_err()
                .pos,
            20
        );
        assert_eq!(string_to_uint8_data_size("16EiB").unwrap_err().pos, 19);
        assert_eq!(
            string_to_uint8_data_size("18446744073709552K")
                .unwrap_err()
                .pos,
            19
        );
    }

    #[test]
    fn error_text_and_position() {
        let e = string_to_uint8_data_size("abc").unwrap_err();
        assert_eq!(e.message, "Cannot convert string 'abc' to Uint8");
        assert_eq!(e.pos, 0);
        // '-' gives m_Pos of the sign; trailing garbage gives m_Pos of the first unconsumed byte
        assert_eq!(string_to_uint8_data_size("-5").unwrap_err().pos, 0);
        assert_eq!(string_to_uint8_data_size("1KBB").unwrap_err().pos, 3);
        assert_eq!(string_to_uint8_data_size("").unwrap_err().pos, 0);
    }
}
