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
    let mut s = String::new();
    s += "Cannot convert string '";
    s += &String::from_utf8_lossy(&printable_string(str_));
    s += "' to ";
    s += "Uint8";
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:569-575
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
    let mut str_ptr: usize = 0;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1572-1579
    const SIGN: ThrowSite = data_size_site(1578);
    if str_ptr < str_end && str_[str_ptr] == b'+' {
        str_ptr += 1;
    } else if str_ptr < str_end && str_[str_ptr] == b'-' {
        // (flags & fMandatorySign) is 0
        return Err(err(SIGN, str_ptr));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1581-1585
    let num_start = str_ptr;
    let mut have_dot = false;
    let mut digs_pre_dot: u32 = 0;
    let mut digs_post_dot: u32 = 0;

    // The `c == ','` branch (lines 1604-1611) needs fAllowCommas and is not ported.
    // allow_dot is true because fDS_ProhibitFractions is 0.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1587-1618
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
    const TRAILING: ThrowSite = data_size_site(1674);
    if str_ptr != str_end {
        return Err(err(TRAILING, str_ptr));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1677-1688
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
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:1708-1708
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
#[path = "data_size_tests.rs"]
mod tests;
