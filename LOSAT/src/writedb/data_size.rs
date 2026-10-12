//! Port of NCBI BLAST+ 2.17.0 `NStr::StringToUInt8_DataSize(str, flags = 0)` (corelib/ncbistr.cpp) and the file-local `Uint8ToString_DataSize(v, minprec = 10)` of makeblastdb.cpp, as used for the `-max_file_sz` option and its "Maximum file size: ..." log line.
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

/// The `CStringException` (`eConvert`) NCBI throws from `StringToUInt8_DataSize`.
///
/// `message` is the text NCBI builds (`CS2N_Guard::Message`), without the
/// `ncbi::NStr::StringToUInt8_DataSize() - ` prefix and without the ` (m_Pos = N)` suffix
/// that the exception printer adds. `pos` is the position given to `NCBI_THROW2`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataSizeError {
    pub message: String,
    pub pos: usize,
}

/// `NStr::PrintableString(str)` with its default mode `fNewLine_Quote | fNonAscii_Passthru` (both 0) and language C,
/// which is what `CS2N_Guard::Message` uses on the input string.
/// Ported branches: `s_PrintableString` with `mode == 0`, `lang == eLanguage_C`
/// (so `fNewLine_Passthru`, `fNonAscii_Quote`, `fPrintable_Full` and `eLanguage_Javascript` are not ported).
fn printable_string(str_: &[u8]) -> String {
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

/// `err_guard.Message(str, "Uint8", kEmptyStr)` as thrown by `S2N_CONVERT_ERROR(Uint8, kEmptyStr, errcode, pos)`:
/// a `CStringException` with this message and position `pos` (the `msg` argument is always empty in the ported paths).
fn convert_error(str_: &[u8], pos: usize) -> DataSizeError {
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
    s += &printable_string(str_);
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
    DataSizeError { message: s, pos }
}

/// `NStr::StringToUInt8_DataSize(text)` with default flags (0).
pub fn string_to_uint8_data_size(text: &str) -> Result<u64, DataSizeError> {
    let str_ = text.as_bytes();
    let str_end = str_.len();
    let err = |pos: usize| convert_error(str_, pos);

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
    if str_ptr < str_end && str_[str_ptr] == b'+' {
        str_ptr += 1;
    } else if str_ptr < str_end && str_[str_ptr] == b'-' {
        // (flags & fMandatorySign) is 0
        return Err(err(str_ptr));
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
    let num_end = str_ptr;
    if num_start == num_end {
        return Err(err(str_ptr));
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
            && str_[str_ptr].to_ascii_uppercase() == b'I'
            && str_[str_ptr + 1].to_ascii_uppercase() == b'B'
        {
            str_ptr += 2;
            binary_suff = true;
        } else if str_ptr < str_end && str_[str_ptr].to_ascii_uppercase() == b'B' {
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
    if str_ptr != str_end {
        return Err(err(str_ptr));
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
    let mut n: u64 = 0;
    for i in 0..digs_to_conv {
        let d: u8 = num_to_conv[i as usize];
        if n >= LIMDIV && (n > LIMDIV || d as i32 > LIMOFF) {
            return Err(err(i as usize));
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
    if digs_post_dot != 0 && num_to_conv[digs_to_conv as usize] >= 5 {
        if n == u64::MAX {
            return Err(err(digs_to_conv as usize));
        }
        n += 1;
    }
    Ok(n)
}

/// makeblastdb's `Uint8ToString_DataSize(v, minprec = 10)`.
///
/// The loop in the C++ has no `break`, so it always runs all 8 iterations and `i` ends as `kMods.size()` (8)
/// whatever `v` is; `rv.append(kMods, 8, 1)` then appends nothing (pos == size is allowed), so only "B" is appended
/// and the K/M/G/... letter never appears. And a value below `minprec * 1024` (10240) is divided by 1024
/// in the first iteration, and the quotient (below 10) stays below 10240 and is divided again in each of the
/// remaining iterations, so it ends as 0; values from 10240 up are never divided.
pub fn uint8_to_string_data_size(value: u64) -> String {
    const MINPREC: u32 = 10;
    let mut v = value;
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:343-350
    // ```c++
    //     static string kMods = "KMGTPEZY";
    //
    //     size_t i(0);
    //     for(i = 0; i < kMods.size(); i++) {
    //         if (v < Uint8(minprec)*1024) {
    //             v /= 1024;
    //         }
    //     }
    // ```
    const K_MODS: &str = "KMGTPEZY";

    let mut i: usize = 0;
    while i < K_MODS.len() {
        if v < (MINPREC as u64) * 1024 {
            v /= 1024;
        }
        i += 1;
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:352-359
    // ```c++
    //     string rv = NStr::UInt8ToString(v);
    //
    //     if (i) {
    //         rv.append(kMods, i, 1);
    //         rv.append("B");
    //     }
    //
    //     return rv;
    // ```
    let mut rv = v.to_string();

    if i != 0 {
        // rv.append(kMods, i, 1): pos == size appends nothing (pos > size would throw, but i <= size)
        rv.push_str(K_MODS.get(i..i + 1).unwrap_or(""));
        rv.push('B');
    }

    rv
}

#[cfg(test)]
mod tests {
    use super::*;

    // Oracle pairs: each input was given to makeblastdb 2.17.0 as `-max_file_sz <input>`
    // (`-blastdb_version 4 -dbtype prot`), and the expected strings are what it printed.

    /// Inputs whose log line is "Maximum file size: <expected>".
    const ORACLE_LOG: &[(&str, &str)] = &[
        ("0", "0B"),
        ("1", "0B"),
        ("5", "0B"),
        ("7", "0B"),
        ("10", "0B"),
        ("9", "0B"),
        ("10239", "0B"),
        ("10240", "10240B"),
        ("10241", "10241B"),
        ("10000", "0B"),
        ("10239B", "0B"),
        ("10240B", "10240B"),
        ("1023", "0B"),
        ("1024", "0B"),
        ("2048", "0B"),
        ("4096", "0B"),
        ("12345", "12345B"),
        ("300B", "0B"),
        ("300", "0B"),
        ("1KB", "0B"),
        ("1000", "0B"),
        ("1K", "0B"),
        ("1 KB", "0B"),
        ("1kb", "0B"),
        ("1KiB", "0B"),
        ("1.5KB", "0B"),
        ("1MiB", "1048576B"),
        ("3GB", "3000000000B"),
        ("4GB", "4000000000B"),
        ("4294967295", "4294967295B"),
        ("+5", "0B"),
        ("5.0", "0B"),
        ("5.5", "0B"),
        ("0.5", "0B"),
        ("0.4", "0B"),
        ("1.4999", "0B"),
        ("2.5", "0B"),
        ("10.5", "0B"),
        ("1.5", "0B"),
        ("1.50", "0B"),
        ("1.4", "0B"),
        ("00010", "0B"),
        ("0000000000000000000000000000000000012", "0B"),
        ("1B", "0B"),
        ("1b", "0B"),
        ("1 B", "0B"),
        ("1 b", "0B"),
        ("9999K", "9999000B"),
        ("10M", "10000000B"),
        ("10m", "10000000B"),
        ("10Mb", "10000000B"),
        ("10MB", "10000000B"),
        ("10MiB", "10485760B"),
        ("10mib", "10485760B"),
        ("10mIb", "10485760B"),
        ("10MIB", "10485760B"),
        ("4096M", "4096000000B"),
        ("4095M", "4095000000B"),
        ("4095MiB", "4293918720B"),
        ("4194303K", "4194303000B"),
        ("4194304K", "4194304000B"),
        ("4194303KiB", "4294966272B"),
        ("0T", "0B"),
        ("0TiB", "0B"),
        ("0EiB", "0B"),
        ("0.5KiB", "0B"),
        ("0.5K", "0B"),
        ("0.0004K", "0B"),
        ("0.0005K", "0B"),
        ("0.0004999K", "0B"),
        ("0.0005KiB", "0B"),
        ("0.0004999KiB", "0B"),
        ("1.0005K", "0B"),
        ("1.2345K", "0B"),
        ("1.2345678K", "0B"),
        ("1.2345678KiB", "0B"),
        ("0.9999M", "999900B"),
        ("0.0009999M", "0B"),
        ("1.5MiB", "1572864B"),
        ("1.5M", "1500000B"),
        ("2.5GiB", "2684354560B"),
        ("2.5G", "2500000000B"),
        ("2.5Gb", "2500000000B"),
        ("1 KiB", "0B"),
        ("10239.5", "10240B"),
        ("10239.4", "0B"),
        ("10240.4", "10240B"),
        ("9.9999K", "0B"),
        ("10.24K", "10240B"),
        ("10.2399K", "10240B"),
        ("10.2394K", "0B"),
        ("10KiB", "10240B"),
        ("10.0KiB", "10240B"),
        ("9.99KiB", "0B"),
        ("9.999KiB", "0B"),
        ("9.9999KiB", "10240B"),
        ("10.2400KiB", "10486B"),
        ("9.7656KiB", "0B"),
        ("9.7657KiB", "0B"),
        ("10239.4999", "0B"),
        ("0.00001KiB", "0B"),
        ("1k", "0B"),
        ("1m", "1000000B"),
        ("1g", "1000000000B"),
        ("1Kib", "0B"),
        ("1kIB", "0B"),
        ("2147483648", "2147483648B"),
        ("2147483648B", "2147483648B"),
        ("4294967295B", "4294967295B"),
        ("4294967295.4", "4294967295B"),
        ("4294967294.5", "4294967295B"),
        ("4294967295 B", "4294967295B"),
    ];

    /// Inputs that parse to a value of 4 GiB or more: makeblastdb 2.17.0 stops with
    /// "BLAST options error: max_file_sz must be < 4 GiB" before the log line, so the parse did not fail.
    const ORACLE_AT_LEAST_4_GIB: &[&str] = &[
        "4294967296",
        "4294967297",
        "5GiB",
        "4GiB",
        "4096MiB",
        "4194304KiB",
        "4194304 KiB",
        "1T",
        "1Tb",
        "1TB",
        "1TiB",
        "1P",
        "1PB",
        "1PiB",
        "1E",
        "1EB",
        "1EiB",
        "1e",
        "18446744073709551615",
        "18446744073709551614",
        "18446744073709551615.4",
        "18446744073709551615.0",
        "1844674407370955161",
        "1844674407370955162",
        "18446744073709551K",
        "18446744073709551.615K",
        "18446744073709551.615KB",
        "15EiB",
        "16E",
        "18E",
        "00000000000000000000018446744073709551615",
        "1t",
        "1p",
        "4294967295.5",
    ];

    /// Inputs for which makeblastdb 2.17.0 ended with the `CStringException` of `StringToUInt8_DataSize`:
    /// (input, text after "StringToUInt8_DataSize() - " and before " (m_Pos = ", the m_Pos).
    const ORACLE_ERRORS: &[(&str, &str, usize)] = &[
        ("-5", "Cannot convert string '-5' to Uint8", 0),
        ("abc", "Cannot convert string 'abc' to Uint8", 0),
        ("", "Cannot convert string '' to Uint8", 0),
        ("+", "Cannot convert string '+' to Uint8", 1),
        ("-", "Cannot convert string '-' to Uint8", 0),
        ("+-5", "Cannot convert string '+-5' to Uint8", 1),
        ("++5", "Cannot convert string '++5' to Uint8", 1),
        (".5", "Cannot convert string '.5' to Uint8", 0),
        ("5.", "Cannot convert string '5.' to Uint8", 1),
        (" 1", "Cannot convert string ' 1' to Uint8", 0),
        (" 1KB", "Cannot convert string ' 1KB' to Uint8", 0),
        ("1 ", "Cannot convert string '1 ' to Uint8", 1),
        ("1KB ", "Cannot convert string '1KB ' to Uint8", 3),
        ("1  KB", "Cannot convert string '1  KB' to Uint8", 1),
        ("1 K B", "Cannot convert string '1 K B' to Uint8", 3),
        ("1K B", "Cannot convert string '1K B' to Uint8", 2),
        ("1KBB", "Cannot convert string '1KBB' to Uint8", 3),
        ("1KiBB", "Cannot convert string '1KiBB' to Uint8", 4),
        ("1KiB ", "Cannot convert string '1KiB ' to Uint8", 4),
        ("1iB", "Cannot convert string '1iB' to Uint8", 1),
        ("1B ", "Cannot convert string '1B ' to Uint8", 2),
        ("1BK", "Cannot convert string '1BK' to Uint8", 2),
        ("1 BK", "Cannot convert string '1 BK' to Uint8", 3),
        ("5Gi", "Cannot convert string '5Gi' to Uint8", 2),
        ("1Z", "Cannot convert string '1Z' to Uint8", 1),
        ("1Y", "Cannot convert string '1Y' to Uint8", 1),
        (
            "18446744073709551616",
            "Cannot convert string '18446744073709551616' to Uint8",
            19,
        ),
        (
            "18446744073709551615.5",
            "Cannot convert string '18446744073709551615.5' to Uint8",
            20,
        ),
        (
            "18446744073709551616.5",
            "Cannot convert string '18446744073709551616.5' to Uint8",
            19,
        ),
        (
            "99999999999999999999",
            "Cannot convert string '99999999999999999999' to Uint8",
            19,
        ),
        (
            "184467440737095516150",
            "Cannot convert string '184467440737095516150' to Uint8",
            20,
        ),
        (
            "18446744073709552K",
            "Cannot convert string '18446744073709552K' to Uint8",
            19,
        ),
        (
            "18446744073709551.616K",
            "Cannot convert string '18446744073709551.616K' to Uint8",
            19,
        ),
        (
            "18446744073709551.6155K",
            "Cannot convert string '18446744073709551.6155K' to Uint8",
            20,
        ),
        (
            "18446744073709551615KiB",
            "Cannot convert string '18446744073709551615KiB' to Uint8",
            19,
        ),
        (
            "17179869184GiB",
            "Cannot convert string '17179869184GiB' to Uint8",
            19,
        ),
        ("16EiB", "Cannot convert string '16EiB' to Uint8", 19),
        (
            "18446744073709551616E",
            "Cannot convert string '18446744073709551616E' to Uint8",
            19,
        ),
        ("1.2.3", "Cannot convert string '1.2.3' to Uint8", 3),
        ("1..2", "Cannot convert string '1..2' to Uint8", 1),
        ("1.K", "Cannot convert string '1.K' to Uint8", 1),
        ("1.KiB", "Cannot convert string '1.KiB' to Uint8", 1),
        ("1. K", "Cannot convert string '1. K' to Uint8", 1),
        ("1 .5", "Cannot convert string '1 .5' to Uint8", 1),
        (".K", "Cannot convert string '.K' to Uint8", 0),
        (".", "Cannot convert string '.' to Uint8", 0),
        ("0x10", "Cannot convert string '0x10' to Uint8", 1),
        ("1e3", "Cannot convert string '1e3' to Uint8", 2),
        ("1,000", "Cannot convert string '1,000' to Uint8", 1),
        ("1_000", "Cannot convert string '1_000' to Uint8", 1),
        ("1KiBs", "Cannot convert string '1KiBs' to Uint8", 4),
        ("1\u{9}KB", "Cannot convert string '1\\tKB' to Uint8", 1),
        ("\u{9}1", "Cannot convert string '\\t1' to Uint8", 0),
        ("1\u{a}", "Cannot convert string '1\\n' to Uint8", 1),
        ("1\u{a}2", "Cannot convert string '1\\n2' to Uint8", 1),
        ("1?", "Cannot convert string '1?' to Uint8", 1),
        ("1??", "Cannot convert string '1\\?\\?' to Uint8", 1),
        ("1???", "Cannot convert string '1\\?\\?\\?' to Uint8", 1),
        ("?1", "Cannot convert string '?1' to Uint8", 0),
        ("??", "Cannot convert string '\\?\\?' to Uint8", 0),
        ("a&b", "Cannot convert string 'a&b' to Uint8", 0),
        ("1&", "Cannot convert string '1&' to Uint8", 1),
        ("\u{e9}1", "Cannot convert string '\\303\\2511' to Uint8", 0),
        ("1\u{e9}", "Cannot convert string '1\\303\\251' to Uint8", 1),
        ("1\"2", "Cannot convert string '1\\\"2' to Uint8", 1),
        ("1'2", "Cannot convert string '1\\'2' to Uint8", 1),
        ("1\\2", "Cannot convert string '1\\\\2' to Uint8", 1),
        ("1\u{7}", "Cannot convert string '1\\a' to Uint8", 1),
        ("1\u{1}", "Cannot convert string '1\\1' to Uint8", 1),
        ("1\u{7f}", "Cannot convert string '1\\177' to Uint8", 1),
        ("1\u{1}2", "Cannot convert string '1\\0012' to Uint8", 1),
        ("1KiB1", "Cannot convert string '1KiB1' to Uint8", 4),
        ("1 KiB ", "Cannot convert string '1 KiB ' to Uint8", 5),
        (
            "4294967295 ",
            "Cannot convert string '4294967295 ' to Uint8",
            10,
        ),
    ];

    fn log_line(input: &str) -> String {
        uint8_to_string_data_size(string_to_uint8_data_size(input).unwrap())
    }

    #[test]
    fn oracle_log_lines() {
        for (input, want) in ORACLE_LOG {
            assert_eq!(&log_line(input), want, "input {:?}", input);
        }
    }

    #[test]
    fn oracle_at_least_4_gib() {
        for input in ORACLE_AT_LEAST_4_GIB {
            let v = string_to_uint8_data_size(input).unwrap();
            assert!(v >= 1 << 32, "input {:?} gave {}", input, v);
            assert_eq!(uint8_to_string_data_size(v), format!("{}B", v));
        }
    }

    #[test]
    fn oracle_errors() {
        for (input, message, pos) in ORACLE_ERRORS {
            let e = string_to_uint8_data_size(input).unwrap_err();
            assert_eq!(&e.message, message, "input {:?}", input);
            assert_eq!(e.pos, *pos, "input {:?}", input);
        }
    }

    #[test]
    fn oracle_pair_count() {
        assert_eq!(
            ORACLE_LOG.len() + ORACLE_AT_LEAST_4_GIB.len() + ORACLE_ERRORS.len(),
            220
        );
    }

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
    fn data_size_formatting_thresholds() {
        // Rule 1: the letter index is always kMods.size(), so only "B" is appended.
        // Rule 2: values below 10 * 1024 are divided down to 0, larger ones are kept.
        assert_eq!(uint8_to_string_data_size(0), "0B");
        assert_eq!(uint8_to_string_data_size(1), "0B");
        assert_eq!(uint8_to_string_data_size(1023), "0B");
        assert_eq!(uint8_to_string_data_size(1024), "0B");
        assert_eq!(uint8_to_string_data_size(10239), "0B");
        assert_eq!(uint8_to_string_data_size(10240), "10240B");
        assert_eq!(uint8_to_string_data_size(10241), "10241B");
        assert_eq!(uint8_to_string_data_size(1 << 20), "1048576B");
        assert_eq!(uint8_to_string_data_size(4294967295), "4294967295B");
        assert_eq!(uint8_to_string_data_size(u64::MAX), "18446744073709551615B");
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
