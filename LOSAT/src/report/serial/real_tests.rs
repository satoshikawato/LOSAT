//! Unit tests of `real` (NCBI's serial real-number digits) and `x87` (the 80-bit `long double`
//! operations behind them).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2447-2654
//! ```c++
//! SIZE_TYPE NStr::DoubleToString_Ecvt(double val, unsigned int precision,
//!                                     char* buffer, SIZE_TYPE bufsize,
//!                                     int* dec, int* sign)
//! ```
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2658-2744
//! ```c++
//! SIZE_TYPE NStr::DoubleToStringPosix(double val, unsigned int precision,
//!                                     char* buffer, SIZE_TYPE bufsize)
//! ```
//!
//! The expected values are NCBI's own outputs. `SAMPLE_ROWS` holds every distinct real that
//! NCBI BLAST+ 2.17.0 printed in the outfmt 5, 8, 12, 13, 14, 15 and 16 samples of
//! outfmt-20261019 (`ncbi-samples/out`); the input is the correctly rounded `strtod` of the
//! printed text, and each row was kept only after NCBI's code printed that same text again
//! (all 293 did). `EDGE_ROWS` and `SWEEP_ROWS` come from the same NCBI code: the harness in
//! outfmt-20261019/S/harness extracts ncbistr.cpp:2422-2744 and the `WriteDouble2` functions
//! of serial/objostrasn.cpp, objostrjson.cpp and objostrxml.cpp verbatim and compiles them with
//! g++ -O2 on x86-64 Linux (80-bit x87 long double, control word 0x037F). Rows tagged
//! `not-correctly-rounded` are values whose 15 NCBI digits differ from the correctly rounded
//! decimal (e.g. 0x3fd613716aefcc27: NCBI 0.34493670886076, correct 0.344936708860759).

use super::real::{
    double_to_string_ecvt, double_to_string_posix, long_double_pow10_table, write_asn_text_double,
    write_json_double, write_xml_double, JsonDoubleError, DBL_DIG,
};
use crate::utils::x87::X87;
use std::cmp::Ordering;

// -------------------------------------------------------------------------------------------
// A tiny unsigned big integer, only for checking the long double literals.

#[derive(Clone, Debug)]
struct Big(Vec<u32>);

impl Big {
    fn from_u64(v: u64) -> Big {
        let mut b = Big(vec![v as u32, (v >> 32) as u32]);
        b.trim();
        b
    }

    fn trim(&mut self) {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
    }

    fn pow10(n: u32) -> Big {
        let mut b = Big::from_u64(1);
        for _ in 0..n {
            let mut carry = 0u64;
            for limb in &mut b.0 {
                let v = (*limb as u64) * 10 + carry;
                *limb = v as u32;
                carry = v >> 32;
            }
            if carry != 0 {
                b.0.push(carry as u32);
            }
        }
        b
    }

    fn bits(&self) -> u32 {
        match self.0.last() {
            None => 0,
            Some(&top) => (self.0.len() as u32 - 1) * 32 + (32 - top.leading_zeros()),
        }
    }

    fn bit(&self, i: u32) -> bool {
        let limb = (i / 32) as usize;
        limb < self.0.len() && (self.0[limb] >> (i % 32)) & 1 == 1
    }

    fn shl(&self, k: u32) -> Big {
        let words = (k / 32) as usize;
        let bits = k % 32;
        let mut out = vec![0u32; words];
        let mut carry = 0u32;
        for &limb in &self.0 {
            if bits == 0 {
                out.push(limb);
            } else {
                out.push((limb << bits) | carry);
                carry = limb >> (32 - bits);
            }
        }
        if carry != 0 {
            out.push(carry);
        }
        let mut b = Big(out);
        b.trim();
        b
    }

    fn cmp(&self, other: &Big) -> Ordering {
        if self.0.len() != other.0.len() {
            return self.0.len().cmp(&other.0.len());
        }
        for i in (0..self.0.len()).rev() {
            if self.0[i] != other.0[i] {
                return self.0[i].cmp(&other.0[i]);
            }
        }
        Ordering::Equal
    }

    fn sub_assign(&mut self, other: &Big) {
        assert!(self.cmp(other) != Ordering::Less);
        let mut borrow = 0i64;
        for i in 0..self.0.len() {
            let o = *other.0.get(i).unwrap_or(&0) as i64;
            let mut v = self.0[i] as i64 - o - borrow;
            borrow = 0;
            if v < 0 {
                v += 1 << 32;
                borrow = 1;
            }
            self.0[i] = v as u32;
        }
        assert_eq!(borrow, 0);
        self.trim();
    }

    /// Long division; the quotient must fit in 128 bits.
    fn divmod(&self, d: &Big) -> (u128, Big) {
        let mut q = 0u128;
        let mut r = Big(Vec::new());
        for i in (0..self.bits()).rev() {
            r = r.shl(1);
            if self.bit(i) {
                if r.0.is_empty() {
                    r.0.push(1);
                } else {
                    r.0[0] |= 1;
                }
            }
            assert!(q >> 127 == 0);
            q <<= 1;
            if r.cmp(d) != Ordering::Less {
                r.sub_assign(d);
                q |= 1;
            }
        }
        (q, r)
    }
}

/// `num / den` rounded to a 64-bit significand, ties to even: `(significand, exponent)` with
/// `num / den ~ significand * 2^(exponent - 63)`.
fn round_ratio_to_x87(num: &Big, den: &Big) -> (u64, i32) {
    let scaled = |e: i32| -> (Big, Big) {
        // num * 2^max(0,-e) and den * 2^max(0,e): compares num/den with 2^e.
        (num.shl((-e).max(0) as u32), den.shl(e.max(0) as u32))
    };
    let mut e = num.bits() as i32 - den.bits() as i32;
    loop {
        let (a, b) = scaled(e);
        if a.cmp(&b) == Ordering::Less {
            e -= 1;
        } else {
            break;
        }
    }
    loop {
        let (a, b) = scaled(e + 1);
        if a.cmp(&b) != Ordering::Less {
            e += 1;
        } else {
            break;
        }
    }
    let k = 63 - e;
    let (n2, d2) = if k >= 0 {
        (num.shl(k as u32), den.clone())
    } else {
        (num.clone(), den.shl((-k) as u32))
    };
    let (q, r) = n2.divmod(&d2);
    assert!(q >> 63 == 1 && q >> 64 == 0, "quotient not normalized");
    let twice_r = r.shl(1);
    let round_up = match twice_r.cmp(&d2) {
        Ordering::Greater => true,
        Ordering::Equal => q & 1 == 1,
        Ordering::Less => false,
    };
    let mut q = q + round_up as u128;
    if q >> 64 != 0 {
        q = 1 << 63;
        e += 1;
    }
    (q as u64, e)
}

/// Every `__NLG(1.e<n>)` literal of DoubleToString_Ecvt is the correctly rounded 80-bit value.
#[test]
fn long_double_literals_are_correctly_rounded() {
    let table = long_double_pow10_table();
    let mut expected_n: Vec<i32> = (0..=24).collect();
    expected_n.extend((-16..=-1).rev());
    expected_n.extend((2..=16).map(|k| 16 * k));
    expected_n.extend((2..=16).map(|k| -16 * k));
    expected_n.sort();
    let table_n: Vec<i32> = table.iter().map(|&(n, _, _)| n).collect();
    assert_eq!(table_n, expected_n);
    for &(n, significand, exponent) in table {
        let (num, den) = if n >= 0 {
            (Big::pow10(n as u32), Big::from_u64(1))
        } else {
            (Big::from_u64(1), Big::pow10((-n) as u32))
        };
        assert_eq!(
            round_ratio_to_x87(&num, &den),
            (significand, exponent),
            "1e{n}L"
        );
    }
    // 10^n for n <= 22 is exact in a double, and the f64 -> x87 conversion is exact.
    for &(n, significand, exponent) in table.iter().filter(|&&(n, _, _)| (0..=22).contains(&n)) {
        let x = X87::from_f64(format!("1e{n}").parse::<f64>().unwrap());
        assert_eq!(x.parts(), (false, exponent, significand), "1e{n}");
    }
}

// -------------------------------------------------------------------------------------------
// x87 operations

#[test]
fn x87_from_f64_is_exact() {
    assert_eq!(X87::from_f64(1.0).parts(), (false, 0, 1 << 63));
    assert_eq!(X87::from_f64(-1.5).parts(), (true, 0, 3 << 62));
    assert_eq!(
        X87::from_f64(f64::MAX).parts(),
        (false, 1023, 0xFFFF_FFFF_FFFF_F800)
    );
    assert_eq!(
        X87::from_f64(f64::MIN_POSITIVE).parts(),
        (false, -1022, 1 << 63)
    );
    // Subnormal doubles become normal 80-bit values.
    assert_eq!(
        X87::from_f64(f64::from_bits(1)).parts(),
        (false, -1074, 1 << 63)
    );
    assert_eq!(
        X87::from_f64(f64::from_bits(3)).parts(),
        (false, -1073, 3 << 62)
    );
    assert_eq!(
        X87::from_f64(f64::from_bits(0x000F_FFFF_FFFF_FFFF)).parts(),
        (false, -1023, 0xFFFF_FFFF_FFFF_F000)
    );
    assert_eq!(X87::from_f64(0.0).parts(), (false, 0, 0));
    assert_eq!(X87::from_f64(-0.0).parts(), (true, 0, 0));
}

#[test]
fn x87_from_u32_and_truncation() {
    assert_eq!(X87::from_u32(1).parts(), (false, 0, 1 << 63));
    assert_eq!(
        X87::from_u32(999_999_999).parts(),
        (false, 29, (999_999_999u64) << 34)
    );
    assert_eq!(X87::from_u32(0).parts(), (false, 0, 0));
    assert_eq!(X87::from_u32(u32::MAX).to_u32_trunc(), u32::MAX);
    assert_eq!(X87::from_f64(123_456_789.999).to_u32_trunc(), 123_456_789);
    assert_eq!(X87::from_f64(0.999).to_u32_trunc(), 0);
    assert_eq!(X87::ZERO.to_u32_trunc(), 0);
}

#[test]
fn x87_mul_rounds_to_nearest_even() {
    let one_and_ulps = |k: u64| X87::from_parts((1 << 63) + k, 0);
    let three_halves = X87::from_parts(3 << 62, 0);
    // 1.5 * (1 + 2^-63) = 1.5 + 1.5 ulp: a tie between odd +1 and even +2 -> +2.
    assert_eq!(
        three_halves.mul(one_and_ulps(1)).parts(),
        (false, 0, (3 << 62) + 2)
    );
    // 1.5 * (1 + 3 * 2^-63) = 1.5 + 4.5 ulp: a tie between even +4 and odd +5 -> +4.
    assert_eq!(
        three_halves.mul(one_and_ulps(3)).parts(),
        (false, 0, (3 << 62) + 4)
    );
    // (1 + 2^-63)^2 = 1 + 2 ulp + 2^-126: below half -> down.
    assert_eq!(
        one_and_ulps(1).mul(one_and_ulps(1)).parts(),
        (false, 0, (1 << 63) + 2)
    );
    // (2 - 2^-63)^2 = 4 - 2^-61 + 2^-126 -> 0xFFFF_FFFF_FFFF_FFFE * 2^-62.
    let max = X87::from_parts(u64::MAX, 0);
    assert_eq!(max.mul(max).parts(), (false, 1, u64::MAX - 1));
    // A product whose rounding carries into the next binade (2^64 * 2^-63 = 2^1).
    let a = X87::from_parts(0x8000_0000_0000_3039, 0);
    let b = X87::from_parts(0xFFFF_FFFF_FFFF_9F8E, 0);
    assert_eq!(a.mul(b).parts(), (false, 1, 1 << 63));
    // Signs and zero.
    assert_eq!(
        X87::from_f64(-2.0).mul(X87::from_f64(3.0)).parts(),
        (true, 2, 3 << 62)
    );
    assert_eq!(X87::from_f64(-2.0).mul(X87::ZERO).parts(), (true, 0, 0));
}

#[test]
fn x87_sub_u32_is_exact() {
    let t1 = X87::from_f64(123_456_789.75);
    assert_eq!(
        t1.sub_u32_exact(123_456_789).parts(),
        X87::from_f64(0.75).parts()
    );
    let t1 = X87::from_f64(999_999_999.5);
    assert_eq!(
        t1.sub_u32_exact(999_999_999).parts(),
        X87::from_f64(0.5).parts()
    );
    // t1 >= 1e9 with first = 999999999.
    let t1 = X87::from_u32(1_000_000_000);
    assert_eq!(
        t1.sub_u32_exact(999_999_999).parts(),
        X87::from_f64(1.0).parts()
    );
    assert_eq!(
        X87::from_u32(100_000_000)
            .sub_u32_exact(100_000_000)
            .parts(),
        (false, 0, 0)
    );
    // The lowest significand bit of a value in [2^26, 2^27) (2^-37) survives.
    let t1 = X87::from_parts((1 << 63) | 1, 26);
    assert_eq!(t1.sub_u32_exact(1 << 26).parts(), (false, -37, 1 << 63));
}

#[test]
fn x87_comparisons() {
    let a = X87::from_f64(1.5);
    let b = X87::from_f64(2.0);
    assert!(a.lt(b) && !b.lt(a) && b.ge(a) && a.le(b) && a.ge(a) && a.le(a) && !a.lt(a));
    assert!(X87::from_f64(-0.0).ge(X87::ZERO) && X87::from_f64(-0.0).le(X87::ZERO));
    assert!(X87::from_f64(-1.0).lt(X87::ZERO) && X87::ZERO.lt(X87::from_f64(f64::from_bits(1))));
    assert!(X87::from_f64(-2.0).lt(X87::from_f64(-1.0)));
}

// -------------------------------------------------------------------------------------------
// DoubleToString_Ecvt, DoubleToStringPosix and the object-stream forms

/// `(f64 bits, DoubleToString_Ecvt digits at DBL_DIG, *dec, *sign, DoubleToStringPosix at
/// DBL_DIG, CObjectOStreamAsn::WriteDouble)`.
type Row = (u64, &'static str, i32, i32, &'static str, &'static str);

fn check_rows(rows: &[Row]) {
    for &(bits, digits, dec, sign, posix, asn) in rows {
        let x = f64::from_bits(bits);
        let e = double_to_string_ecvt(x, DBL_DIG);
        assert_eq!(
            (std::str::from_utf8(e.digits()).unwrap(), e.dec(), e.sign()),
            (digits, dec, sign),
            "DoubleToString_Ecvt {bits:#018x}"
        );
        let mut out = Vec::new();
        double_to_string_posix(x, DBL_DIG, &mut out);
        assert_eq!(out, posix.as_bytes(), "DoubleToStringPosix {bits:#018x}");
        out.clear();
        write_asn_text_double(x, &mut out);
        assert_eq!(out, asn.as_bytes(), "ASN.1 text {bits:#018x}");
        out.clear();
        write_xml_double(x, &mut out);
        assert_eq!(out, posix.as_bytes(), "XML {bits:#018x}");
        out.clear();
        write_json_double(x, &mut out).unwrap();
        assert_eq!(out, posix.as_bytes(), "JSON {bits:#018x}");
    }
}

#[test]
fn ncbi_sample_reals() {
    assert_eq!(SAMPLE_ROWS.len(), 293);
    check_rows(SAMPLE_ROWS);
}

#[test]
fn ncbi_edge_reals() {
    check_rows(EDGE_ROWS);
}

/// `(f64 bits, precision, Ecvt digits, *dec, *sign, DoubleToStringPosix)` for precision
/// 0..=17 (Ecvt clamps to DBL_DIG; DoubleToStringPosix keeps the caller's precision for its
/// `e` format decision).
type SweepRow = (u64, u32, &'static str, i32, i32, &'static str);

#[test]
fn ncbi_precision_sweep() {
    for &(bits, precision, digits, dec, sign, posix) in SWEEP_ROWS {
        let x = f64::from_bits(bits);
        let e = double_to_string_ecvt(x, precision);
        assert_eq!(
            (std::str::from_utf8(e.digits()).unwrap(), e.dec(), e.sign()),
            (digits, dec, sign),
            "DoubleToString_Ecvt {bits:#018x} precision {precision}"
        );
        let mut out = Vec::new();
        double_to_string_posix(x, precision, &mut out);
        assert_eq!(
            out,
            posix.as_bytes(),
            "DoubleToStringPosix {bits:#018x} precision {precision}"
        );
    }
}

/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:185-196
/// ```c++
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
/// ```
#[test]
fn non_finite_reals() {
    let cases: [(f64, &str, &str, Option<JsonDoubleError>); 4] = [
        (
            f64::NAN,
            "NOT-A-NUMBER",
            "NaN",
            Some(JsonDoubleError::NotANumber),
        ),
        (
            -f64::NAN,
            "NOT-A-NUMBER",
            "NaN",
            Some(JsonDoubleError::NotANumber),
        ),
        (
            f64::INFINITY,
            "PLUS-INFINITY",
            "INF",
            Some(JsonDoubleError::Infinite),
        ),
        (
            f64::NEG_INFINITY,
            "MINUS-INFINITY",
            "-INF",
            Some(JsonDoubleError::Infinite),
        ),
    ];
    for (x, asn, xml, json) in cases {
        let mut out = Vec::new();
        write_asn_text_double(x, &mut out);
        assert_eq!(out, asn.as_bytes());
        out.clear();
        write_xml_double(x, &mut out);
        assert_eq!(out, xml.as_bytes());
        out.clear();
        assert_eq!(write_json_double(x, &mut out).err(), json);
        assert!(out.is_empty());
    }
    assert_eq!(
        JsonDoubleError::NotANumber.message(),
        "invalid double: not a number"
    );
    assert_eq!(
        JsonDoubleError::Infinite.message(),
        "invalid double: infinite"
    );
}

// -------------------------------------------------------------------------------------------
// Tables (generated by outfmt-20261019/S/harness/assemble_tests.py from NCBI's outputs)

#[rustfmt::skip]
const SAMPLE_ROWS: &[Row] = &[
    (0x0000000000000000, "0", 0, 0, "0", "{ 0, 10, 0 }"), // json blastn-blastn_db_bf_fmt12 (+2)
    (0x408382cd42bfa376, "624350224969078", 2, 1, "624.350224969078", "{ 624350224969078, 10, -12 }"), // json blastn-blastn_db_bf_fmt12 (+1)
    (0x4056280000000000, "88625", 1, 1, "88.625", "{ 88625, 10, -3 }"), // json blastn-blastn_db_bf_fmt12 (+1)
    (0x22975ddd6754cd3d, "479047620391651", -142, 1, "4.79047620391651e-142", "{ 479047620391651, 10, -156 }"), // json blastn-blastn_db_bf_fmt12 (+1)
    (0x407e832253935c80, "488195880485182", 2, 1, "488.195880485182", "{ 488195880485182, 10, -12 }"), // json blastn-blastn_db_bf_fmt12 (+1)
    (0x4055bd5555555553, "869583333333333", 1, 1, "86.9583333333333", "{ 869583333333333, 10, -13 }"), // json blastn-blastn_db_bf_fmt12 (+1)
    (0x1b10eab42ac922cb, "260916226603389", -178, 1, "2.60916226603389e-178", "{ 260916226603389, 10, -192 }"), // json blastn-megablast_db_bf_fmt12 (+1)
    (0x40830558200916e7, "608668029852878", 2, 1, "608.668029852878", "{ 608668029852878, 10, -12 }"), // json blastn-megablast_db_bf_fmt12 (+1)
    (0x251ab3da375dcf16, "601914206611342", -130, 1, "6.01914206611342e-130", "{ 601914206611342, 10, -144 }"), // json blastn-megablast_db_bf_fmt12 (+1)
    (0x407c0026f40fe42e, "448009510099483", 2, 1, "448.009510099483", "{ 448009510099483, 10, -12 }"), // json blastn-megablast_db_bf_fmt12 (+1)
    (0x2c91b3132a20ac13, "530325651018409", -94, 1, "5.30325651018409e-94", "{ 530325651018409, 10, -108 }"), // json blastp_db_bf_fmt12 (+1)
    (0x406ff84473563d7b, "255758355778158", 2, 1, "255.758355778158", "{ 255758355778158, 10, -12 }"), // json blastp_db_bf_fmt12 (+1)
    (0x4059000000000000, "1", 2, 1, "100", "{ 1, 10, 2 }"), // json blastp_db_bf_fmt12 (+1)
    (0x321e6ed8a5474331, "28220561436493", -67, 1, "2.8220561436493e-67", "{ 28220561436493, 10, -80 }"), // json blastp_db_bf_fmt12 (+1)
    (0x40678b2656a7e839, "188348429992622", 2, 1, "188.348429992622", "{ 188348429992622, 10, -12 }"), // json blastp_db_bf_fmt12 (+1)
    (0x4058d24924924926, "992857142857143", 1, 1, "99.2857142857143", "{ 992857142857143, 10, -13 }"), // json blastp_db_bf_fmt12 (+1)
    (0x3ff83e70790499e9, "151524398110052", 0, 1, "1.51524398110052", "{ 151524398110052, 10, -14 }"), // json blastp_db_bf_fmt12 (+1)
    (0x40296514ad303db8, "126974233743084", 1, 1, "12.6974233743084", "{ 126974233743084, 10, -13 }"), // json blastp_db_bf_fmt12 (+1)
    (0x402baaaaaaaaaa98, "138333333333333", 1, 1, "13.8333333333333", "{ 138333333333333, 10, -13 }"), // json blastp_db_bf_fmt12 (+1)
    (0x3ffe29db6c812812, "188521902451112", 0, 1, "1.88521902451112", "{ 188521902451112, 10, -14 }"), // json blastp_db_bf_fmt12 (+1)
    (0x40289fdbcc363ecf, "123122237983911", 1, 1, "12.3122237983911", "{ 123122237983911, 10, -13 }"), // json blastp_db_bf_fmt12 (+1)
    (0x4034800000000000, "205", 1, 1, "20.5", "{ 205, 10, -1 }"), // json blastp_db_bf_fmt12 (+1)
    (0x2c8161df70066518, "260410803546674", -94, 1, "2.60410803546674e-94", "{ 260410803546674, 10, -108 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x320efa7798a44bc1, "14363153132552", -67, 1, "1.4363153132552e-67", "{ 14363153132552, 10, -80 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x3fe7cf367c5ac51d, "74404453552248", -1, 1, "0.74404453552248", "{ 74404453552248, 10, -14 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x3fed9f78fb88ef5b, "92571686866675", -1, 1, "0.92571686866675", "{ 92571686866675, 10, -14 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x40185b9590865c8d, "608943773098066", 0, 1, "6.08943773098066", "{ 608943773098066, 10, -14 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x40240086865a4477, "100010263428869", 1, 1, "10.0010263428869", "{ 100010263428869, 10, -13 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x4022c71c71c71c72, "938888888888889", 0, 1, "9.38888888888889", "{ 938888888888889, 10, -14 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x402239050a0f4eec, "9111366571778", 0, 1, "9.111366571778", "{ 9111366571778, 10, -12 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x40233b4da560458d, "96158267669696", 0, 1, "9.6158267669696", "{ 96158267669696, 10, -13 }"), // asn blastp_subject_fmt8_bl2seqlegacy
    (0x2bd0dc09f7b30ad6, "123329452098788", -97, 1, "1.23329452098788e-97", "{ 123329452098788, 10, -111 }"), // json blastx_db_bf_fmt12 (+1)
    (0x40709c6070763e02, "26577354475201", 2, 1, "265.77354475201", "{ 26577354475201, 10, -11 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4056000000000000, "88", 1, 1, "88", "{ 88, 10, 0 }"), // json blastx_db_bf_fmt12 (+1)
    (0x3fd56a47d320bc77, "334611850909773", -1, 1, "0.334611850909773", "{ 334611850909773, 10, -15 }"), // json blastx_db_bf_fmt12 (+1)
    (0x40308ca6bb7a1a13, "16549419133482", 1, 1, "16.549419133482", "{ 16549419133482, 10, -12 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4020000000000000, "8", 0, 1, "8", "{ 8, 10, 0 }"), // json blastx_db_bf_fmt12 (+1)
    (0x40137eb2896ce7f2, "487372793891153", 0, 1, "4.87372793891153", "{ 487372793891153, 10, -14 }"), // json blastx_db_bf_fmt12 (+1)
    (0x402a2a4d8e2a3cd9, "130826229502258", 1, 1, "13.0826229502258", "{ 130826229502258, 10, -13 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4023c00000000000, "9875", 0, 1, "9.875", "{ 9875, 10, -3 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4018c6c3e97565e6, "619410671977082", 0, 1, "6.19410671977082", "{ 619410671977082, 10, -14 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4022800000000000, "925", 0, 1, "9.25", "{ 925, 10, -2 }"), // json blastx_db_bf_fmt12 (+1)
    (0x31ef51534bcdbdaf, "363011628351791", -68, 1, "3.63011628351791e-68", "{ 363011628351791, 10, -82 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4067edc2c724e7a0, "19143002659996", 2, 1, "191.43002659996", "{ 19143002659996, 10, -11 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4055d80000000000, "87375", 1, 1, "87.375", "{ 87375, 10, -3 }"), // json blastx_db_bf_fmt12 (+1)
    (0x3faa2d3490c7003d, "511261393289719", -2, 1, "0.0511261393289719", "{ 511261393289719, 10, -16 }"), // json blastx_db_bf_fmt12 (+1)
    (0x4032dc515e681723, "188606165889861", 1, 1, "18.8606165889861", "{ 188606165889861, 10, -13 }"), // json blastx_db_bf_fmt12 (+1)
    (0x402dc00000000000, "14875", 1, 1, "14.875", "{ 14875, 10, -3 }"), // json blastx_db_bf_fmt12 (+1)
    (0x40109356d283317d, "414388588834515", 0, 1, "4.14388588834515", "{ 414388588834515, 10, -14 }"), // json blastx_db_bf_fmt12 (+1)
    (0x402f000000000000, "155", 1, 1, "15.5", "{ 155, 10, -1 }"), // json blastx_db_bf_fmt12 (+1)
    (0x3fe567809bf39199, "668884567825944", -1, 1, "0.668884567825944", "{ 668884567825944, 10, -15 }"), // json blastx_db_bf_fmt12 (+1)
    (0x40584aaaaaaaaaad, "971666666666667", 1, 1, "97.1666666666667", "{ 971666666666667, 10, -13 }"), // json blastx_db_bf_fmt12 (+1)
    (0x32890f7262b9601d, "297452685001389", -65, 1, "2.97452685001389e-65", "{ 297452685001389, 10, -79 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4067038f3bfc08e1, "184111234657531", 2, 1, "184.111234657531", "{ 184111234657531, 10, -12 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4058c4924924924b, "990714285714286", 1, 1, "99.0714285714286", "{ 990714285714286, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x3fd4de20b6c1d975, "32605760428796", -1, 1, "0.32605760428796", "{ 32605760428796, 10, -14 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x403d12492492492d, "290714285714286", 1, 1, "29.0714285714286", "{ 290714285714286, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x3ff06f63a44fe634, "102719463524262", 0, 1, "1.02719463524262", "{ 102719463524262, 10, -14 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x402e0469f30c37d8, "150086208298125", 1, 1, "15.0086208298125", "{ 150086208298125, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4027db6db6db6da7, "119285714285714", 1, 1, "11.9285714285714", "{ 119285714285714, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x35ebccb54ece49d0, "594415206848111", -49, 1, "5.94415206848111e-49", "{ 594415206848111, 10, -63 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4060148f5332119c, "128642495725432", 2, 1, "128.642495725432", "{ 128642495725432, 10, -12 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4047d2492492491f, "476428571428571", 1, 1, "47.6428571428571", "{ 476428571428571, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x40510aa26739f63a, "681661623064073", 1, 1, "68.1661623064073", "{ 681661623064073, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x403e800000000000, "305", 1, 1, "30.5", "{ 305, 10, -1 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x3d704a62826f9e80, "926011445922103", -13, 1, "9.26011445922103e-13", "{ 926011445922103, 10, -27 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4048d69a42bbf818, "496765826623743", 1, 1, "49.6765826623743", "{ 496765826623743, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4048892492492496, "490714285714286", 1, 1, "49.0714285714286", "{ 490714285714286, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x40019af70018e56f, "220066642833688", 0, 1, "2.20066642833688", "{ 220066642833688, 10, -14 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x403b2aaaaaaaaab4, "271666666666667", 1, 1, "27.1666666666667", "{ 271666666666667, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x40068d5f31df36cb, "281902922594495", 0, 1, "2.81902922594495", "{ 281902922594495, 10, -14 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4027daa2eb3c3fad, "119270242224737", 1, 1, "11.9270242224737", "{ 119270242224737, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x40300e38e38e38f0, "160555555555556", 1, 1, "16.0555555555556", "{ 160555555555556, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x401ad3dde424cc49, "670690113520238", 0, 1, "6.70690113520238", "{ 670690113520238, 10, -14 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x40265031294841a3, "11156625070639", 1, 1, "11.156625070639", "{ 11156625070639, 10, -12 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x4049ce38e38e38e2, "516111111111111", 1, 1, "51.6111111111111", "{ 516111111111111, 10, -13 }"), // json tblastn_db_bf_fmt12 (+1)
    (0x301b49403b0fca51, "589119181469459", -77, 1, "5.89119181469459e-77", "{ 589119181469459, 10, -91 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x406334d32a4e6001, "153650777009083", 2, 1, "153.650777009083", "{ 153650777009083, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402c800000000000, "1425", 1, 1, "14.25", "{ 1425, 10, -2 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40588d486538e1f6, "982075436645026", 1, 1, "98.2075436645026", "{ 982075436645026, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4022eaaaaaaaaaa9, "945833333333333", 0, 1, "9.45833333333333", "{ 945833333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404cb6a7660d65c2, "574269836011331", 1, 1, "57.4269836011331", "{ 574269836011331, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3079dac60867320f, "357255275884826", -75, 1, "3.57255275884826e-75", "{ 357255275884826, 10, -89 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x406410c3d4a60511, "160523905109651", 2, 1, "160.523905109651", "{ 160523905109651, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4055003271f0c04a, "840030789233289", 1, 1, "84.0030789233289", "{ 840030789233289, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404d669b2186b69b, "588016092212467", 1, 1, "58.8016092212467", "{ 588016092212467, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401e555555555552, "758333333333333", 0, 1, "7.58333333333333", "{ 758333333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3216c0b8765cca38, "210987200318381", -67, 1, "2.10987200318381e-67", "{ 210987200318381, 10, -81 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4063177fe064d2a2, "152734359929008", 2, 1, "152.734359929008", "{ 152734359929008, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4050e07c0d18db47, "67507571481966", 1, 1, "67.507571481966", "{ 67507571481966, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404c7c00d23a4ad8, "569687750610953", 1, 1, "56.9687750610953", "{ 569687750610953, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401d800000000000, "7375", 0, 1, "7.375", "{ 7375, 10, -3 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x321e7e06475b3f7a, "282755411966957", -67, 1, "2.82755411966957e-67", "{ 282755411966957, 10, -81 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40621e3bec23a00f, "144944814748364", 2, 1, "144.944814748364", "{ 144944814748364, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402ceaaaaaaaaa98, "144583333333333", 1, 1, "14.4583333333333", "{ 144583333333333, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405592d2e38083a7, "862941216235182", 1, 1, "86.2941216235182", "{ 862941216235182, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4023555555555557, "966666666666667", 0, 1, "9.66666666666667", "{ 966666666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4046c1bc629ca932, "455135615601488", 1, 1, "45.5135615601488", "{ 455135615601488, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4016d55555555552, "570833333333333", 0, 1, "5.70833333333333", "{ 570833333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x32f0eddf5382d0e0, "257203616700464", -63, 1, "2.57203616700464e-63", "{ 257203616700464, 10, -77 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x406308d63b700be0, "15227615138897", 2, 1, "152.27615138897", "{ 15227615138897, 10, -11 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404cf14df9e080b9, "57885192141171", 1, 1, "57.885192141171", "{ 57885192141171, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4024955555555568, "102916666666667", 1, 1, "10.2916666666667", "{ 102916666666667, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404aa6cc33a17344, "533031067407924", 1, 1, "53.3031067407924", "{ 533031067407924, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x34f6025ca31cae8e, "143616900956157", -53, 1, "1.43616900956157e-53", "{ 143616900956157, 10, -67 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4053a04afafe1e9e, "785045764428746", 1, 1, "78.5045764428746", "{ 785045764428746, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4052f0573f84cdcc, "757553252026475", 1, 1, "75.7553252026475", "{ 757553252026475, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401f2aaaaaaaaaae, "779166666666667", 0, 1, "7.79166666666667", "{ 779166666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404e8bdc04a63d56, "61092651921436", 1, 1, "61.092651921436", "{ 61092651921436, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40401cdda3b29bc9, "322255138990509", 1, 1, "32.2255138990509", "{ 322255138990509, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4008000000000000, "3", 0, 1, "3", "{ 3, 10, 0 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3f7b6026418ca169, "668349212427475", -3, 1, "0.00668349212427475", "{ 668349212427475, 10, -17 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40378501560f37aa, "235195516383316", 1, 1, "23.5195516383316", "{ 235195516383316, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3ff15804f0875863, "108398908573647", 0, 1, "1.08398908573647", "{ 108398908573647, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4030302edbabd95a, "161882149977258", 1, 1, "16.1882149977258", "{ 161882149977258, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4009aaaaaaaaaaa3, "320833333333333", 0, 1, "3.20833333333333", "{ 320833333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3ff1b7e18ccc6624, "110739283561123", 0, 1, "1.10739283561123", "{ 110739283561123, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40102aaaaaaaaaae, "404166666666667", 0, 1, "4.04166666666667", "{ 404166666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402da08ec9726f87, "148135893776123", 1, 1, "14.8135893776123", "{ 148135893776123, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4003000000000000, "2375", 0, 1, "2.375", "{ 2375, 10, -3 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3ff7d3da612097b3, "148922193470342", 0, 1, "1.48922193470342", "{ 148922193470342, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402f75c3680b470d, "15730006457688", 1, 1, "15.730006457688", "{ 15730006457688, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40153878a6cf5d98, "53051477493349", 0, 1, "5.3051477493349", "{ 53051477493349, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402bcb5a2ad997c8, "138971722975365", 1, 1, "13.8971722975365", "{ 138971722975365, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4011d55555555552, "445833333333333", 0, 1, "4.45833333333333", "{ 445833333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x33046fc430a4b03b, "620982399626209", -63, 1, "6.20982399626209e-63", "{ 620982399626209, 10, -77 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405cacfeca10c6d6, "114703051105865", 2, 1, "114.703051105865", "{ 114703051105865, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4027eaaaaaaaaa98, "119583333333333", 1, 1, "11.9583333333333", "{ 119583333333333, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4057faa7f3a91e99, "959165009643133", 1, 1, "95.9165009643133", "{ 959165009643133, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4025000000000000, "105", 1, 1, "10.5", "{ 105, 10, -1 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4049bc31e4550781, "51470272580641", 1, 1, "51.470272580641", "{ 51470272580641, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401a2aaaaaaaaaae, "654166666666667", 0, 1, "6.54166666666667", "{ 654166666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3369b0ec7f78fba3, "499613005537967", -61, 1, "4.99613005537967e-61", "{ 499613005537967, 10, -75 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405ad7ca2b77ef7a, "10737171446526", 2, 1, "107.37171446526", "{ 10737171446526, 10, -11 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4057680782195b3c, "93625458264124", 1, 1, "93.625458264124", "{ 93625458264124, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4027155555555568, "115416666666667", 1, 1, "11.5416666666667", "{ 115416666666667, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404b56bfef1ac41d, "54677732360906", 1, 1, "54.677732360906", "{ 54677732360906, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3529b1a3d8a290b3, "13412844949906", -52, 1, "1.3412844949906e-52", "{ 13412844949906, 10, -65 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405a4529b9e82beb, "10508067176507", 2, 1, "105.08067176507", "{ 10508067176507, 10, -11 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4027800000000000, "1175", 1, 1, "11.75", "{ 1175, 10, -2 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404737098a42df13, "464299786402245", 1, 1, "46.4299786402245", "{ 464299786402245, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4017aaaaaaaaaaae, "591666666666667", 0, 1, "5.91666666666667", "{ 591666666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3563962b5d15c021, "163594909590429", -51, 1, "1.63594909590429e-51", "{ 163594909590429, 10, -65 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4056f2ba5a732553, "917926241039725", 1, 1, "91.7926241039725", "{ 917926241039725, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4054c58bde1da559, "830866618432532", 1, 1, "83.0866618432532", "{ 830866618432532, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4025d55555555568, "109166666666667", 1, 1, "10.9166666666667", "{ 109166666666667, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x40485c4a6d6265ce, "487210213404138", 1, 1, "48.7210213404138", "{ 487210213404138, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401bd55555555552, "695833333333333", 0, 1, "6.95833333333333", "{ 695833333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3771dbe9d92aa299, "128132306489295", -41, 1, "1.28132306489295e-41", "{ 128132306489295, 10, -55 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405432eb6c8de1fc, "807956191430639", 1, 1, "80.7956191430639", "{ 807956191430639, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404f01292c4c7337, "620090690015117", 1, 1, "62.0090690015117", "{ 620090690015117, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4047ac56b1e914f5, "473463957203002", 1, 1, "47.3463957203002", "{ 473463957203002, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x37d4b357682296d5, "950522753486266", -40, 1, "9.50522753486266e-40", "{ 950522753486266, 10, -54 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x405a0a832615112c, "104164254684995", 2, 1, "104.164254684995", "{ 104164254684995, 10, -12 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x404a6c259fce584c, "528448982007545", 1, 1, "52.8448982007545", "{ 528448982007545, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x401b000000000000, "675", 0, 1, "6.75", "{ 675, 10, -2 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4004aaaaaaaaaaa3, "258333333333333", 0, 1, "2.58333333333333", "{ 258333333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4033da9818dd8882, "198538833180287", 1, 1, "19.8538833180287", "{ 198538833180287, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x400d000000000000, "3625", 0, 1, "3.625", "{ 3625, 10, -3 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x3ffb3bdd15787a3c, "17021151388832", 0, 1, "1.7021151388832", "{ 17021151388832, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x400155555555555d, "216666666666667", 0, 1, "2.16666666666667", "{ 216666666666667, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402cb5f47a2603a7, "143553808375744", 1, 1, "14.3553808375744", "{ 143553808375744, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4030a57c03520f49, "166464235377637", 1, 1, "16.6464235377637", "{ 166464235377637, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x402ae0bfdb8d2c21, "134389637574987", 1, 1, "13.4389637574987", "{ 134389637574987, 10, -13 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x4010194be618c7cd, "402470359350441", 0, 1, "4.02470359350441", "{ 402470359350441, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x400eaaaaaaaaaaa3, "383333333333333", 0, 1, "3.83333333333333", "{ 383333333333333, 10, -14 }"), // json tblastx_db_bf_fmt12 (+1)
    (0x408382cccccccccd, "62435", 2, 1, "624.35", "{ 62435, 10, -2 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x407e8322d0e56042, "488196", 2, 1, "488.196", "{ 488196, 10, -3 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x22975dde9dfb9c91, "479048", -142, 1, "4.79048e-142", "{ 479048, 10, -147 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x3fda3d70a3d70a3d, "41", -1, 1, "0.41", "{ 41, 10, -2 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x3fe4000000000000, "625", -1, 1, "0.625", "{ 625, 10, -3 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x3fe8f5c28f5c28f6, "78", -1, 1, "0.78", "{ 78, 10, -2 }"), // json blastn-blastn_db_bf_fmt15 (+1)
    (0x4083055810624dd3, "608668", 2, 1, "608.668", "{ 608668, 10, -3 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x1b10eab3344b0de3, "260916", -178, 1, "2.60916e-178", "{ 260916, 10, -183 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x407c0028f5c28f5c, "44801", 2, 1, "448.01", "{ 44801, 10, -2 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x251ab3d99d96d351, "601914", -130, 1, "6.01914e-130", "{ 601914, 10, -135 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x3fdd70a3d70a3d71, "46", -1, 1, "0.46", "{ 46, 10, -2 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x3ff47ae147ae147b, "128", 0, 1, "1.28", "{ 128, 10, -2 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x3feb333333333333, "85", -1, 1, "0.85", "{ 85, 10, -2 }"), // json blastn-megablast_db_bf_fmt15 (+1)
    (0x406ff84189374bc7, "255758", 2, 1, "255.758", "{ 255758, 10, -3 }"), // json blastp_db_bf_fmt15 (+1)
    (0x2c91b313ed88f8d0, "530326", -94, 1, "5.30326e-94", "{ 530326, 10, -99 }"), // json blastp_db_bf_fmt15 (+1)
    (0x40678b22d0e56042, "188348", 2, 1, "188.348", "{ 188348, 10, -3 }"), // json blastp_db_bf_fmt15 (+1)
    (0x321e6edb5efd3966, "282206", -67, 1, "2.82206e-67", "{ 282206, 10, -72 }"), // json blastp_db_bf_fmt15 (+1)
    (0x3fa4fdf3b645a1cb, "41", -2, 1, "0.041", "{ 41, 10, -3 }"), // json blastp_db_bf_fmt15 (+1)
    (0x3fd116872b020c4a, "267", -1, 1, "0.267", "{ 267, 10, -3 }"), // json blastp_db_bf_fmt15 (+1)
    (0x3fc1eb851eb851ec, "14", -1, 1, "0.14", "{ 14, 10, -2 }"), // json blastp_db_bf_fmt15 (+1)
    (0x402965119ce075f7, "126974", 1, 1, "12.6974", "{ 126974, 10, -4 }"), // json blastp_db_bf_fmt15 (+1)
    (0x3ff83e6c4c5974e6, "151524", 0, 1, "1.51524", "{ 151524, 10, -5 }"), // json blastp_db_bf_fmt15 (+1)
    (0x40289fd8adab9f56, "123122", 1, 1, "12.3122", "{ 123122, 10, -4 }"), // json blastp_db_bf_fmt15 (+1)
    (0x3ffe29dc725c3dee, "188522", 0, 1, "1.88522", "{ 188522, 10, -5 }"), // json blastp_db_bf_fmt15 (+1)
    (0x40709c624dd2f1aa, "265774", 2, 1, "265.774", "{ 265774, 10, -3 }"), // json blastx_db_bf_fmt15 (+1)
    (0x2bd0dc05eacfba80, "123329", -97, 1, "1.23329e-97", "{ 123329, 10, -102 }"), // json blastx_db_bf_fmt15 (+1)
    (0x40308ca57a786c22, "165494", 1, 1, "16.5494", "{ 165494, 10, -4 }"), // json blastx_db_bf_fmt15 (+1)
    (0x3fd56a4873365882, "334612", -1, 1, "0.334612", "{ 334612, 10, -6 }"), // json blastx_db_bf_fmt15 (+1)
    (0x402a2a4a8c154c98, "130826", 1, 1, "13.0826", "{ 130826, 10, -4 }"), // json blastx_db_bf_fmt15 (+1)
    (0x40137eb313be22e6, "487373", 0, 1, "4.87373", "{ 487373, 10, -5 }"), // json blastx_db_bf_fmt15 (+1)
    (0x4018c6c4c5974e66, "619411", 0, 1, "6.19411", "{ 619411, 10, -5 }"), // json blastx_db_bf_fmt15 (+1)
    (0x4067edc28f5c28f6, "19143", 2, 1, "191.43", "{ 19143, 10, -2 }"), // json blastx_db_bf_fmt15 (+1)
    (0x31ef515565ba4b10, "363012", -68, 1, "3.63012e-68", "{ 363012, 10, -73 }"), // json blastx_db_bf_fmt15 (+1)
    (0x4032dc504816f007, "188606", 1, 1, "18.8606", "{ 188606, 10, -4 }"), // json blastx_db_bf_fmt15 (+1)
    (0x3faa2d333ef1ad57, "511261", -2, 1, "0.0511261", "{ 511261, 10, -7 }"), // json blastx_db_bf_fmt15 (+1)
    (0x40109357e670e2c1, "414389", 0, 1, "4.14389", "{ 414389, 10, -5 }"), // json blastx_db_bf_fmt15 (+1)
    (0x3fe5678183f91e64, "668885", -1, 1, "0.668885", "{ 668885, 10, -6 }"), // json blastx_db_bf_fmt15 (+1)
    (0x4067038d4fdf3b64, "184111", 2, 1, "184.111", "{ 184111, 10, -3 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x32890f741ff7b5ee, "297453", -65, 1, "2.97453e-65", "{ 297453, 10, -70 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x3fd4de225fa658c5, "326058", -1, 1, "0.326058", "{ 326058, 10, -6 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x402e0467381d7dbf, "150086", 1, 1, "15.0086", "{ 150086, 10, -4 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x3ff06f5ec80c73ac, "102719", 0, 1, "1.02719", "{ 102719, 10, -5 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x4060148b43958106, "128642", 2, 1, "128.642", "{ 128642, 10, -3 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x35ebccb4ac8146a6, "594415", -49, 1, "5.94415e-49", "{ 594415, 10, -54 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x40510aa305532618, "681662", 1, 1, "68.1662", "{ 681662, 10, -4 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x4048d69ad42c3c9f, "496766", 1, 1, "49.6766", "{ 496766, 10, -4 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x3d704a61fed2b299, "926011", -13, 1, "9.26011e-13", "{ 926011, 10, -18 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x40019af8df7a4e7b, "220067", 0, 1, "2.20067", "{ 220067, 10, -5 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x4027da9fbe76c8b4, "11927", 1, 1, "11.927", "{ 11927, 10, -3 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x40068d5f99c38b05, "281903", 0, 1, "2.81903", "{ 281903, 10, -5 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x4026502de00d1b71, "111566", 1, 1, "11.1566", "{ 111566, 10, -4 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x401ad3dd97f62b6b, "67069", 0, 1, "6.7069", "{ 67069, 10, -4 }"), // json tblastn_db_bf_fmt15 (+1)
    (0x406334d4fdf3b646, "153651", 2, 1, "153.651", "{ 153651, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x301b493fae0c32e7, "589119", -77, 1, "5.89119e-77", "{ 589119, 10, -82 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40588d47ae147ae1, "982075", 1, 1, "98.2075", "{ 982075, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404cb6a7ef9db22d, "57427", 1, 1, "57.427", "{ 57427, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x406410c49ba5e354, "160524", 2, 1, "160.524", "{ 160524, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3079dac4b96ebc88, "357255", -75, 1, "3.57255e-75", "{ 357255, 10, -80 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40550032ca57a787, "840031", 1, 1, "84.0031", "{ 840031, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404d669ad42c3c9f, "588016", 1, 1, "58.8016", "{ 588016, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4063177ced916873, "152734", 2, 1, "152.734", "{ 152734, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3216c0b70befa78f, "210987", -67, 1, "2.10987e-67", "{ 210987, 10, -72 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4050e07c84b5dcc6, "675076", 1, 1, "67.5076", "{ 675076, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404c7c01a36e2eb2, "569688", 1, 1, "56.9688", "{ 569688, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40621e3d70a3d70a, "144945", 2, 1, "144.945", "{ 144945, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x321e7e035e01307a, "282755", -67, 1, "2.82755e-67", "{ 282755, 10, -72 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405592d288ce703b, "862941", 1, 1, "86.2941", "{ 862941, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4046c1bda5119ce0, "455136", 1, 1, "45.5136", "{ 455136, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x406308d4fdf3b646, "152276", 2, 1, "152.276", "{ 152276, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x32f0ede0fac7e1ac, "257204", -63, 1, "2.57204e-63", "{ 257204, 10, -68 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404cf14e3bcd35a8, "578852", 1, 1, "57.8852", "{ 578852, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404aa6cbfb15b574, "533031", 1, 1, "53.3031", "{ 533031, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4053a04b5dcc63f1, "785046", 1, 1, "78.5046", "{ 785046, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x34f6025da1c38515, "143617", -53, 1, "1.43617e-53", "{ 143617, 10, -58 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4052f056d5cfaace, "757553", 1, 1, "75.7553", "{ 757553, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404e8bdd97f62b6b, "610927", 1, 1, "61.0927", "{ 610927, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40401cdd2f1a9fbe, "322255", 1, 1, "32.2255", "{ 322255, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40378504816f0069, "235196", 1, 1, "23.5196", "{ 235196, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3f7b6025af91fe70, "668349", -3, 1, "0.00668349", "{ 668349, 10, -8 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4030302de00d1b71, "161882", 1, 1, "16.1882", "{ 161882, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3ff15805e5f30e80, "108399", 0, 1, "1.08399", "{ 108399, 10, -5 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3ff1b7de939eadd6, "110739", 0, 1, "1.10739", "{ 110739, 10, -5 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x402da0902de00d1b, "148136", 1, 1, "14.8136", "{ 148136, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x402f75c28f5c28f6, "1573", 1, 1, "15.73", "{ 1573, 10, -2 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3ff7d3d859c8c932, "148922", 0, 1, "1.48922", "{ 148922, 10, -5 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x402bcb5dcc63f141, "138972", 1, 1, "13.8972", "{ 138972, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x401538793dd97f63, "530515", 0, 1, "5.30515", "{ 530515, 10, -5 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405cacfdf3b645a2, "114703", 2, 1, "114.703", "{ 114703, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x33046fc353fe6ab7, "620982", -63, 1, "6.20982e-63", "{ 620982, 10, -68 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4057faa7ef9db22d, "959165", 1, 1, "95.9165", "{ 959165, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4049bc32ca57a787, "514703", 1, 1, "51.4703", "{ 514703, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405ad7ced916872b, "107372", 2, 1, "107.372", "{ 107372, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3369b0ec7ab1e3ef, "499613", -61, 1, "4.99613e-61", "{ 499613, 10, -66 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405768083126e979, "936255", 1, 1, "93.6255", "{ 936255, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404b56bedfa43fe6, "546777", 1, 1, "54.6777", "{ 546777, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405a452f1a9fbe77, "105081", 2, 1, "105.081", "{ 105081, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3529b19e3400307a, "134128", -52, 1, "1.34128e-52", "{ 134128, 10, -57 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4047370a3d70a3d7, "4643", 1, 1, "46.43", "{ 4643, 10, -2 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4056f2b9f559b3d0, "917926", 1, 1, "91.7926", "{ 917926, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3563962c12b02d04, "163595", -51, 1, "1.63595e-51", "{ 163595, 10, -56 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4054c58c7e28240b, "830867", 1, 1, "83.0867", "{ 830867, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x40485c49ba5e353f, "48721", 1, 1, "48.721", "{ 48721, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405432eb1c432ca5, "807956", 1, 1, "80.7956", "{ 807956, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3771dbe70c789ecf, "128132", -41, 1, "1.28132e-41", "{ 128132, 10, -46 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404f012a30553261, "620091", 1, 1, "62.0091", "{ 620091, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4047ac56d5cfaace, "473464", 1, 1, "47.3464", "{ 473464, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x405a0a7ef9db22d1, "104164", 2, 1, "104.164", "{ 104164, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x37d4b357c234860d, "950523", -40, 1, "9.50523e-40", "{ 950523, 10, -45 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x404a6c25aee631f9, "528449", 1, 1, "52.8449", "{ 528449, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4033da9930be0ded, "198539", 1, 1, "19.8539", "{ 198539, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3ffb3be22e5de15d, "170212", 0, 1, "1.70212", "{ 170212, 10, -5 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x402cb5f6fd21ff2e, "143554", 1, 1, "14.3554", "{ 143554, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4030a57a786c2268, "166464", 1, 1, "16.6464", "{ 166464, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x402ae0c49ba5e354, "13439", 1, 1, "13.439", "{ 13439, 10, -3 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x4010194af4f0d845, "40247", 0, 1, "4.0247", "{ 40247, 10, -4 }"), // json tblastx_db_bf_fmt15 (+1)
    (0xbff0000000000000, "1", 0, -1, "-1", "{ -1, 10, 0 }"), // json tblastx_db_bf_fmt15 (+1)
    (0x3fc1257995d6905f, "133956144488482", -1, 1, "0.133956144488482", "{ 133956144488482, 10, -15 }"), // json tblastx_db_fmt13 (+1)
    (0x3fd453a7f043dabf, "317605957635731", -1, 1, "0.317605957635731", "{ 317605957635731, 10, -15 }"), // json tblastx_db_fmt13 (+1)
    (0x3fd9ad7faf594cdf, "401214524497119", -1, 1, "0.401214524497119", "{ 401214524497119, 10, -15 }"), // json tblastx_db_fmt13 (+1)
    (0x2c8161e04c0712e5, "260411", -94, 1, "2.60411e-94", "{ 260411, 10, -99 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x320efa7e388aac71, "143632", -67, 1, "1.43632e-67", "{ 143632, 10, -72 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x3fe7cf3775b81301, "744045", -1, 1, "0.744045", "{ 744045, 10, -6 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x3fed9f79420b3d4b, "925717", -1, 1, "0.925717", "{ 925717, 10, -6 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x40240083126e978d, "10001", 1, 1, "10.001", "{ 10001, 10, -3 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x40185b9628cbd124, "608944", 0, 1, "6.08944", "{ 608944, 10, -5 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x40233b4e11dbca97, "961583", 0, 1, "9.61583", "{ 961583, 10, -5 }"), // xml blastp_subject_fmt16_bl2seqlegacy
    (0x402239057d1782d4, "911137", 0, 1, "9.11137", "{ 911137, 10, -5 }"), // xml blastp_subject_fmt16_bl2seqlegacy
];

#[rustfmt::skip]
const EDGE_ROWS: &[Row] = &[
    (0x8000000000000000, "-0", 0, -1, "-0", "{ -0, 10, 0 }"), // negative-zero
    (0x0010000000000000, "22250738585072", -308, 1, "2.2250738585072e-308", "{ 22250738585072, 10, -321 }"), // DBL_MIN
    (0x8010000000000000, "22250738585072", -308, -1, "-2.2250738585072e-308", "{ -22250738585072, 10, -321 }"), // -DBL_MIN
    (0x7fefffffffffffff, "179769313486232", 308, 1, "1.79769313486232e+308", "{ 179769313486232, 10, 294 }"), // DBL_MAX
    (0xffefffffffffffff, "179769313486232", 308, -1, "-1.79769313486232e+308", "{ -179769313486232, 10, 294 }"), // -DBL_MAX
    (0x0000000000000001, "494065645841247", -324, 1, "4.94065645841247e-324", "{ 494065645841247, 10, -338 }"), // DBL_TRUE_MIN
    (0x8000000000000001, "494065645841247", -324, -1, "-4.94065645841247e-324", "{ -494065645841247, 10, -338 }"), // -DBL_TRUE_MIN
    (0x000fffffffffffff, "22250738585072", -308, 1, "2.2250738585072e-308", "{ 22250738585072, 10, -321 }"), // largest-subnormal
    (0x0000000000000003, "148219693752374", -323, 1, "1.48219693752374e-323", "{ 148219693752374, 10, -337 }"), // subnormal-3
    (0x0000000000003039, "609924039791019", -320, 1, "6.09924039791019e-320", "{ 609924039791019, 10, -334 }"), // subnormal-12345
    (0x0000000000004000, "809477154146298", -320, 1, "8.09477154146298e-320", "{ 809477154146298, 10, -334 }"), // subnormal-2^-1060
    (0x7feffffffffffffe, "179769313486232", 308, 1, "1.79769313486232e+308", "{ 179769313486232, 10, 294 }"), // DBL_MAX-prev
    (0x3ff0000000000000, "1", 0, 1, "1", "{ 1, 10, 0 }"), // value:1
    (0x4000000000000000, "2", 0, 1, "2", "{ 2, 10, 0 }"), // value:2
    (0x4014000000000000, "5", 0, 1, "5", "{ 5, 10, 0 }"), // value:5
    (0xc014000000000000, "5", 0, -1, "-5", "{ -5, 10, 0 }"), // value:-5
    (0x4024000000000000, "1", 1, 1, "10", "{ 1, 10, 1 }"), // value:10
    (0x3fb999999999999a, "1", -1, 1, "0.1", "{ 1, 10, -1 }"), // value:0.1
    (0x3fe0000000000000, "5", -1, 1, "0.5", "{ 5, 10, -1 }"), // value:0.5
    (0x408f400000000000, "1", 3, 1, "1000", "{ 1, 10, 3 }"), // value:1000
    (0x40c81c8000000000, "12345", 4, 1, "12345", "{ 12345, 10, 0 }"), // value:12345
    (0x3f1a36e2eb1c432d, "1", -4, 1, "0.0001", "{ 1, 10, -4 }"), // value:0.0001
    (0x3ee4f8b588e368f1, "1", -5, 1, "1e-05", "{ 1, 10, -5 }"), // value:0.00001
    (0x3f50624dd2f1a9fc, "1", -3, 1, "0.001", "{ 1, 10, -3 }"), // value:0.001
    (0xbf1a36e2eb1c432d, "1", -4, -1, "-0.0001", "{ -1, 10, -4 }"), // value:-0.0001
    (0x42d6bcc41e900000, "1", 14, 1, "100000000000000", "{ 1, 10, 14 }"), // value:1e14
    (0x430c6bf526340000, "1", 15, 1, "1e+15", "{ 1, 10, 15 }"), // value:1e15
    (0x4341c37937e08000, "1", 16, 1, "1e+16", "{ 1, 10, 16 }"), // value:1e16
    (0x42dc12218377de40, "123456789012345", 14, 1, "123456789012345", "{ 123456789012345, 10, 0 }"), // value:123456789012345
    (0x43118b54f22aeb00, "123456789012346", 15, 1, "1.23456789012346e+15", "{ 123456789012346, 10, 1 }"), // value:1234567890123456
    (0x4345ee2a2eb5a5c4, "123456789012346", 16, 1, "1.23456789012346e+16", "{ 123456789012346, 10, 2 }"), // value:12345678901234567
    (0x42d6bcc41e8ffffd, "1", 14, 1, "100000000000000", "{ 1, 10, 14 }"), // value:99999999999999.95
    (0x430c6bf52633fff8, "999999999999999", 14, 1, "999999999999999", "{ 999999999999999, 10, 0 }"), // value:999999999999999
    (0x3fd3333333333334, "3", -1, 1, "0.3", "{ 3, 10, -1 }"), // value:0.30000000000000004
    (0x7e37e43c8800759c, "1", 300, 1, "1e+300", "{ 1, 10, 300 }"), // value:1e300
    (0x01a56e1fc2f8f359, "1", -300, 1, "1e-300", "{ 1, 10, -300 }"), // value:1e-300
    (0x7fe1ccf385ebc8a0, "1", 308, 1, "1e+308", "{ 1, 10, 308 }"), // value:1e308
    (0x000730d67819e8d2, "1", -308, 1, "1e-308", "{ 1, 10, -308 }"), // value:1e-308
    (0x0000000000000002, "988131291682493", -324, 1, "9.88131291682493e-324", "{ 988131291682493, 10, -338 }"), // value:1e-323
    (0x4480f0cf064dd592, "1", 22, 1, "1e+22", "{ 1, 10, 22 }"), // value:1e22
    (0x44b52d02c7e14af6, "1", 23, 1, "1e+23", "{ 1, 10, 23 }"), // value:1e23
    (0x44b52d02c7e14af7, "1", 23, 1, "1e+23", "{ 1, 10, 23 }"), // value:1.0000000000000001e23
    (0x3f0a36e2eb1c432d, "5", -5, 1, "5e-05", "{ 5, 10, -5 }"), // value:5e-5
    (0x3f0a36e2eb1c432b, "5", -5, 1, "5e-05", "{ 5, 10, -5 }"), // value:4.999999999999999e-5
    (0x3efa36e2eb1c432d, "25", -5, 1, "2.5e-05", "{ 25, 10, -6 }"), // value:2.5e-5
    (0x3f1a36d1bd105b06, "99999", -5, 1, "9.9999e-05", "{ 99999, 10, -9 }"), // value:9.9999e-5
    (0x3ff0000000000017, "100000000000001", 0, 1, "1.00000000000001", "{ 100000000000001, 10, -14 }"), // value:1.000000000000005
    (0x3ff0000000000016, "1", 0, 1, "1", "{ 1, 10, 0 }"), // value:1.0000000000000049
    (0x3ff3c0ca428c59f3, "123456789012345", 0, 1, "1.23456789012345", "{ 123456789012345, 10, -14 }"), // value:1.234567890123455
    (0x4023fffffffffffd, "999999999999999", 0, 1, "9.99999999999999", "{ 999999999999999, 10, -14 }"), // value:9.999999999999995
    (0x4023fffffffffffa, "999999999999999", 0, 1, "9.99999999999999", "{ 999999999999999, 10, -14 }"), // value:9.99999999999999
    (0x42d6bcc41e8fffe0, "999999999999995", 13, 1, "99999999999999.5", "{ 999999999999995, 10, -1 }"), // value:99999999999999.5
    (0x3feffffffffffffb, "999999999999999", -1, 1, "0.999999999999999", "{ 999999999999999, 10, -15 }"), // value:0.9999999999999995
    (0x41cdcd64ffffffff, "1", 9, 1, "1000000000", "{ 1, 10, 9 }"), // value:999999999.9999999
    (0x4011c71c71c71c72, "444444444444444", 0, 1, "4.44444444444444", "{ 444444444444444, 10, -14 }"), // value:4.444444444444445
    (0x401638e38e38e38e, "555555555555556", 0, 1, "5.55555555555556", "{ 555555555555556, 10, -14 }"), // value:5.555555555555555
    (0x3fd4538ef34d6a16, "3176", -1, 1, "0.3176", "{ 3176, 10, -4 }"), // value:0.3176
    (0x3fc126e978d4fdf4, "134", -1, 1, "0.134", "{ 134, 10, -3 }"), // value:0.134
    (0x3fd9ad42c3c9eecc, "4012", -1, 1, "0.4012", "{ 4012, 10, -4 }"), // value:0.4012
    (0x3fca9fbe76c8b439, "208", -1, 1, "0.208", "{ 208, 10, -3 }"), // value:0.208
    (0x3ff5fc2656abde40, "137406", 0, 1, "1.37406", "{ 137406, 10, -5 }"), // value:1.37406
    (0x3fe6bd4280ae1049, "710603", -1, 1, "0.710603", "{ 710603, 10, -6 }"), // value:0.710603
    (0x3fb3404ea4a8c155, "752", -2, 1, "0.0752", "{ 752, 10, -4 }"), // value:0.0752
    (0x1a90991a9bfa58c8, "1", -180, 1, "1e-180", "{ 1, 10, -180 }"), // value:1e-180
    (0x1ac4bf6142f8eefa, "1", -179, 1, "1e-179", "{ 1, 10, -179 }"), // value:1e-179
    (0x1aa0991a9bfa58c8, "2", -180, 1, "2e-180", "{ 2, 10, -180 }"), // value:2e-180
    (0x2b2bff2ee48e0530, "1", -100, 1, "1e-100", "{ 1, 10, -100 }"), // value:1e-100
    (0x01b01297d23ab683, "15", -300, 1, "1.5e-300", "{ 15, 10, -301 }"), // value:1.5e-300
    (0x358dee7a4ad4b81f, "1", -50, 1, "1e-50", "{ 1, 10, -50 }"), // value:1e-50
    (0x3f4d7dbf487fcb92, "9", -4, 1, "0.0009", "{ 9, 10, -4 }"), // value:0.0009
    (0x40c3880000000000, "1", 4, 1, "10000", "{ 1, 10, 4 }"), // value:10000
    (0xc056280000000000, "88625", 1, -1, "-88.625", "{ -88625, 10, -3 }"), // value:-88.625
    (0xab2bff2ee48e0530, "1", -100, -1, "-1e-100", "{ -1, 10, -100 }"), // value:-1e-100
    (0x9b10eab42ac922cb, "260916226603389", -178, -1, "-2.60916226603389e-178", "{ -260916226603389, 10, -192 }"), // value:-2.60916226603389e-178
    (0x4058600000000000, "975", 1, 1, "97.5", "{ 975, 10, -1 }"), // value:97.5
    (0x4058e00000000000, "995", 1, 1, "99.5", "{ 995, 10, -1 }"), // value:99.5
    (0x4058a66666666666, "986", 1, 1, "98.6", "{ 986, 10, -1 }"), // value:98.6
    (0x4040aaaaaaaaaaa6, "333333333333333", 1, 1, "33.3333333333333", "{ 333333333333333, 10, -13 }"), // value:33.3333333333333
    (0x4050aaaaaaaaaaad, "666666666666667", 1, 1, "66.6666666666667", "{ 666666666666667, 10, -13 }"), // value:66.6666666666667
    (0x4341c37937e08001, "1", 16, 1, "1e+16", "{ 1, 10, 16 }"), // 1e16-next
    (0x4341c37937e07fff, "1", 16, 1, "1e+16", "{ 1, 10, 16 }"), // 1e16-prev
    (0x3c9cd2b297d889bc, "1", -16, 1, "1e-16", "{ 1, 10, -16 }"), // 1e-16
    (0x3c9cd2b297d889bd, "1", -16, 1, "1e-16", "{ 1, 10, -16 }"), // 1e-16-next
    (0x3c9cd2b297d889bb, "1", -16, 1, "1e-16", "{ 1, 10, -16 }"), // 1e-16-prev
    (0x4693b8b5b5056e17, "1", 32, 1, "1e+32", "{ 1, 10, 32 }"), // 1e32
    (0x4693b8b5b5056e18, "1", 32, 1, "1e+32", "{ 1, 10, 32 }"), // 1e32-next
    (0x4693b8b5b5056e16, "1", 32, 1, "1e+32", "{ 1, 10, 32 }"), // 1e32-prev
    (0x3949f623d5a8a733, "1", -32, 1, "1e-32", "{ 1, 10, -32 }"), // 1e-32
    (0x3949f623d5a8a734, "1", -32, 1, "1e-32", "{ 1, 10, -32 }"), // 1e-32-next
    (0x3949f623d5a8a732, "1", -32, 1, "1e-32", "{ 1, 10, -32 }"), // 1e-32-prev
    (0x49e5e531a0a1c873, "1", 48, 1, "1e+48", "{ 1, 10, 48 }"), // 1e48
    (0x49e5e531a0a1c874, "1", 48, 1, "1e+48", "{ 1, 10, 48 }"), // 1e48-next
    (0x49e5e531a0a1c872, "1", 48, 1, "1e+48", "{ 1, 10, 48 }"), // 1e48-prev
    (0x35f7624f8a762fd8, "1", -48, 1, "1e-48", "{ 1, 10, -48 }"), // 1e-48
    (0x35f7624f8a762fd9, "1", -48, 1, "1e-48", "{ 1, 10, -48 }"), // 1e-48-next
    (0x35f7624f8a762fd7, "1", -48, 1, "1e-48", "{ 1, 10, -48 }"), // 1e-48-prev
    (0x4d384f03e93ff9f5, "1", 64, 1, "1e+64", "{ 1, 10, 64 }"), // 1e64
    (0x4d384f03e93ff9f6, "1", 64, 1, "1e+64", "{ 1, 10, 64 }"), // 1e64-next
    (0x4d384f03e93ff9f4, "1", 64, 1, "1e+64", "{ 1, 10, 64 }"), // 1e64-prev
    (0x32a50ffd44f4a73d, "1", -64, 1, "1e-64", "{ 1, 10, -64 }"), // 1e-64
    (0x32a50ffd44f4a73e, "1", -64, 1, "1e-64", "{ 1, 10, -64 }"), // 1e-64-next
    (0x32a50ffd44f4a73c, "1", -64, 1, "1e-64", "{ 1, 10, -64 }"), // 1e-64-prev
    (0x508afcef51f0fb5f, "1", 80, 1, "1e+80", "{ 1, 10, 80 }"), // 1e80
    (0x508afcef51f0fb60, "1", 80, 1, "1e+80", "{ 1, 10, 80 }"), // 1e80-next
    (0x508afcef51f0fb5e, "1", 80, 1, "1e+80", "{ 1, 10, 80 }"), // 1e80-prev
    (0x2f52f8ac174d6123, "1", -80, 1, "1e-80", "{ 1, 10, -80 }"), // 1e-80
    (0x2f52f8ac174d6124, "1", -80, 1, "1e-80", "{ 1, 10, -80 }"), // 1e-80-next
    (0x2f52f8ac174d6122, "1", -80, 1, "1e-80", "{ 1, 10, -80 }"), // 1e-80-prev
    (0x53ddf67562d8b363, "1", 96, 1, "1e+96", "{ 1, 10, 96 }"), // 1e96
    (0x53ddf67562d8b364, "1", 96, 1, "1e+96", "{ 1, 10, 96 }"), // 1e96-next
    (0x53ddf67562d8b362, "1", 96, 1, "1e+96", "{ 1, 10, 96 }"), // 1e96-prev
    (0x2c0116805effaeaa, "1", -96, 1, "1e-96", "{ 1, 10, -96 }"), // 1e-96
    (0x2c0116805effaeab, "1", -96, 1, "1e-96", "{ 1, 10, -96 }"), // 1e-96-next
    (0x2c0116805effaea9, "1", -96, 1, "1e-96", "{ 1, 10, -96 }"), // 1e-96-prev
    (0x5730a1f5b8132466, "1", 112, 1, "1e+112", "{ 1, 10, 112 }"), // 1e112
    (0x5730a1f5b8132467, "1", 112, 1, "1e+112", "{ 1, 10, 112 }"), // 1e112-next
    (0x5730a1f5b8132465, "1", 112, 1, "1e+112", "{ 1, 10, 112 }"), // 1e112-prev
    (0x28aec866b79e0cba, "1", -112, 1, "1e-112", "{ 1, 10, -112 }"), // 1e-112
    (0x28aec866b79e0cbb, "1", -112, 1, "1e-112", "{ 1, 10, -112 }"), // 1e-112-next
    (0x28aec866b79e0cb9, "1", -112, 1, "1e-112", "{ 1, 10, -112 }"), // 1e-112-prev
    (0x5a827748f9301d32, "1", 128, 1, "1e+128", "{ 1, 10, 128 }"), // 1e128
    (0x5a827748f9301d33, "1", 128, 1, "1e+128", "{ 1, 10, 128 }"), // 1e128-next
    (0x5a827748f9301d31, "1", 128, 1, "1e+128", "{ 1, 10, 128 }"), // 1e128-prev
    (0x255bba08cf8c979d, "1", -128, 1, "1e-128", "{ 1, 10, -128 }"), // 1e-128
    (0x255bba08cf8c979e, "1", -128, 1, "1e-128", "{ 1, 10, -128 }"), // 1e-128-next
    (0x255bba08cf8c979c, "1", -128, 1, "1e-128", "{ 1, 10, -128 }"), // 1e-128-prev
    (0x5dd4805738b51a75, "1", 144, 1, "1e+144", "{ 1, 10, 144 }"), // 1e144
    (0x5dd4805738b51a76, "1", 144, 1, "1e+144", "{ 1, 10, 144 }"), // 1e144-next
    (0x5dd4805738b51a74, "1", 144, 1, "1e+144", "{ 1, 10, 144 }"), // 1e144-prev
    (0x2208f9574dcf8a70, "1", -144, 1, "1e-144", "{ 1, 10, -144 }"), // 1e-144
    (0x2208f9574dcf8a71, "1", -144, 1, "1e-144", "{ 1, 10, -144 }"), // 1e-144-next
    (0x2208f9574dcf8a6f, "1", -144, 1, "1e-144", "{ 1, 10, -144 }"), // 1e-144-prev
    (0x6126c2d4256ffcc3, "1", 160, 1, "1e+160", "{ 1, 10, 160 }"), // 1e160
    (0x6126c2d4256ffcc4, "1", 160, 1, "1e+160", "{ 1, 10, 160 }"), // 1e160-next
    (0x6126c2d4256ffcc2, "1", 160, 1, "1e+160", "{ 1, 10, 160 }"), // 1e160-prev
    (0x1eb67e9c127b6e74, "1", -160, 1, "1e-160", "{ 1, 10, -160 }"), // 1e-160
    (0x1eb67e9c127b6e75, "1", -160, 1, "1e-160", "{ 1, 10, -160 }"), // 1e-160-next
    (0x1eb67e9c127b6e73, "1", -160, 1, "1e-160", "{ 1, 10, -160 }"), // 1e-160-prev
    (0x647945145230b378, "1", 176, 1, "1e+176", "{ 1, 10, 176 }"), // 1e176
    (0x647945145230b379, "1", 176, 1, "1e+176", "{ 1, 10, 176 }"), // 1e176-next
    (0x647945145230b377, "1", 176, 1, "1e+176", "{ 1, 10, 176 }"), // 1e176-prev
    (0x1b6442e4fb671960, "1", -176, 1, "1e-176", "{ 1, 10, -176 }"), // 1e-176
    (0x1b6442e4fb671961, "1", -176, 1, "1e-176", "{ 1, 10, -176 }"), // 1e-176-next
    (0x1b6442e4fb67195f, "1", -176, 1, "1e-176", "{ 1, 10, -176 }"), // 1e-176-prev
    (0x67cc0e1ef1a724eb, "1", 192, 1, "1e+192", "{ 1, 10, 192 }"), // 1e192
    (0x67cc0e1ef1a724ec, "1", 192, 1, "1e+192", "{ 1, 10, 192 }"), // 1e192-next
    (0x67cc0e1ef1a724ea, "1", 192, 1, "1e+192", "{ 1, 10, 192 }"), // 1e192-prev
    (0x18123ff06eea847a, "1", -192, 1, "1e-192", "{ 1, 10, -192 }"), // 1e-192
    (0x18123ff06eea847b, "1", -192, 1, "1e-192", "{ 1, 10, -192 }"), // 1e-192-next
    (0x18123ff06eea8479, "1", -192, 1, "1e-192", "{ 1, 10, -192 }"), // 1e-192-prev
    (0x6b1f25c186a6f04c, "1", 208, 1, "1e+208", "{ 1, 10, 208 }"), // 1e208
    (0x6b1f25c186a6f04d, "1", 208, 1, "1e+208", "{ 1, 10, 208 }"), // 1e208-next
    (0x6b1f25c186a6f04b, "1", 208, 1, "1e+208", "{ 1, 10, 208 }"), // 1e208-prev
    (0x14c0701bd527b498, "1", -208, 1, "1e-208", "{ 1, 10, -208 }"), // 1e-208
    (0x14c0701bd527b499, "1", -208, 1, "1e-208", "{ 1, 10, -208 }"), // 1e-208-next
    (0x14c0701bd527b497, "1", -208, 1, "1e-208", "{ 1, 10, -208 }"), // 1e-208-prev
    (0x6e714a52dffc6799, "1", 224, 1, "1e+224", "{ 1, 10, 224 }"), // 1e224
    (0x6e714a52dffc679a, "1", 224, 1, "1e+224", "{ 1, 10, 224 }"), // 1e224-next
    (0x6e714a52dffc6798, "1", 224, 1, "1e+224", "{ 1, 10, 224 }"), // 1e224-prev
    (0x116d9ca79d89462a, "1", -224, 1, "1e-224", "{ 1, 10, -224 }"), // 1e-224
    (0x116d9ca79d89462b, "1", -224, 1, "1e-224", "{ 1, 10, -224 }"), // 1e-224-next
    (0x116d9ca79d894629, "1", -224, 1, "1e-224", "{ 1, 10, -224 }"), // 1e-224-prev
    (0x71c33234de7ad7e3, "1", 240, 1, "1e+240", "{ 1, 10, 240 }"), // 1e240
    (0x71c33234de7ad7e4, "1", 240, 1, "1e+240", "{ 1, 10, 240 }"), // 1e240-next
    (0x71c33234de7ad7e2, "1", 240, 1, "1e+240", "{ 1, 10, 240 }"), // 1e240-prev
    (0x0e1aac0bf9b9e65c, "1", -240, 1, "1e-240", "{ 1, 10, -240 }"), // 1e-240
    (0x0e1aac0bf9b9e65d, "1", -240, 1, "1e-240", "{ 1, 10, -240 }"), // 1e-240-next
    (0x0e1aac0bf9b9e65b, "1", -240, 1, "1e-240", "{ 1, 10, -240 }"), // 1e-240-prev
    (0x75154fdd7f73bf3c, "1", 256, 1, "1e+256", "{ 1, 10, 256 }"), // 1e256
    (0x75154fdd7f73bf3d, "1", 256, 1, "1e+256", "{ 1, 10, 256 }"), // 1e256-next
    (0x75154fdd7f73bf3b, "1", 256, 1, "1e+256", "{ 1, 10, 256 }"), // 1e256-prev
    (0x0ac8062864ac6f43, "1", -256, 1, "1e-256", "{ 1, 10, -256 }"), // 1e-256
    (0x0ac8062864ac6f44, "1", -256, 1, "1e-256", "{ 1, 10, -256 }"), // 1e-256-next
    (0x0ac8062864ac6f42, "1", -256, 1, "1e-256", "{ 1, 10, -256 }"), // 1e-256-prev
    (0x4002aaaaaaaaaaab, "233333333333333", 0, 1, "2.33333333333333", "{ 233333333333333, 10, -14 }"), // computed
    (0x3fe5555555555555, "666666666666667", -1, 1, "0.666666666666667", "{ 666666666666667, 10, -15 }"), // computed
    (0x3fd5555555555555, "333333333333333", -1, 1, "0.333333333333333", "{ 333333333333333, 10, -15 }"), // computed
    (0x3fc2492492492492, "142857142857143", -1, 1, "0.142857142857143", "{ 142857142857143, 10, -15 }"), // computed
    (0x4051b55555555555, "708333333333333", 1, 1, "70.8333333333333", "{ 708333333333333, 10, -13 }"), // computed
    (0x4024955555555555, "102916666666667", 1, 1, "10.2916666666667", "{ 102916666666667, 10, -13 }"), // computed
    (0x400921fb54442d18, "314159265358979", 0, 1, "3.14159265358979", "{ 314159265358979, 10, -14 }"), // computed
    (0x4005bf0a8b145769, "271828182845905", 0, 1, "2.71828182845905", "{ 271828182845905, 10, -14 }"), // computed
    (0x3fe62e42fefa39ef, "693147180559945", -1, 1, "0.693147180559945", "{ 693147180559945, 10, -15 }"), // computed
    (0x3ff71547652b82fe, "144269504088896", 0, 1, "1.44269504088896", "{ 144269504088896, 10, -14 }"), // computed
    (0x3d3a56e0c2ac7f75, "935762296884017", -14, 1, "9.35762296884017e-14", "{ 935762296884017, 10, -28 }"), // computed
    (0x0028000000000000, "66752215755216", -308, 1, "6.6752215755216e-308", "{ 66752215755216, 10, -321 }"), // computed
    (0x7fe0000000000000, "898846567431158", 307, 1, "8.98846567431158e+307", "{ 898846567431158, 10, 293 }"), // computed
    (0x3feb5c09e32d2a5a, "854985183433495", -1, 1, "0.854985183433495", "{ 854985183433495, 10, -15 }"), // not-correctly-rounded:blast
    (0x40bd57a18b9539bd, "751163103611621", 3, 1, "7511.63103611621", "{ 751163103611621, 10, -11 }"), // not-correctly-rounded:blast
    (0x3fd613716aefcc27, "34493670886076", -1, 1, "0.34493670886076", "{ 34493670886076, 10, -14 }"), // not-correctly-rounded:blast
    (0x40b990e13c913975, "654487983043339", 3, 1, "6544.87983043339", "{ 654487983043339, 10, -11 }"), // not-correctly-rounded:blast
    (0x3fd9ec8e951033d9, "40506329113924", -1, 1, "0.40506329113924", "{ 40506329113924, 10, -14 }"), // not-correctly-rounded:blast
    (0x3fdb99ba1b65c73b, "431257750277939", -1, 1, "0.431257750277939", "{ 431257750277939, 10, -15 }"), // not-correctly-rounded:blast
    (0x3fe6f6f6f7b02b1b, "717647060171063", -1, 1, "0.717647060171063", "{ 717647060171063, 10, -15 }"), // not-correctly-rounded:blast
    (0x3fe3892d3892d389, "61049519586105", -1, 1, "0.61049519586105", "{ 61049519586105, 10, -14 }"), // not-correctly-rounded:blast
    (0x4041e7bd9b21d5de, "358104738154613", 1, 1, "35.8104738154613", "{ 358104738154613, 10, -13 }"), // not-correctly-rounded:blast
    (0x3ff9aa48817effaf, "160407305321131", 0, 1, "1.60407305321131", "{ 160407305321131, 10, -14 }"), // not-correctly-rounded:blast
    (0x10c2246226f8543b, "59830328936817", -228, 1, "5.9830328936817e-228", "{ 59830328936817, 10, -241 }"), // not-correctly-rounded:blast
    (0x3fef062acb56a6c0, "969502827783081", -1, 1, "0.969502827783081", "{ 969502827783081, 10, -15 }"), // not-correctly-rounded:blast
    (0x3fce17b3c7eb8d0d, "235098335854766", -1, 1, "0.235098335854766", "{ 235098335854766, 10, -15 }"), // not-correctly-rounded:blast
    (0x40b72052d0dad4ad, "592032349937145", 3, 1, "5920.32349937145", "{ 592032349937145, 10, -11 }"), // not-correctly-rounded:blast
    (0x3fd386ee031c3308, "305110457455883", -1, 1, "0.305110457455883", "{ 305110457455883, 10, -15 }"), // not-correctly-rounded:blast
    (0x3fd1ec8e951033d9, "28006329113924", -1, 1, "0.28006329113924", "{ 28006329113924, 10, -14 }"), // not-correctly-rounded:blast
    (0x3fe2a49fb6ae4593, "582595688627565", -1, 1, "0.582595688627565", "{ 582595688627565, 10, -15 }"), // not-correctly-rounded:blast
    (0x409e21ba6103814f, "192843201070287", 3, 1, "1928.43201070287", "{ 192843201070287, 10, -11 }"), // not-correctly-rounded:blast
    (0x0919d9bdf6afac69, "801703341494089", -265, 1, "8.01703341494089e-265", "{ 801703341494089, 10, -279 }"), // not-correctly-rounded:blast
    (0x42d6bcc41e900020, "100000000000001", 14, 1, "100000000000001", "{ 100000000000001, 10, 0 }"), // not-correctly-rounded:edge
    (0x43118b54f22aeafc, "123456789012345", 15, 1, "1.23456789012345e+15", "{ 123456789012345, 10, 1 }"), // not-correctly-rounded:edge
    (0x4341c37937e08019, "100000000000001", 16, 1, "1.00000000000001e+16", "{ 100000000000001, 10, 2 }"), // not-correctly-rounded:edge
    (0x4345ee2a2eb5a5bb, "123456789012345", 16, 1, "1.23456789012345e+16", "{ 123456789012345, 10, 2 }"), // not-correctly-rounded:edge
    (0x3e90000000000000, "238418579101563", -7, 1, "2.38418579101563e-07", "{ 238418579101563, 10, -21 }"), // not-correctly-rounded:edge
    (0x0000126812ccc426, "9999", -311, 1, "9.999e-311", "{ 9999, 10, -314 }"), // not-correctly-rounded:dec
    (0x4303d4da74e41dd4, "697757410821051", 14, 1, "697757410821051", "{ 697757410821051, 10, 0 }"), // not-correctly-rounded:dec
    (0x0000b810bbffa97c, "9999", -310, 1, "9.999e-310", "{ 9999, 10, -313 }"), // not-correctly-rounded:dec
    (0x0000b815665893e1, "999998999999999", -310, 1, "9.99998999999999e-310", "{ 999998999999999, 10, -324 }"), // not-correctly-rounded:dec
    (0x8000b7e6524bb3b3, "998999999999999", -310, -1, "-9.98999999999999e-310", "{ -998999999999999, 10, -324 }"), // not-correctly-rounded:dec
    (0xc2e6148486308750, "194219028350011", 14, -1, "-194219028350011", "{ -194219028350011, 10, 0 }"), // not-correctly-rounded:dec
    (0xc2d466981b6f2030, "897234846178567", 13, -1, "-89723484617856.7", "{ -897234846178567, 10, -1 }"), // not-correctly-rounded:dec
    (0x00004424c1ae142c, "370176999999999", -310, 1, "3.70176999999999e-310", "{ 370176999999999, 10, -324 }"), // not-correctly-rounded:dec
    (0x000012d4635285a3, "10228843", -310, 1, "1.0228843e-310", "{ 10228843, 10, -317 }"), // not-correctly-rounded:dec
    (0x4322fe97a43cfc56, "267323841488439", 15, 1, "2.67323841488439e+15", "{ 267323841488439, 10, 1 }"), // not-correctly-rounded:dec
    (0xc338ed30ff78b3b1, "701619414136927", 15, -1, "-7.01619414136927e+15", "{ -701619414136927, 10, 1 }"), // not-correctly-rounded:dec
    (0x4309c544b44c71d4, "906721661390395", 14, 1, "906721661390395", "{ 906721661390395, 10, 0 }"), // not-correctly-rounded:dec
    (0xc2d6827af3fe2b10, "989986574808763", 13, -1, "-98998657480876.3", "{ -989986574808763, 10, -1 }"), // not-correctly-rounded:dec
    (0x433ae7877907ae83, "757291843162483", 15, 1, "7.57291843162483e+15", "{ 757291843162483, 10, 1 }"), // not-correctly-rounded:dec
    (0x8000126812ccc426, "9999", -311, -1, "-9.999e-311", "{ -9999, 10, -314 }"), // not-correctly-rounded:dec
    (0x8000b810bbffa97c, "9999", -310, -1, "-9.999e-310", "{ -9999, 10, -313 }"), // not-correctly-rounded:dec
    (0x00000c2350ef780e, "659371179999999", -311, 1, "6.59371179999999e-311", "{ 659371179999999, 10, -325 }"), // not-correctly-rounded:dec
    (0x00000b6e7d813f10, "621", -311, 1, "6.21e-311", "{ 621, 10, -313 }"), // not-correctly-rounded:dec
    (0xcc2e2a5f51569a76, "946760093448025", 58, -1, "-9.46760093448025e+58", "{ -946760093448025, 10, 44 }"), // not-correctly-rounded:dec
    (0xa664a061d492d6fa, "975073109668934", -124, -1, "-9.75073109668934e-124", "{ -975073109668934, 10, -138 }"), // not-correctly-rounded:dec
    (0x000079adeba005e4, "660999999999999", -310, 1, "6.60999999999999e-310", "{ 660999999999999, 10, -324 }"), // not-correctly-rounded:dec
    (0x00007b560ca79fba, "67", -310, 1, "6.7e-310", "{ 67, 10, -311 }"), // not-correctly-rounded:dec
    (0xc31342bd84b83e1c, "135535157510541", 15, -1, "-1.35535157510541e+15", "{ -135535157510541, 10, 1 }"), // not-correctly-rounded:dec
    (0x430787a821d72564, "827885082043565", 14, 1, "827885082043565", "{ 827885082043565, 10, 0 }"), // not-correctly-rounded:dec
    (0x00009247c3dff6b9, "794639999999999", -310, 1, "7.94639999999999e-310", "{ 794639999999999, 10, -324 }"), // not-correctly-rounded:dec
    (0xc2bcceef66f019c0, "316751053414657", 13, -1, "-31675105341465.7", "{ -316751053414657, 10, -1 }"), // not-correctly-rounded:dec
    (0x4317dceec05e006c, "167921061301455", 15, 1, "1.67921061301455e+15", "{ 167921061301455, 10, 1 }"), // not-correctly-rounded:bits
    (0xc3085000b037c0d4, "855420415965211", 14, -1, "-855420415965211", "{ -855420415965211, 10, 0 }"), // not-correctly-rounded:bits
    (0x42f9e4f9a64312e8, "455539706573103", 14, 1, "455539706573103", "{ 455539706573103, 10, 0 }"), // not-correctly-rounded:bits
    (0x433f637fc1d8da1f, "883512464224105", 15, 1, "8.83512464224105e+15", "{ 883512464224105, 10, 1 }"), // not-correctly-rounded:bits
    (0xc2a49be35b266480, "113298834440183", 13, -1, "-11329883444018.3", "{ -113298834440183, 10, -1 }"), // not-correctly-rounded:bits
    (0xc2fe96176a061608, "538079788032353", 14, -1, "-538079788032353", "{ -538079788032353, 10, 0 }"), // not-correctly-rounded:bits
    (0x433b6560307c17ab, "771128817589239", 15, 1, "7.71128817589239e+15", "{ 771128817589239, 10, 1 }"), // not-correctly-rounded:bits
    (0x7d4ad6f950f2c4ad, "342834166456578", 295, 1, "3.42834166456578e+295", "{ 342834166456578, 10, 281 }"), // not-correctly-rounded:bits
    (0xc309e802b2373854, "911496586913547", 14, -1, "-911496586913547", "{ -911496586913547, 10, 0 }"), // not-correctly-rounded:bits
    (0xc3270c1cae8fd3ae, "324362089581205", 15, -1, "-3.24362089581205e+15", "{ -324362089581205, 10, 1 }"), // not-correctly-rounded:bits
    (0x433478123acea481, "576151922557863", 15, 1, "5.76151922557863e+15", "{ 576151922557863, 10, 1 }"), // not-correctly-rounded:bits
    (0x27c4267f5c88743a, "399535114666255", -117, 1, "3.99535114666255e-117", "{ 399535114666255, 10, -131 }"), // not-correctly-rounded:bits
    (0xc306652b515c9524, "787960776331941", 14, -1, "-787960776331941", "{ -787960776331941, 10, 0 }"), // not-correctly-rounded:bits
    (0xc2d140370a97a890, "758699960808023", 13, -1, "-75869996080802.3", "{ -758699960808023, 10, -1 }"), // not-correctly-rounded:bits
    (0x43081aa70ee812b4, "848088031625815", 14, 1, "848088031625815", "{ 848088031625815, 10, 0 }"), // not-correctly-rounded:bits
    (0x42f19a0d8db9d768, "309653600443767", 14, 1, "309653600443767", "{ 309653600443767, 10, 0 }"), // not-correctly-rounded:bits
    (0x42e59aca8f437650, "190036540267443", 14, 1, "190036540267443", "{ 190036540267443, 10, 0 }"), // not-correctly-rounded:bits
    (0xc2f2bc465801ff68, "329597493190647", 14, -1, "-329597493190647", "{ -329597493190647, 10, 0 }"), // not-correctly-rounded:bits
    (0x42e960114fba8ed0, "223203183940727", 14, 1, "223203183940727", "{ 223203183940727, 10, 0 }"), // not-correctly-rounded:bits
    (0xc2cdff663012d060, "659655365851207", 13, -1, "-65965536585120.7", "{ -659655365851207, 10, -1 }"), // not-correctly-rounded:bits
    (0x430b66abc04d9f04, "964089027867617", 14, 1, "964089027867617", "{ 964089027867617, 10, 0 }"), // not-correctly-rounded:bits
    (0x42a8bacd5e923580, "135952942553867", 13, 1, "13595294255386.7", "{ 135952942553867, 10, -1 }"), // not-correctly-rounded:bits
    (0xcfb2f44aa56f2654, "857325714763101", 75, -1, "-8.57325714763101e+75", "{ -857325714763101, 10, 61 }"), // not-correctly-rounded:bits
    (0x42dfd68d85f2dea0, "140025431313275", 14, 1, "140025431313275", "{ 140025431313275, 10, 0 }"), // not-correctly-rounded:bits
    (0x42ef6699e84afbd0, "276204234102751", 14, 1, "276204234102751", "{ 276204234102751, 10, 0 }"), // not-correctly-rounded:bits
    (0xc330d323eeda5a4d, "473575091197397", 15, -1, "-4.73575091197397e+15", "{ -473575091197397, 10, 1 }"), // not-correctly-rounded:bits
    (0xc3296094d5a3cc06, "357153338673715", 15, -1, "-3.57153338673715e+15", "{ -357153338673715, 10, 1 }"), // not-correctly-rounded:bits
    (0x4306f359bf471214, "807502034297411", 14, 1, "807502034297411", "{ 807502034297411, 10, 0 }"), // not-correctly-rounded:bits
    (0xc31c7e1e54d34cdc, "200499202128773", 15, -1, "-2.00499202128773e+15", "{ -200499202128773, 10, 1 }"), // not-correctly-rounded:bits
    (0xc2f8c54e53097668, "435771227019111", 14, -1, "-435771227019111", "{ -435771227019111, 10, 0 }"), // not-correctly-rounded:bits
];

#[rustfmt::skip]
const SWEEP_ROWS: &[SweepRow] = &[
    (0x0000000000000000, 0, "", 0, 0, ""),
    (0x0000000000000000, 1, "0", 0, 0, "0"),
    (0x0000000000000000, 2, "0", 0, 0, "0"),
    (0x0000000000000000, 3, "0", 0, 0, "0"),
    (0x0000000000000000, 4, "0", 0, 0, "0"),
    (0x0000000000000000, 5, "0", 0, 0, "0"),
    (0x0000000000000000, 6, "0", 0, 0, "0"),
    (0x0000000000000000, 7, "0", 0, 0, "0"),
    (0x0000000000000000, 8, "0", 0, 0, "0"),
    (0x0000000000000000, 9, "0", 0, 0, "0"),
    (0x0000000000000000, 10, "0", 0, 0, "0"),
    (0x0000000000000000, 11, "0", 0, 0, "0"),
    (0x0000000000000000, 12, "0", 0, 0, "0"),
    (0x0000000000000000, 13, "0", 0, 0, "0"),
    (0x0000000000000000, 14, "0", 0, 0, "0"),
    (0x0000000000000000, 15, "0", 0, 0, "0"),
    (0x0000000000000000, 16, "0", 0, 0, "0"),
    (0x0000000000000000, 17, "0", 0, 0, "0"),
    (0x8000000000000000, 0, "", 0, 0, ""),
    (0x8000000000000000, 1, "-0", 0, -1, "-0"),
    (0x8000000000000000, 2, "-0", 0, -1, "-0"),
    (0x8000000000000000, 3, "-0", 0, -1, "-0"),
    (0x8000000000000000, 4, "-0", 0, -1, "-0"),
    (0x8000000000000000, 5, "-0", 0, -1, "-0"),
    (0x8000000000000000, 6, "-0", 0, -1, "-0"),
    (0x8000000000000000, 7, "-0", 0, -1, "-0"),
    (0x8000000000000000, 8, "-0", 0, -1, "-0"),
    (0x8000000000000000, 9, "-0", 0, -1, "-0"),
    (0x8000000000000000, 10, "-0", 0, -1, "-0"),
    (0x8000000000000000, 11, "-0", 0, -1, "-0"),
    (0x8000000000000000, 12, "-0", 0, -1, "-0"),
    (0x8000000000000000, 13, "-0", 0, -1, "-0"),
    (0x8000000000000000, 14, "-0", 0, -1, "-0"),
    (0x8000000000000000, 15, "-0", 0, -1, "-0"),
    (0x8000000000000000, 16, "-0", 0, -1, "-0"),
    (0x8000000000000000, 17, "-0", 0, -1, "-0"),
    (0x3ff0000000000000, 0, "", 0, 0, ""),
    (0x3ff0000000000000, 1, "1", 0, 1, "1"),
    (0x3ff0000000000000, 2, "1", 0, 1, "1"),
    (0x3ff0000000000000, 3, "1", 0, 1, "1"),
    (0x3ff0000000000000, 4, "1", 0, 1, "1"),
    (0x3ff0000000000000, 5, "1", 0, 1, "1"),
    (0x3ff0000000000000, 6, "1", 0, 1, "1"),
    (0x3ff0000000000000, 7, "1", 0, 1, "1"),
    (0x3ff0000000000000, 8, "1", 0, 1, "1"),
    (0x3ff0000000000000, 9, "1", 0, 1, "1"),
    (0x3ff0000000000000, 10, "1", 0, 1, "1"),
    (0x3ff0000000000000, 11, "1", 0, 1, "1"),
    (0x3ff0000000000000, 12, "1", 0, 1, "1"),
    (0x3ff0000000000000, 13, "1", 0, 1, "1"),
    (0x3ff0000000000000, 14, "1", 0, 1, "1"),
    (0x3ff0000000000000, 15, "1", 0, 1, "1"),
    (0x3ff0000000000000, 16, "1", 0, 1, "1"),
    (0x3ff0000000000000, 17, "1", 0, 1, "1"),
    (0x3fe0000000000000, 0, "", 0, 0, ""),
    (0x3fe0000000000000, 1, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 2, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 3, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 4, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 5, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 6, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 7, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 8, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 9, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 10, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 11, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 12, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 13, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 14, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 15, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 16, "5", -1, 1, "0.5"),
    (0x3fe0000000000000, 17, "5", -1, 1, "0.5"),
    (0x4056280000000000, 0, "", 0, 0, ""),
    (0x4056280000000000, 1, "9", 1, 1, "9e+01"),
    (0x4056280000000000, 2, "89", 1, 1, "89"),
    (0x4056280000000000, 3, "886", 1, 1, "88.6"),
    (0x4056280000000000, 4, "8863", 1, 1, "88.63"),
    (0x4056280000000000, 5, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 6, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 7, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 8, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 9, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 10, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 11, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 12, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 13, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 14, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 15, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 16, "88625", 1, 1, "88.625"),
    (0x4056280000000000, 17, "88625", 1, 1, "88.625"),
    (0x408382cd42bfa376, 0, "", 0, 0, ""),
    (0x408382cd42bfa376, 1, "6", 2, 1, "6e+02"),
    (0x408382cd42bfa376, 2, "62", 2, 1, "6.2e+02"),
    (0x408382cd42bfa376, 3, "624", 2, 1, "624"),
    (0x408382cd42bfa376, 4, "6244", 2, 1, "624.4"),
    (0x408382cd42bfa376, 5, "62435", 2, 1, "624.35"),
    (0x408382cd42bfa376, 6, "62435", 2, 1, "624.35"),
    (0x408382cd42bfa376, 7, "6243502", 2, 1, "624.3502"),
    (0x408382cd42bfa376, 8, "62435022", 2, 1, "624.35022"),
    (0x408382cd42bfa376, 9, "624350224", 2, 1, "624.350224"),
    (0x408382cd42bfa376, 10, "624350225", 2, 1, "624.350225"),
    (0x408382cd42bfa376, 11, "62435022497", 2, 1, "624.35022497"),
    (0x408382cd42bfa376, 12, "624350224969", 2, 1, "624.350224969"),
    (0x408382cd42bfa376, 13, "6243502249691", 2, 1, "624.3502249691"),
    (0x408382cd42bfa376, 14, "62435022496908", 2, 1, "624.35022496908"),
    (0x408382cd42bfa376, 15, "624350224969078", 2, 1, "624.350224969078"),
    (0x408382cd42bfa376, 16, "624350224969078", 2, 1, "624.350224969078"),
    (0x408382cd42bfa376, 17, "624350224969078", 2, 1, "624.350224969078"),
    (0x1b10eab42ac922cb, 0, "", 0, 0, ""),
    (0x1b10eab42ac922cb, 1, "3", -178, 1, "3e-178"),
    (0x1b10eab42ac922cb, 2, "26", -178, 1, "2.6e-178"),
    (0x1b10eab42ac922cb, 3, "261", -178, 1, "2.61e-178"),
    (0x1b10eab42ac922cb, 4, "2609", -178, 1, "2.609e-178"),
    (0x1b10eab42ac922cb, 5, "26092", -178, 1, "2.6092e-178"),
    (0x1b10eab42ac922cb, 6, "260916", -178, 1, "2.60916e-178"),
    (0x1b10eab42ac922cb, 7, "2609162", -178, 1, "2.609162e-178"),
    (0x1b10eab42ac922cb, 8, "26091623", -178, 1, "2.6091623e-178"),
    (0x1b10eab42ac922cb, 9, "260916226", -178, 1, "2.60916226e-178"),
    (0x1b10eab42ac922cb, 10, "2609162266", -178, 1, "2.609162266e-178"),
    (0x1b10eab42ac922cb, 11, "2609162266", -178, 1, "2.609162266e-178"),
    (0x1b10eab42ac922cb, 12, "260916226603", -178, 1, "2.60916226603e-178"),
    (0x1b10eab42ac922cb, 13, "2609162266034", -178, 1, "2.609162266034e-178"),
    (0x1b10eab42ac922cb, 14, "26091622660339", -178, 1, "2.6091622660339e-178"),
    (0x1b10eab42ac922cb, 15, "260916226603389", -178, 1, "2.60916226603389e-178"),
    (0x1b10eab42ac922cb, 16, "260916226603389", -178, 1, "2.60916226603389e-178"),
    (0x1b10eab42ac922cb, 17, "260916226603389", -178, 1, "2.60916226603389e-178"),
    (0x3ee4f8b588e368f1, 0, "", 0, 0, ""),
    (0x3ee4f8b588e368f1, 1, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 2, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 3, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 4, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 5, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 6, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 7, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 8, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 9, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 10, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 11, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 12, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 13, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 14, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 15, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 16, "1", -5, 1, "1e-05"),
    (0x3ee4f8b588e368f1, 17, "1", -5, 1, "1e-05"),
    (0x3f1a36e2eb1c432d, 0, "", 0, 0, ""),
    (0x3f1a36e2eb1c432d, 1, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 2, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 3, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 4, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 5, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 6, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 7, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 8, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 9, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 10, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 11, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 12, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 13, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 14, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 15, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 16, "1", -4, 1, "0.0001"),
    (0x3f1a36e2eb1c432d, 17, "1", -4, 1, "0.0001"),
    (0x7fefffffffffffff, 0, "", 0, 0, ""),
    (0x7fefffffffffffff, 1, "2", 308, 1, "2e+308"),
    (0x7fefffffffffffff, 2, "18", 308, 1, "1.8e+308"),
    (0x7fefffffffffffff, 3, "18", 308, 1, "1.8e+308"),
    (0x7fefffffffffffff, 4, "1798", 308, 1, "1.798e+308"),
    (0x7fefffffffffffff, 5, "17977", 308, 1, "1.7977e+308"),
    (0x7fefffffffffffff, 6, "179769", 308, 1, "1.79769e+308"),
    (0x7fefffffffffffff, 7, "1797693", 308, 1, "1.797693e+308"),
    (0x7fefffffffffffff, 8, "17976931", 308, 1, "1.7976931e+308"),
    (0x7fefffffffffffff, 9, "179769313", 308, 1, "1.79769313e+308"),
    (0x7fefffffffffffff, 10, "1797693135", 308, 1, "1.797693135e+308"),
    (0x7fefffffffffffff, 11, "17976931349", 308, 1, "1.7976931349e+308"),
    (0x7fefffffffffffff, 12, "179769313486", 308, 1, "1.79769313486e+308"),
    (0x7fefffffffffffff, 13, "1797693134862", 308, 1, "1.797693134862e+308"),
    (0x7fefffffffffffff, 14, "17976931348623", 308, 1, "1.7976931348623e+308"),
    (0x7fefffffffffffff, 15, "179769313486232", 308, 1, "1.79769313486232e+308"),
    (0x7fefffffffffffff, 16, "179769313486232", 308, 1, "1.79769313486232e+308"),
    (0x7fefffffffffffff, 17, "179769313486232", 308, 1, "1.79769313486232e+308"),
    (0x0000000000000001, 0, "", 0, 0, ""),
    (0x0000000000000001, 1, "5", -324, 1, "5e-324"),
    (0x0000000000000001, 2, "49", -324, 1, "4.9e-324"),
    (0x0000000000000001, 3, "494", -324, 1, "4.94e-324"),
    (0x0000000000000001, 4, "4941", -324, 1, "4.941e-324"),
    (0x0000000000000001, 5, "49407", -324, 1, "4.9407e-324"),
    (0x0000000000000001, 6, "494066", -324, 1, "4.94066e-324"),
    (0x0000000000000001, 7, "4940656", -324, 1, "4.940656e-324"),
    (0x0000000000000001, 8, "49406565", -324, 1, "4.9406565e-324"),
    (0x0000000000000001, 9, "494065645", -324, 1, "4.94065645e-324"),
    (0x0000000000000001, 10, "4940656458", -324, 1, "4.940656458e-324"),
    (0x0000000000000001, 11, "49406564584", -324, 1, "4.9406564584e-324"),
    (0x0000000000000001, 12, "494065645841", -324, 1, "4.94065645841e-324"),
    (0x0000000000000001, 13, "4940656458412", -324, 1, "4.940656458412e-324"),
    (0x0000000000000001, 14, "49406564584125", -324, 1, "4.9406564584125e-324"),
    (0x0000000000000001, 15, "494065645841247", -324, 1, "4.94065645841247e-324"),
    (0x0000000000000001, 16, "494065645841247", -324, 1, "4.94065645841247e-324"),
    (0x0000000000000001, 17, "494065645841247", -324, 1, "4.94065645841247e-324"),
    (0x4023fffffffffffd, 0, "", 0, 0, ""),
    (0x4023fffffffffffd, 1, "1", 1, 1, "1e+01"),
    (0x4023fffffffffffd, 2, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 3, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 4, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 5, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 6, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 7, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 8, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 9, "999999999", 0, 1, "9.99999999"),
    (0x4023fffffffffffd, 10, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 11, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 12, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 13, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 14, "1", 1, 1, "10"),
    (0x4023fffffffffffd, 15, "999999999999999", 0, 1, "9.99999999999999"),
    (0x4023fffffffffffd, 16, "999999999999999", 0, 1, "9.99999999999999"),
    (0x4023fffffffffffd, 17, "999999999999999", 0, 1, "9.99999999999999"),
    (0x42dc12218377de40, 0, "", 0, 0, ""),
    (0x42dc12218377de40, 1, "1", 14, 1, "1e+14"),
    (0x42dc12218377de40, 2, "12", 14, 1, "1.2e+14"),
    (0x42dc12218377de40, 3, "123", 14, 1, "1.23e+14"),
    (0x42dc12218377de40, 4, "1235", 14, 1, "1.235e+14"),
    (0x42dc12218377de40, 5, "12346", 14, 1, "1.2346e+14"),
    (0x42dc12218377de40, 6, "123457", 14, 1, "1.23457e+14"),
    (0x42dc12218377de40, 7, "1234568", 14, 1, "1.234568e+14"),
    (0x42dc12218377de40, 8, "12345679", 14, 1, "1.2345679e+14"),
    (0x42dc12218377de40, 9, "123456789", 14, 1, "1.23456789e+14"),
    (0x42dc12218377de40, 10, "123456789", 14, 1, "1.23456789e+14"),
    (0x42dc12218377de40, 11, "12345678901", 14, 1, "1.2345678901e+14"),
    (0x42dc12218377de40, 12, "123456789012", 14, 1, "1.23456789012e+14"),
    (0x42dc12218377de40, 13, "1234567890123", 14, 1, "1.234567890123e+14"),
    (0x42dc12218377de40, 14, "12345678901235", 14, 1, "1.2345678901235e+14"),
    (0x42dc12218377de40, 15, "123456789012345", 14, 1, "123456789012345"),
    (0x42dc12218377de40, 16, "123456789012345", 14, 1, "123456789012345"),
    (0x42dc12218377de40, 17, "123456789012345", 14, 1, "123456789012345"),
    (0x43118b54f22aeb00, 0, "", 0, 0, ""),
    (0x43118b54f22aeb00, 1, "1", 15, 1, "1e+15"),
    (0x43118b54f22aeb00, 2, "12", 15, 1, "1.2e+15"),
    (0x43118b54f22aeb00, 3, "123", 15, 1, "1.23e+15"),
    (0x43118b54f22aeb00, 4, "1235", 15, 1, "1.235e+15"),
    (0x43118b54f22aeb00, 5, "12346", 15, 1, "1.2346e+15"),
    (0x43118b54f22aeb00, 6, "123457", 15, 1, "1.23457e+15"),
    (0x43118b54f22aeb00, 7, "1234568", 15, 1, "1.234568e+15"),
    (0x43118b54f22aeb00, 8, "12345679", 15, 1, "1.2345679e+15"),
    (0x43118b54f22aeb00, 9, "123456789", 15, 1, "1.23456789e+15"),
    (0x43118b54f22aeb00, 10, "123456789", 15, 1, "1.23456789e+15"),
    (0x43118b54f22aeb00, 11, "12345678901", 15, 1, "1.2345678901e+15"),
    (0x43118b54f22aeb00, 12, "123456789012", 15, 1, "1.23456789012e+15"),
    (0x43118b54f22aeb00, 13, "1234567890123", 15, 1, "1.234567890123e+15"),
    (0x43118b54f22aeb00, 14, "12345678901235", 15, 1, "1.2345678901235e+15"),
    (0x43118b54f22aeb00, 15, "123456789012346", 15, 1, "1.23456789012346e+15"),
    (0x43118b54f22aeb00, 16, "123456789012346", 15, 1, "1234567890123460"),
    (0x43118b54f22aeb00, 17, "123456789012346", 15, 1, "1234567890123460"),
    (0x42d6bcc41e8ffffd, 0, "", 0, 0, ""),
    (0x42d6bcc41e8ffffd, 1, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 2, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 3, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 4, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 5, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 6, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 7, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 8, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 9, "999999999", 13, 1, "9.99999999e+13"),
    (0x42d6bcc41e8ffffd, 10, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 11, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 12, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 13, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 14, "1", 14, 1, "1e+14"),
    (0x42d6bcc41e8ffffd, 15, "1", 14, 1, "100000000000000"),
    (0x42d6bcc41e8ffffd, 16, "1", 14, 1, "100000000000000"),
    (0x42d6bcc41e8ffffd, 17, "1", 14, 1, "100000000000000"),
    (0x4341c37937e08000, 0, "", 0, 0, ""),
    (0x4341c37937e08000, 1, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 2, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 3, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 4, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 5, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 6, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 7, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 8, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 9, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 10, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 11, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 12, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 13, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 14, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 15, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 16, "1", 16, 1, "1e+16"),
    (0x4341c37937e08000, 17, "1", 16, 1, "10000000000000000"),
    (0xbe8091b5aeffdb8e, 0, "", 0, 0, ""),
    (0xbe8091b5aeffdb8e, 1, "1", -7, -1, "-1e-07"),
    (0xbe8091b5aeffdb8e, 2, "12", -7, -1, "-1.2e-07"),
    (0xbe8091b5aeffdb8e, 3, "123", -7, -1, "-1.23e-07"),
    (0xbe8091b5aeffdb8e, 4, "1235", -7, -1, "-1.235e-07"),
    (0xbe8091b5aeffdb8e, 5, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 6, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 7, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 8, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 9, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 10, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 11, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 12, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 13, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 14, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 15, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 16, "12345", -7, -1, "-1.2345e-07"),
    (0xbe8091b5aeffdb8e, 17, "12345", -7, -1, "-1.2345e-07"),
    (0x3fd613716aefcc27, 0, "", 0, 0, ""),
    (0x3fd613716aefcc27, 1, "3", -1, 1, "0.3"),
    (0x3fd613716aefcc27, 2, "34", -1, 1, "0.34"),
    (0x3fd613716aefcc27, 3, "345", -1, 1, "0.345"),
    (0x3fd613716aefcc27, 4, "3449", -1, 1, "0.3449"),
    (0x3fd613716aefcc27, 5, "34494", -1, 1, "0.34494"),
    (0x3fd613716aefcc27, 6, "344937", -1, 1, "0.344937"),
    (0x3fd613716aefcc27, 7, "3449367", -1, 1, "0.3449367"),
    (0x3fd613716aefcc27, 8, "34493671", -1, 1, "0.34493671"),
    (0x3fd613716aefcc27, 9, "344936708", -1, 1, "0.344936708"),
    (0x3fd613716aefcc27, 10, "3449367089", -1, 1, "0.3449367089"),
    (0x3fd613716aefcc27, 11, "34493670886", -1, 1, "0.34493670886"),
    (0x3fd613716aefcc27, 12, "344936708861", -1, 1, "0.344936708861"),
    (0x3fd613716aefcc27, 13, "3449367088608", -1, 1, "0.3449367088608"),
    (0x3fd613716aefcc27, 14, "34493670886076", -1, 1, "0.34493670886076"),
    (0x3fd613716aefcc27, 15, "34493670886076", -1, 1, "0.34493670886076"),
    (0x3fd613716aefcc27, 16, "34493670886076", -1, 1, "0.34493670886076"),
    (0x3fd613716aefcc27, 17, "34493670886076", -1, 1, "0.34493670886076"),
    (0x406ff84473563d7b, 0, "", 0, 0, ""),
    (0x406ff84473563d7b, 1, "3", 2, 1, "3e+02"),
    (0x406ff84473563d7b, 2, "26", 2, 1, "2.6e+02"),
    (0x406ff84473563d7b, 3, "256", 2, 1, "256"),
    (0x406ff84473563d7b, 4, "2558", 2, 1, "255.8"),
    (0x406ff84473563d7b, 5, "25576", 2, 1, "255.76"),
    (0x406ff84473563d7b, 6, "255758", 2, 1, "255.758"),
    (0x406ff84473563d7b, 7, "2557584", 2, 1, "255.7584"),
    (0x406ff84473563d7b, 8, "25575836", 2, 1, "255.75836"),
    (0x406ff84473563d7b, 9, "255758355", 2, 1, "255.758355"),
    (0x406ff84473563d7b, 10, "2557583558", 2, 1, "255.7583558"),
    (0x406ff84473563d7b, 11, "25575835578", 2, 1, "255.75835578"),
    (0x406ff84473563d7b, 12, "255758355778", 2, 1, "255.758355778"),
    (0x406ff84473563d7b, 13, "2557583557782", 2, 1, "255.7583557782"),
    (0x406ff84473563d7b, 14, "25575835577816", 2, 1, "255.75835577816"),
    (0x406ff84473563d7b, 15, "255758355778158", 2, 1, "255.758355778158"),
    (0x406ff84473563d7b, 16, "255758355778158", 2, 1, "255.758355778158"),
    (0x406ff84473563d7b, 17, "255758355778158", 2, 1, "255.758355778158"),
    (0x430c6bf526340000, 0, "", 0, 0, ""),
    (0x430c6bf526340000, 1, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 2, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 3, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 4, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 5, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 6, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 7, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 8, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 9, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 10, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 11, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 12, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 13, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 14, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 15, "1", 15, 1, "1e+15"),
    (0x430c6bf526340000, 16, "1", 15, 1, "1000000000000000"),
    (0x430c6bf526340000, 17, "1", 15, 1, "1000000000000000"),
];
