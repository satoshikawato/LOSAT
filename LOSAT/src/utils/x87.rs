//! The x87 80-bit extended-precision (`long double`) operations that NCBI's x86-64 build
//! evaluates, in integer arithmetic: those of the serial real-number digit generator
//! (`report::serial::real`) and of the line `result *= expl(-z * z) / z;` of ErfC
//! (`stats::spouge`).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2479-2481
//! ```c++
//!     unsigned int first, second=0;
//!     long double mult = __NLG(1.);
//!     long double value = val;
//! ```
//! NCBI reference (598d8ae6): c++/include/corelib/ncbitype.h:243
//! ```c++
//! #  define NCBI_CONST_LONGDOUBLE(v)   v##L
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
//! ```c
//!             result *= expl(-z * z) / z;
//! ```
//!
//! The oracle (NCBI BLAST+ 2.17.0 for x86-64 Linux, built with gcc) evaluates every `long double`
//! expression of `NStr::DoubleToString_Ecvt` and of `ErfImpl` on the x87 unit: a sign, a 15-bit
//! exponent and a 64-bit significand with an explicit integer bit, precision control "64 bits"
//! (the Linux control word 0x037F) and round to nearest, ties to even. LOSAT must print the same
//! digits on every platform, including those whose `long double` is the 64-bit double (MSVC,
//! wasm32) or a 128-bit quad (aarch64 Linux), so the operations are modelled here with integers.
//!
//! Only the operations that the two ports perform exist. `DoubleToString_Ecvt`: exact conversion
//! from `f64` (`fld` of a double) and from `u32` (`t1 = first`), multiplication rounded to 64 bits
//! (`value *= ...`, `value * mult`, `(t1-first) * 1e8`), the exact subtraction `t1 - first`,
//! comparisons and the truncating conversion to `unsigned int`. All its operands are finite and
//! nonnegative, and every product stays between 2^-1100 and 2^1100, far inside the x87 normal
//! range (2^-16382 .. 2^16384), so no denormal, overflow or NaN case arises. `ErfImpl`: the
//! 64-bit `expl` result as the value nearest a double-double exponential, division and
//! multiplication rounded to 64 bits (`/ z`, `result *= ...`) and the store to a double
//! (`fstpl`), which rounds the 64-bit product a second time, to 53 bits. Its operands are finite
//! and nonnegative (LOSAT's double-double stand-in for `expl` is zero where it underflows; NCBI's
//! 80-bit `expl` is still normal there), and the store can give a subnormal double or zero.

use std::cmp::Ordering;

/// An x87 extended-precision value: `(-1)^negative * significand * 2^(exponent - 63)`.
///
/// A nonzero value always has bit 63 of `significand` set (the explicit integer bit of the
/// 80-bit format). Zero has `significand == 0`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct X87 {
    negative: bool,
    /// Unbiased binary exponent of the integer bit.
    exponent: i32,
    significand: u64,
}

const INTEGER_BIT: u64 = 1 << 63;

impl X87 {
    pub(crate) const ZERO: X87 = X87 {
        negative: false,
        exponent: 0,
        significand: 0,
    };

    /// A positive normal value from its 64-bit significand (integer bit set) and unbiased
    /// exponent, as the 80-bit memory format stores them.
    pub(crate) const fn from_parts(significand: u64, exponent: i32) -> X87 {
        assert!(
            significand & INTEGER_BIT != 0,
            "x87 significand is not normalized"
        );
        X87 {
            negative: false,
            exponent,
            significand,
        }
    }

    /// `(negative, unbiased exponent, significand)`.
    #[cfg(test)]
    pub(crate) const fn parts(self) -> (bool, i32, u64) {
        (self.negative, self.exponent, self.significand)
    }

    /// `long double value = val;` (`fldl`): exact, subnormal doubles included, because the
    /// 80-bit format has 11 more significand bits and a wider exponent range than a double.
    pub(crate) fn from_f64(x: f64) -> X87 {
        debug_assert!(x.is_finite());
        let bits = x.to_bits();
        let negative = (bits >> 63) != 0;
        let biased = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1u64 << 52) - 1);
        if biased == 0 {
            if fraction == 0 {
                return X87 {
                    negative,
                    ..X87::ZERO
                };
            }
            // A subnormal double: fraction * 2^-1074, normalized here.
            let lz = fraction.leading_zeros() as i32;
            return X87 {
                negative,
                exponent: -1011 - lz,
                significand: fraction << lz,
            };
        }
        X87 {
            negative,
            exponent: biased - 1023,
            significand: INTEGER_BIT | (fraction << 11),
        }
    }

    /// `t1 = first;` (`fild` of an unsigned int): exact.
    pub(crate) fn from_u32(v: u32) -> X87 {
        if v == 0 {
            return X87::ZERO;
        }
        let lz = v.leading_zeros() as i32;
        X87 {
            negative: false,
            exponent: 31 - lz,
            significand: (v as u64) << (32 + lz),
        }
    }

    /// `fmul` with precision control 64 bits and round to nearest, ties to even: the exact
    /// 128-bit product of the significands rounded once to 64 bits.
    pub(crate) fn mul(self, rhs: X87) -> X87 {
        let negative = self.negative != rhs.negative;
        if self.significand == 0 || rhs.significand == 0 {
            return X87 {
                negative,
                ..X87::ZERO
            };
        }
        // Both significands lie in [2^63, 2^64), so the product lies in [2^126, 2^128).
        let product = (self.significand as u128) * (rhs.significand as u128);
        let mut exponent = self.exponent + rhs.exponent;
        let (mut significand, rest, half) = if (product >> 127) != 0 {
            exponent += 1;
            ((product >> 64) as u64, product as u64, 1u64 << 63)
        } else {
            (
                (product >> 63) as u64,
                (product as u64) & ((1u64 << 63) - 1),
                1u64 << 62,
            )
        };
        if rest > half || (rest == half && (significand & 1) != 0) {
            significand = significand.wrapping_add(1);
            if significand == 0 {
                // 0xFFFF_FFFF_FFFF_FFFF rounded up to 2^64.
                significand = INTEGER_BIT;
                exponent += 1;
            }
        }
        debug_assert!((-16382..=16383).contains(&exponent));
        X87 {
            negative,
            exponent,
            significand,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
    // ```c
    //             result *= expl(-z * z) / z;
    // ```
    // No NCBI counterpart: the x87 `fdiv` of `expl(...) / z`.
    /// `fdiv` with precision control 64 bits and round to nearest, ties to even: the exact
    /// quotient of the significands rounded once to 64 bits. `rhs` is nonzero. (A quotient of
    /// two 64-bit significands is never exactly halfway, so the ties-to-even terms never decide.)
    pub(crate) fn div(self, rhs: X87) -> X87 {
        debug_assert!(rhs.significand != 0, "x87 division by zero");
        let negative = self.negative != rhs.negative;
        if self.significand == 0 {
            return X87 {
                negative,
                ..X87::ZERO
            };
        }
        // The numerator lies in [2^127, 2^128) and the divisor in [2^63, 2^64), so the quotient
        // lies in (2^63, 2^65).
        let divisor = rhs.significand as u128;
        let numerator = (self.significand as u128) << 64;
        let (quotient, remainder) = (numerator / divisor, numerator % divisor);
        let exponent = self.exponent - rhs.exponent;
        let (exponent, significand, up) = if (quotient >> 64) != 0 {
            // The dropped part is the quotient's lowest bit and the remainder below it.
            let significand = (quotient >> 1) as u64;
            let up = (quotient & 1) != 0 && (remainder != 0 || (significand & 1) != 0);
            (exponent, significand, up)
        } else {
            // The dropped part is remainder / divisor, compared with one half.
            let significand = quotient as u64;
            let twice = remainder << 1;
            let up = twice > divisor || (twice == divisor && (significand & 1) != 0);
            (exponent - 1, significand, up)
        };
        X87::round_up_if(negative, exponent, significand, up)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
    // ```c
    //             result *= expl(-z * z) / z;
    // ```
    // No NCBI counterpart: the 64-bit result of `expl` from LOSAT's double-double stand-in.
    /// The value nearest the exact sum `hi + lo` (round to nearest, ties to even) of a
    /// normalized double-double: `hi` finite and nonnegative, and `hi + lo == hi` (the sum of
    /// the double-double's own operations). This is the 64-bit `expl(-z * z)` that the
    /// double-double exponential of `stats::spouge` holds to more than 64 bits.
    pub(crate) fn nearest_to_sum(hi: f64, lo: f64) -> X87 {
        debug_assert!(hi.is_finite() && hi >= 0.0 && hi + lo == hi);
        let h = X87::from_f64(hi);
        if h.significand == 0 {
            return X87::ZERO;
        }
        // hi = (h.significand << 63) * 2^base exactly, in [2^126, 2^127) units of 2^base; |lo|
        // is below 2^75 such units.
        let base = h.exponent - 126;
        let mut sum = (h.significand as u128) << 63;
        // Set when bits of lo below 2^base are dropped: the exact sum then lies strictly between
        // `sum` and `sum + 1` units. (For a double lo that happens only when |lo| is below 2^53
        // units, far from the halfway point, so it never decides a result; it keeps the rounding
        // exact for any input.)
        let mut sticky = false;
        if lo != 0.0 {
            let l = X87::from_f64(lo);
            // |lo| = l.significand * 2^shift units.
            let shift = l.exponent - 63 - base;
            let v = l.significand as u128;
            // floor(|lo|) in units, and whether a fraction of a unit remains.
            let (part, inexact) = if shift >= 0 {
                (v << shift, false)
            } else if shift > -128 {
                let s = (-shift) as u32;
                (v >> s, (v >> s) << s != v)
            } else {
                (0, true)
            };
            if l.negative {
                sum -= part + inexact as u128;
            } else {
                sum += part;
            }
            sticky = inexact;
        }
        // The sum is 126 to 128 bits wide; its top 64 bits are kept.
        let drop = 64 - sum.leading_zeros() as i32;
        let significand = (sum >> drop) as u64;
        let rest = sum & ((1u128 << drop) - 1);
        let half = 1u128 << (drop - 1);
        let up = rest > half || (rest == half && (sticky || (significand & 1) != 0));
        X87::round_up_if(false, base + drop + 63, significand, up)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
    // ```c
    //             result *= expl(-z * z) / z;
    // ```
    // No NCBI counterpart: the x87 `fstpl` of `result *= ...` to the double `result`.
    /// `fstpl` (the store to a double) with round to nearest, ties to even: the nearest double,
    /// subnormal results and the underflow to zero included. Values from 2^1024 up give
    /// infinity.
    pub(crate) fn to_f64(self) -> f64 {
        let sign = (self.negative as u64) << 63;
        if self.significand == 0 {
            return f64::from_bits(sign);
        }
        if self.exponent > 1023 {
            return f64::from_bits(sign | f64::INFINITY.to_bits());
        }
        // Significand bits below the last place of the double: 11 for a normal result, more for
        // a subnormal one.
        let (drop, biased) = if self.exponent >= -1022 {
            (11u32, (self.exponent + 1023) as u64)
        } else {
            ((-1011 - self.exponent) as u32, 0)
        };
        if drop > 64 {
            // Below half the smallest subnormal (2^-1075).
            return f64::from_bits(sign);
        }
        let significand = self.significand as u128;
        let kept = (significand >> drop) as u64;
        let rest = significand & ((1u128 << drop) - 1);
        let half = 1u128 << (drop - 1);
        let up = rest > half || (rest == half && (kept & 1) != 0);
        // A normal result drops its integer bit for the biased exponent. Rounding up carries
        // into the exponent field: past the largest subnormal to the smallest normal, past the
        // largest finite value to infinity.
        let magnitude = if biased == 0 {
            kept
        } else {
            (biased << 52) | (kept & ((1u64 << 52) - 1))
        };
        f64::from_bits(sign | (magnitude + up as u64))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
    // ```c
    //             result *= expl(-z * z) / z;
    // ```
    // No NCBI counterpart: the rounding step of `fdiv`, `fmul`'s result and `expl`'s result.
    /// `significand` (bit 63 set) plus one unit in the last place when `up`; the carry out of
    /// 0xFFFF_FFFF_FFFF_FFFF gives 2^63 and increments the exponent.
    fn round_up_if(negative: bool, exponent: i32, significand: u64, up: bool) -> X87 {
        let (significand, exponent) = match significand.checked_add(up as u64) {
            Some(significand) => (significand, exponent),
            None => (INTEGER_BIT, exponent + 1),
        };
        X87 {
            negative,
            exponent,
            significand,
        }
    }

    /// `t1 - first` (`fsub` of an unsigned int converted exactly), for `0 <= first <= t1` and
    /// `1 <= t1 < 2^63`. The difference is always exact, so rounding never applies:
    /// `t1 = m * 2^(e-63)` with `0 <= e <= 62`, and `first` is an integer, hence a multiple of
    /// `2^(e-63) <= 1/2`; the difference is `d * 2^(e-63)` with the integer
    /// `d = m - first * 2^(63-e)`, and `0 <= d <= m < 2^64` fits the 64-bit significand.
    /// (`DoubleToString_Ecvt` subtracts `first = trunc(t1)` for t1 in [1e8, 1e9),
    /// `first = 999999999 <= t1` for t1 >= 1e9, and `first == t1` after `t1 = first`.)
    pub(crate) fn sub_u32_exact(self, first: u32) -> X87 {
        if first == 0 {
            return self;
        }
        debug_assert!(!self.negative && (0..=62).contains(&self.exponent));
        let shift = (63 - self.exponent) as u32;
        let minuend = self.significand as u128;
        let subtrahend = (first as u128) << shift;
        debug_assert!(subtrahend <= minuend);
        let d = (minuend - subtrahend) as u64;
        if d == 0 {
            return X87::ZERO;
        }
        let lz = d.leading_zeros() as i32;
        X87 {
            negative: false,
            exponent: self.exponent - lz,
            significand: d << lz,
        }
    }

    /// `(unsigned int)x` for `0 <= x < 2^32`: gcc stores the value as a 64-bit integer with
    /// truncation toward zero (`fistp` under a chopping control word) and keeps the low 32 bits.
    pub(crate) fn to_u32_trunc(self) -> u32 {
        debug_assert!(!self.negative || self.significand == 0);
        if self.significand == 0 || self.exponent < 0 {
            return 0;
        }
        debug_assert!(self.exponent <= 31);
        (self.significand >> (63 - self.exponent.min(63))) as u32
    }

    /// The ordering of the represented values (`fcomp`; no NaN exists here, and `-0 == +0`).
    fn compare(self, other: X87) -> Ordering {
        let a_zero = self.significand == 0;
        let b_zero = other.significand == 0;
        if a_zero && b_zero {
            return Ordering::Equal;
        }
        let a_neg = self.negative && !a_zero;
        let b_neg = other.negative && !b_zero;
        if a_neg != b_neg {
            return if a_neg {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
        let magnitude = if a_zero {
            Ordering::Less
        } else if b_zero {
            Ordering::Greater
        } else {
            (self.exponent, self.significand).cmp(&(other.exponent, other.significand))
        };
        if a_neg {
            magnitude.reverse()
        } else {
            magnitude
        }
    }

    /// `a >= b`
    pub(crate) fn ge(self, other: X87) -> bool {
        self.compare(other) != Ordering::Less
    }

    /// `a < b`
    pub(crate) fn lt(self, other: X87) -> bool {
        self.compare(other) == Ordering::Less
    }

    /// `a <= b`
    pub(crate) fn le(self, other: X87) -> bool {
        self.compare(other) != Ordering::Greater
    }
}

#[cfg(test)]
mod tests {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143,171,199,229
    // ```c
    //             result *= expl(-z * z) / z;
    // ```
    // The expected values are the results of the x87 unit itself (gcc `long double` on x86-64
    // Linux, control word 0x037F): `fmul`, `fdiv`, `fadd` of two doubles and `fstpl`.
    use super::X87;

    fn ld(exponent: i32, significand: u64) -> X87 {
        X87::from_parts(significand, exponent)
    }

    #[test]
    fn mul_rounds_ties_to_even_and_carries() {
        // 1.5 * (1 + 2^-63): the product lies halfway; the kept significand is odd, so up.
        assert_eq!(
            ld(0, 0xC000_0000_0000_0000)
                .mul(ld(0, 0x8000_0000_0000_0001))
                .parts(),
            (false, 0, 0xC000_0000_0000_0002)
        );
        // 1.5 * (1 + 3 * 2^-63): halfway with an even kept significand, so it stays.
        assert_eq!(
            ld(0, 0xC000_0000_0000_0000)
                .mul(ld(0, 0x8000_0000_0000_0003))
                .parts(),
            (false, 0, 0xC000_0000_0000_0004)
        );
        // The product rounds up to 2^64, which carries into the exponent.
        assert_eq!(
            ld(0, 0xFFFF_FFFF_FFFF_FFFE)
                .mul(ld(0, 0x8000_0000_0000_0001))
                .parts(),
            (false, 1, 0x8000_0000_0000_0000)
        );
    }

    #[test]
    fn div_rounds_to_nearest() {
        // 1 / 1.5 rounds up, 1 / 7 down, 3 / 1.75 up (quotient below 2^64 of the significands).
        assert_eq!(
            ld(0, 0x8000_0000_0000_0000)
                .div(ld(0, 0xC000_0000_0000_0000))
                .parts(),
            (false, -1, 0xAAAA_AAAA_AAAA_AAAB)
        );
        assert_eq!(
            ld(0, 0x8000_0000_0000_0000)
                .div(ld(2, 0xE000_0000_0000_0000))
                .parts(),
            (false, -3, 0x9249_2492_4924_9249)
        );
        assert_eq!(
            ld(1, 0xC000_0000_0000_0000)
                .div(ld(0, 0xE000_0000_0000_0000))
                .parts(),
            (false, 0, 0xDB6D_B6DB_6DB6_DB6E)
        );
        // Quotients of 2^64 and above of the significands drop one more bit.
        assert_eq!(
            ld(0, 0xFFFF_FFFF_FFFF_FFFF)
                .div(ld(0, 0xAAAA_AAAA_AAAA_AAAB))
                .parts(),
            (false, 0, 0xBFFF_FFFF_FFFF_FFFF)
        );
        assert_eq!(
            ld(0, 0xFFFF_FFFF_FFFF_FFFF)
                .div(ld(0, 0xC000_0000_0000_0001))
                .parts(),
            (false, 0, 0xAAAA_AAAA_AAAA_AAA9)
        );
        // -36 / 3 is exact.
        assert_eq!(
            X87::from_f64(-36.0).div(X87::from_f64(3.0)).parts(),
            (true, 3, 0xC000_0000_0000_0000)
        );
        assert_eq!(X87::ZERO.div(X87::from_f64(3.0)).parts(), (false, 0, 0));
    }

    #[test]
    fn nearest_to_sum_rounds_ties_to_even_with_lower_bits_and_borrows() {
        let sum = |hi: f64, lo: f64| X87::nearest_to_sum(hi, lo).parts();
        let p = |e: i32| 2f64.powi(e);
        // 1 + 2^-64: halfway, even kept significand.
        assert_eq!(sum(1.0, p(-64)), (false, 0, 0x8000_0000_0000_0000));
        // 1 + 2^-63 + 2^-64 and (1 + 2^-52) + 2^-63 + 2^-64: halfway, odd, so up.
        assert_eq!(sum(1.0, 1.5 * p(-63)), (false, 0, 0x8000_0000_0000_0002));
        assert_eq!(
            sum(1.0 + p(-52), 1.5 * p(-63)),
            (false, 0, 0x8000_0000_0000_0802)
        );
        // Bits of lo far below the last place (above 2^-126, exact) break the tie.
        assert_eq!(
            sum(1.0, p(-64) + p(-100)),
            (false, 0, 0x8000_0000_0000_0001)
        );
        // A negative lo below a power of two: halfway (to even = 1), exact, below halfway.
        assert_eq!(sum(1.0, -p(-65)), (false, 0, 0x8000_0000_0000_0000));
        assert_eq!(sum(1.0, -p(-64)), (false, -1, 0xFFFF_FFFF_FFFF_FFFF));
        assert_eq!(
            sum(1.0, -(p(-65) + p(-117))),
            (false, -1, 0xFFFF_FFFF_FFFF_FFFF)
        );
        assert_eq!(
            sum(1.0 + p(-52), -(p(-64) + p(-110))),
            (false, 0, 0x8000_0000_0000_07FF)
        );
        // A subnormal hi, and zero (the exponential below the double range).
        assert_eq!(
            sum(f64::from_bits(0x10), 0.0),
            (false, -1070, 0x8000_0000_0000_0000)
        );
        assert_eq!(sum(0.0, 0.0), (false, 0, 0));
    }

    #[test]
    fn to_f64_rounds_ties_to_even_carries_and_gives_subnormals() {
        let store = |exponent: i32, significand: u64| ld(exponent, significand).to_f64().to_bits();
        assert_eq!(store(0, 0x8000_0000_0000_0400), 0x3FF0_0000_0000_0000);
        assert_eq!(store(0, 0x8000_0000_0000_0C00), 0x3FF0_0000_0000_0002);
        assert_eq!(store(0, 0x8000_0000_0000_0401), 0x3FF0_0000_0000_0001);
        assert_eq!(store(0, 0xFFFF_FFFF_FFFF_FC00), 0x4000_0000_0000_0000);
        assert_eq!(store(1023, 0xFFFF_FFFF_FFFF_FC00), 0x7FF0_0000_0000_0000);
        assert_eq!(store(-1074, 0x8000_0000_0000_0000), 0x0000_0000_0000_0001);
        assert_eq!(store(-1075, 0x8000_0000_0000_0000), 0x0000_0000_0000_0000);
        assert_eq!(store(-1074, 0xC000_0000_0000_0000), 0x0000_0000_0000_0002);
        assert_eq!(store(-1075, 0x8000_0000_0000_0001), 0x0000_0000_0000_0001);
        assert_eq!(store(-1023, 0xFFFF_FFFF_FFFF_FFFF), 0x0010_0000_0000_0000);
        assert_eq!(store(-1023, 0xFFFF_FFFF_FFFF_E800), 0x000F_FFFF_FFFF_FFFE);
        assert_eq!(store(-1100, 0x8000_0000_0000_0000), 0x0000_0000_0000_0000);
        assert_eq!(
            X87::from_f64(-1.0)
                .mul(ld(0, 0xC000_0000_0000_0400))
                .to_f64()
                .to_bits(),
            0xBFF8_0000_0000_0000
        );
        assert_eq!(X87::from_f64(-0.0).to_f64().to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn f64_round_trips_exactly() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for _ in 0..100_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            // Every finite double, subnormals included (the top exponent field, 0x7FF, is skipped).
            let bits = (state & !(0x7FF << 52)) | (((state >> 52) % 0x7FF) << 52);
            let x = f64::from_bits(bits);
            assert_eq!(X87::from_f64(x).to_f64().to_bits(), bits, "{x:e}");
            if x > 0.0 {
                assert_eq!(
                    X87::nearest_to_sum(x, 0.0).parts(),
                    X87::from_f64(x).parts()
                );
            }
        }
    }
}
