//! The x87 80-bit extended-precision (`long double`) operations of NCBI's serial real-number
//! digit generator, in integer arithmetic.
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
//!
//! The oracle (NCBI BLAST+ 2.17.0 for x86-64 Linux, built with gcc) evaluates every `long double`
//! expression of `NStr::DoubleToString_Ecvt` on the x87 unit: a sign, a 15-bit exponent and a
//! 64-bit significand with an explicit integer bit, precision control "64 bits" (the Linux
//! control word 0x037F) and round to nearest, ties to even. LOSAT must print the same digits on
//! every platform, including those whose `long double` is the 64-bit double (MSVC, wasm32) or a
//! 128-bit quad (aarch64 Linux), so the operations are modelled here with integers.
//!
//! Only the operations that `DoubleToString_Ecvt` performs exist: exact conversion from `f64`
//! (`fld` of a double) and from `u32` (`t1 = first`), multiplication rounded to 64 bits
//! (`value *= ...`, `value * mult`, `(t1-first) * 1e8`), the exact subtraction `t1 - first`,
//! comparisons and the truncating conversion to `unsigned int`. All operands of the port are
//! finite and nonnegative, and every product stays between 2^-1100 and 2^1100, far inside the
//! x87 normal range (2^-16382 .. 2^16384), so no denormal, overflow or NaN case arises.

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
