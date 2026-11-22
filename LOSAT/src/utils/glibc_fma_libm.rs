// The algorithms, constants and tables of `exp`, `log` and `pow` here are the code of Arm
// Optimized Routines (math/exp.c, exp_data.c, log.c, log_data.c, pow.c, pow_log_data.c,
// math_config.h, math_err.c at commit 11253b0b9d6b), which glibc 2.35 also carries
// (sysdeps/ieee754/dbl-64; musl carries the same tables). The notices of those files
// (math_config.h has the first, the others the second):
//
// Copyright (c) 2017-2018, Arm Limited.
// SPDX-License-Identifier: MIT
//
// Copyright (c) 2018, Arm Limited.
// SPDX-License-Identifier: MIT
//
// The MIT License text is in THIRD_PARTY_NOTICES, section 4 (Arm Optimized Routines).
// `log10` follows sysdeps/ieee754/dbl-64/e_log10.c of glibc 2.35, which is Sun's fdlibm code
// and carries only this notice (THIRD_PARTY_NOTICES, section 5 (Sun's fdlibm)):
//
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.

//! No NCBI counterpart: the C library that NCBI BLAST+ calls. `exp`, `log`, `log10` and `pow` with
//! the bits of glibc 2.35 on x86-64 with FMA, where NCBI BLAST+ 2.17.0 wrote the frozen bytes of the
//! tests, on every target LOSAT builds for.
//!
//! On x86-64 glibc picks FMA variants of `exp`, `log` and `pow` at run time (`__ieee754_exp_fma`,
//! `__ieee754_log_fma`, `__ieee754_pow_fma`: sysdeps/ieee754/dbl-64/e_exp.c, e_log.c and e_pow.c
//! built with -mfma -mavx2); `log10` (e_log10.c, built without FMA) calls that `log`. Another libm
//! (glibc on a CPU without FMA, musl and so wasi-libc, Apple's, the Windows CRT) gives a result 1 ULP
//! apart on a small share of the arguments (exp 0.07 %, log 0.01 %), and a TBLASTX sum E value
//! printed with 15 digits (`-outfmt 8`/`11`) multiplies such a difference by ~50.
//!
//! `exp`, `log` and `pow` are the code of Arm Optimized Routines (MIT, see the notice above), which
//! glibc 2.35 carries, built to give glibc's bits; `log10` is Sun's fdlibm code as in glibc's
//! e_log10.c. Every fused multiply-add of the compiled glibc function (read from the disassembly of
//! the library, libm.so.6 of Ubuntu 2.35-0ubuntu3.18) is one call of an [`Fma`] method here, every
//! other operation is the same IEEE operation in the same order, and the tables are the same bits
//! as glibc's (`tables.rs`, read from that library). An FMA rounds once, so its result is the same
//! on every target; how it is computed is the only thing that differs:
//! - [`Fused`] (`f64::mul_add`): the FMA instruction on aarch64, in the x86-64 copies compiled with
//!   the `fma` feature (taken when the CPU has it, `x86_fma`), and on x86-64 builds that enable
//!   `fma`; elsewhere the C library's correctly rounded `fma`.
//! - [`Emulated`]: the same result from IEEE additions and multiplications only, for WebAssembly
//!   (no FMA instruction; wasi-libc's software `fma` costs ~16 ns per `exp` call) and x86-64 CPUs
//!   without FMA. Most of glibc's FMAs add a small product to a larger sum, where a cheap test
//!   proves the unfused result right; the others say what makes them cheaper (an exact product, an
//!   exact sum, an integer factor). The polynomials of `exp` and `log` run unfused, and one test of
//!   their last rounding, with a proved bound on the unfused minus the fused value, keeps the bits;
//!   the arguments that fail a test (about 1 in 50) take every FMA exactly.
//!
//! `utils::glibc_fma_libm` is the only place where LOSAT calls a libm function on `f64`:
//! `clippy.toml` refuses `f64::exp`, `f64::ln`, `f64::log10`, `f64::powf` and the other
//! transcendental methods everywhere else.
#![allow(
    clippy::unreadable_literal,
    clippy::eq_op,
    clippy::approx_constant,
    reason = "the upstream code: its bit patterns, its `(x - x) / (x - x)` for NaN and its digits of 1/ln 10"
)]

mod tables;

use tables::{EXP_HEAD, EXP_TAB, LOG_HEAD, LOG_TAB, POW_LOG_HEAD, POW_LOG_TAB};

#[inline(always)]
fn d(u: u64) -> f64 {
    f64::from_bits(u)
}

/// One fused multiply-add of glibc's code: `a * b + c` rounded once to nearest. Each variant
/// names what the operands of that FMA make true; [`Fused`] computes them all with one
/// instruction, [`Emulated`] uses what they say to do less work.
pub(crate) trait Fma: Sized {
    /// Finite operands.
    fn fma(a: f64, b: f64, c: f64) -> f64;
    /// `a * b` is a double: `b` is a power of 2, or `a` is an integer and `b` has the trailing zero
    /// bits for it.
    #[inline(always)]
    fn fma_exact_product(a: f64, b: f64, c: f64) -> f64 {
        Self::fma(a, b, c)
    }
    /// `RN(a * b) + c` is a double: `c` is `-RN(a * b)`, or `-1` with `a * b` in `[0.5, 2]`.
    #[inline(always)]
    fn fma_exact_sum(a: f64, b: f64, c: f64) -> f64 {
        Self::fma(a, b, c)
    }
    /// [`Fma::fma_exact_sum`] with `b` from a table and `b_split` its Veltkamp split
    /// ([`Emulated::split`], computed with the table: [`invc_split`]).
    #[inline(always)]
    fn fma_exact_sum_split(a: f64, b: f64, _b_split: (f64, f64), c: f64) -> f64 {
        Self::fma_exact_sum(a, b, c)
    }
    /// `a` is an integer below 2^26 in magnitude; the product's rounding error reaches the result.
    #[inline(always)]
    fn fma_small_int(a: f64, b: f64, c: f64) -> f64 {
        Self::fma(a, b, c)
    }
    /// `|a * b|` is near or above `|c|`, so the product's rounding error reaches the result.
    #[inline(always)]
    fn fma_product_led(a: f64, b: f64, c: f64) -> f64 {
        Self::fma(a, b, c)
    }
    /// `b` is `+-(2^m + 2^n)` with `n < m` (two bits); the product's rounding error reaches the
    /// result.
    #[inline(always)]
    fn fma_two_bit(a: f64, b: f64, c: f64) -> f64 {
        Self::fma(a, b, c)
    }
    /// `exp`'s last FMA, `scale * tmp + scale` with [`exp_tmp`], when the result is normal.
    #[inline(always)]
    fn exp_tail(r: f64, rt: f64, scale: f64) -> f64 {
        Self::fma(scale, exp_tmp::<Self>(r, rt), scale)
    }
    /// `log`'s result `y + hi` with `y` of [`ln_poly`] (outside the interval near 1).
    #[inline(always)]
    fn ln_tail(r: f64, kd: f64, w: f64, hi: f64) -> f64 {
        ln_poly::<Self>(r, kd, w, hi) + hi
    }
    /// `log`'s result near 1, `hi + y` with `y` of [`ln1_poly`].
    #[inline(always)]
    fn ln1_tail(r: f64, r2: f64, r3: f64, parts: [f64; 3], hi: f64) -> f64 {
        hi + ln1_poly::<Self>(r, r2, r3, parts)
    }
}

/// `f64::mul_add`.
#[cfg_attr(
    all(target_arch = "wasm32", not(test)),
    expect(dead_code, reason = "WebAssembly has no FMA instruction")
)]
pub(crate) struct Fused;

impl Fma for Fused {
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        a.mul_add(b, c)
    }
}

/// Each FMA as a rounded product and a rounded sum: not glibc's bits. Only [`Emulated`]'s
/// `exp_tail` and `ln_tail` use it, as an estimate whose distance from the fused value they bound.
struct Unfused;

impl Fma for Unfused {
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        a * b + c
    }
}

/// `a * b + c` rounded once, from IEEE additions and multiplications. The general way
/// ([`Emulated::full`]) is Boldo and Melquiond's ("Emulation of FMA and correctly rounded sums:
/// proved algorithms using rounding to odd", IEEE Trans. Computers 57(4), 2008, algorithm 5.4):
/// the exact product `p + e` (Dekker, with Veltkamp's split), the exact sum `c + p = sh + sl`
/// (Knuth's TwoSum), `v = sl + e` rounded to odd, and `sh + v` rounded to nearest. The variants of
/// [`Fma`] skip what their operands make unneeded, and [`Fma::fma`] first tries a test that
/// proves the two-rounding `c + RN(a * b)` right.
///
/// Exact for finite operands whose split and product do not overflow (`|a|`, `|b|` below
/// 2^995) and whose product's exponents sum to -969 or more (no bit of `e` below 2^-1074), or
/// that have a zero factor. Every call in this module meets that: the special inputs (zero,
/// subnormal results, infinity, NaN) leave the functions before their first `fma`, and the
/// smallest operand products are `tmp * scale` of `exp` (2^-200 * 2^-739, glibc's own bound) and
/// the powers of the reduced argument `r` (nonzero `|r|` >= 2^-110).
pub(crate) struct Emulated;

/// The bit of the sign and the bits of the exponent of an `f64`.
const SIGN_BIT: u64 = 0x8000_0000_0000_0000;
const EXP_BITS: u64 = 0x7ff0_0000_0000_0000;
/// 2^-53 + 2^-105: `|p| * E_BOUND` rounded to nearest is at least `|p| * 2^-53` when normal.
const E_BOUND: u64 = 0x3ca0_0000_0000_0001;
/// `exp_tail`: `|r|` at most 0.00271 (ln 2 / 256 = 0.0027076 plus `pow`'s tail, < 2^-41), and
/// twice the bound on `|tmp|` unfused minus fused there, 2^-59 (see `exp_tail`).
const EXP_R_MAX: f64 = 0.00271;
const EXP_TMP_SLACK: u64 = 0x3c40_0000_0000_0000; // 2^-59
/// `ln_tail`: `|r|` at most 0.0041 (glibc: `|r| < 1/256` and a little) and the half width of the
/// bracket around the unfused `y`, 2^-66: three times the bound on `|y|` unfused minus fused there
/// (see `ln_tail`).
const LN_R_MAX: f64 = 0.0041;
const LN_Y_SLACK: u64 = 0x3bd0_0000_0000_0000; // 2^-66

impl Emulated {
    /// Veltkamp's split: `a = hi + lo` with 26 and 27 significant bits.
    #[inline(always)]
    const fn split(a: f64) -> (f64, f64) {
        let t = 134217729.0 * a; // 2^27 + 1
        let hi = t - (t - a);
        (hi, a - hi)
    }

    /// `a * b - p` for `p = RN(a * b)`, exact (Dekker).
    #[inline(always)]
    fn product_error(a: f64, b: f64, p: f64) -> f64 {
        Self::product_error_split(a, Self::split(b), p)
    }

    /// [`Emulated::product_error`] with `b` given by its split.
    #[inline(always)]
    fn product_error_split(a: f64, (bh, bl): (f64, f64), p: f64) -> f64 {
        let (ah, al) = Self::split(a);
        (((ah * bh - p) + ah * bl) + al * bh) + al * bl
    }

    /// `c + p + e` rounded once to nearest, for `p + e` an exact product with `p` rounded to
    /// nearest (Boldo and Melquiond's algorithm 5.4 after the product).
    #[inline(always)]
    fn round_sum(c: f64, p: f64, e: f64) -> f64 {
        let sh = c + p;
        let cb = sh - p;
        let sl = (c - cb) + (p - (sh - cb));
        let v = sl + e;
        let vb = v - sl;
        let ve = (sl - (v - vb)) + (e - vb);
        // Round to odd: an inexact sum with an even last bit moves one step towards `ve`.
        let vbits = v.to_bits();
        let v = if ve != 0.0 && vbits & 1 == 0 {
            d(if (vbits ^ ve.to_bits()) >> 63 == 0 {
                vbits + 1
            } else {
                vbits - 1
            })
        } else {
            v
        };
        // `v` is zero only when `c + p + e` is `sh` exactly; `sh + 0.0` would lose the sign of -0.
        if v == 0.0 {
            sh
        } else {
            sh + v
        }
    }

    /// The general way: any operands within the bounds above.
    #[inline(always)]
    fn full(a: f64, b: f64, c: f64) -> f64 {
        let p = a * b;
        Self::round_sum(c, p, Self::product_error(a, b, p))
    }

    #[cold]
    #[inline(never)]
    fn full_cold(a: f64, b: f64, c: f64) -> f64 {
        Self::full(a, b, c)
    }

    /// Whether every real within `bound` of `u + z` rounds to `u`: `|z| + bound` below half the
    /// gap from `u` to its nearer neighbour (2^(E-53) for `u` in [2^E, 2^(E+1)), 2^(E-54) when `u`
    /// is 2^E). `RN(|z| + bound)` below that half gap (a power of 2) means `|z| + bound` is below
    /// it. False for `u` below 2^-968, infinite or NaN.
    #[inline(always)]
    fn rounds_to(u: f64, z: f64, bound: f64) -> bool {
        let ub = u.to_bits() & !SIGN_BIT;
        let eb = ub & EXP_BITS;
        let half_gap = d(eb.wrapping_sub(if ub == eb { 54 << 52 } else { 53 << 52 }));
        eb.wrapping_sub(55 << 52) < EXP_BITS - (55 << 52) && z.abs() + bound < half_gap
    }

    /// `u = c + RN(a * b)` when it is provably `c + a * b'` rounded once for every `a * b'` within
    /// `slack / 1.01` of `a * b` (`slack` 0: `b' = b`). With `u + z = c + p` exactly (TwoSum) and
    /// `a * b = p + e`, the exact sum is `u + z + e`; `|e|` is at most `|p| * 2^-53` for a normal
    /// `p` and 2^-1075 below, so at most `RN(|p| * E_BOUND) + 2^-1022`; the 1 % covers the
    /// rounding of the bound. Fails for `|p|` near `|u|` (about `2 |p| / |u|` of the arguments).
    #[inline(always)]
    fn unfused_if_exact(a: f64, b: f64, c: f64, slack: f64) -> Option<f64> {
        let p = a * b;
        let u = c + p;
        let pb = u - c;
        let z = (c - (u - pb)) + (p - pb);
        let bound = p.abs() * d(E_BOUND) + (slack + f64::MIN_POSITIVE);
        Self::rounds_to(u, z, bound).then_some(u)
    }

    /// `r`, which a variant computed for `a * b + c`; debug builds compare it with the general way.
    #[inline(always)]
    fn checked(r: f64, a: f64, b: f64, c: f64) -> f64 {
        debug_assert_eq!(
            r.to_bits(),
            Self::full(a, b, c).to_bits(),
            "{a:e} * {b:e} + {c:e}"
        );
        r
    }
}

impl Fma for Emulated {
    /// `c + RN(a * b)` when [`Emulated::unfused_if_exact`] proves it, else [`Emulated::full`].
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        match Self::unfused_if_exact(a, b, c, 0.0) {
            Some(u) => u,
            None => Self::full_cold(a, b, c),
        }
    }

    #[inline(always)]
    fn fma_exact_product(a: f64, b: f64, c: f64) -> f64 {
        let r = a * b + c;
        Self::checked(r, a, b, c)
    }

    /// `RN(a * b) + c` exact plus the exact rounding error of the product: one rounding.
    #[inline(always)]
    fn fma_exact_sum(a: f64, b: f64, c: f64) -> f64 {
        Self::fma_exact_sum_split(a, b, Self::split(b), c)
    }

    /// The same, with the split of `b` read from its table instead of computed.
    #[inline(always)]
    fn fma_exact_sum_split(a: f64, b: f64, b_split: (f64, f64), c: f64) -> f64 {
        debug_assert_eq!(b_split, Self::split(b));
        let p = a * b;
        let r = (p + c) + Self::product_error_split(a, b_split, p);
        Self::checked(r, a, b, c)
    }

    /// Dekker's product with `a` as its own high part (an integer below 2^26 splits into itself
    /// and 0), so only `b` is split; for a constant `b` the compiler folds the split.
    #[inline(always)]
    fn fma_small_int(a: f64, b: f64, c: f64) -> f64 {
        let p = a * b;
        let (bh, bl) = Self::split(b);
        let r = Self::round_sum(c, p, (a * bh - p) + a * bl);
        Self::checked(r, a, b, c)
    }

    #[inline(always)]
    fn fma_product_led(a: f64, b: f64, c: f64) -> f64 {
        Self::full(a, b, c)
    }

    /// `a * 2^m` and `a * 2^n` are exact and `a * 2^m - p` too (Sterbenz), so the product's
    /// rounding error needs no split.
    #[inline(always)]
    fn fma_two_bit(a: f64, b: f64, c: f64) -> f64 {
        let p = a * b;
        let bh = d(b.to_bits() & 0xfff0_0000_0000_0000);
        let r = Self::round_sum(c, p, (a * bh - p) + a * (b - bh));
        Self::checked(r, a, b, c)
    }

    /// The polynomial unfused and one test of the last FMA. For `|r| <= 0.00271`, unfused minus
    /// fused (`u - f`, each step `|RN(x) - RN(y)| <= |x - y| + ulp`): `C2 + r C3` in
    /// [0.4995, 0.5005] differs by at most 2^-53 + |r C3| 2^-53 = 1.00046 * 2^-53, `C4 + r C5`
    /// in [2^-5, 2^-4) by 1.0004 * 2^-57; `rt + p23 r2` (`|r2| <= 7.35e-6`, result below 2^-8)
    /// by 2^-61 + (1.00046 + 0.5005) * 7.35e-6 * 2^-53 = 1.0029 * 2^-61; `tmp` (below 2^-8) by
    /// that, `r2^2 (2^-57 + p45 2^-53)` (< 1e-9 * 2^-61) and 2^-61: below 2.003 * 2^-61. So
    /// `scale * tmp + scale` is within `|scale| * 2.003 * 2^-61` of glibc's exact last sum (`pow`'s
    /// `scale` is negative for an odd power of a negative `x`), and
    /// `EXP_TMP_SLACK` (2^-59) is twice that. Fails for about 2^-6 of the arguments, which take
    /// every FMA exactly.
    #[inline(always)]
    fn exp_tail(r: f64, rt: f64, scale: f64) -> f64 {
        if r.abs() <= EXP_R_MAX {
            let tmp = exp_tmp::<Unfused>(r, rt);
            if let Some(v) =
                Self::unfused_if_exact(scale, tmp, scale, scale.abs() * d(EXP_TMP_SLACK))
            {
                debug_assert_eq!(
                    v.to_bits(),
                    Self::fma(scale, exp_tmp::<Self>(r, rt), scale).to_bits()
                );
                return v;
            }
        }
        Self::fma(scale, exp_tmp::<Self>(r, rt), scale)
    }

    /// `y` unfused and one test of `y + hi`. For `|r| <= 0.0041`, `|kd| <= 1075`, `|hi| < 1024`,
    /// unfused minus fused, steps as in `exp_tail`: `lo` (below 2^-33; `|kd Ln2lo| <= 5.91e-11`)
    /// by 2^-86 + 5.91e-11 2^-53 = 1.51 2^-86; `r2 A0 + lo` (below 2^-16; `|r2 A0| <= 8.41e-6`)
    /// by 2^-69 + 8.41e-6 2^-53 + that = 1.551 2^-69; `A1 + r A2` in [0.332, 0.335) by 1.0021
    /// 2^-54 and `A3 + r A4` in [0.199, 0.201) by 1.003 2^-55, so `p` in [0.332, 0.335) by
    /// 2.0021 2^-54; `y` (below 2^-16; `|r3| <= 6.9e-8`) by 2^-69 + 1.551 2^-69 + 6.9e-8 (2.0021
    /// 2^-54 + 0.335 2^-53) = 2.557 2^-69.
    ///
    /// The test brackets the fused `y`: with `s` = `LN_Y_SLACK` (2^-66), `RN(y - s)` is at most
    /// `y - s` plus half an ulp of `|y| + s` (below 2^-16 + 2^-66, so at most 2^-69), so at most
    /// `y - 7 2^-69`, below the fused `y` (at least `y - 2.557 2^-69`); `RN(y + s)` is above it the
    /// same way. Rounding to nearest is monotone, so glibc's `RN(y + hi)` lies between
    /// `RN(RN(y - s) + hi)` and `RN(RN(y + s) + hi)`, and is both when they are equal. Fails when
    /// a rounding boundary of `y + hi` lies within about `s` of the unfused sum: at most 2^-9 of
    /// the arguments (the result is at least 0.0625, an ulp at least 2^-56), which take every FMA
    /// exactly.
    #[inline(always)]
    fn ln_tail(r: f64, kd: f64, w: f64, hi: f64) -> f64 {
        if r.abs() <= LN_R_MAX {
            let y = ln_poly::<Unfused>(r, kd, w, hi);
            let below = (y - d(LN_Y_SLACK)) + hi;
            let above = (y + d(LN_Y_SLACK)) + hi;
            if below.to_bits() == above.to_bits() {
                debug_assert_eq!(
                    below.to_bits(),
                    (ln_poly::<Self>(r, kd, w, hi) + hi).to_bits()
                );
                return below;
            }
        }
        ln_poly::<Self>(r, kd, w, hi) + hi
    }

    /// `y` unfused and one test of `hi + y`, with a bound that scales with `r = x - 1` (exact;
    /// `|r| <= 0.0648`). Unfused minus fused, steps as in `exp_tail`: `B1 + r B2` in [0.317, 0.35)
    /// by 1.0324 * 2^-54, `q1` by 2.0341 * 2^-54; `q2` (in [0.157, 0.177]) by 2.0391 * 2^-55;
    /// `q3` (in [0.104, 0.118]) by 3.0551 * 2^-56; `q3 r3 + q2` by 3.0397 * 2^-55; `p` (in
    /// [0.317, 0.351)) by 3.035 * 2^-54. `lo` (below 1.04 |r| 2^-53) plus `B0 rlo (r + rhi)`
    /// (below `r^2 2^-26`) by `r^2 2^-79` and the ulp of the sum; `y` by that, `|r3| (3.035 *
    /// 2^-54 + 0.3506 * 2^-53)` and its ulp: below 2.569 |r3| 2^-53 + 1.25 r^2 2^-77 + 2.08 |r|
    /// 2^-105. The slack is about twice each term.
    #[inline(always)]
    fn ln1_tail(r: f64, r2: f64, r3: f64, parts: [f64; 3], hi: f64) -> f64 {
        let y = ln1_poly::<Unfused>(r, r2, r3, parts);
        let u = hi + y;
        let yb = u - hi;
        let z = (hi - (u - yb)) + (y - yb);
        let slack = r3.abs() * d(0x3cc0_0000_0000_0000) // 2^-51
            + (r2 * d(0x3b30_0000_0000_0000) + r.abs() * d(0x3990_0000_0000_0000)); // 2^-76, 2^-102
        if Self::rounds_to(u, z, slack) {
            debug_assert_eq!(
                u.to_bits(),
                (hi + ln1_poly::<Self>(r, r2, r3, parts)).to_bits()
            );
            return u;
        }
        hi + ln1_poly::<Self>(r, r2, r3, parts)
    }
}

/// [`Emulated::split`] of every `invc` of a table of `N` doubles in rows of `STRIDE`, `invc` first
/// (`LOG_TAB`, `POW_LOG_TAB`), evaluated by the compiler, for [`Fma::fma_exact_sum_split`].
const fn invc_split<const N: usize, const STRIDE: usize, const ROWS: usize>(
    tab: &[u64; N],
) -> [(f64, f64); ROWS] {
    assert!(N == STRIDE * ROWS);
    let mut out = [(0.0, 0.0); ROWS];
    let mut i = 0;
    while i < ROWS {
        out[i] = Emulated::split(f64::from_bits(tab[STRIDE * i]));
        i += 1;
    }
    out
}

static LOG_INVC_SPLIT: [(f64, f64); 128] = invc_split::<256, 2, 128>(&LOG_TAB);
static POW_LOG_INVC_SPLIT: [(f64, f64); 128] = invc_split::<384, 3, 128>(&POW_LOG_TAB);

#[inline(always)]
fn top12(x: f64) -> u32 {
    (x.to_bits() >> 52) as u32
}

/// `__math_oflow` / `__math_uflow` of math_err.c (`xflow`): `+-0x1p769 * 0x1p769` or `+-0x1p-767 * 0x1p-767`.
#[cold]
fn xflow(sign: u64, y: f64) -> f64 {
    let s = if sign != 0 { -y } else { y };
    s * y
}

fn math_oflow(sign: u64) -> f64 {
    xflow(sign, d(0x7000_0000_0000_0000)) // 0x1p769
}

fn math_uflow(sign: u64) -> f64 {
    xflow(sign, d(0x1000_0000_0000_0000)) // 0x1p-767
}

const EXP_TABLE_BITS: u32 = 7;
const EXP_N_MASK: u64 = (1 << EXP_TABLE_BITS) - 1;
/// `SIGN_BIAS` of pow.c: added to `ki` it sets the sign bit of the scale.
const SIGN_BIAS: u64 = 0x800 << EXP_TABLE_BITS;

/// `specialcase` of pow.c (the one of exp.c is the same for its positive results): the result
/// `scale * (1 + tmp)` near overflow (`k > 0`) or in the subnormal range (`k < 0`).
#[cold]
#[inline(never)]
fn exp_specialcase<F: Fma>(tmp: f64, sbits: u64, ki: u64) -> f64 {
    if ki & 0x8000_0000 == 0 {
        // k > 0: the exponent of scale might have overflowed by <= 460.
        let scale = d(sbits.wrapping_sub(1009u64 << 52));
        return d(0x7f00_0000_0000_0000) * F::fma(scale, tmp, scale); // 0x1p1009
    }
    // k < 0: round in the normal range before scaling into the subnormal range.
    let sbits = sbits.wrapping_add(1022u64 << 52);
    let scale = d(sbits);
    let st = scale * tmp;
    let mut y = scale + st;
    if y.abs() < 1.0 {
        let one = if y < 0.0 { -1.0 } else { 1.0 };
        let lo = (scale - y) + st;
        let hi = one + y;
        let lo = ((one - hi) + y) + lo;
        y = (hi + lo) - one;
        if y == 0.0 {
            y = d(sbits & 0x8000_0000_0000_0000);
        }
    }
    y * d(0x0010_0000_0000_0000) // 0x1p-1022
}

/// The body of exp.c's `exp` and of pow.c's `exp_inline` once the input is in range:
/// `sign * exp(x + xtail)` (`TAIL`: `pow` adds its `xtail`); `special` when `|x| >= 512`.
#[inline(always)]
fn exp_core<F: Fma, const TAIL: bool>(x: f64, xtail: f64, sign_bias: u64, special: bool) -> f64 {
    let kd = F::fma(x, d(EXP_HEAD[0]), d(EXP_HEAD[1])); // fma(x, InvLn2N, Shift)
    let ki = kd.to_bits();
    let kd = kd - d(EXP_HEAD[1]);
    // kd * NegLn2hiN is a double for |kd| <= 2^17 (NegLn2hiN has 36 significant bits): |x| < 512.
    let t1 = if special {
        F::fma(kd, d(EXP_HEAD[2]), x)
    } else {
        F::fma_exact_product(kd, d(EXP_HEAD[2]), x)
    };
    let mut r = F::fma(kd, d(EXP_HEAD[3]), t1);
    if TAIL {
        r += xtail;
    }
    let idx = 2 * (ki & EXP_N_MASK) as usize;
    let top = ki.wrapping_add(sign_bias) << (52 - EXP_TABLE_BITS);
    let tail = d(EXP_TAB[idx]);
    let sbits = EXP_TAB[idx + 1].wrapping_add(top);
    let rt = r + tail;
    if special {
        return exp_specialcase::<F>(exp_tmp::<F>(r, rt), sbits, ki);
    }
    F::exp_tail(r, rt, d(sbits))
}

/// `tmp` of exp.c and pow.c: `tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5)`, with
/// glibc's three FMAs and `rt = r + tail`.
#[inline(always)]
fn exp_tmp<F: Fma>(r: f64, rt: f64) -> f64 {
    let r2 = r * r;
    let p23 = F::fma(r, d(EXP_HEAD[5]), d(EXP_HEAD[4])); // C2 + r * C3
    let p45 = F::fma(r, d(EXP_HEAD[7]), d(EXP_HEAD[6])); // C4 + r * C5
    F::fma(r2 * r2, p45, F::fma(p23, r2, rt))
}

/// `exp` of exp.c, with glibc's bits.
#[inline(always)]
fn exp_impl<F: Fma>(x: f64) -> f64 {
    let abstop = top12(x) & 0x7ff;
    // top12(0x1p-54) = 0x3c9, top12(512.0) = 0x408, top12(1024.0) = 0x409
    if abstop.wrapping_sub(0x3c9) >= 0x408 - 0x3c9 {
        if abstop.wrapping_sub(0x3c9) >= 0x8000_0000 {
            return 1.0 + x; // tiny x
        }
        if abstop >= 0x409 {
            if x.to_bits() == f64::NEG_INFINITY.to_bits() {
                return 0.0;
            }
            if abstop >= 0x7ff {
                return 1.0 + x;
            }
            return if x.to_bits() >> 63 != 0 {
                math_uflow(0)
            } else {
                math_oflow(0)
            };
        }
        return exp_core::<F, false>(x, 0.0, 0, true); // large |x|
    }
    exp_core::<F, false>(x, 0.0, 0, false)
}

/// `log` of log.c, with glibc's bits.
#[inline(always)]
fn ln_impl<F: Fma>(x: f64) -> f64 {
    let mut ix = x.to_bits();
    let top = (ix >> 48) as u32;
    const LO: u64 = 0x3fee_0000_0000_0000; // 1.0 - 0x1p-4
    const HI: u64 = 0x3ff1_0900_0000_0000; // 1.0 + 0x1.09p-4
    if ix.wrapping_sub(LO) < HI - LO {
        if ix == 0x3ff0_0000_0000_0000 {
            return 0.0;
        }
        let b0 = d(LOG_HEAD[7]); // B[0] = -0.5
        let r = x - 1.0;
        let r2 = r * r;
        let r3 = r * r2;
        let c27 = d(0x41a0_0000_0000_0000); // 0x1p27
        let t = F::fma_exact_product(r, c27, r);
        let rhi = F::fma_exact_product(-c27, r, t);
        let rlo = r - rhi;
        let rr = rhi * rhi;
        let hi = F::fma_exact_product(rr, b0, r);
        let lo = F::fma_exact_product(rr, b0, r - hi);
        return F::ln1_tail(r, r2, r3, [rhi, rlo, lo], hi);
    }
    if top.wrapping_sub(0x0010) >= 0x7ff0 - 0x0010 {
        if ix << 1 == 0 {
            return -1.0 / 0.0_f64; // __math_divzero(1)
        }
        if ix == f64::INFINITY.to_bits() {
            return x;
        }
        if (top & 0x8000) != 0 || (top & 0x7ff0) == 0x7ff0 {
            return (x - x) / (x - x); // __math_invalid
        }
        // subnormal x: normalize it
        ix = (x * d(0x4330_0000_0000_0000))
            .to_bits()
            .wrapping_sub(52u64 << 52); // x * 0x1p52
    }
    const OFF: u64 = 0x3fe6_0000_0000_0000;
    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> 45) % 128) as usize;
    let k = (tmp as i64) >> 52;
    let iz = ix.wrapping_sub(tmp & (0xfffu64 << 52));
    let invc = d(LOG_TAB[2 * i]);
    let logc = d(LOG_TAB[2 * i + 1]);
    let z = d(iz);
    // z * invc is within 2^-7 of 1; kd * Ln2hi is a double (|k| < 2^11, Ln2hi has 42 bits).
    let r = F::fma_exact_sum_split(z, invc, LOG_INVC_SPLIT[i], -1.0);
    let kd = k as f64;
    let w = F::fma_exact_product(kd, d(LOG_HEAD[0]), logc);
    let hi = r + w;
    F::ln_tail(r, kd, w, hi)
}

/// `y` of log.c near 1, whose result is `hi + y`: `lo + B[0] * rlo * (r + rhi) + r3 * p` with
/// `p` log.c's polynomial in `r` (`B[1]` to `B[10]`), with its FMAs; `parts` is `[rhi, rlo, lo]`.
#[inline(always)]
fn ln1_poly<F: Fma>(r: f64, r2: f64, r3: f64, parts: [f64; 3]) -> f64 {
    let [rhi, rlo, lo] = parts;
    let b = |i: usize| d(LOG_HEAD[7 + i]);
    let q1 = F::fma(r2, b(3), F::fma(r, b(2), b(1)));
    let q2 = F::fma(r2, b(6), F::fma(r, b(5), b(4)));
    let q3 = F::fma(r3, b(10), F::fma(r2, b(9), F::fma(r, b(8), b(7))));
    let p = F::fma(F::fma(q3, r3, q2), r3, q1);
    let lo = F::fma_product_led(b(0) * rlo, r + rhi, lo);
    F::fma_product_led(p, r3, lo)
}

/// `y` of log.c, whose result is `y + hi`: `lo + r2 * A[0] + r * r2 * (A[1] + r * A[2] + r2 *
/// (A[3] + r * A[4]))` with `lo = w - hi + r + kd * Ln2lo`, with glibc's FMAs.
#[inline(always)]
fn ln_poly<F: Fma>(r: f64, kd: f64, w: f64, hi: f64) -> f64 {
    let lo = F::fma_small_int(kd, d(LOG_HEAD[1]), (w - hi) + r);
    let r2 = r * r;
    let a = |j: usize| d(LOG_HEAD[2 + j]);
    let p = F::fma(F::fma(r, a(4), a(3)), r2, F::fma(r, a(2), a(1)));
    F::fma(r * r2, p, F::fma_two_bit(r2, a(0), lo)) // A[0] = -0x1.0000000000001p-1
}

/// glibc's `log10` (e_log10.c `__ieee754_log10`, built without FMA; its `__ieee754_log` is the
/// FMA `log` above).
#[inline(always)]
fn log10_impl<F: Fma>(x: f64) -> f64 {
    const TWO54: f64 = 1.80143985094819840000e+16; // 0x4350000000000000
    const IVLN10: f64 = 4.34294481903251816668e-01; // 0x3FDBCB7B1526E50E
    const LOG10_2HI: f64 = 3.01029995663611771306e-01; // 0x3FD34413509F6000
    const LOG10_2LO: f64 = 3.69423907715893078616e-13; // 0x3D59FEF311F12B36
    let mut x = x;
    let mut hx = x.to_bits() as i64;
    let mut k: i32 = 0;
    if hx < 0x0010_0000_0000_0000 {
        if hx & 0x7fff_ffff_ffff_ffff == 0 {
            return -TWO54 / x.abs(); // log(+-0) = -inf
        }
        if hx < 0 {
            return (x - x) / (x - x); // log(-#) = NaN
        }
        k -= 54;
        x *= TWO54; // subnormal number, scale up x
        hx = x.to_bits() as i64;
    }
    if hx >= 0x7ff0_0000_0000_0000 {
        return x + x;
    }
    k += ((hx >> 52) - 1023) as i32;
    let i = i64::from(k < 0);
    let hx = (hx & 0x000f_ffff_ffff_ffff) | ((0x3ff - i) << 52);
    let y = f64::from(k + i as i32);
    let z = y * LOG10_2LO + IVLN10 * ln_impl::<F>(d(hx as u64));
    z + y * LOG10_2HI
}

/// `checkint` of pow.c: 0 if `iy` is not an integer, 1 if odd, 2 if even.
#[inline(always)]
fn checkint(iy: u64) -> u32 {
    let e = (iy >> 52 & 0x7ff) as u32;
    if e < 0x3ff {
        return 0;
    }
    if e > 0x3ff + 52 {
        return 2;
    }
    if iy & ((1u64 << (0x3ff + 52 - e)) - 1) != 0 {
        return 0;
    }
    if iy & (1u64 << (0x3ff + 52 - e)) != 0 {
        return 1;
    }
    2
}

/// `zeroinfnan` of pow.c: the bits of 0, infinity or NaN.
#[inline(always)]
fn zeroinfnan(i: u64) -> bool {
    (i << 1).wrapping_sub(1) >= (f64::INFINITY.to_bits() << 1) - 1
}

/// `issignaling_inline` of math_config.h.
fn issignaling(x: f64) -> bool {
    (x.to_bits() ^ 0x0008_0000_0000_0000) << 1 > 0x7ff8_0000_0000_0000 << 1
}

/// `log_inline` of pow.c: `log(x)` as `y + tail` with ~15 more bits in `tail`.
#[inline(always)]
fn pow_log_inline<F: Fma>(ix: u64) -> (f64, f64) {
    const OFF: u64 = 0x3fe6_9555_0000_0000;
    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> 45) % 128) as usize;
    let k = (tmp as i64) >> 52;
    let iz = ix.wrapping_sub(tmp & (0xfffu64 << 52));
    let z = d(iz);
    let kd = k as f64;
    let invc = d(POW_LOG_TAB[3 * i]);
    let logc = d(POW_LOG_TAB[3 * i + 1]);
    let logctail = d(POW_LOG_TAB[3 * i + 2]);
    let a = |j: usize| d(POW_LOG_HEAD[2 + j]);
    let r = F::fma_exact_sum_split(z, invc, POW_LOG_INVC_SPLIT[i], -1.0);
    // k*Ln2 + log(c) + r.
    let t1 = F::fma_exact_product(kd, d(POW_LOG_HEAD[0]), logc);
    let t2 = t1 + r;
    let lo1 = F::fma_small_int(kd, d(POW_LOG_HEAD[1]), logctail);
    let lo2 = (t1 - t2) + r;
    let ar = a(0) * r; // A[0] = -0.5
    let ar2 = r * ar;
    let ar3 = r * ar2;
    // k*Ln2 + log(c) + r + A[0]*r*r.
    let hi = t2 + ar2;
    let lo3 = F::fma_exact_sum(ar, r, -ar2);
    let lo4 = (t2 - hi) + ar2;
    // p = log1p(r) - r - A[0]*r*r, added to the low parts by the last fma.
    let q = F::fma(F::fma(r, a(6), a(5)), ar2, F::fma(r, a(4), a(3)));
    let p = F::fma(ar2, q, F::fma(r, a(2), a(1)));
    let lo = F::fma_product_led(ar3, p, ((lo1 + lo2) + lo3) + lo4);
    let y = hi + lo;
    (y, (hi - y) + lo)
}

/// The special inputs of pow.c's `pow`: the result, or the bits of `|x|` (normalized when
/// subnormal) and the sign bias for the general path.
#[cold]
#[inline(never)]
fn pow_special(x: f64, y: f64, topx: u32, topy: u32) -> Result<f64, (u64, u64)> {
    const ONE: u64 = 0x3ff0_0000_0000_0000;
    const INF: u64 = 0x7ff0_0000_0000_0000;
    let mut ix = x.to_bits();
    let iy = y.to_bits();
    let mut topx = topx;
    let mut sign_bias = 0;
    if zeroinfnan(iy) {
        if iy << 1 == 0 {
            return Ok(if issignaling(x) { x + y } else { 1.0 });
        }
        if ix == ONE {
            return Ok(if issignaling(y) { x + y } else { 1.0 });
        }
        if ix << 1 > INF << 1 || iy << 1 > INF << 1 {
            return Ok(x + y);
        }
        if ix << 1 == ONE << 1 {
            return Ok(1.0);
        }
        if (ix << 1 < ONE << 1) == (iy >> 63 == 0) {
            return Ok(0.0); // |x|<1 && y==inf or |x|>1 && y==-inf.
        }
        return Ok(y * y);
    }
    if zeroinfnan(ix) {
        let mut x2 = x * x;
        if ix >> 63 != 0 && checkint(iy) == 1 {
            x2 = -x2;
            sign_bias = 1;
        }
        if ix << 1 == 0 && iy >> 63 != 0 {
            return Ok(if sign_bias != 0 { -1.0 } else { 1.0 } / 0.0); // __math_divzero
        }
        return Ok(if iy >> 63 != 0 { 1.0 / x2 } else { x2 });
    }
    // Here x and y are non-zero finite.
    if ix >> 63 != 0 {
        // Finite x < 0.
        let yint = checkint(iy);
        if yint == 0 {
            return Ok((x - x) / (x - x)); // __math_invalid
        }
        if yint == 1 {
            sign_bias = SIGN_BIAS;
        }
        ix &= 0x7fff_ffff_ffff_ffff;
        topx &= 0x7ff;
    }
    if (topy & 0x7ff).wrapping_sub(0x3be) >= 0x43e - 0x3be {
        // Note: sign_bias == 0 here because y is not odd.
        if ix == ONE {
            return Ok(1.0);
        }
        if (topy & 0x7ff) < 0x3be {
            // |y| < 2^-65, x^y ~= 1 + y*log(x).
            return Ok(if ix > ONE { 1.0 + y } else { 1.0 - y });
        }
        return Ok(if (ix > ONE) == (topy < 0x800) {
            math_oflow(0)
        } else {
            math_uflow(0)
        });
    }
    if topx == 0 {
        // Normalize subnormal x so exponent becomes negative.
        ix = (x * d(0x4330_0000_0000_0000)).to_bits(); // x * 0x1p52
        ix &= 0x7fff_ffff_ffff_ffff;
        ix = ix.wrapping_sub(52u64 << 52);
    }
    Err((ix, sign_bias))
}

/// `pow` of pow.c, with glibc's bits.
#[inline(always)]
fn pow_impl<F: Fma>(x: f64, y: f64) -> f64 {
    let mut ix = x.to_bits();
    let topx = top12(x);
    let topy = top12(y);
    let mut sign_bias = 0;
    if topx.wrapping_sub(0x001) >= 0x7ff - 0x001
        || (topy & 0x7ff).wrapping_sub(0x3be) >= 0x43e - 0x3be
    {
        match pow_special(x, y, topx, topy) {
            Ok(v) => return v,
            Err((jx, bias)) => {
                ix = jx;
                sign_bias = bias;
            }
        }
    }
    let (hi, lo) = pow_log_inline::<F>(ix);
    let ehi = y * hi;
    let elo = F::fma_product_led(y, lo, F::fma_exact_sum(hi, y, -ehi));
    // exp_inline (ehi, elo, sign_bias)
    let abstop = top12(ehi) & 0x7ff;
    if abstop.wrapping_sub(0x3c9) >= 0x408 - 0x3c9 {
        if abstop.wrapping_sub(0x3c9) >= 0x8000_0000 {
            // Avoid spurious underflow for tiny x.
            let one = 1.0 + ehi;
            return if sign_bias != 0 { -one } else { one };
        }
        if abstop >= 0x409 {
            return if ehi.to_bits() >> 63 != 0 {
                math_uflow(sign_bias)
            } else {
                math_oflow(sign_bias)
            };
        }
        return exp_core::<F, true>(ehi, elo, sign_bias, true); // large |x|
    }
    exp_core::<F, true>(ehi, elo, sign_bias, false)
}

/// The library's [`Fma`] where the CPU is not asked: [`Emulated`] on WebAssembly, [`Fused`] on
/// aarch64 and on x86-64 builds with the `fma` feature. x86-64 builds without it ask the CPU
/// (`x86_fma`).
#[cfg(not(all(target_arch = "x86_64", not(target_feature = "fma"))))]
#[cfg(target_arch = "wasm32")]
type Portable = Emulated;
#[cfg(not(all(target_arch = "x86_64", not(target_feature = "fma"))))]
#[cfg(not(target_arch = "wasm32"))]
type Portable = Fused;

/// Defines the entry points `$name`, each `$imp::<Portable>`. Where [`Portable`] is [`Fused`] they
/// are never inlined, so their FMA instructions stay in this module, the only one besides
/// `stats::spouge` where `LOSAT/tests/check_fp_contract.py binary` allows them; on WebAssembly
/// they hold none and the compiler decides.
#[cfg(not(all(target_arch = "x86_64", not(target_feature = "fma"))))]
macro_rules! entry_points {
    ($(fn $name:ident($($arg:ident),+) = $imp:ident;)+) => {$(
        #[cfg_attr(not(target_arch = "wasm32"), inline(never))]
        fn $name($($arg: f64),+) -> f64 {
            $imp::<Portable>($($arg),+)
        }
    )+};
}

#[cfg(not(all(target_arch = "x86_64", not(target_feature = "fma"))))]
entry_points! {
    fn exp_portable(x) = exp_impl;
    fn ln_portable(x) = ln_impl;
    fn log10_portable(x) = log10_impl;
    fn pow_portable(x, y) = pow_impl;
}

/// x86-64 builds without the `fma` feature: each entry point is chosen on its first call, a copy
/// compiled with `fma` ([`Fused`]) when the CPU has FMA, else [`Emulated`] (the same bits), and
/// then called through a pointer, as glibc's ifunc does. A direct call would let a caller compiled
/// with a feature that implies `fma` (the `avx512f` helper of `x_newton_exact.rs`) inline the
/// copy and take its FMA instructions out of this module, where
/// `LOSAT/tests/check_fp_contract.py binary` allows them: rustc 1.92 drops `#[inline(never)]` on a
/// `#[target_feature]` function.
#[cfg(all(target_arch = "x86_64", not(target_feature = "fma")))]
mod x86_fma {
    use super::{exp_impl, ln_impl, log10_impl, pow_impl, Emulated, Fused};
    use std::sync::atomic::{AtomicPtr, Ordering};

    macro_rules! ifunc {
        (@f64 $arg:ident) => { f64 };
        ($(fn $name:ident($($arg:ident),+) = $imp:ident;)+) => {$(
            pub(super) fn $name($($arg: f64),+) -> f64 {
                type Entry = unsafe fn($(ifunc!(@f64 $arg)),+) -> f64;
                #[target_feature(enable = "fma")]
                unsafe fn fused($($arg: f64),+) -> f64 {
                    $imp::<Fused>($($arg),+)
                }
                fn emulated($($arg: f64),+) -> f64 {
                    $imp::<Emulated>($($arg),+)
                }
                fn resolve($($arg: f64),+) -> f64 {
                    let entry: Entry = if std::arch::is_x86_feature_detected!("fma") {
                        fused
                    } else {
                        emulated
                    };
                    SLOT.store(entry as *mut (), Ordering::Relaxed);
                    // SAFETY: `fused` is chosen only when the CPU has FMA.
                    unsafe { entry($($arg),+) }
                }
                static SLOT: AtomicPtr<()> = AtomicPtr::new(resolve as Entry as *mut ());
                // SAFETY: `SLOT` holds `resolve`, `emulated` or `fused` (stored only when the CPU
                // has FMA), each an `Entry`.
                unsafe {
                    let entry = std::mem::transmute::<*mut (), Entry>(SLOT.load(Ordering::Relaxed));
                    entry($($arg),+)
                }
            }
        )+};
    }

    ifunc! {
        fn exp(x) = exp_impl;
        fn ln(x) = ln_impl;
        fn log10(x) = log10_impl;
        fn pow(x, y) = pow_impl;
    }
}
#[cfg(all(target_arch = "x86_64", not(target_feature = "fma")))]
use x86_fma::{
    exp as exp_portable, ln as ln_portable, log10 as log10_portable, pow as pow_portable,
};

/// glibc's `exp` on x86-64 with FMA, bit for bit.
#[inline]
pub fn exp(x: f64) -> f64 {
    exp_portable(x)
}

/// glibc's `log` on x86-64 with FMA, bit for bit.
#[inline]
pub fn ln(x: f64) -> f64 {
    ln_portable(x)
}

/// glibc's `log10` on x86-64 with FMA, bit for bit.
#[inline]
pub fn log10(x: f64) -> f64 {
    log10_portable(x)
}

/// glibc's `pow` on x86-64 with FMA, bit for bit.
#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    pow_portable(x, y)
}

/// The functions above as methods: `x.glibc_exp()` in place of `x.exp()`.
pub trait GlibcLibm {
    fn glibc_exp(self) -> f64;
    fn glibc_ln(self) -> f64;
    fn glibc_log10(self) -> f64;
    fn glibc_pow(self, y: f64) -> f64;
}

impl GlibcLibm for f64 {
    #[inline]
    fn glibc_exp(self) -> f64 {
        exp(self)
    }
    #[inline]
    fn glibc_ln(self) -> f64 {
        ln(self)
    }
    #[inline]
    fn glibc_log10(self) -> f64 {
        log10(self)
    }
    #[inline]
    fn glibc_pow(self, y: f64) -> f64 {
        pow(self, y)
    }
}

/// The two ways to compute each function, for the tests that compare them.
#[cfg(test)]
#[path = "glibc_fma_libm_variants_tests.rs"]
pub(crate) mod variants;

#[cfg(test)]
#[path = "glibc_fma_libm_tests.rs"]
mod tests;
