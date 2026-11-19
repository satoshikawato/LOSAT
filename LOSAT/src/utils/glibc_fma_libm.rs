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
//! Each function below is glibc's code: every fused multiply-add of the compiled glibc function
//! (read from the disassembly of the library, libm.so.6 of Ubuntu 2.35-0ubuntu3.18) is one
//! [`Fma::fma`] here, every other operation is the same IEEE operation in the same order, and the
//! tables are glibc's (`tables.rs`). An FMA rounds once, so its result is the same on every target;
//! how it is computed is the only thing that differs:
//! - [`Fused`] (`f64::mul_add`): the FMA instruction on aarch64, in the x86-64 copies compiled with
//!   the `fma` feature (taken when the CPU has it, `x86_fma`), and on x86-64 builds that enable
//!   `fma`; elsewhere the C library's correctly rounded `fma`.
//! - [`Emulated`]: the same result from IEEE additions and multiplications only, for WebAssembly
//!   (no FMA instruction; wasi-libc's software `fma` costs ~16 ns per `exp` call) and x86-64 CPUs
//!   without FMA.
//!
//! `utils::glibc_fma_libm` is the only place where LOSAT calls a libm function on `f64`:
//! `clippy.toml` refuses `f64::exp`, `f64::ln`, `f64::log10`, `f64::powf` and the other
//! transcendental methods everywhere else.
#![allow(
    clippy::unreadable_literal,
    clippy::eq_op,
    clippy::approx_constant,
    reason = "glibc's code: its bit patterns, its `(x - x) / (x - x)` for NaN and its digits of 1/ln 10"
)]

mod tables;

use tables::{EXP_HEAD, EXP_TAB, LOG_HEAD, LOG_TAB, POW_LOG_HEAD, POW_LOG_TAB};

#[inline(always)]
fn d(u: u64) -> f64 {
    f64::from_bits(u)
}

/// One fused multiply-add: `a * b + c` rounded once to nearest.
pub(crate) trait Fma {
    fn fma(a: f64, b: f64, c: f64) -> f64;
}

/// `f64::mul_add`.
pub(crate) struct Fused;

impl Fma for Fused {
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        a.mul_add(b, c)
    }
}

/// `a * b + c` rounded once, from IEEE additions and multiplications (Boldo and Melquiond,
/// "Emulation of FMA and correctly rounded sums: proved algorithms using rounding to odd", IEEE
/// Trans. Computers 57(4), 2008, algorithm 5.4): the exact product `ph + pl` (Dekker, with
/// Veltkamp's split), the exact sum `c + ph = sh + sl` (Knuth's TwoSum), `v = sl + pl` rounded to
/// odd, and `sh + v` rounded to nearest.
///
/// Exact for finite operands whose split and product do not overflow (`|a|`, `|b|` below
/// 2^995) and whose product's exponents sum to -969 or more (no bit of `pl` below 2^-1074), or
/// that have a zero factor. Every call in this module meets that: the special inputs (zero,
/// subnormal results, infinity, NaN) leave the functions before their first `fma`, and the
/// smallest operand products are `tmp * scale` of `exp` (2^-200 * 2^-739, glibc's own bound) and
/// the powers of the reduced argument `r` (nonzero `|r|` >= 2^-110).
pub(crate) struct Emulated;

impl Emulated {
    /// Veltkamp's split: `a = hi + lo` with 26 and 27 significant bits.
    #[inline(always)]
    fn split(a: f64) -> (f64, f64) {
        let t = 134217729.0 * a; // 2^27 + 1
        let hi = t - (t - a);
        (hi, a - hi)
    }
}

impl Fma for Emulated {
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        let ph = a * b;
        let (ah, al) = Self::split(a);
        let (bh, bl) = Self::split(b);
        let pl = (((ah * bh - ph) + ah * bl) + al * bh) + al * bl;
        let sh = c + ph;
        let cb = sh - ph;
        let sl = (c - cb) + (ph - (sh - cb));
        let v = sl + pl;
        let vb = v - sl;
        let ve = (sl - (v - vb)) + (pl - vb);
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
        // `v` is zero only when `c + a * b` is `sh` exactly; `sh + 0.0` would lose the sign of -0.
        if v == 0.0 {
            sh
        } else {
            sh + v
        }
    }
}

#[inline(always)]
fn top12(x: f64) -> u32 {
    (x.to_bits() >> 52) as u32
}

/// glibc's `__math_oflow` / `__math_uflow` (`xflow`): `+-0x1p769 * 0x1p769` or `+-0x1p-767 * 0x1p-767`.
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
/// `SIGN_BIAS` of e_pow.c: added to `ki` it sets the sign bit of the scale.
const SIGN_BIAS: u64 = 0x800 << EXP_TABLE_BITS;

/// `specialcase` of e_pow.c (the one of e_exp.c is the same for its positive results): the result
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

/// The body of e_exp.c's `__exp` and of e_pow.c's `exp_inline` once the input is in range:
/// `sign * exp(x + xtail)` (`TAIL`: `pow` adds its `xtail`); `special` when `|x| >= 512`.
#[inline(always)]
fn exp_core<F: Fma, const TAIL: bool>(x: f64, xtail: f64, sign_bias: u64, special: bool) -> f64 {
    let kd = F::fma(x, d(EXP_HEAD[0]), d(EXP_HEAD[1])); // fma(x, InvLn2N, Shift)
    let ki = kd.to_bits();
    let kd = kd - d(EXP_HEAD[1]);
    let mut r = F::fma(kd, d(EXP_HEAD[3]), F::fma(kd, d(EXP_HEAD[2]), x));
    if TAIL {
        r += xtail;
    }
    let idx = 2 * (ki & EXP_N_MASK) as usize;
    let top = ki.wrapping_add(sign_bias) << (52 - EXP_TABLE_BITS);
    let tail = d(EXP_TAB[idx]);
    let sbits = EXP_TAB[idx + 1].wrapping_add(top);
    let r2 = r * r;
    let p23 = F::fma(r, d(EXP_HEAD[5]), d(EXP_HEAD[4])); // C2 + r * C3
    let p45 = F::fma(r, d(EXP_HEAD[7]), d(EXP_HEAD[6])); // C4 + r * C5
    let tmp = F::fma(r2 * r2, p45, F::fma(p23, r2, r + tail));
    if special {
        return exp_specialcase::<F>(tmp, sbits, ki);
    }
    let scale = d(sbits);
    F::fma(scale, tmp, scale)
}

/// glibc's `exp` (e_exp.c `__exp`).
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

/// glibc's `log` (e_log.c `__log`).
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
        let b = |i: usize| d(LOG_HEAD[7 + i]);
        let r = x - 1.0;
        let r2 = r * r;
        let r3 = r * r2;
        let q1 = F::fma(r2, b(3), F::fma(r, b(2), b(1)));
        let q2 = F::fma(r2, b(6), F::fma(r, b(5), b(4)));
        let q3 = F::fma(r3, b(10), F::fma(r2, b(9), F::fma(r, b(8), b(7))));
        let p = F::fma(F::fma(q3, r3, q2), r3, q1);
        let c27 = d(0x41a0_0000_0000_0000); // 0x1p27
        let t = F::fma(r, c27, r);
        let rhi = F::fma(-c27, r, t);
        let rlo = r - rhi;
        let rr = rhi * rhi;
        let hi = F::fma(rr, b(0), r);
        let lo = F::fma(rr, b(0), r - hi);
        let lo = F::fma(b(0) * rlo, r + rhi, lo);
        let y = F::fma(p, r3, lo);
        return hi + y;
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
    let r = F::fma(z, invc, -1.0);
    let kd = k as f64;
    let w = F::fma(kd, d(LOG_HEAD[0]), logc);
    let hi = r + w;
    let lo = F::fma(kd, d(LOG_HEAD[1]), (w - hi) + r);
    let r2 = r * r;
    let a = |j: usize| d(LOG_HEAD[2 + j]);
    let p = F::fma(F::fma(r, a(4), a(3)), r2, F::fma(r, a(2), a(1)));
    let y = F::fma(r * r2, p, F::fma(r2, a(0), lo));
    y + hi
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

/// `checkint` of e_pow.c: 0 if `iy` is not an integer, 1 if odd, 2 if even.
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

/// `zeroinfnan` of e_pow.c: the bits of 0, infinity or NaN.
#[inline(always)]
fn zeroinfnan(i: u64) -> bool {
    (i << 1).wrapping_sub(1) >= (f64::INFINITY.to_bits() << 1) - 1
}

/// `issignaling_inline` of glibc's math_config.h.
fn issignaling(x: f64) -> bool {
    (x.to_bits() ^ 0x0008_0000_0000_0000) << 1 > 0x7ff8_0000_0000_0000 << 1
}

/// `log_inline` of e_pow.c: `log(x)` as `y + tail` with ~15 more bits in `tail`.
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
    let r = F::fma(z, invc, -1.0);
    // k*Ln2 + log(c) + r.
    let t1 = F::fma(kd, d(POW_LOG_HEAD[0]), logc);
    let t2 = t1 + r;
    let lo1 = F::fma(kd, d(POW_LOG_HEAD[1]), logctail);
    let lo2 = (t1 - t2) + r;
    let ar = a(0) * r; // A[0] = -0.5
    let ar2 = r * ar;
    let ar3 = r * ar2;
    // k*Ln2 + log(c) + r + A[0]*r*r.
    let hi = t2 + ar2;
    let lo3 = F::fma(ar, r, -ar2);
    let lo4 = (t2 - hi) + ar2;
    // p = log1p(r) - r - A[0]*r*r, added to the low parts by the last fma.
    let q = F::fma(F::fma(r, a(6), a(5)), ar2, F::fma(r, a(4), a(3)));
    let p = F::fma(ar2, q, F::fma(r, a(2), a(1)));
    let lo = F::fma(ar3, p, ((lo1 + lo2) + lo3) + lo4);
    let y = hi + lo;
    (y, (hi - y) + lo)
}

/// The special inputs of e_pow.c's `__pow`: the result, or the bits of `|x|` (normalized when
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

/// glibc's `pow` (e_pow.c `__pow`).
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
    let elo = F::fma(y, lo, F::fma(hi, y, -ehi));
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

/// The library's [`Fma`] where the CPU is not asked: [`Emulated`] on WebAssembly and on x86-64
/// builds without the `fma` feature (whose CPUs with FMA take `x86_fma`), [`Fused`] elsewhere.
#[cfg(any(
    target_arch = "wasm32",
    all(target_arch = "x86_64", not(target_feature = "fma"))
))]
type Portable = Emulated;
#[cfg(not(any(
    target_arch = "wasm32",
    all(target_arch = "x86_64", not(target_feature = "fma"))
)))]
type Portable = Fused;

/// x86-64 builds without the `fma` feature: copies compiled with it, called when the CPU has FMA
/// (the same bits as [`Emulated`]; the instruction is the fast way to get them).
#[cfg(all(target_arch = "x86_64", not(target_feature = "fma")))]
mod x86_fma {
    use super::{exp_impl, ln_impl, log10_impl, pow_impl, Fused};

    #[inline]
    pub(super) fn available() -> bool {
        std::arch::is_x86_feature_detected!("fma")
    }

    #[target_feature(enable = "fma")]
    pub(super) fn exp(x: f64) -> f64 {
        exp_impl::<Fused>(x)
    }

    #[target_feature(enable = "fma")]
    pub(super) fn ln(x: f64) -> f64 {
        ln_impl::<Fused>(x)
    }

    #[target_feature(enable = "fma")]
    pub(super) fn log10(x: f64) -> f64 {
        log10_impl::<Fused>(x)
    }

    #[target_feature(enable = "fma")]
    pub(super) fn pow(x: f64, y: f64) -> f64 {
        pow_impl::<Fused>(x, y)
    }
}

/// Calls the `x86_fma` copy (`$fast`) when the CPU has FMA, else `$portable`.
macro_rules! dispatch {
    ($fast:expr, $portable:expr) => {{
        #[cfg(all(target_arch = "x86_64", not(target_feature = "fma")))]
        {
            if x86_fma::available() {
                // SAFETY: the CPU has the `fma` feature that the copy is compiled with.
                return unsafe { $fast };
            }
        }
        $portable
    }};
}

fn exp_portable(x: f64) -> f64 {
    exp_impl::<Portable>(x)
}

fn ln_portable(x: f64) -> f64 {
    ln_impl::<Portable>(x)
}

fn log10_portable(x: f64) -> f64 {
    log10_impl::<Portable>(x)
}

fn pow_portable(x: f64, y: f64) -> f64 {
    pow_impl::<Portable>(x, y)
}

/// glibc's `exp` on x86-64 with FMA, bit for bit.
#[inline]
pub fn exp(x: f64) -> f64 {
    dispatch!(x86_fma::exp(x), exp_portable(x))
}

/// glibc's `log` on x86-64 with FMA, bit for bit.
#[inline]
pub fn ln(x: f64) -> f64 {
    dispatch!(x86_fma::ln(x), ln_portable(x))
}

/// glibc's `log10` on x86-64 with FMA, bit for bit.
#[inline]
pub fn log10(x: f64) -> f64 {
    dispatch!(x86_fma::log10(x), log10_portable(x))
}

/// glibc's `pow` on x86-64 with FMA, bit for bit.
#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    dispatch!(x86_fma::pow(x, y), pow_portable(x, y))
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
pub(crate) mod variants {
    use super::{exp_impl, ln_impl, log10_impl, pow_impl, Emulated, Fma, Fused};

    pub fn fma_emulated(a: f64, b: f64, c: f64) -> f64 {
        Emulated::fma(a, b, c)
    }
    pub fn exp_emulated(x: f64) -> f64 {
        exp_impl::<Emulated>(x)
    }
    pub fn ln_emulated(x: f64) -> f64 {
        ln_impl::<Emulated>(x)
    }
    pub fn log10_emulated(x: f64) -> f64 {
        log10_impl::<Emulated>(x)
    }
    pub fn pow_emulated(x: f64, y: f64) -> f64 {
        pow_impl::<Emulated>(x, y)
    }
    pub fn exp_fused(x: f64) -> f64 {
        exp_impl::<Fused>(x)
    }
    pub fn ln_fused(x: f64) -> f64 {
        ln_impl::<Fused>(x)
    }
    pub fn log10_fused(x: f64) -> f64 {
        log10_impl::<Fused>(x)
    }
    pub fn pow_fused(x: f64, y: f64) -> f64 {
        pow_impl::<Fused>(x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::variants as v;
    use super::{exp, ln, log10, pow};

    /// (function, x, y, result) as bits, function 0 `exp`, 1 `ln`, 2 `log10`, 3 `pow`: arguments
    /// where glibc 2.35 and 2.39 on x86-64 with FMA give `result` and glibc with its FMA variants
    /// off (the same algorithm as musl's, so WASI's) gives other bits.
    #[rustfmt::skip]
    const PINNED: [(u8, u64, u64, u64); 80] = [
        (0, 0xbfb7346e152d2b50, 0x0000000000000000, 0x3fed3a1b212e1efd), // exp[-1,1]
        (0, 0x3fe95f2fc2688c1e, 0x0000000000000000, 0x4001ad85e637996d), // exp[-1,1]
        (0, 0xbfe9b99b5117812a, 0x0000000000000000, 0x3fdca51949f00589), // exp[-1,1]
        (0, 0xc008608684fe4de0, 0x0000000000000000, 0x3fa8514291f81ddb), // exp[-60,0]
        (0, 0xc03ea6da0f9bb927, 0x0000000000000000, 0x3d2b73c49f291119), // exp[-60,0]
        (0, 0xc04a250d7859b4ac, 0x0000000000000000, 0x3b379ff8df28f1c1), // exp[-60,0]
        (0, 0xc04b7466a52670f0, 0x0000000000000000, 0x3afb8598c6df27a7), // exp[-745,710]
        (0, 0xc04e5001b3d27720, 0x0000000000000000, 0x3a773532ba6abab1), // exp[-745,710]
        (0, 0x408292cbb12e1cd6, 0x0000000000000000, 0x758615d8d0597efd), // exp[-745,710]
        (0, 0xc086253befc3fa92, 0x0000000000000000, 0x000c5d08786e7dd3), // exp[-745,-708]
        (0, 0xc086247de02b3c29, 0x0000000000000000, 0x000d90d1b2f6a67f), // exp[-745,-708]
        (0, 0xc08630542322e779, 0x0000000000000000, 0x000316d8552aba69), // exp[-745,-708]
        (0, 0xc005af8431352386, 0x0000000000000000, 0x3fb1057fba5e6c13), // exp[bits]
        (0, 0xc04ec2ddd59eb197, 0x0000000000000000, 0x3a62ebe01ba5b68b), // exp[bits]
        (0, 0x404456c5103f2292, 0x0000000000000000, 0x4399bcbf6fa78141), // exp[bits]
        (1, 0x3ff173bfb5d375f0, 0x0000000000000000, 0x3fb63d5c42031633), // ln[0.5,2]
        (1, 0x3ff2203196002bec, 0x0000000000000000, 0x3fbfef4d89eeeb90), // ln[0.5,2]
        (1, 0x3ff50394db9b6778, 0x0000000000000000, 0x3fd17245f941fc78), // ln[0.5,2]
        (1, 0x3feee9f1015f2ff7, 0x0000000000000000, 0xbfa1ae31e42e66fb), // ln[near1]
        (1, 0x3ff069c9c9d1c7db, 0x0000000000000000, 0x3f9a1c7e3882fb17), // ln[near1]
        (1, 0x3ff09fa57208a91e, 0x0000000000000000, 0x3fa393a2cd7cabbf), // ln[near1]
        (1, 0x413fc35fc9b9d5f0, 0x0000000000000000, 0x402d18ea4fd58f10), // ln[1,1e7]
        (1, 0x40f9633d3eacb7e9, 0x0000000000000000, 0x40271aa38fef55bb), // ln[1,1e7]
        (1, 0x40f206ccd46ef0cf, 0x0000000000000000, 0x40266b52564646c8), // ln[1,1e7]
        (1, 0x729b79abeb172e68, 0x0000000000000000, 0x40818feb5823828d), // ln[posbits]
        (1, 0x3d61dd35777df42b, 0x0000000000000000, 0xc03c4f0f7e9132e0), // ln[posbits]
        (1, 0x3bb38ab5519c1979, 0x0000000000000000, 0xc047778e9d465d6d), // ln[posbits]
        (1, 0x000976060b0abca6, 0x0000000000000000, 0xc086275fe65325c5), // ln[subnormal]
        (1, 0x3fece4e6c6ce0e67, 0x0000000000000000, 0xbfba2305570b4413), // ln[bits]
        (1, 0x3b89ce2fc16490f6, 0x0000000000000000, 0xc0485e235c8054cc), // ln[bits]
        (1, 0x54fc61b6b5fcf89e, 0x0000000000000000, 0x406d2f0f392a0aa9), // ln[bits]
        (2, 0x401e4458d91a5290, 0x0000000000000000, 0x3fec2005febfa766), // log10[posbits]
        (2, 0x3ffa3fe2cec11f1b, 0x0000000000000000, 0x3fcb852f1f8b31de), // log10[posbits]
        (2, 0x38083b5e0cf0cd5c, 0x0000000000000000, 0xc04306782824ebdc), // log10[posbits]
        (2, 0x3fe7684d532711f3, 0x0000000000000000, 0xbfc161c57c4cf001), // log10[1e-300,1e3]
        (2, 0x3ff97f340663b4e1, 0x0000000000000000, 0x3fc9e72b000ac001), // log10[1e-300,1e3]
        (2, 0x3fe6835ae0cdc81e, 0x0000000000000000, 0xbfc38c26b0345000), // log10[1e-300,1e3]
        (2, 0x3fb894eb91f41dc5, 0x0000000000000000, 0xbff04833fd669cf2), // log10[bits]
        (2, 0x44a334e05ee3f129, 0x0000000000000000, 0x4036a8159013bc5f), // log10[bits]
        (2, 0x40727763e100434c, 0x0000000000000000, 0x4003c39645736ba2), // log10[bits]
        (3, 0x40135012c00000dc, 0x3ff8bce027903aa4, 0x4026d0e1007ec7a7), // pow[0,10]x[-5,5]
        (3, 0x401e613a32570e30, 0x4005d8a80d052564, 0x406fba46c387441f), // pow[0,10]x[-5,5]
        (3, 0x4023da2f44fa623a, 0x400c131266523fbc, 0x40a897b395667ccf), // pow[0,10]x[-5,5]
        (3, 0x400ea5eb0693542b, 0xbfa1b636d3a50080, 0x3fee8c14724b6659), // pow[0,10]x[-5,5]
        (3, 0x3ff2bd6cd5900238, 0x3fe0f1193b7fb1d6, 0x3ff16587e3f39be7), // pow[0,4]x[0,1]
        (3, 0x400d2b0e306f1031, 0x3fe4e31006bda2ac, 0x40029cbdec3e1a9d), // pow[0,4]x[0,1]
        (3, 0x3ffc9ab45a50f7fa, 0x3fc3c0db53373ff8, 0x3ff180334679ac79), // pow[0,4]x[0,1]
        (3, 0x4000372b47465044, 0x3fae92b556245c80, 0x3ff0b080a4bcc283), // pow[0,4]x[0,1]
        (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
        (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
        (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
        (3, 0x4024000000000000, 0x405f800000000000, 0x5a17a2ecc414a03f), // pow10^int
        (3, 0x3feffda9549186e8, 0xc0f15ccfe8c6db48, 0x41c39b28f0000f03), // pow[near1]x[big]
        (3, 0x3feffb0c1bcfdce7, 0x4117e74593989f50, 0x2a93c6d03463cbcf), // pow[near1]x[big]
        (3, 0x3fefff84d6b72c19, 0x40e2d94b6f0aed90, 0x3fba86534d13fc21), // pow[near1]x[big]
        (3, 0x3ff00242ac054991, 0x412adc2e0866e670, 0x6bb7909b8b62dc17), // pow[near1]x[big]
        (3, 0xbffd5aafcb95bcf2, 0xc058c00000000000, 0xba840da9e7914d69), // pow[neg]x[int]
        (3, 0xbff1a261f339bae2, 0x4024000000000000, 0x400528713c1e5757), // pow[neg]x[int]
        (3, 0xbff60ac87808c33a, 0xc046800000000000, 0xbea265a0852bf0db), // pow[neg]x[int]
        (3, 0xbffac9c1474afc0a, 0x4008000000000000, 0xc012c5cdb12e4433), // pow[neg]x[int]
        (3, 0x300a8d6a7a6a0dd1, 0xbfdc074dc018c9e0, 0x46e47adabdcedd37), // pow[posbits]x[-2,2]
        (3, 0x51ff461ae6bf320e, 0x3fe089adcdf84630, 0x4944406bad06ff6d), // pow[posbits]x[-2,2]
        (3, 0x54b56930780fc9c6, 0xbff18c4d7623a760, 0x29255c51b09dad81), // pow[posbits]x[-2,2]
        (3, 0x58814ada472ea9ea, 0xbfe20860cdac6170, 0x321632ed94d182f7), // pow[posbits]x[-2,2]
        (3, 0x0002a8b242943bb4, 0x3fd3b506dc3c7710, 0x2c36b26193f158e9), // pow[subn]x[-1,1]
        (3, 0x00080d8a74ecc080, 0xbf9688e2935d3240, 0x4156d2fdd889844b), // pow[subn]x[-1,1]
        (3, 0x000ec724df16f69c, 0xbfd2decfcea072a0, 0x52c4acde4ad4d013), // pow[subn]x[-1,1]
        (3, 0x000ff1f02d7edff4, 0xbfe9921dbff5bee2, 0x72f97678a909b4f9), // pow[subn]x[-1,1]
        (3, 0x46698bef4eb624f5, 0x3f9cf6d49d358f15, 0x401e89750613d523), // pow[bits]x[bits]
        (3, 0x24df72154c46d6ec, 0x3ff136ff4a783591, 0x22d11b37cb48b749), // pow[bits]x[bits]
        (3, 0x2801c12f0e55700f, 0x3f34b1cbaa38d9c1, 0x3fed6d8242cc45a3), // pow[bits]x[bits]
        (3, 0x1821408f368c1705, 0xbf4a699f2074e54e, 0x3ff6d687cd32677d), // pow[bits]x[bits]
        (3, 0x3ff9f95f9d1b139c, 0x40921fcc90e0ec26, 0x729c035093405175), // pow[ovf]
        (3, 0x3ffd969a00f66d59, 0x4090d874d00a9f83, 0x7bb2de524f617177), // pow[ovf]
        (3, 0x3ffbc0a346964152, 0x4091bbfc0a276b1c, 0x784bf8f8a0a60da1), // pow[ovf]
        (3, 0x3ffb5c4386785a90, 0x4093fc0d93a96f72, 0x7dcf890f34c08f57), // pow[ovf]
        (3, 0x3fda8ee00ea9a966, 0x407b161b1bebc578, 0x1d90d6437bd90ae9), // pow[unf]
        (3, 0x3fd87ccca77b5bc0, 0x407ea46956b476f2, 0x1576152909f85fb6), // pow[unf]
        (3, 0x3fdeab6e8bc8116e, 0x407d49fc164ceaee, 0x20d96c56eab1fe33), // pow[unf]
        (3, 0x3fe96ebcd599258e, 0x407b1665889bbd2c, 0x36f4bab57ab5f09f), // pow[unf]
    ];

    fn call(f: u8, x: f64, y: f64) -> [f64; 3] {
        match f {
            0 => [exp(x), v::exp_emulated(x), v::exp_fused(x)],
            1 => [ln(x), v::ln_emulated(x), v::ln_fused(x)],
            2 => [log10(x), v::log10_emulated(x), v::log10_fused(x)],
            _ => [pow(x, y), v::pow_emulated(x, y), v::pow_fused(x, y)],
        }
    }

    fn bits(x: f64) -> u64 {
        if x.is_nan() {
            0x7ff8_0000_0000_0000
        } else {
            x.to_bits()
        }
    }

    #[test]
    fn every_target_writes_glibcs_fma_bits_where_other_libms_differ() {
        for (f, x, y, want) in PINNED {
            let (x, y) = (f64::from_bits(x), f64::from_bits(y));
            for (i, got) in call(f, x, y).into_iter().enumerate() {
                assert_eq!(got.to_bits(), want, "function {f} x={x:e} y={y:e} way {i}");
            }
        }
    }

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }

        fn unif(&mut self, a: f64, b: f64) -> f64 {
            a + (b - a) * ((self.next() >> 11) as f64 / (1u64 << 53) as f64)
        }
    }

    /// Arguments of every range the functions treat apart: (function, x, y).
    fn arguments(n: usize) -> Vec<(u8, f64, f64)> {
        let mut r = Rng(12345);
        let mut out = Vec::new();
        for _ in 0..n {
            let any = f64::from_bits(r.next());
            out.extend([
                (0, r.unif(-1.0, 1.0), 0.0),
                (0, r.unif(-60.0, 0.0), 0.0),
                (0, r.unif(-745.2, 709.8), 0.0),
                (0, r.unif(-745.2, -708.0), 0.0),
                (0, any, 0.0),
                (0, r.unif(-1.0, 1.0) * 1e-17, 0.0),
                (1, r.unif(0.5, 2.0), 0.0),
                (1, r.unif(0.937, 1.0647), 0.0),
                (1, r.unif(1.0, 1e7), 0.0),
                (1, f64::from_bits(r.next() >> 1), 0.0),
                (1, f64::from_bits(r.next() >> 12), 0.0),
                (1, any, 0.0),
                (2, f64::from_bits(r.next() >> 1), 0.0),
                (2, any, 0.0),
                (3, r.unif(0.0, 10.0), r.unif(-5.0, 5.0)),
                (3, 10.0, ((r.next() % 700) as i32 - 400) as f64),
                (3, r.unif(0.999, 1.001), r.unif(-1e6, 1e6)),
                (3, -r.unif(0.0, 3.0), ((r.next() % 200) as i64 - 100) as f64),
                (3, f64::from_bits(r.next() >> 1), r.unif(-2.0, 2.0)),
                (3, f64::from_bits(r.next() >> 12), r.unif(-1.0, 1.0)),
                (3, any, f64::from_bits(r.next())),
                (3, r.unif(1.5, 3.0), r.unif(900.0, 1900.0)),
                (3, r.unif(0.1, 0.9), r.unif(300.0, 800.0)),
            ]);
        }
        let special = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            2.0,
            -2.0,
            0.5,
            10.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            709.782712893384,
            709.7827128933841,
            -745.1332191019411,
            -745.1332191019412,
            -708.3964185322641,
            1e-300,
            5e-324,
            f64::MAX,
            f64::MIN_POSITIVE,
            3.0,
            -3.0,
            0.1,
            1e300,
            1024.0,
            -1075.0,
            1e-20,
            f64::from_bits(0x3be0_0000_0000_0000),
            f64::from_bits(0x43e0_0000_0000_0000),
        ];
        for x in special {
            out.extend([(0, x, 0.0), (1, x, 0.0), (2, x, 0.0)]);
            for y in special {
                out.push((3, x, y));
            }
        }
        out
    }

    #[test]
    fn the_emulated_fma_way_gives_the_fused_way_bits() {
        for (f, x, y) in arguments(20_000) {
            let [public, emulated, fused] = call(f, x, y);
            assert_eq!(bits(emulated), bits(fused), "function {f} x={x:e} y={y:e}");
            assert_eq!(bits(public), bits(fused), "function {f} x={x:e} y={y:e}");
        }
    }

    #[test]
    fn the_emulated_fma_rounds_once() {
        let mut r = Rng(99);
        for i in 0..1_000_000 {
            let ea = (r.next() % 120) as i32 - 60;
            let eb = (r.next() % 120) as i32 - 60;
            let sign = |r: &mut Rng| if r.next() & 1 == 0 { 1.0 } else { -1.0 };
            let a = r.unif(1.0, 2.0) * 2f64.powi(ea) * sign(&mut r);
            let b = r.unif(1.0, 2.0) * 2f64.powi(eb) * sign(&mut r);
            let ph = a * b;
            let c = match i % 5 {
                0 => -ph,
                1 => -ph + r.unif(-1.0, 1.0) * ph.abs() * 2f64.powi(-((r.next() % 60) as i32)),
                2 => r.unif(-4.0, 4.0) * 2f64.powi(ea + eb + (r.next() % 140) as i32 - 70),
                3 => {
                    // ph + c on a midpoint of the sum's precision
                    let s = ph * r.unif(1.0, 3.0);
                    let u = f64::from_bits(s.to_bits() + 1) - s;
                    s - ph + u * 0.5
                }
                _ => a.mul_add(b, -ph) * r.unif(-3.0, 3.0) - ph * 0.5,
            };
            assert_eq!(
                bits(v::fma_emulated(a, b, c)),
                bits(a.mul_add(b, c)),
                "a={a:e} b={b:e} c={c:e}"
            );
        }
    }

    #[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
    #[test]
    #[expect(
        clippy::disallowed_methods,
        reason = "compares the port with glibc itself"
    )]
    fn glibc_with_fma_gives_the_same_bits() {
        if !std::arch::is_x86_feature_detected!("fma") {
            return; // glibc then takes its variants without FMA
        }
        for (f, x, y) in arguments(20_000) {
            let glibc = match f {
                0 => x.exp(),
                1 => x.ln(),
                2 => x.log10(),
                _ => x.powf(y),
            };
            assert_eq!(
                bits(call(f, x, y)[0]),
                bits(glibc),
                "function {f} x={x:e} y={y:e}"
            );
        }
    }
}
