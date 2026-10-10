//! `log_slice`: the natural logarithm of every element, through libm.
//!
//! The strict set calls libm's `log` here.  The optional round-4 patch
//! replaces this file by a bit-identical SIMD clone of this machine's
//! `__log_fma` behind `LOSAT_X_LOGCLONE` (see that patch's README).

/// No shadow mode in the strict version (nothing to compare).
pub fn print_shadow_stats() {}

/// `out[k] = ln(src[k])`.
#[inline]
pub(crate) fn log_slice(src: &[f64], out: &mut [f64]) {
    debug_assert_eq!(src.len(), out.len());
    for (o, &x) in out.iter_mut().zip(src) {
        *o = x.ln();
    }
}
