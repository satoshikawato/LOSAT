//! `log_slice`: the natural logarithm of every element, through libm.
//!
//! The strict set calls libm's `log` here.  The optional round-4 patch
//! replaces this file by a bit-identical SIMD clone of this machine's
//! `__log_fma` behind `LOSAT_X_LOGCLONE` (see that patch's README).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:636,679
//! ```c
//! temp = log(x[k] / q[k]);
//! ...
//!     scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
//! ```
//! NCBI calls libm `log` element by element in EvaluateReFunctions and ComputeScoresFromProbs.
//! In this branch `log_slice` calls the same libm function (`f64::ln`) on the same arguments in
//! the same order, so each result is the one NCBI gets on this machine.

// No NCBI counterpart: an empty function so that main.rs builds with or without the optional patch.
/// No shadow mode in the strict version (nothing to compare).
pub fn print_shadow_stats() {}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:636,679
// ```c
// temp = log(x[k] / q[k]);
// ...
//     scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
// ```
// Same libm call as the C `log` above, one element at a time.
/// `out[k] = ln(src[k])`.
#[inline]
pub(crate) fn log_slice(src: &[f64], out: &mut [f64]) {
    debug_assert_eq!(src.len(), out.len());
    for (o, &x) in out.iter_mut().zip(src) {
        *o = x.ln();
    }
}
