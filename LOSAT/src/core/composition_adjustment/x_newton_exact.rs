//! EXPERIMENT (LOSAT_X_NEWTONEXACT / LOSAT_X_NEWTONEXACTSHADOW):
//! `Blast_OptimizeTargetFrequencies` for the relative-entropy-constrained
//! problem, bit for bit.
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:602-713
//! (`Blast_OptimizeTargetFrequencies`) and nlm_linear_algebra.c:120-171
//! (`Nlm_FactorLtriangPosDef`, `Nlm_SolveLtriangPosDef`).
//!
//! `optimize_target_frequencies_reference` in `adjust_scores.rs` is a
//! line-by-line port of the reference.  This module computes **the same
//! IEEE-754 operations on every scalar value, in the same order per value**,
//! so that each value it produces is the same bit pattern.  What changes:
//!
//! * storage: fixed-size arrays, no allocation per call;
//! * which values are advanced next: independent values (different array
//!   elements, different columns of the Cholesky factor) are computed side by
//!   side so that the compiler puts them in SIMD lanes.  Every lane performs
//!   the reference's own operation on the reference's own operands.  Rust
//!   never fuses `a * b + c` and never re-associates sums; SIMD add / sub /
//!   mul / div / sqrt are IEEE-754 per element;
//! * the Cholesky factorisation of the 40 x 40 matrix `W = J D^-1 J^T` uses
//!   its structure (optimize_target_freq.c:95-126, 408-470): rows 0..19 are
//!   diagonal, rows 20..38 are dense only in columns 0..19 and on the
//!   diagonal, row 39 is dense.  Element `(i, j)` still receives
//!   `temp -= L[i][k] * L[j][k]` for `k = 0, 1, ..., j - 1` in that order and
//!   is then divided by `L[j][j]` (square-rooted on the diagonal), exactly as
//!   `Nlm_FactorLtriangPosDef` does.  Terms that the reference computes as
//!   `temp -= L[i][k] * 0.0` are skipped where `temp` is finite and not
//!   `-0.0` and `L[i][k]` is finite (then `L[i][k] * 0.0` is `+0.0` or
//!   `-0.0`, and `temp - (+0.0) == temp`, `temp - (-0.0) == temp`); the
//!   conditions are checked (`structure_ok`), otherwise the literal dense
//!   loops run;
//! * the triangular solves are done column by column (forward) and as
//!   descending dot products (backward); each element receives the same
//!   terms in the same order as in `Nlm_SolveLtriangPosDef`.  The zero
//!   products of the diagonal block are performed literally;
//! * the first iteration starts from `x == q`, where every quotient
//!   `x[k] / q[k]` is exactly `1.0` (IEEE-754: `v / v == 1` for finite
//!   non-zero `v`), so the 400 calls `ln(1.0)` of that iteration are replaced
//!   by one call `ln(1.0)` of the same libm function, made once per process
//!   on an input the optimiser cannot see;
//! * `values[0]` of `EvaluateReFunctions` (the objective) is never read by
//!   the reference after it is computed, so it is not computed here;
//! * sequential sums (`values[1]`, `W[39][39]`, `z[39]`) are accumulated in the
//!   same order, in the same loop as independent element-wise work.
//!
//! What this module does **not** change: the libm `ln` (same function, same
//! arguments), the Euclidean norms, the step bound, the stopping test.
//!
//! `LOSAT_X_NEWTONEXACTSHADOW=1` runs both this code and the reference on every
//! call and aborts on the first value that differs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:762-768,774-777,781-785
//! ```c
//!         if ( !(rnorm > tol) ) {
//!             /* We converged at the current iterate */
//!             break;
//!         } else {
//!             /* we did not converge, so increment the iteration counter
//!                and start a new iteration */
//!             if (++its <= maxits) {
//! ...
//!                 FactorReNewtonSystem(newton_system, x, z, grads,
//!                                      constrain_rel_entropy, workspace);
//!                 SolveReNewtonSystem(resids_x, resids_z, newton_system,
//!                                     workspace);
//! ...
//!                 alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
//!                 alpha *= 0.95;
//!                 Nlm_AddVectors(x, n, alpha, resids_x);
//!                 Nlm_AddVectors(z, m, alpha, resids_z);
//! ```
//! The Newton loop of `optimize` below is this loop. Each call here is a step of one
//! iteration: residuals, convergence test, factor, solve, step bound, update of x and z.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,171-177
//! ```c
//!     for (i = 0;  i < n;  i++) {
//!         for (j = 0;  j < i;  j++) {
//!             temp = A[i][j];
//!             for (k = 0;  k < j;  k++) {
//!                 temp -= A[i][k] * A[j][k];
//!             }
//!             A[i][j] = temp/A[j][j];
//! ...
//!     for (i = 0;  i < n;  i++) {
//!         temp = x[i];
//!         for (j = 0;  j < i;  j++) {
//!             temp -= L[i][j] * x[j];
//!         }
//!         x[i] = temp/L[i][i];
//!     }
//! ```
//! These are the Cholesky factorisation and the forward solve whose per-element operations
//! and order are kept (see the list above). The rest of this module changes only the data layout
//! and which independent elements are advanced together.

use super::adjust_scores::{
    COMPO_NUM_TRUE_AA, K_COMPO_ADJUST_ERR_TOLERANCE, K_COMPO_ADJUST_ITERATION_LIMIT,
};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// No NCBI counterpart: cycle counters for the profile output; it does not change any value NCBI computes.
/// Cycle counts per phase (LOSAT_X_NEWTONEXACT_PROF=1): logs, residuals,
/// factor preparation, factor, solve, step, calls, iterations.
pub(crate) static PROF: [AtomicU64; 12] = [
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
    AtomicU64::new(0),
];
// No NCBI counterpart: reads LOSAT_X_NEWTONEXACT_PROF once; it does not change any value NCBI computes.
pub(crate) fn prof_enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_NEWTONEXACT_PROF").is_some())
}
// No NCBI counterpart: reads the time-stamp counter for the profile output; it does not change any value NCBI computes.
#[inline(always)]
fn tsc() -> u64 {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: rdtsc has no preconditions.
    unsafe {
        core::arch::x86_64::_rdtsc()
    }
    #[cfg(not(target_arch = "x86_64"))]
    0
}
// No NCBI counterpart: adds elapsed cycles to a profile counter; it does not change any value NCBI computes.
#[inline(always)]
fn prof_add(i: usize, t0: u64, on: bool) -> u64 {
    if on {
        let t = tsc();
        PROF[i].fetch_add(t - t0, Ordering::Relaxed);
        t
    } else {
        0
    }
}
// No NCBI counterpart: prints the profile counters; it does not change any value NCBI computes.
pub(crate) fn prof_print() {
    if !prof_enabled() {
        return;
    }
    let names = [
        "logs", "resid", "fprep", "factor", "solve", "step", "calls", "iters", "f_sparse",
        "f_terms", "f_block", "f_check",
    ];
    for (i, n) in names.iter().enumerate() {
        eprintln!(
            "[NEWTONEXACT_PROF] {}={}",
            n,
            PROF[i].load(Ordering::Relaxed)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:725-727
// ```c
//     n  = alphsize * alphsize;
//     mA = 2 * alphsize - 1;
//     m  = constrain_rel_entropy ? mA + 1 : mA;
// ```
// `N` is `alphsize`, `NN` is `n`, and `M` is `m` for `constrain_rel_entropy` set. This module
// handles only that case; other cases run the reference.
const N: usize = COMPO_NUM_TRUE_AA;
const NN: usize = N * N;
/// Number of constraints with the relative-entropy constraint (`m`).
const M: usize = 2 * N;
/// Lane padding for the 20-row vectors of the trailing block.
const P: usize = 24;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// ```c
// int
// Blast_OptimizeTargetFrequencies(double x[],
// ...
//                                 int constrain_rel_entropy,
//                                 double relative_entropy,
//                                 double tol,
//                                 int maxits)
// ```
// No NCBI counterpart: reads LOSAT_X_NEWTONEXACT and LOSAT_X_NEWTONEXACTSHADOW once. It chooses
// between the port of this function (reference), this module, or both with a comparison.
/// 0 = reference, 1 = this module, 2 = both (checked bit for bit on every call).
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_NEWTONEXACTSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_NEWTONEXACT").is_some() {
            1
        } else {
            0
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:747-748,635-636
// ```c
//     /* Use q as the initial value for x */
//     memcpy(x, q, n * sizeof(double));
// ...
//     for (k = 0;  k < alphsize * alphsize;  k++) {
//         temp = log(x[k] / q[k]);
// ```
// NCBI starts from `x = q`, so in the first iteration every `log(x[k] / q[k])` is `log(1.0)`.
// This is the value of one such call, made once per process.
/// `ln(1.0)` as the libm of this process computes it (not constant-folded).
fn ln_one() -> f64 {
    static LN_ONE: OnceLock<f64> = OnceLock::new();
    *LN_ONE.get_or_init(|| std::hint::black_box(1.0f64).ln())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:193-208
// ```c
//     double sum   = 1.0;   /* sum of squares of elements in v */
//     double scale = 0.0;   /* a scale factor for the elements in v */
//     int i;                /* iteration index */
//     for (i = 0;  i < n;  i++) {
//         if (v[i] != 0.0) {
//             double absvi = fabs(v[i]);
//             if (scale < absvi) {
//                 sum = 1.0 + sum * (scale/absvi) * (scale/absvi);
//                 scale = absvi;
//             } else {
//                 sum += (absvi/scale) * (absvi/scale);
//             }
//         }
//     }
//     return scale * sqrt(sum);
// ```
/// The reference's `Nlm_EuclideanNorm`, copied verbatim.
#[inline(always)]
fn euclidean_norm(v: &[f64]) -> f64 {
    let mut sum = 1.0;
    let mut scale = 0.0;
    for &value in v {
        if value != 0.0 {
            let abs_value = value.abs();
            if scale < abs_value {
                sum = 1.0 + sum * (scale / abs_value) * (scale / abs_value);
                scale = abs_value;
            } else {
                sum += (abs_value / scale) * (abs_value / scale);
            }
        }
    }
    scale * sum.sqrt()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:132-135,525-531
// ```c
//             W[j][j] += dd;
//             if (i > 0) {
//                 W[i + alphsize - 1][j] += dd;
//                 W[i + alphsize - 1][i + alphsize - 1] += dd;
// ...
//         W[m - 1][m - 1] = 0.0;
//         for (i = 0;  i < n;  i++) {
//             workspace[i] = Dinv[i] * grad_re[i];
//             W[m - 1][m - 1] += grad_re[i] * workspace[i];
//         }
//         MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
// ```
// The parts are the sums that `ScaledSymmetricProductA` and `FactorReNewtonSystem` put into
// the lower triangle of `W` (column sums `W[j][j]`, row sums `W[i+19][i+19]`, the entries
// `W[i+19][j]`, and the last row).
/// The lower triangle of `W` as the reference holds it, from its parts.
struct Wparts<'a> {
    /// `W[j][j]`, `j < 20`: column sums of `dinv`.
    wdiag: &'a [f64; N],
    /// `W[i + 19][i + 19]`, `1 <= i < 20`: row sums of `dinv` (index `i`).
    wrow: &'a [f64; N],
    /// `W[i + 19][j] = 0.0 + dinv[i * 20 + j]`.
    dinv: &'a [f64; NN],
    /// `W[39][c]`, `c < 40`.
    w39: &'a [f64; M],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
// ```c
//     for (i = 0;  i < n;  i++) {
//         for (j = 0;  j < i;  j++) {
//             temp = A[i][j];
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
//             A[i][j] = temp/A[j][j];
// ...
//         temp = A[i][i];
//         for (k = 0;  k < i;  k++) {
//             temp -= A[i][k] * A[i][k];
//         }
//         A[i][i] = sqrt(temp);
// ```
// The entries of `A` after this function, kept in blocks. Only the storage differs.
/// The Cholesky factor in structured storage.
struct Factor {
    /// `L[j][j]` for `j < 20`.
    d: [f64; N],
    /// `lsp[k][r] = L[20 + r][k]` for `k < 20`, `r < 20` (rows 20..39).
    lsp: [[f64; P]; N],
    /// `c[jj][ii] = L[20 + ii][20 + jj]` for `ii >= jj` (the trailing block,
    /// column-major).  Entries with `ii < jj` and the padding are scratch.
    c: [[f64; P]; N],
}

// No NCBI counterpart: constructor and test accessor of the structured storage; it does not change any value NCBI computes.
impl Factor {
    #[inline(always)]
    fn new() -> Self {
        Factor {
            d: [0.0; N],
            lsp: [[0.0; P]; N],
            c: [[0.0; P]; N],
        }
    }
    /// `L[i][k]` for any `i >= k`.
    #[cfg(test)]
    fn get(&self, i: usize, k: usize) -> f64 {
        if i < N {
            if i == k {
                self.d[i]
            } else {
                0.0
            }
        } else if k < N {
            self.lsp[k][i - N]
        } else {
            self.c[k - N][i - N]
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
// ```c
//     for (i = 0;  i < n;  i++) {
//         for (j = 0;  j < i;  j++) {
//             temp = A[i][j];
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
//             A[i][j] = temp/A[j][j];
// ...
//         temp = A[i][i];
//         for (k = 0;  k < i;  k++) {
//             temp -= A[i][k] * A[i][k];
//         }
//         A[i][i] = sqrt(temp);
// ```
// The loops below are this function on a dense 40 x 40 copy of `W`, in the same order.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:121-125,132-135,525-531
// ```c
//     for (rowW = 0;  rowW < m;  rowW++) {
//         for (colW = 0;  colW <= rowW;  colW++) {
//             W[rowW][colW] = 0.0;
//         }
//     }
// ...
//             W[j][j] += dd;
//             if (i > 0) {
//                 W[i + alphsize - 1][j] += dd;
//                 W[i + alphsize - 1][i + alphsize - 1] += dd;
// ...
//         W[m - 1][m - 1] = 0.0;
//         for (i = 0;  i < n;  i++) {
//             workspace[i] = Dinv[i] * grad_re[i];
//             W[m - 1][m - 1] += grad_re[i] * workspace[i];
//         }
//         MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
// ```
// `W` is filled as `ScaledSymmetricProductA` and `FactorReNewtonSystem` fill it.
/// `Nlm_FactorLtriangPosDef` on `W`, literally (dense lower triangle).
#[inline(never)]
fn factor_literal(w: &Wparts, out: &mut Factor) {
    let mut a = [[0.0f64; M]; M];
    for j in 0..N {
        a[j][j] = w.wdiag[j];
    }
    for i in 1..N {
        for j in 0..N {
            a[i + N - 1][j] = 0.0 + w.dinv[i * N + j];
        }
        a[i + N - 1][i + N - 1] = w.wrow[i];
    }
    a[M - 1] = *w.w39;
    for i in 0..M {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
    for j in 0..N {
        out.d[j] = a[j][j];
    }
    for k in 0..N {
        for r in 0..N {
            out.lsp[k][r] = a[N + r][k];
        }
    }
    for jj in 0..N {
        for ii in jj..N {
            out.c[jj][ii] = a[N + ii][N + jj];
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:146-149
// ```c
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
//             A[i][j] = temp/A[j][j];
// ```
// The terms `temp -= A[i][k] * A[j][k]` that this check allows to skip are those with
// `A[j][k] == +0.0`, as the comment below explains.
/// Whether the structured factorisation reproduces the literal one.
///
/// The skipped terms are `temp -= L[i][k] * L[j][k]` with `L[j][k] == +0.0`
/// (rows 0..19 off the diagonal).  `L[j][k]` is `+0.0` exactly: it is
/// `(0.0 - 0.0 * 0.0 - ...) / L[k][k]` with `L[k][k] = sqrt(wdiag[k]) > 0`.
/// With `L[i][k]` finite the product is `+0.0` or `-0.0`, and subtracting
/// it leaves `temp` unchanged unless `temp` is `-0.0` (which `-0.0 - (-0.0)`
/// turns into `+0.0`) or NaN.  `temp` is `W[i][j]`: `0.0 + dinv > 0` for rows
/// 20..38, and `W[39][j]` for row 39, which is required finite and not `-0.0`.
/// `L[i][k]` finite needs `dinv` finite and `wdiag > 0` (rows 20..38) and
/// `W[39][k]` finite (row 39).  Everything else is computed literally.
#[inline(always)]
fn structure_ok(w: &Wparts, dinv_ok: bool) -> bool {
    let mut ok = dinv_ok;
    for &v in w.w39.iter() {
        ok &= v.is_finite() & (v.to_bits() != (-0.0f64).to_bits());
    }
    for &v in w.wdiag.iter() {
        ok &= v.is_finite() & (v > 0.0);
    }
    ok
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
// ```c
//     for (i = 0;  i < n;  i++) {
//         for (j = 0;  j < i;  j++) {
//             temp = A[i][j];
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
//             A[i][j] = temp/A[j][j];
// ...
//         temp = A[i][i];
//         for (k = 0;  k < i;  k++) {
//             temp -= A[i][k] * A[i][k];
//         }
//         A[i][i] = sqrt(temp);
// ```
// Element `(i, j)` receives the same terms in the same order and the same final division or
// square root; only the terms known to be zero are left out (checked by `structure_ok`).
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
// ```c
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             double dd;     /* an individual diagonal element */
//             dd = diagonal[i * alphsize + j];
//             W[j][j] += dd;
//             if (i > 0) {
//                 W[i + alphsize - 1][j] += dd;
//                 W[i + alphsize - 1][i + alphsize - 1] += dd;
//             }
//         }
// ```
// This loop fills `W`, so rows 0..19 are diagonal and rows 20..38 are dense only in columns 0..19.
/// `Nlm_FactorLtriangPosDef` on `W` using its structure (see module comment).
#[inline(always)]
fn factor_structured(w: &Wparts, f: &mut Factor) {
    let on = prof_enabled();
    let mut t = if on { tsc() } else { 0 };
    // Rows 0..19: L[j][j] = sqrt(W[j][j]); the rest of each row is
    // (0.0 - 0.0 * 0.0 - ...) / L[j][j] = +0.0, never stored.
    for j in 0..N {
        f.d[j] = w.wdiag[j].sqrt();
    }
    // Rows 20..39, columns 0..19: temp = W[i][j]; the k < j terms are
    // L[i][k] * L[j][k] with L[j][k] == +0.0 (skipped, see structure_ok);
    // L[i][j] = temp / L[j][j].
    for r in 0..N - 1 {
        let dr = &w.dinv[(r + 1) * N..(r + 2) * N];
        let mut row = [0.0f64; N];
        for j in 0..N {
            row[j] = (0.0 + dr[j]) / f.d[j];
        }
        for j in 0..N {
            f.lsp[j][r] = row[j];
        }
    }
    {
        let mut row = [0.0f64; N];
        for j in 0..N {
            row[j] = w.w39[j] / f.d[j];
        }
        for j in 0..N {
            f.lsp[j][N - 1] = row[j];
        }
    }
    // The trailing block T[ii][jj] = W[20 + ii][20 + jj] (ii >= jj), column-major.
    for jj in 0..N {
        f.c[jj] = [0.0; P];
    }
    for jj in 0..N - 1 {
        f.c[jj][jj] = w.wrow[jj + 1];
        f.c[jj][N - 1] = w.w39[N + jj];
    }
    f.c[N - 1][N - 1] = w.w39[M - 1];
    t = prof_add(8, t, on);
    // Terms k = 0..19: T[ii][jj] -= L[20+ii][k] * L[20+jj][k], in k order,
    // column jj in lanes from the 4-aligned row below jj (rows ii < jj and
    // the padding are scratch).
    for jj in 0..N {
        match jj / 4 {
            0 => column_terms::<5>(&mut f.c[jj], &f.lsp, jj),
            1 => column_terms::<4>(&mut f.c[jj], &f.lsp, jj),
            2 => column_terms::<3>(&mut f.c[jj], &f.lsp, jj),
            3 => column_terms::<2>(&mut f.c[jj], &f.lsp, jj),
            _ => column_terms::<1>(&mut f.c[jj], &f.lsp, jj),
        }
    }
    t = prof_add(9, t, on);
    // The block's own factorisation, right-looking: column kk finalised
    // (square root, then the divisions), then subtracted from the later
    // columns, so that element (ii, jj) receives its terms k = 20 + kk in
    // increasing kk and is divided after the last one.
    for kk in 0..N {
        let dkk = f.c[kk][kk].sqrt();
        f.c[kk][kk] = dkk;
        // lanes kk+1..20 (from the 4-aligned lane at or below kk+1; lanes
        // at or above the diagonal are scratch)
        let start = 4 * ((kk + 1) / 4);
        match (kk + 1) / 4 {
            0 => column_divide::<5>(&mut f.c[kk], dkk, start),
            1 => column_divide::<4>(&mut f.c[kk], dkk, start),
            2 => column_divide::<3>(&mut f.c[kk], dkk, start),
            3 => column_divide::<2>(&mut f.c[kk], dkk, start),
            4 => column_divide::<1>(&mut f.c[kk], dkk, start),
            _ => {}
        }
        // lane kk may have been in the tile: it is the diagonal, restore it.
        f.c[kk][kk] = dkk;
        let (done, rest) = f.c.split_at_mut(kk + 1);
        let ck = &done[kk];
        for (offset, cj) in rest.iter_mut().enumerate() {
            let jj = kk + 1 + offset;
            let ljk = ck[jj];
            match jj / 4 {
                0 => column_update::<5>(cj, ck, ljk, jj),
                1 => column_update::<4>(cj, ck, ljk, jj),
                2 => column_update::<3>(cj, ck, ljk, jj),
                3 => column_update::<2>(cj, ck, ljk, jj),
                _ => column_update::<1>(cj, ck, ljk, jj),
            }
        }
    }
    prof_add(10, t, on);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:146-148
// ```c
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
// ```
// Several elements of one column are advanced together; each gets its own
// `temp -= A[i][k] * A[j][k]` for k in increasing order.
/// `col[ii] -= lsp[k][ii] * lsp[k][jj]` for `k = 0..20` in order, over the
/// `4 * V` lanes from `4 * (jj / 4)` (the lanes above the diagonal and the
/// padding are scratch).
#[inline(always)]
fn column_terms<const V: usize>(col: &mut [f64; P], lsp: &[[f64; P]; N], jj: usize) {
    let start = 4 * (jj / 4);
    let mut acc = [0.0f64; 24];
    acc[..4 * V].copy_from_slice(&col[start..start + 4 * V]);
    for lk in lsp.iter() {
        let ljk = lk[jj];
        for l in 0..4 * V {
            acc[l] -= lk[start + l] * ljk;
        }
    }
    col[start..start + 4 * V].copy_from_slice(&acc[..4 * V]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:149-149
// ```c
//             A[i][j] = temp/A[j][j];
// ```
/// `col[l] /= d` over the `4 * V` lanes from `start`.
#[inline(always)]
fn column_divide<const V: usize>(col: &mut [f64; P], d: f64, start: usize) {
    for l in 0..4 * V {
        col[start + l] /= d;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:146-148
// ```c
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
// ```
// Right-looking form of the same subtraction: after column `kk` is final, term `k = kk` is
// subtracted from each later column. Each element still receives its terms in increasing k.
/// `cj[ii] -= ck[ii] * ljk` over the `4 * V` lanes from `4 * (jj / 4)`.
#[inline(always)]
fn column_update<const V: usize>(cj: &mut [f64; P], ck: &[f64; P], ljk: f64, jj: usize) {
    let start = 4 * (jj / 4);
    for l in 0..4 * V {
        cj[start + l] -= ck[start + l] * ljk;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-177,179-183
// ```c
//     for (i = 0;  i < n;  i++) {
//         temp = x[i];
//         for (j = 0;  j < i;  j++) {
//             temp -= L[i][j] * x[j];
//         }
//         x[i] = temp/L[i][i];
//     }
// ...
//     for (j = n - 1;  j >= 0;  j--) {
//         x[j] /= L[j][j];
//         for (i = 0;  i < j;  i++) {
//             x[i] -= L[j][i] * x[j];
//         }
// ```
/// `Nlm_SolveLtriangPosDef` with the structured factor.
///
/// Forward substitution by columns: `x[i]` receives `-= L[i][j] * x[j]` for
/// `j = 0, 1, ..., i - 1` in that order (the zero products of the diagonal
/// block included, literally) and is then divided by `L[i][i]`.
/// Back substitution: `x[i]` receives `-= L[j][i] * x[j]` for
/// `j = 39, 38, ..., i + 1` (each `x[j]` already final) and is then divided
/// by `L[i][i]`, as in the reference.
#[inline(always)]
fn solve_structured(x: &mut [f64; M], f: &Factor) {
    for j in 0..N {
        let xj = x[j] / f.d[j];
        x[j] = xj;
        for i in j + 1..N {
            x[i] -= 0.0 * xj;
        }
        let lj = &f.lsp[j];
        for r in 0..N {
            x[N + r] -= lj[r] * xj;
        }
    }
    for jj in 0..N {
        let cj = &f.c[jj];
        let xj = x[N + jj] / cj[jj];
        x[N + jj] = xj;
        for ii in jj + 1..N {
            x[N + ii] -= cj[ii] * xj;
        }
    }
    // Back substitution, descending dot products over column i of L.
    for ii in (0..N).rev() {
        let ci = &f.c[ii];
        let mut temp = x[N + ii];
        for jj in (ii + 1..N).rev() {
            temp -= ci[jj] * x[N + jj];
        }
        x[N + ii] = temp / ci[ii];
    }
    for i in (0..N).rev() {
        let li = &f.lsp[i];
        let mut temp = x[i];
        for r in (0..N).rev() {
            temp -= li[r] * x[N + r];
        }
        for j in (i + 1..N).rev() {
            temp -= 0.0 * x[j];
        }
        x[i] = temp / f.d[i];
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-758,762-768
// ```c
//     while (its <= maxits) {
//         /* Compute the residuals */
//         EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
//                             constrain_rel_entropy);
//         CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
//                            grads, row_sums, col_sums, x, z,
//                            constrain_rel_entropy, relative_entropy);
// ...
//         if ( !(rnorm > tol) ) {
//             /* We converged at the current iterate */
//             break;
//         } else {
//             /* we did not converge, so increment the iteration counter
//                and start a new iteration */
//             if (++its <= maxits) {
// ```
// The function body ports this loop for `alphsize == 20` and `constrain_rel_entropy`. The
// stopping test, the iteration limit and the status are the reference's.
/// `Blast_OptimizeTargetFrequencies` with `constrain_rel_entropy` set.
/// Returns 0 when the iteration converged (the reference's status), 1
/// otherwise; `x` holds the final iterate either way.
#[inline(always)]
fn optimize(
    x: &mut [f64; NN],
    q: &[f64; NN],
    row_sums: &[f64; N],
    col_sums: &[f64; N],
    relative_entropy: f64,
    ln_one: f64,
    first_iteration_shortcut: bool,
) -> i32 {
    let on = prof_enabled();
    let mut t = if on { tsc() } else { 0 };
    if on {
        PROF[6].fetch_add(1, Ordering::Relaxed);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:675-679
    // ```c
    //     for (i = 0;  i < alphsize;  i++) {
    //         for (j = 0;  j < alphsize;  j++) {
    //             k = i * alphsize + j;
    //             scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
    // ```
    // ComputeScoresFromProbs
    let mut scores = [0.0f64; NN];
    {
        let mut ratios = [0.0f64; NN];
        for i in 0..N {
            for j in 0..N {
                ratios[i * N + j] = q[i * N + j] / (row_sums[i] * col_sums[j]);
            }
        }
        crate::utils::x_logclone::log_slice(&ratios, &mut scores);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:747-748
    // ```c
    //     /* Use q as the initial value for x */
    //     memcpy(x, q, n * sizeof(double));
    // ```
    *x = *q;
    let mut z = [0.0f64; M];
    let mut resids_x = [0.0f64; NN];
    let mut resids_z = [0.0f64; M];
    let mut grads1 = [0.0f64; NN];
    let mut dinv = [0.0f64; NN];
    let mut workspace = [0.0f64; NN];
    let mut factor = Factor::new();
    let mut rnorm = 0.0f64;
    let mut its = 0usize;
    while its <= K_COMPO_ADJUST_ITERATION_LIMIT {
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:634-645
        // ```c
        //     values[0] = 0.0; values[1] = 0.0;
        //     for (k = 0;  k < alphsize * alphsize;  k++) {
        //         temp = log(x[k] / q[k]);
        //         values[0]   += x[k] * temp;
        //         grads[0][k]  = temp + 1;
        //         if (constrain_rel_entropy) {
        //             temp += scores[k];
        //             values[1]   += x[k] * temp;
        //             grads[1][k]  = temp + 1;
        // ```
        // `values[0]` is not used by the reference after this call, so it is not computed here.
        // ---- EvaluateReFunctions (values[0] is never read by the reference)
        let mut logs = [0.0f64; NN];
        if its == 0 && first_iteration_shortcut {
            logs = [ln_one; NN];
        } else {
            let mut ratios = [0.0f64; NN];
            for k in 0..NN {
                ratios[k] = x[k] / q[k];
            }
            crate::utils::x_logclone::log_slice(&ratios, &mut logs);
        }
        t = prof_add(0, t, on);
        if on {
            PROF[7].fetch_add(1, Ordering::Relaxed);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:277-279,286-286
        // ```c
        //         eta = z[2 * alphsize - 1];
        //         for (i = 0;  i < n;  i++) {
        //             resids_x[i] = -grads[0][i] + eta * grads[1][i];
        // ...
        //     MultiplyByAtranspose(1.0, resids_x, alphsize, 1.0, z);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-221
        // ```c
        //     for (i = 0;  i < alphsize;  i++) {
        //         for (j = 0;  j < alphsize;  j++) {
        //             k = i * alphsize + j;
        //             y[k] += alpha * x[j];
        //             if (i > 0) {
        //                 y[k] += alpha * x[i + alphsize - 1];
        //             }
        // ```
        // `DualResiduals` and `MultiplyByAtranspose`: the per-element operations below are these, in this order.
        // grads[0] = temp + 1, grads[1] = (temp + scores) + 1, values[1] += x * (temp + scores);
        // DualResiduals: resids_x = -grads[0] + eta * grads[1], then
        // MultiplyByAtranspose(1.0, resids_x, 1.0, z): += 1.0 * z[j], += 1.0 * z[i + 19].
        let eta = z[M - 1];
        let mut values1 = 0.0f64;
        for i in 0..N {
            let base = i * N;
            let zr = z[i + N - 1];
            for j in 0..N {
                let k = base + j;
                let temp = logs[k];
                let g0 = temp + 1.0;
                let temp = temp + scores[k];
                let g1 = temp + 1.0;
                grads1[k] = g1;
                let mut r = -g0 + eta * g1;
                r += 1.0 * z[j];
                if i > 0 {
                    r += 1.0 * zr;
                }
                resids_x[k] = r;
                values1 += x[k] * temp;
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:331-332
        // ```c
        //     DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
        //     norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);
        // ```
        let norm_resids_x = euclidean_norm(&resids_x);
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:243-249,172-172,177-177
        // ```c
        //     for (i = 0;  i < alphsize;  i++) {
        //         rA[i] = col_sums[i];
        //     }
        //     for (i = 1;  i < alphsize;  i++) {
        //         rA[i + alphsize - 1] = row_sums[i];
        //     }
        //     MultiplyByA(1.0, rA, alphsize, -1.0, x);
        // ...
        //             y[j] += alpha * x[i * alphsize + j];
        // ...
        //             y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        // ```
        // `MultiplyByA(1.0, rA, alphsize, -1.0, x)` adds `-1.0 * x` to each entry, first the column
        // sums over i, then the row sums over j, in the order of the loops.
        // ResidualsLinearConstraints: (col_sums, row_sums[1..]) - A x
        resids_z[..N].copy_from_slice(col_sums);
        for i in 1..N {
            resids_z[i + N - 1] = row_sums[i];
        }
        for i in 0..N {
            let xr = &x[i * N..i * N + N];
            for j in 0..N {
                resids_z[j] += -1.0 * xr[j];
            }
        }
        for i in 1..N {
            let xr = &x[i * N..i * N + N];
            let mut s = resids_z[i + N - 1];
            for &v in xr {
                s += -1.0 * v;
            }
            resids_z[i + N - 1] = s;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:337-344
        // ```c
        //         resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        //         norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize);
        //     } else {
        //         norm_resids_z = Nlm_EuclideanNorm(resids_z, 2 * alphsize - 1);
        //     }
        //     *rnorm =
        //         sqrt(norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z);
        // ```
        resids_z[M - 1] = relative_entropy - values1;
        let norm_resids_z = euclidean_norm(&resids_z);
        rnorm = (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt();
        t = prof_add(1, t, on);
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:762-768
        // ```c
        //         if ( !(rnorm > tol) ) {
        //             /* We converged at the current iterate */
        //             break;
        //         } else {
        //             /* we did not converge, so increment the iteration counter
        //                and start a new iteration */
        //             if (++its <= maxits) {
        // ```
        if !(rnorm > K_COMPO_ADJUST_ERR_TOLERANCE) {
            break;
        }
        its += 1;
        if its <= K_COMPO_ADJUST_ITERATION_LIMIT {
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:504-510,525-531,534-534
            // ```c
            //      if (constrain_rel_entropy) {
            //         double eta;             /* dual variable for the relative
            //                                    entropy constraint */
            //         eta = z[m - 1];
            //         for (i = 0;  i < n;  i++) {
            //             Dinv[i] = x[i] / (1 - eta);
            //         }
            // ...
            //         W[m - 1][m - 1] = 0.0;
            //         for (i = 0;  i < n;  i++) {
            //             workspace[i] = Dinv[i] * grad_re[i];
            //             W[m - 1][m - 1] += grad_re[i] * workspace[i];
            //         }
            //         MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
            // ...
            //     Nlm_FactorLtriangPosDef(W, m);
            // ```
            // The loops below compute `Dinv`, the sums that `ScaledSymmetricProductA` adds into `W`, and
            // the last row of `W`, in the same order per value.
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
            // ```c
            //     for (i = 0;  i < alphsize;  i++) {
            //         for (j = 0;  j < alphsize;  j++) {
            //             double dd;     /* an individual diagonal element */
            //             dd = diagonal[i * alphsize + j];
            //             W[j][j] += dd;
            //             if (i > 0) {
            //                 W[i + alphsize - 1][j] += dd;
            //                 W[i + alphsize - 1][i + alphsize - 1] += dd;
            //             }
            //         }
            // ```
            // ---- FactorReNewtonSystem
            // dinv = x / (1 - eta); W[j][j] += dinv over i in order;
            // W[i+19][i+19] += dinv over j in order; workspace = dinv * grad_re;
            // W[39][39] += grad_re * workspace in order; W[39][c] column / row sums.
            let mut wdiag = [0.0f64; N];
            let mut wrow = [0.0f64; N];
            let mut w39 = [0.0f64; M];
            let mut w39_39 = 0.0f64;
            // every dinv finite and positive (for structure_ok)
            let mut dinv_pos = true;
            let mut dinv_fin = true;
            for i in 0..N {
                let base = i * N;
                let mut srow = 0.0f64;
                let mut s39 = 0.0f64;
                for j in 0..N {
                    let k = base + j;
                    let dd = x[k] / (1.0 - eta);
                    dinv[k] = dd;
                    dinv_pos &= dd > 0.0;
                    dinv_fin &= dd < f64::INFINITY;
                    wdiag[j] += dd;
                    srow += dd;
                    let ws = dd * grads1[k];
                    workspace[k] = ws;
                    w39_39 += grads1[k] * ws;
                    w39[j] += 1.0 * ws;
                    s39 += 1.0 * ws;
                }
                if i > 0 {
                    wrow[i] = srow;
                    w39[i + N - 1] = s39;
                }
            }
            w39[M - 1] = w39_39;
            let w = Wparts {
                wdiag: &wdiag,
                wrow: &wrow,
                dinv: &dinv,
                w39: &w39,
            };
            t = prof_add(2, t, on);
            let ok = structure_ok(&w, dinv_pos & dinv_fin);
            t = prof_add(11, t, on);
            if ok {
                factor_structured(&w, &mut factor);
            } else {
                factor_literal(&w, &mut factor);
            }
            t = prof_add(3, t, on);

            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:575-578,581-583,587-587
            // ```c
            //     for (i = 0;  i < n;  i++) {
            //         workspace[i] = x[i] * Dinv[i];
            //     }
            //     MultiplyByA(1.0, z, alphsize, -1.0, workspace);
            // ...
            //         for (i = 0;  i < n;  i++) {
            //             z[m - 1] -= grad_re[i] * workspace[i];
            //         }
            // ...
            //     Nlm_SolveLtriangPosDef(z, m, W);
            // ```
            // ---- SolveReNewtonSystem(resids_x, resids_z)
            // workspace = resids_x * dinv; resids_z -= A workspace (column sums
            // over i, row sums over j); resids_z[39] -= grad_re . workspace.
            let mut s39 = resids_z[M - 1];
            for i in 0..N {
                let base = i * N;
                let mut srow = resids_z[i + N - 1];
                for j in 0..N {
                    let k = base + j;
                    let ws = resids_x[k] * dinv[k];
                    workspace[k] = ws;
                    resids_z[j] += -1.0 * ws;
                    srow += -1.0 * ws;
                    s39 -= grads1[k] * ws;
                }
                if i > 0 {
                    resids_z[i + N - 1] = srow;
                }
            }
            resids_z[M - 1] = s39;
            solve_structured(&mut resids_z, &factor);
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:593-602
            // ```c
            //     if (constrain_rel_entropy) {
            //         for(i = 0; i < n; i++) {
            //             x[i] += grad_re[i] * z[m - 1];
            //         }
            //     }
            //     MultiplyByAtranspose(1.0, x, alphsize, 1.0, z);
            //     for (i = 0;  i < n;  i++) {
            //         x[i] *= Dinv[i];
            //     }
            // ```
            // resids_x += grad_re * z[39]; MultiplyByAtranspose(1.0, resids_x, 1.0, z); *= dinv
            let z_re = resids_z[M - 1];
            for i in 0..N {
                let base = i * N;
                let zr = resids_z[i + N - 1];
                for j in 0..N {
                    let k = base + j;
                    let mut r = resids_x[k];
                    r += grads1[k] * z_re;
                    r += 1.0 * resids_z[j];
                    if i > 0 {
                        r += 1.0 * zr;
                    }
                    resids_x[k] = r * dinv[k];
                }
            }
            t = prof_add(4, t, on);

            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:781-785
            // ```c
            //                 alpha = Nlm_StepBound(x, n, resids_x, 1.0 / .95);
            //                 alpha *= 0.95;
            //                 Nlm_AddVectors(x, n, alpha, resids_x);
            //                 Nlm_AddVectors(z, m, alpha, resids_z);
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:230-238
            // ```c
            //     for (i = 0; i < n; i++) {
            //         double alpha_i;    /* a step to the boundary for the current i */
            //         alpha_i = -x[i] / step_x[i];
            //         if (alpha_i >= 0 && alpha_i < alpha) {
            //             alpha = alpha_i;
            //         }
            //     }
            //     return alpha;
            // ```
            // ---- Nlm_StepBound(x, 400, resids_x, 1 / 0.95) * 0.95, then Nlm_AddVectors
            let mut ratios = [0.0f64; NN];
            for k in 0..NN {
                ratios[k] = -x[k] / resids_x[k];
            }
            let mut alpha = 1.0 / 0.95;
            for &alpha_i in &ratios {
                if alpha_i >= 0.0 && alpha_i < alpha {
                    alpha = alpha_i;
                }
            }
            alpha *= 0.95;
            for k in 0..NN {
                x[k] += alpha * resids_x[k];
            }
            for k in 0..M {
                z[k] += alpha * resids_z[k];
            }
            t = prof_add(5, t, on);
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:789-797
    // ```c
    //     converged = 0;
    //     if (its <= maxits && rnorm <= tol) {
    //         /* Newton's iteration converged */
    //         if ( !constrain_rel_entropy || z[m - 1] < 1 ) {
    //             /* and the final iterate is a minimizer */
    //             converged = 1;
    //         }
    //     }
    //     status = converged ? 0 : 1;
    // ```
    let converged = its <= K_COMPO_ADJUST_ITERATION_LIMIT
        && rnorm <= K_COMPO_ADJUST_ERR_TOLERANCE
        && z[M - 1] < 1.0;
    // No NCBI counterpart: iteration counter; it does not change any value NCBI computes.
    crate::utils::xstats::add(&crate::utils::xstats::NEWTON_ITERS, its as u64);
    if converged {
        0
    } else {
        1
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the same `optimize` body compiled with AVX2 enabled so that the compiler may
// use wider vectors. Rust does not fuse `a * b + c` or reorder sums, so it computes the same IEEE-754
// operations as the body without the feature; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn optimize_avx2(
    x: &mut [f64; NN],
    q: &[f64; NN],
    row_sums: &[f64; N],
    col_sums: &[f64; N],
    relative_entropy: f64,
    ln_one: f64,
    first_iteration_shortcut: bool,
) -> i32 {
    optimize(
        x,
        q,
        row_sums,
        col_sums,
        relative_entropy,
        ln_one,
        first_iteration_shortcut,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the same `optimize` body compiled with AVX-512F enabled so that the compiler may
// use wider vectors. Rust does not fuse `a * b + c` or reorder sums, so it computes the same IEEE-754
// operations as the body without the feature; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn optimize_avx512(
    x: &mut [f64; NN],
    q: &[f64; NN],
    row_sums: &[f64; N],
    col_sums: &[f64; N],
    relative_entropy: f64,
    ln_one: f64,
    first_iteration_shortcut: bool,
) -> i32 {
    optimize(
        x,
        q,
        row_sums,
        col_sums,
        relative_entropy,
        ln_one,
        first_iteration_shortcut,
    )
}

// No NCBI counterpart: run-time choice of the vector width (LOSAT_X_NEWTONEXACT_LEVEL); it only selects which compilation of `optimize` runs; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn x86_level() -> u8 {
    static LEVEL: OnceLock<u8> = OnceLock::new();
    *LEVEL.get_or_init(|| {
        let forced = std::env::var("LOSAT_X_NEWTONEXACT_LEVEL").ok();
        match forced.as_deref() {
            Some("0") => 0,
            Some("2") => {
                if std::is_x86_feature_detected!("avx512f") {
                    2
                } else {
                    u8::from(std::is_x86_feature_detected!("avx2"))
                }
            }
            _ => u8::from(std::is_x86_feature_detected!("avx2")),
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// ```c
// int
// Blast_OptimizeTargetFrequencies(double x[],
// ...
//                                 int constrain_rel_entropy,
//                                 double relative_entropy,
//                                 double tol,
//                                 int maxits)
// ```
// The caller chooses between this function and the port of the reference (see `mode`).
/// The entry point: the same contract as the reference
/// `optimize_target_frequencies` with `constrain_rel_entropy == true`.
pub(crate) fn optimize_target_frequencies_exact(
    x: &mut [f64],
    q: &[f64],
    row_sums: &[f64; N],
    col_sums: &[f64; N],
    relative_entropy: f64,
) -> i32 {
    let x: &mut [f64; NN] = x.try_into().expect("400 target frequencies");
    let q: &[f64; NN] = q.try_into().expect("400 standard frequencies");
    // The shortcut needs x[k] / q[k] == 1.0 exactly, i.e. finite non-zero q.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:747-748,636-636
    // ```c
    //     /* Use q as the initial value for x */
    //     memcpy(x, q, n * sizeof(double));
    // ...
    //         temp = log(x[k] / q[k]);
    // ```
    let shortcut = q.iter().all(|&v| v.is_finite() && v != 0.0);
    let ln_one = ln_one();
    #[cfg(target_arch = "x86_64")]
    {
        match x86_level() {
            // SAFETY: the feature was detected at run time.
            2 => {
                return unsafe {
                    optimize_avx512(x, q, row_sums, col_sums, relative_entropy, ln_one, shortcut)
                }
            }
            1 => {
                return unsafe {
                    optimize_avx2(x, q, row_sums, col_sums, relative_entropy, ln_one, shortcut)
                }
            }
            _ => {}
        }
    }
    optimize(x, q, row_sums, col_sums, relative_entropy, ln_one, shortcut)
}

#[cfg(test)]
mod tests {
    use super::*;

    // No NCBI counterpart: random number generator of the test; it does not change any value NCBI computes.
    fn xorshift(state: &mut u64) -> f64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
    // ```c
    //     for (i = 0;  i < n;  i++) {
    //         for (j = 0;  j < i;  j++) {
    //             temp = A[i][j];
    //             for (k = 0;  k < j;  k++) {
    //                 temp -= A[i][k] * A[j][k];
    //             }
    //             A[i][j] = temp/A[j][j];
    // ...
    //         temp = A[i][i];
    //         for (k = 0;  k < i;  k++) {
    //             temp -= A[i][k] * A[i][k];
    //         }
    //         A[i][i] = sqrt(temp);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-177,179-183
    // ```c
    //     for (i = 0;  i < n;  i++) {
    //         temp = x[i];
    //         for (j = 0;  j < i;  j++) {
    //             temp -= L[i][j] * x[j];
    //         }
    //         x[i] = temp/L[i][i];
    //     }
    // ...
    //     for (j = n - 1;  j >= 0;  j--) {
    //         x[j] /= L[j][j];
    //         for (i = 0;  i < j;  i++) {
    //             x[i] -= L[j][i] * x[j];
    //         }
    // ```
    /// Reference `Nlm_FactorLtriangPosDef` + `Nlm_SolveLtriangPosDef` on a dense copy.
    fn reference_factor_solve(w: &[[f64; M]; M], b: &[f64; M]) -> ([[f64; M]; M], [f64; M]) {
        let mut a = *w;
        for i in 0..M {
            for j in 0..i {
                let mut temp = a[i][j];
                for k in 0..j {
                    temp -= a[i][k] * a[j][k];
                }
                a[i][j] = temp / a[j][j];
            }
            let mut temp = a[i][i];
            for k in 0..i {
                temp -= a[i][k] * a[i][k];
            }
            a[i][i] = temp.sqrt();
        }
        let mut x = *b;
        for i in 0..M {
            let mut temp = x[i];
            for j in 0..i {
                temp -= a[i][j] * x[j];
            }
            x[i] = temp / a[i][i];
        }
        for j in (0..M).rev() {
            x[j] /= a[j][j];
            for i in 0..j {
                x[i] -= a[j][i] * x[j];
            }
        }
        (a, x)
    }

    // No NCBI counterpart: random test input for `dinv` and the gradient; it does not change any value NCBI computes.
    fn random_parts(state: &mut u64, wide: bool) -> ([f64; NN], [f64; NN]) {
        let mut dinv = [0.0f64; NN];
        for v in dinv.iter_mut() {
            *v = if wide {
                2f64.powf(xorshift(state) * 40.0 - 30.0)
            } else {
                1.0e-4 + xorshift(state) * 0.05
            };
        }
        let mut grad = [0.0f64; NN];
        for v in grad.iter_mut() {
            *v = xorshift(state) * 4.0 - 2.0;
        }
        (dinv, grad)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:132-135,527-531
    // ```c
    //             W[j][j] += dd;
    //             if (i > 0) {
    //                 W[i + alphsize - 1][j] += dd;
    //                 W[i + alphsize - 1][i + alphsize - 1] += dd;
    // ...
    //             workspace[i] = Dinv[i] * grad_re[i];
    //             W[m - 1][m - 1] += grad_re[i] * workspace[i];
    //         }
    //         MultiplyByA(0.0, &W[m - 1][0], alphsize, 1.0, workspace);
    // ```
    // Builds `W` the way `ScaledSymmetricProductA` and `FactorReNewtonSystem` do, for the test.
    fn build_w(
        dinv: &[f64; NN],
        grad: &[f64; NN],
    ) -> ([f64; N], [f64; N], [f64; M], [[f64; M]; M]) {
        let mut wdiag = [0.0f64; N];
        let mut wrow = [0.0f64; N];
        for i in 0..N {
            for j in 0..N {
                wdiag[j] += dinv[i * N + j];
            }
        }
        for i in 1..N {
            let mut s = 0.0;
            for j in 0..N {
                s += dinv[i * N + j];
            }
            wrow[i] = s;
        }
        let mut ws = [0.0f64; NN];
        for k in 0..NN {
            ws[k] = dinv[k] * grad[k];
        }
        let mut w39 = [0.0f64; M];
        for i in 0..N {
            for j in 0..N {
                w39[j] += 1.0 * ws[i * N + j];
            }
        }
        for i in 1..N {
            for j in 0..N {
                w39[i + N - 1] += 1.0 * ws[i * N + j];
            }
        }
        let mut s = 0.0;
        for k in 0..NN {
            s += grad[k] * ws[k];
        }
        w39[M - 1] = s;
        let mut w = [[0.0f64; M]; M];
        for j in 0..N {
            w[j][j] = wdiag[j];
        }
        for i in 1..N {
            for j in 0..N {
                w[i + N - 1][j] = 0.0 + dinv[i * N + j];
            }
            w[i + N - 1][i + N - 1] = wrow[i];
        }
        w[M - 1] = w39;
        (wdiag, wrow, w39, w)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155,171-177
    // ```c
    //     for (i = 0;  i < n;  i++) {
    //         for (j = 0;  j < i;  j++) {
    //             temp = A[i][j];
    //             for (k = 0;  k < j;  k++) {
    //                 temp -= A[i][k] * A[j][k];
    //             }
    //             A[i][j] = temp/A[j][j];
    // ...
    //         temp = A[i][i];
    //         for (k = 0;  k < i;  k++) {
    //             temp -= A[i][k] * A[i][k];
    //         }
    //         A[i][i] = sqrt(temp);
    // ...
    //     for (i = 0;  i < n;  i++) {
    //         temp = x[i];
    //         for (j = 0;  j < i;  j++) {
    //             temp -= L[i][j] * x[j];
    //         }
    //         x[i] = temp/L[i][i];
    //     }
    // ```
    // The test factors and solves 400 random Newton matrices with the structured code and with a
    // copy of these loops, and compares every entry by its bits.
    #[test]
    fn x_structured_factor_and_solve_match_reference_bitwise_on_random_newton_matrices() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut structured = 0usize;
        for round in 0..400 {
            let (dinv, grad) = random_parts(&mut state, round % 2 == 1);
            let (wdiag, wrow, w39, w) = build_w(&dinv, &grad);
            let parts = Wparts {
                wdiag: &wdiag,
                wrow: &wrow,
                dinv: &dinv,
                w39: &w39,
            };
            let mut b = [0.0f64; M];
            for v in b.iter_mut() {
                *v = xorshift(&mut state) * 2.0 - 1.0;
            }
            let (l_ref, x_ref) = reference_factor_solve(&w, &b);
            let mut f = Factor::new();
            let dinv_ok = dinv.iter().all(|&v| v > 0.0 && v < f64::INFINITY);
            if structure_ok(&parts, dinv_ok) {
                structured += 1;
                factor_structured(&parts, &mut f);
            } else {
                factor_literal(&parts, &mut f);
            }
            for i in 0..M {
                for j in 0..=i {
                    assert_eq!(
                        f.get(i, j).to_bits(),
                        l_ref[i][j].to_bits(),
                        "round {round}: L[{i}][{j}] differs: {} vs {}",
                        f.get(i, j),
                        l_ref[i][j]
                    );
                }
            }
            let mut xs = b;
            solve_structured(&mut xs, &f);
            for i in 0..M {
                assert_eq!(
                    xs[i].to_bits(),
                    x_ref[i].to_bits(),
                    "round {round}: x[{i}] differs"
                );
            }
            // The literal path must agree too (it is the fallback).
            let mut g = Factor::new();
            factor_literal(&parts, &mut g);
            for i in 0..M {
                for j in 0..=i {
                    assert_eq!(g.get(i, j).to_bits(), l_ref[i][j].to_bits());
                }
            }
        }
        assert!(
            structured > 300,
            "structured path exercised {structured} times"
        );
    }
}
