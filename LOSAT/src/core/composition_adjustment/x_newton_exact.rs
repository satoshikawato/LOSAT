//! `Blast_OptimizeTargetFrequencies` for the relative-entropy-constrained
//! problem, bit for bit.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:602-713
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:762-768,774-777,781-785
//! The Newton loop of `optimize` below is this loop. Each call here is a step of one
//! iteration: residuals, convergence test, factor, solve, step bound, update of x and z.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,171-177
//! These are the Cholesky factorisation and the forward solve whose per-element operations
//! and order are kept (see the list above). The rest of this module changes only the data layout
//! and which independent elements are advanced together.

use super::adjust_scores::{
    COMPO_NUM_TRUE_AA, K_COMPO_ADJUST_ERR_TOLERANCE, K_COMPO_ADJUST_ITERATION_LIMIT,
};
use crate::utils::glibc_fma_libm::GlibcLibm;
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
    if !crate::utils::xdiag::ON {
        return false;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xdiag::var_os("LOSAT_X_NEWTONEXACT_PROF").is_some())
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
#[cfg(test)]
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
// `N` is `alphsize`, `NN` is `n`, and `M` is `m` for `constrain_rel_entropy` set. This module
// handles only that case; other cases run the reference.
const N: usize = COMPO_NUM_TRUE_AA;
const NN: usize = N * N;
/// Number of constraints with the relative-entropy constraint (`m`).
const M: usize = 2 * N;
/// Lane padding for the 20-row vectors of the trailing block.
const P: usize = 24;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:747-748,635-636
// NCBI starts from `x = q`, so in the first iteration every `log(x[k] / q[k])` is `log(1.0)`.
// This is the value of one such call, made once per process.
/// `ln(1.0)` as the libm of this process computes it (not constant-folded).
pub(super) fn ln_one() -> f64 {
    static LN_ONE: OnceLock<f64> = OnceLock::new();
    *LN_ONE.get_or_init(|| std::hint::black_box(1.0f64).glibc_ln())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:193-208
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
// The parts are the sums that `ScaledSymmetricProductA` and `FactorReNewtonSystem` put into
// the lower triangle of `W` (column sums `W[j][j]`, row sums `W[i+19][i+19]`, the entries
// `W[i+19][j]`, and the last row).
/// The lower triangle of `W` as the reference holds it, from its parts.
pub(super) struct Wparts<'a> {
    /// `W[j][j]`, `j < 20`: column sums of `dinv`.
    pub(super) wdiag: &'a [f64; N],
    /// `W[i + 19][i + 19]`, `1 <= i < 20`: row sums of `dinv` (index `i`).
    pub(super) wrow: &'a [f64; N],
    /// `W[i + 19][j] = 0.0 + dinv[i * 20 + j]`.
    pub(super) dinv: &'a [f64; NN],
    /// `W[39][c]`, `c < 40`.
    pub(super) w39: &'a [f64; M],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
// The entries of `A` after this function, kept in blocks. Only the storage differs.
/// The Cholesky factor in structured storage.
pub(super) struct Factor {
    /// `L[j][j]` for `j < 20`.
    pub(super) d: [f64; N],
    /// `lsp[k][r] = L[20 + r][k]` for `k < 20`, `r < 20` (rows 20..39).
    pub(super) lsp: [[f64; P]; N],
    /// `c[jj][ii] = L[20 + ii][20 + jj]` for `ii >= jj` (the trailing block,
    /// column-major).  Entries with `ii < jj` and the padding are scratch.
    pub(super) c: [[f64; P]; N],
}

// No NCBI counterpart: constructor and test accessor of the structured storage; it does not change any value NCBI computes.
impl Factor {
    #[inline(always)]
    pub(super) fn new() -> Self {
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
// The loops below are this function on a dense 40 x 40 copy of `W`, in the same order.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:121-125,132-135,525-531
// `W` is filled as `ScaledSymmetricProductA` and `FactorReNewtonSystem` fill it.
/// `Nlm_FactorLtriangPosDef` on `W`, literally (dense lower triangle).
#[inline(never)]
pub(super) fn factor_literal(w: &Wparts, out: &mut Factor) {
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
pub(super) fn structure_ok(w: &Wparts, dinv_ok: bool) -> bool {
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
// Element `(i, j)` receives the same terms in the same order and the same final division or
// square root; only the terms known to be zero are left out (checked by `structure_ok`).
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
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
/// `col[l] /= d` over the `4 * V` lanes from `start`.
#[inline(always)]
fn column_divide<const V: usize>(col: &mut [f64; P], d: f64, start: usize) {
    for l in 0..4 * V {
        col[start + l] /= d;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:146-148
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-221
        // `DualResiduals` and `MultiplyByAtranspose`: the per-element operations below are these, in this order.
        // `grads[0]` is `temp` plus 1, `grads[1]` is the sum of `temp` and `scores` plus 1, and
        // `values[1]` grows by `x` times the sum of `temp` and `scores`.
        // `DualResiduals` sets each entry of `resids_x` to minus `grads[0]` plus `eta` times
        // `grads[1]`, then `MultiplyByAtranspose` with the scalars 1.0 and 1.0, `resids_x` and `z`
        // adds 1.0 times `z[j]` and, in a row `i` above 0, 1.0 times `z` at `i` plus 19.
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
        let norm_resids_x = euclidean_norm(&resids_x);
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:243-249,172-172,177-177
        // `MultiplyByA` with the scalars 1.0 and -1.0, `rA`, `alphsize` and `x` adds -1.0 times
        // `x` to each entry, first the column sums over i, then the row sums over j, in the order of
        // the loops.
        // ResidualsLinearConstraints: (col_sums, row_sums[1..]) - A x
        resids_z[..N].copy_from_slice(col_sums);
        for i in 1..N {
            resids_z[i + N - 1] = row_sums[i];
        }
        for i in 0..N {
            let xr = &x[i * N..i * N + N];
            for j in 0..N {
                resids_z[j] += -xr[j];
            }
        }
        for i in 1..N {
            let xr = &x[i * N..i * N + N];
            let mut s = resids_z[i + N - 1];
            for &v in xr {
                s += -v;
            }
            resids_z[i + N - 1] = s;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:337-344
        resids_z[M - 1] = relative_entropy - values1;
        let norm_resids_z = euclidean_norm(&resids_z);
        rnorm = (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt();
        t = prof_add(1, t, on);
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:762-768
        #[expect(
            clippy::neg_cmp_op_on_partial_ord,
            reason = "NCBI's convergence test `!(rnorm > tol)` stops on a NaN rnorm (optimize_target_freq.c:760-762); `<=` would not"
        )]
        if !(rnorm > K_COMPO_ADJUST_ERR_TOLERANCE) {
            break;
        }
        its += 1;
        if its <= K_COMPO_ADJUST_ITERATION_LIMIT {
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:504-510,525-531,534-534
            // The loops below compute `Dinv`, the sums that `ScaledSymmetricProductA` adds into `W`, and
            // the last row of `W`, in the same order per value.
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
            // ---- FactorReNewtonSystem
            // `dinv = x / (1 - eta)`; `W[j][j] += dinv` over i in order;
            // `W[i+19][i+19] += dinv` over j in order; `workspace` is set to `dinv` times `grad_re`;
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
                    resids_z[j] += -ws;
                    srow += -ws;
                    s39 -= grads1[k] * ws;
                }
                if i > 0 {
                    resids_z[i + N - 1] = srow;
                }
            }
            resids_z[M - 1] = s39;
            solve_structured(&mut resids_z, &factor);
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:593-602
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:230-238
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
    // No NCBI counterpart (memory placement only): see `x_frame_align`.
    crate::utils::x_frame::x_frame_align();
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
pub(super) fn x86_level() -> u8 {
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
            // SAFETY: `x86_level` returns 1 only when `is_x86_feature_detected!("avx2")` was true,
            // which the `#[target_feature(enable = "avx2")]` of `optimize_avx2` needs.
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
#[path = "x_newton_exact_tests.rs"]
mod tests;
