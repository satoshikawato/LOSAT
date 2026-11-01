//! `Blast_OptimizeTargetFrequencies` for four independent problems at once,
//! one problem per vector lane, bit for bit.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-758,762-768,774-777,781-785
//! ```c
//!     while (its <= maxits) {
//!         /* Compute the residuals */
//!         EvaluateReFunctions(values, grads, alphsize, x, q, old_scores,
//!                             constrain_rel_entropy);
//!         CalculateResiduals(&rnorm, resids_x, alphsize, resids_z, values,
//!                            grads, row_sums, col_sums, x, z,
//!                            constrain_rel_entropy, relative_entropy);
//! ...
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
//! Each lane runs this loop for its own problem with exactly the operation
//! sequence of `x_newton_exact::optimize` (which is the reference's sequence
//! of IEEE-754 operations per value): the same operations on the same
//! operands in the same order for every value, including the literal
//! `1.0 *`, `-1.0 *` and `0.0 *` products, the sequential sums, the
//! Euclidean norms (both branches evaluated, then selected per lane), the
//! sequential step bound, the stopping test, the iteration count and the
//! status rule.  The libm `log` is called per value through
//! `x_logclone::log_slice` on a contiguous buffer per lane.  What changes is
//! only which independent values are computed together: value `k` of four
//! different problems sits in one `[f64; 4]`.  Rust never fuses `a * b + c`
//! and never re-associates sums; SIMD add / sub / mul / div / sqrt are
//! IEEE-754 per element, so a lane gives the scalar result.
//!
//! On aarch64 (no variable; `LOSAT_X_NEWTONLANES8` graduated 2026-10-22): eight
//! problems per round, value `k` of eight problems in one `[f64; 8]` (`F8`), so
//! each operation is four 128-bit NEON operations on independent values.  The
//! lane body is the same generic code; each lane still runs the operation
//! sequence above.  When the queue is empty and at most four lanes hold a
//! problem, those problems move to a four-lane state (bit copies of the values
//! that persist between rounds) and finish there.  On x86-64 and Wasm the lanes
//! stay four wide (on x86-64 the two 256-bit
//! operations per value gave no gain; Wasm is not measured); the unit tests run
//! the eight-lane body on every target.
//!
//! No NCBI counterpart: batching of independent problems, the result store
//! and its lookup by bitwise-equal input. NCBI solves one problem per call;
//! a result is used only for a call whose input (q, row and column
//! probabilities, relative entropy) has the same bits, so it does not change
//! any value NCBI computes.

use super::adjust_scores::{
    COMPO_NUM_TRUE_AA, K_COMPO_ADJUST_ERR_TOLERANCE, K_COMPO_ADJUST_ITERATION_LIMIT,
};
use super::x_newton_exact::{
    factor_literal, ln_one, optimize_target_frequencies_exact, structure_ok, Factor, Wparts,
};
use std::cell::RefCell;
use std::collections::VecDeque;
use std::ops::{Index, IndexMut};
use std::sync::atomic::{AtomicU64, Ordering};

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:725-727
// ```c
//     n  = alphsize * alphsize;
//     mA = 2 * alphsize - 1;
//     m  = constrain_rel_entropy ? mA + 1 : mA;
// ```
// `N` is `alphsize`, `NN` is `n`, `M` is `m` with `constrain_rel_entropy` set (the only case here).
const N: usize = COMPO_NUM_TRUE_AA;
const NN: usize = N * N;
const M: usize = 2 * N;
// No NCBI counterpart: batching constants (problems per vector, smallest batch for the vector
// kernel, results kept per thread); they do not change any value NCBI computes.
const LANES: usize = 4;
const MIN_LANE_BATCH: usize = 3;
const STORE_CAPACITY: usize = 512;
const NO_PROBLEM: usize = usize::MAX;
// No NCBI counterpart: the eight-lane width (aarch64), the smallest batch that uses
// it, and the number of lanes with a problem at or below which an eight-lane batch whose queue is
// empty moves to four lanes; they do not change any value NCBI computes.
const LANES8: usize = 8;
const MAX_LANES: usize = LANES8;
const MIN_LANE8_BATCH: usize = LANES + 1;
const HAND_OVER_LANES: usize = LANES;

// No NCBI counterpart: on aarch64 the lanes share a round among eight problems; on x86-64 and
// Wasm the lanes stay four wide (x86-64 AVX2 measured no gain). It does not change any value NCBI
// computes.
/// Whether batches of `MIN_LANE8_BATCH` or more problems use eight lanes.
fn eight_lanes() -> bool {
    cfg!(target_arch = "aarch64")
}

// No NCBI counterpart: counters of the LOSAT_X_STATS report; they do not change any value NCBI
// computes.
static PREFETCHED: AtomicU64 = AtomicU64::new(0);
static HITS: AtomicU64 = AtomicU64::new(0);
static MISSES: AtomicU64 = AtomicU64::new(0);
static BATCHES: AtomicU64 = AtomicU64::new(0);
static LANE_ROUNDS: AtomicU64 = AtomicU64::new(0);
static BUSY_LANE_ROUNDS: AtomicU64 = AtomicU64::new(0);
static BATCHES8: AtomicU64 = AtomicU64::new(0);
static ROUNDS8: AtomicU64 = AtomicU64::new(0);
static TAIL_ROUNDS: AtomicU64 = AtomicU64::new(0);
static DROPPED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: prints the counters at exit; it does not change any value NCBI computes.
pub(crate) fn print_stats() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let prefetched = PREFETCHED.load(Ordering::Relaxed);
        let hits = HITS.load(Ordering::Relaxed);
        eprintln!(
            "[X_STATS] NEWTONLANES prefetched={} hits={} misses={} unused={} batches={} lane_rounds={} busy_lane_rounds={}",
            prefetched,
            hits,
            MISSES.load(Ordering::Relaxed),
            prefetched.saturating_sub(hits),
            BATCHES.load(Ordering::Relaxed),
            LANE_ROUNDS.load(Ordering::Relaxed),
            BUSY_LANE_ROUNDS.load(Ordering::Relaxed),
        );
        if prefetched > 0 {
            eprintln!(
                "[X_STATS] NEWTONDROP dropped={}",
                DROPPED.load(Ordering::Relaxed)
            );
        }
        if eight_lanes() {
            eprintln!(
                "[X_STATS] NEWTONLANES8 batches8={} rounds8={} tail_rounds={}",
                BATCHES8.load(Ordering::Relaxed),
                ROUNDS8.load(Ordering::Relaxed),
                TAIL_ROUNDS.load(Ordering::Relaxed),
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1385-1399
// ```c
//     Blast_ApplyPseudocounts(row_probs, length1,
//                             NRrecord->first_standard_freq, pseudocounts);
//     Blast_ApplyPseudocounts(col_probs, length2,
//                             NRrecord->second_standard_freq, pseudocounts);
//
//     status =
//         Blast_OptimizeTargetFrequencies(&NRrecord->mat_final[0][0],
//                                         COMPO_NUM_TRUE_AA,
//                                         &iteration_count,
//                                         &NRrecord->mat_b[0][0],
//                                         row_probs, col_probs,
//                                         (desired_re > 0.0),
//                                         desired_re,
// ```
// The arguments of this call that determine its result (`q = mat_b`, `row_probs`, `col_probs`,
// `desired_re` with `desired_re > 0.0`).
/// Input of one relative-entropy-constrained problem (alphabet of 20).
#[derive(Clone)]
pub(crate) struct NewtonInput {
    pub(crate) q: [f64; NN],
    pub(crate) row: [f64; N],
    pub(crate) col: [f64; N],
    pub(crate) relative_entropy: f64,
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
/// Result of one problem: the final iterate `x` and the status (0 converged).
#[derive(Clone)]
pub(crate) struct NewtonOutput {
    pub(crate) x: [f64; NN],
    pub(crate) status: i32,
}

impl NewtonOutput {
    // No NCBI counterpart: empty output slot; it does not change any value NCBI computes.
    pub(crate) fn new() -> Self {
        Self {
            x: [0.0; NN],
            status: 0,
        }
    }
}

// No NCBI counterpart: a vector of f64 values, one per problem, with element-wise IEEE-754
// operations written as plain scalar Rust operations (the compiler may put them in vector
// registers), and its per-lane comparison result (all bits set or clear). Each lane gets exactly
// the scalar operation. It does not change any value NCBI computes.
trait LaneVec: Copy + Index<usize, Output = f64> + IndexMut<usize> {
    /// Problems per value.
    const L: usize;
    type Mask: LaneMask;
    fn splat(v: f64) -> Self;
    fn add(self, o: Self) -> Self;
    fn sub(self, o: Self) -> Self;
    fn mul(self, o: Self) -> Self;
    fn div(self, o: Self) -> Self;
    fn lt(self, o: Self) -> Self::Mask;
    fn gt(self, o: Self) -> Self::Mask;
    fn ge(self, o: Self) -> Self::Mask;
    fn ne(self, o: Self) -> Self::Mask;
    fn neg(self) -> Self;
    fn abs(self) -> Self;
    fn sqrt(self) -> Self;
    /// Per lane: `m ? a : b` (bit select; no arithmetic on either value).
    fn select(m: Self::Mask, a: Self, b: Self) -> Self;
}

// No NCBI counterpart: per-lane flags of a `LaneVec` comparison; it does not change any value NCBI
// computes.
trait LaneMask: Copy {
    fn all() -> Self;
    fn and(self, o: Self) -> Self;
    fn any(self) -> bool;
    fn lane(self, l: usize) -> bool;
}

// No NCBI counterpart: a comparison result as all bits set / clear.
#[inline(always)]
fn lane_mask(b: bool) -> u64 {
    (b as u64).wrapping_neg()
}

// No NCBI counterpart: defines a `LaneVec` (`$f`, `$n` values) and its mask (`$m`) with every
// element operation written out per lane index; it does not change any value NCBI computes.
macro_rules! lane_vec {
    ($f:ident, $m:ident, $n:literal, $align:literal, [$($i:literal),+]) => {
        #[derive(Clone, Copy)]
        #[repr(C, align($align))]
        struct $f([f64; $n]);

        #[derive(Clone, Copy)]
        #[repr(C, align($align))]
        struct $m([u64; $n]);

        impl Index<usize> for $f {
            type Output = f64;
            #[inline(always)]
            fn index(&self, l: usize) -> &f64 {
                &self.0[l]
            }
        }

        impl IndexMut<usize> for $f {
            #[inline(always)]
            fn index_mut(&mut self, l: usize) -> &mut f64 {
                &mut self.0[l]
            }
        }

        impl LaneVec for $f {
            const L: usize = $n;
            type Mask = $m;
            #[inline(always)]
            fn splat(v: f64) -> $f {
                $f([v; $n])
            }
            #[inline(always)]
            fn add(self, o: $f) -> $f {
                let (a, b) = (self.0, o.0);
                $f([$(a[$i] + b[$i]),+])
            }
            #[inline(always)]
            fn sub(self, o: $f) -> $f {
                let (a, b) = (self.0, o.0);
                $f([$(a[$i] - b[$i]),+])
            }
            #[inline(always)]
            fn mul(self, o: $f) -> $f {
                let (a, b) = (self.0, o.0);
                $f([$(a[$i] * b[$i]),+])
            }
            #[inline(always)]
            fn div(self, o: $f) -> $f {
                let (a, b) = (self.0, o.0);
                $f([$(a[$i] / b[$i]),+])
            }
            #[inline(always)]
            fn lt(self, o: $f) -> $m {
                let (a, b) = (self.0, o.0);
                $m([$(lane_mask(a[$i] < b[$i])),+])
            }
            #[inline(always)]
            fn gt(self, o: $f) -> $m {
                let (a, b) = (self.0, o.0);
                $m([$(lane_mask(a[$i] > b[$i])),+])
            }
            #[inline(always)]
            fn ge(self, o: $f) -> $m {
                let (a, b) = (self.0, o.0);
                $m([$(lane_mask(a[$i] >= b[$i])),+])
            }
            #[inline(always)]
            fn ne(self, o: $f) -> $m {
                let (a, b) = (self.0, o.0);
                $m([$(lane_mask(a[$i] != b[$i])),+])
            }
            #[inline(always)]
            fn neg(self) -> $f {
                let a = self.0;
                $f([$(-a[$i]),+])
            }
            #[inline(always)]
            fn abs(self) -> $f {
                let a = self.0;
                $f([$(a[$i].abs()),+])
            }
            #[inline(always)]
            fn sqrt(self) -> $f {
                let a = self.0;
                $f([$(a[$i].sqrt()),+])
            }
            #[inline(always)]
            fn select(m: $m, a: $f, b: $f) -> $f {
                let (m, a, b) = (m.0, a.0, b.0);
                let pick =
                    |l: usize| f64::from_bits((a[l].to_bits() & m[l]) | (b[l].to_bits() & !m[l]));
                $f([$(pick($i)),+])
            }
        }

        impl LaneMask for $m {
            #[inline(always)]
            fn all() -> $m {
                $m([u64::MAX; $n])
            }
            #[inline(always)]
            fn and(self, o: $m) -> $m {
                let (a, b) = (self.0, o.0);
                $m([$(a[$i] & b[$i]),+])
            }
            #[inline(always)]
            fn any(self) -> bool {
                (0 $(| self.0[$i])+) != 0
            }
            #[inline(always)]
            fn lane(self, l: usize) -> bool {
                self.0[l] != 0
            }
        }
    };
}

// No NCBI counterpart: four problems per value (one 256-bit or two 128-bit vectors) and, for
// eight on aarch64; it does not change any value NCBI computes.
lane_vec!(F4, M4, 4, 32, [0, 1, 2, 3]);
lane_vec!(F8, M8, 8, 64, [0, 1, 2, 3, 4, 5, 6, 7]);

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:706-722
// ```c
//     double * scores;       /* the scores of the old matrix */
// ...
//     double * resids_x;     /* the residuals in the x variables */
//     double * resids_z;     /* the residuals in the z variables */
// ...
//     double *z;             /* the dual variables */
// ```
// The working arrays of `Blast_OptimizeTargetFrequencies` (and of `x_newton_exact::optimize`) for
// `V::L` problems, value `k` of lane `l` at `[k][l]`. `d`, `lsp`, `c` are the Cholesky factor in the
// structured storage of `x_newton_exact::Factor` (lower triangle of `c` only).
/// Lane-innermost working storage of four (`F4`) or eight (`F8`) problems.
struct LaneState<V> {
    x: [V; NN],
    q: [V; NN],
    scores: [V; NN],
    logs: [V; NN],
    grads1: [V; NN],
    dinv: [V; NN],
    resids_x: [V; NN],
    z: [V; M],
    resids_z: [V; M],
    row: [V; N],
    col: [V; N],
    re: V,
    d: [V; N],
    lsp: [[V; N]; N],
    c: [[V; N]; N],
}

// No NCBI counterpart: zero-filled heap allocation of the working storage; it does not change any
// value NCBI computes.
fn new_lane_state<V: LaneVec>() -> Box<LaneState<V>> {
    let mut state = Box::<LaneState<V>>::new_uninit();
    // SAFETY: `LaneState` holds only f64 arrays (through `F4` or `F8`, the only `LaneVec` types),
    // for which the all-zero byte pattern is the valid value 0.0; every byte is written before
    // `assume_init`.
    unsafe {
        std::ptr::write_bytes(state.as_mut_ptr(), 0, 1);
        state.assume_init()
    }
}

// No NCBI counterpart: test-only switch that makes one lane use the literal dense factorisation in
// one round (the fallback path); it does not change any value NCBI computes.
#[cfg(test)]
thread_local! {
    static FORCE_LITERAL: std::cell::Cell<Option<(usize, usize)>> = const { std::cell::Cell::new(None) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:675-679,747-748
// ```c
//     for (i = 0;  i < alphsize;  i++) {
//         for (j = 0;  j < alphsize;  j++) {
//             k = i * alphsize + j;
//             scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
// ...
//     /* Use q as the initial value for x */
//     memcpy(x, q, n * sizeof(double));
// ```
// `ComputeScoresFromProbs`, `x = q` and `z = 0` for one problem in lane `l`, as at the start of
// `x_newton_exact::optimize`. Returns whether every `q` is finite and non-zero (then the first
// iteration's quotients `x[k] / q[k]` are exactly 1.0).
#[inline(always)]
fn load_lane<V: LaneVec>(st: &mut LaneState<V>, l: usize, input: &NewtonInput) -> bool {
    let mut ratios = [0.0f64; NN];
    let mut scores = [0.0f64; NN];
    for i in 0..N {
        for j in 0..N {
            ratios[i * N + j] = input.q[i * N + j] / (input.row[i] * input.col[j]);
        }
    }
    crate::utils::x_logclone::log_slice(&ratios, &mut scores);
    for k in 0..NN {
        st.scores[k][l] = scores[k];
        st.q[k][l] = input.q[k];
        st.x[k][l] = input.q[k];
    }
    for k in 0..M {
        st.z[k][l] = 0.0;
    }
    for i in 0..N {
        st.row[i][l] = input.row[i];
        st.col[i][l] = input.col[i];
    }
    st.re[l] = input.relative_entropy;
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:747-748,636-636
    // ```c
    //     /* Use q as the initial value for x */
    //     memcpy(x, q, n * sizeof(double));
    // ...
    //         temp = log(x[k] / q[k]);
    // ```
    input.q.iter().all(|&v| v.is_finite() && v != 0.0)
}

// No NCBI counterpart: a lane without a problem repeats the state of a running lane (so that it
// computes the same finite values and is discarded); it does not change any value NCBI computes.
#[inline(always)]
fn mirror_lane<V: LaneVec>(st: &mut LaneState<V>, l: usize, src: usize) {
    for k in 0..NN {
        st.x[k][l] = st.x[k][src];
        st.q[k][l] = st.q[k][src];
        st.scores[k][l] = st.scores[k][src];
    }
    for k in 0..M {
        st.z[k][l] = st.z[k][src];
    }
    for i in 0..N {
        st.row[i][l] = st.row[i][src];
        st.col[i][l] = st.col[i][src];
    }
    st.re[l] = st.re[src];
}

// No NCBI counterpart: moves the problem of lane `sl` of an eight-lane state to lane `dl` of a
// four-lane state between rounds: bit copies of every value that persists from one round to the
// next (`x`, `q`, `scores`, `z`, the probabilities and the relative entropy; the other arrays are
// recomputed in each round before they are read). The problem continues with the same operations;
// it does not change any value NCBI computes.
#[inline(always)]
fn move_lane<A: LaneVec, B: LaneVec>(
    dst: &mut LaneState<A>,
    dl: usize,
    src: &LaneState<B>,
    sl: usize,
) {
    for k in 0..NN {
        dst.x[k][dl] = src.x[k][sl];
        dst.q[k][dl] = src.q[k][sl];
        dst.scores[k][dl] = src.scores[k][sl];
    }
    for k in 0..M {
        dst.z[k][dl] = src.z[k][sl];
    }
    for i in 0..N {
        dst.row[i][dl] = src.row[i][sl];
        dst.col[i][dl] = src.col[i][sl];
    }
    dst.re[dl] = src.re[sl];
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
// `Nlm_EuclideanNorm` per lane. The branch taken by each lane is chosen by a select: when no lane
// has a new maximum, only the `else` branch is evaluated (lanes with `v[i] == 0.0` keep `sum`);
// otherwise both branches are evaluated and selected per lane. `sum * (scale/absvi) * (scale/absvi)`
// is `(sum * r) * r`; `v[i] != 0.0` is true for NaN, and `scale < NaN` is false, as in C.
/// `Nlm_EuclideanNorm` of `V::L` vectors.
#[inline(always)]
fn euclidean_norm<V: LaneVec>(v: &[V]) -> V {
    let mut sum = V::splat(1.0);
    let mut scale = V::splat(0.0);
    for &value in v {
        euclidean_norm_step(&mut sum, &mut scale, value);
    }
    scale.mul(sum.sqrt())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:196-206
// ```c
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
// ```
/// One element of `Nlm_EuclideanNorm` (the body of the loop), `V::L` vectors.
#[inline(always)]
fn euclidean_norm_step<V: LaneVec>(sum: &mut V, scale: &mut V, value: V) {
    let zero = V::splat(0.0);
    let one = V::splat(1.0);
    let abs_value = value.abs();
    let nonzero = value.ne(zero);
    let larger = nonzero.and(scale.lt(abs_value));
    if larger.any() {
        let r_new = scale.div(abs_value);
        let sum_new = one.add(sum.mul(r_new).mul(r_new));
        let r_old = abs_value.div(*scale);
        let sum_old = sum.add(r_old.mul(r_old));
        *sum = V::select(nonzero, V::select(larger, sum_new, sum_old), *sum);
        *scale = V::select(larger, abs_value, *scale);
    } else {
        let r_old = abs_value.div(*scale);
        *sum = V::select(nonzero, sum.add(r_old.mul(r_old)), *sum);
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
// `x_newton_exact::solve_structured` per lane: forward substitution by columns (the zero products
// of the diagonal block included), back substitution as descending dot products; each element
// receives the same terms in the same order and the same final division.
/// `Nlm_SolveLtriangPosDef` with the structured factor, `V::L` problems.
#[inline(always)]
fn solve_factored<V: LaneVec>(x: &mut [V; M], d: &[V; N], lsp: &[[V; N]; N], c: &[[V; N]; N]) {
    let zero = V::splat(0.0);
    for j in 0..N {
        let xj = x[j].div(d[j]);
        x[j] = xj;
        for i in j + 1..N {
            x[i] = x[i].sub(zero.mul(xj));
        }
        let lj = &lsp[j];
        for r in 0..N {
            x[N + r] = x[N + r].sub(lj[r].mul(xj));
        }
    }
    for jj in 0..N {
        let cj = &c[jj];
        let xj = x[N + jj].div(cj[jj]);
        x[N + jj] = xj;
        for ii in jj + 1..N {
            x[N + ii] = x[N + ii].sub(cj[ii].mul(xj));
        }
    }
    for ii in (0..N).rev() {
        let ci = &c[ii];
        let mut temp = x[N + ii];
        for jj in (ii + 1..N).rev() {
            temp = temp.sub(ci[jj].mul(x[N + jj]));
        }
        x[N + ii] = temp.div(ci[ii]);
    }
    for i in (0..N).rev() {
        let li = &lsp[i];
        let mut temp = x[i];
        for r in (0..N).rev() {
            temp = temp.sub(li[r].mul(x[N + r]));
        }
        for j in (i + 1..N).rev() {
            temp = temp.sub(zero.mul(x[j]));
        }
        x[i] = temp.div(d[i]);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:146-148
// ```c
//             for (k = 0;  k < j;  k++) {
//                 temp -= A[i][k] * A[j][k];
//             }
// ```
// `x_newton_exact::column_terms` per lane: `B` elements of column `jj` from row `start`, each receiving
// `temp -= L[20+ii][k] * L[20+jj][k]` for k = 0..19 in increasing order (kept in registers).
/// `col[start + t] -= lsp[k][start + t] * lsp[k][jj]` for `k = 0..20`, `t < B`.
#[inline(always)]
fn column_terms_block<const B: usize, V: LaneVec>(
    col: &mut [V; N],
    lsp: &[[V; N]; N],
    jj: usize,
    start: usize,
) {
    let mut acc: [V; B] = col[start..start + B].try_into().expect("column block");
    for lk in lsp.iter() {
        let ljk = lk[jj];
        let segment: &[V; B] = lk[start..start + B].try_into().expect("column block");
        for t in 0..B {
            acc[t] = acc[t].sub(segment[t].mul(ljk));
        }
    }
    col[start..start + B].copy_from_slice(&acc);
}

// No NCBI counterpart: `structure_ok` reads only `wdiag` and `w39` of its argument (the `dinv`
// condition arrives as its flag), so the per-lane check passes this array as `dinv`; it does not
// change any value NCBI computes.
static UNREAD_DINV: [f64; NN] = [0.0; NN];

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
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:575-578,581-583,587-587,593-602
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
// ...
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
// One Newton step (factor, solve, step bound, update of x and z) for every lane, in the operation
// order of `x_newton_exact::optimize`; the Cholesky factor of a lane whose `structure_ok` fails (or,
// in tests, a forced lane) comes from `factor_literal` on that lane's values. Lanes with `step[l]`
// false compute values that are discarded.
#[inline(always)]
fn newton_step<V: LaneVec>(st: &mut LaneState<V>, step: &[bool; MAX_LANES], _round: usize) {
    let zero = V::splat(0.0);
    let one = V::splat(1.0);
    let minus_one = V::splat(-1.0);
    let infinity = V::splat(f64::INFINITY);
    // ---- FactorReNewtonSystem: dinv, the sums of ScaledSymmetricProductA, the last row of W.
    let eta = st.z[M - 1];
    let one_minus_eta = one.sub(eta);
    let mut wdiag = [zero; N];
    let mut wrow = [zero; N];
    let mut w39 = [zero; M];
    let mut w39_39 = zero;
    // every dinv finite and positive (for structure_ok), two flags as in x_newton_exact
    let mut dinv_pos = V::Mask::all();
    let mut dinv_fin = V::Mask::all();
    for i in 0..N {
        let base = i * N;
        let mut srow = zero;
        let mut s39 = zero;
        for j in 0..N {
            let k = base + j;
            let dd = st.x[k].div(one_minus_eta);
            st.dinv[k] = dd;
            dinv_pos = dinv_pos.and(dd.gt(zero));
            dinv_fin = dinv_fin.and(dd.lt(infinity));
            wdiag[j] = wdiag[j].add(dd);
            srow = srow.add(dd);
            let g = st.grads1[k];
            let ws = dd.mul(g);
            w39_39 = w39_39.add(g.mul(ws));
            w39[j] = w39[j].add(one.mul(ws));
            s39 = s39.add(one.mul(ws));
        }
        if i > 0 {
            wrow[i] = srow;
            w39[i + N - 1] = s39;
        }
    }
    w39[M - 1] = w39_39;
    // Per-lane choice between the structured and the literal factorisation.
    let mut literal = [false; MAX_LANES];
    for l in 0..V::L {
        if !step[l] {
            continue;
        }
        let wd: [f64; N] = std::array::from_fn(|j| wdiag[j][l]);
        let wr: [f64; N] = std::array::from_fn(|j| wrow[j][l]);
        let w3: [f64; M] = std::array::from_fn(|j| w39[j][l]);
        let parts = Wparts {
            wdiag: &wd,
            wrow: &wr,
            dinv: &UNREAD_DINV,
            w39: &w3,
        };
        literal[l] = !structure_ok(&parts, dinv_pos.lane(l) & dinv_fin.lane(l));
        #[cfg(test)]
        if FORCE_LITERAL.with(|f| f.get()) == Some((_round, l)) {
            literal[l] = true;
        }
    }
    // ---- Structured factorisation (x_newton_exact::factor_structured), every lane.
    for j in 0..N {
        st.d[j] = wdiag[j].sqrt();
    }
    for r in 0..N - 1 {
        let base = (r + 1) * N;
        for j in 0..N {
            st.lsp[j][r] = zero.add(st.dinv[base + j]).div(st.d[j]);
        }
    }
    for j in 0..N {
        st.lsp[j][N - 1] = w39[j].div(st.d[j]);
    }
    for jj in 0..N {
        st.c[jj] = [zero; N];
    }
    for jj in 0..N - 1 {
        st.c[jj][jj] = wrow[jj + 1];
        st.c[jj][N - 1] = w39[N + jj];
    }
    st.c[N - 1][N - 1] = w39[M - 1];
    // Terms k = 0..19: T[ii][jj] -= L[20+ii][k] * L[20+jj][k], in k order, column jj in blocks
    // from the 4-aligned row at or below jj (rows ii < jj are scratch, never read). Eight lanes
    // use blocks of four rows (eight accumulators of eight lanes do not fit in the vector
    // registers); each element receives the same terms in the same order with either blocking.
    for jj in 0..N {
        let start = 4 * (jj / 4);
        let col = &mut st.c[jj];
        let lsp = &st.lsp;
        if V::L > LANES {
            for block in (start..N).step_by(4) {
                column_terms_block::<4, V>(col, lsp, jj, block);
            }
            continue;
        }
        match N - start {
            20 => {
                column_terms_block::<8, V>(col, lsp, jj, start);
                column_terms_block::<8, V>(col, lsp, jj, start + 8);
                column_terms_block::<4, V>(col, lsp, jj, start + 16);
            }
            16 => {
                column_terms_block::<8, V>(col, lsp, jj, start);
                column_terms_block::<8, V>(col, lsp, jj, start + 8);
            }
            12 => {
                column_terms_block::<8, V>(col, lsp, jj, start);
                column_terms_block::<4, V>(col, lsp, jj, start + 8);
            }
            8 => column_terms_block::<8, V>(col, lsp, jj, start),
            _ => column_terms_block::<4, V>(col, lsp, jj, start),
        }
    }
    // The block's own factorisation, right-looking: column kk finalised (square root, then the
    // divisions), then subtracted from the later columns in increasing kk.
    for kk in 0..N {
        let dkk = st.c[kk][kk].sqrt();
        st.c[kk][kk] = dkk;
        for ii in kk + 1..N {
            st.c[kk][ii] = st.c[kk][ii].div(dkk);
        }
        let (done, rest) = st.c.split_at_mut(kk + 1);
        let ck = &done[kk];
        for (offset, cj) in rest.iter_mut().enumerate() {
            let jj = kk + 1 + offset;
            let ljk = ck[jj];
            for ii in jj..N {
                cj[ii] = cj[ii].sub(ck[ii].mul(ljk));
            }
        }
    }
    // ---- Literal factorisation for the lanes that need it (scalar, then put into the lane).
    for l in 0..V::L {
        if !literal[l] {
            continue;
        }
        let wd: [f64; N] = std::array::from_fn(|j| wdiag[j][l]);
        let wr: [f64; N] = std::array::from_fn(|j| wrow[j][l]);
        let w3: [f64; M] = std::array::from_fn(|j| w39[j][l]);
        let dv: [f64; NN] = std::array::from_fn(|k| st.dinv[k][l]);
        let parts = Wparts {
            wdiag: &wd,
            wrow: &wr,
            dinv: &dv,
            w39: &w3,
        };
        let mut f = Factor::new();
        factor_literal(&parts, &mut f);
        for j in 0..N {
            st.d[j][l] = f.d[j];
        }
        for k in 0..N {
            for r in 0..N {
                st.lsp[k][r][l] = f.lsp[k][r];
            }
        }
        for jj in 0..N {
            for ii in jj..N {
                st.c[jj][ii][l] = f.c[jj][ii];
            }
        }
    }
    // ---- SolveReNewtonSystem(resids_x, resids_z)
    let mut s39 = st.resids_z[M - 1];
    for i in 0..N {
        let base = i * N;
        let mut srow = st.resids_z[i + N - 1];
        for j in 0..N {
            let k = base + j;
            let ws = st.resids_x[k].mul(st.dinv[k]);
            st.resids_z[j] = st.resids_z[j].add(minus_one.mul(ws));
            srow = srow.add(minus_one.mul(ws));
            s39 = s39.sub(st.grads1[k].mul(ws));
        }
        if i > 0 {
            st.resids_z[i + N - 1] = srow;
        }
    }
    st.resids_z[M - 1] = s39;
    solve_factored(&mut st.resids_z, &st.d, &st.lsp, &st.c);
    // The final resids_x[k] is followed at once by term k of Nlm_StepBound (same values, k in
    // increasing order).
    let z_re = st.resids_z[M - 1];
    let mut alpha = V::splat(1.0 / 0.95);
    for i in 0..N {
        let base = i * N;
        let zr = st.resids_z[i + N - 1];
        for j in 0..N {
            let k = base + j;
            let mut r = st.resids_x[k];
            r = r.add(st.grads1[k].mul(z_re));
            r = r.add(one.mul(st.resids_z[j]));
            if i > 0 {
                r = r.add(one.mul(zr));
            }
            let step_k = r.mul(st.dinv[k]);
            st.resids_x[k] = step_k;
            // ---- Nlm_StepBound(x, 400, resids_x, 1 / 0.95), element k
            let alpha_k = st.x[k].neg().div(step_k);
            let smaller = alpha_k.ge(zero).and(alpha_k.lt(alpha));
            if smaller.any() {
                alpha = V::select(smaller, alpha_k, alpha);
            }
        }
    }
    // ---- alpha *= 0.95, then Nlm_AddVectors
    alpha = alpha.mul(V::splat(0.95));
    for k in 0..NN {
        st.x[k] = st.x[k].add(alpha.mul(st.resids_x[k]));
    }
    for k in 0..M {
        st.z[k] = st.z[k].add(alpha.mul(st.resids_z[k]));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:749-751
// ```c
//     its = 0;        /* Initialize the iteration count. Note that we may
//                        converge in zero iterations if the initial x is
//                        optimal. */
// ```
// Per lane: the problem it holds (`NO_PROBLEM` for none), that problem's `its`, and whether its first
// round may take log(1) for every quotient (`load_lane`).
struct LaneSlots {
    slot: [usize; MAX_LANES],
    its: [usize; MAX_LANES],
    shortcut: [bool; MAX_LANES],
}

impl LaneSlots {
    // No NCBI counterpart: no lane holds a problem; it does not change any value NCBI computes.
    fn empty() -> Self {
        Self {
            slot: [NO_PROBLEM; MAX_LANES],
            its: [0; MAX_LANES],
            shortcut: [false; MAX_LANES],
        }
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
// The Newton loop for a queue of problems, `V::L` lanes at a time. In each round every lane with a
// problem performs one pass of the loop body of `x_newton_exact::optimize` (residuals, stopping
// test, iteration count, Newton step); a lane whose problem ends (converged or `its > maxits`)
// writes its result and takes the next queued problem in the next round. The loop starts from the
// lanes and the queue position it is given, and ends when no lane holds a problem or, with
// `hand_over_at` above zero, when the queue is empty and at most `hand_over_at` lanes hold one.
/// Returns (rounds, rounds × lanes that held a problem).
#[inline(always)]
fn solve_lanes_body<V: LaneVec>(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st: &mut LaneState<V>,
    lanes: &mut LaneSlots,
    next: &mut usize,
    hand_over_at: usize,
) -> (u64, u64) {
    let ln1 = ln_one();
    let zero = V::splat(0.0);
    let one = V::splat(1.0);
    let minus_one = V::splat(-1.0);
    let LaneSlots {
        slot,
        its,
        shortcut,
    } = lanes;
    let mut rounds = 0u64;
    let mut busy = 0u64;
    let mut buf = [0.0f64; NN];
    let mut out = [0.0f64; NN];
    loop {
        // ---- Fill the lanes without a problem from the queue.
        for l in 0..V::L {
            if slot[l] == NO_PROBLEM && *next < inputs.len() {
                shortcut[l] = load_lane(st, l, &inputs[*next]);
                its[l] = 0;
                slot[l] = *next;
                *next += 1;
            }
        }
        let Some(src) = (0..V::L).find(|&l| slot[l] != NO_PROBLEM) else {
            break;
        };
        if hand_over_at > 0
            && *next == inputs.len()
            && (0..V::L).filter(|&l| slot[l] != NO_PROBLEM).count() <= hand_over_at
        {
            break;
        }
        for l in 0..V::L {
            if slot[l] == NO_PROBLEM {
                mirror_lane(st, l, src);
            } else {
                busy += 1;
            }
        }
        let round = rounds as usize;
        rounds += 1;

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
        // ---- EvaluateReFunctions: log(x / q) per lane (values[0] is never read by the reference).
        for k in 0..NN {
            st.logs[k] = st.x[k].div(st.q[k]);
        }
        for l in 0..V::L {
            if slot[l] == NO_PROBLEM {
                continue;
            }
            if its[l] == 0 && shortcut[l] {
                for k in 0..NN {
                    st.logs[k][l] = ln1;
                }
            } else {
                for k in 0..NN {
                    buf[k] = st.logs[k][l];
                }
                crate::utils::x_logclone::log_slice(&buf, &mut out);
                for k in 0..NN {
                    st.logs[k][l] = out[k];
                }
            }
        }
        for l in 0..V::L {
            if slot[l] == NO_PROBLEM {
                for k in 0..NN {
                    st.logs[k][l] = st.logs[k][src];
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:277-279,286-286
        // ```c
        //         eta = z[2 * alphsize - 1];
        //         for (i = 0;  i < n;  i++) {
        //             resids_x[i] = -grads[0][i] + eta * grads[1][i];
        // ...
        //     MultiplyByAtranspose(1.0, resids_x, alphsize, 1.0, z);
        // ```
        // ---- grads, values[1], DualResiduals and MultiplyByAtranspose, as in x_newton_exact.
        let eta = st.z[M - 1];
        let mut values1 = zero;
        // Nlm_EuclideanNorm(resids_x) runs element by element as resids_x[k] is formed (k in order).
        let mut norm_sum = one;
        let mut norm_scale = zero;
        for i in 0..N {
            let base = i * N;
            let zr = st.z[i + N - 1];
            for j in 0..N {
                let k = base + j;
                let temp = st.logs[k];
                let g0 = temp.add(one);
                let temp = temp.add(st.scores[k]);
                let g1 = temp.add(one);
                st.grads1[k] = g1;
                let mut r = g0.neg().add(eta.mul(g1));
                r = r.add(one.mul(st.z[j]));
                if i > 0 {
                    r = r.add(one.mul(zr));
                }
                st.resids_x[k] = r;
                values1 = values1.add(st.x[k].mul(temp));
                euclidean_norm_step(&mut norm_sum, &mut norm_scale, r);
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:331-332
        // ```c
        //     DualResiduals(resids_x, alphsize, grads, z, constrain_rel_entropy);
        //     norm_resids_x = Nlm_EuclideanNorm(resids_x, alphsize * alphsize);
        // ```
        let norm_resids_x = norm_scale.mul(norm_sum.sqrt());
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:243-249
        // ```c
        //     for (i = 0;  i < alphsize;  i++) {
        //         rA[i] = col_sums[i];
        //     }
        //     for (i = 1;  i < alphsize;  i++) {
        //         rA[i + alphsize - 1] = row_sums[i];
        //     }
        //     MultiplyByA(1.0, rA, alphsize, -1.0, x);
        // ```
        // ---- ResidualsLinearConstraints: column sums over i, then row sums over j.
        st.resids_z[..N].copy_from_slice(&st.col);
        for i in 1..N {
            st.resids_z[i + N - 1] = st.row[i];
        }
        for i in 0..N {
            let base = i * N;
            for j in 0..N {
                st.resids_z[j] = st.resids_z[j].add(minus_one.mul(st.x[base + j]));
            }
        }
        for i in 1..N {
            let base = i * N;
            let mut s = st.resids_z[i + N - 1];
            for j in 0..N {
                s = s.add(minus_one.mul(st.x[base + j]));
            }
            st.resids_z[i + N - 1] = s;
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
        st.resids_z[M - 1] = st.re.sub(values1);
        let norm_resids_z = euclidean_norm(&st.resids_z);
        let rnorm = norm_resids_x
            .mul(norm_resids_x)
            .add(norm_resids_z.mul(norm_resids_z))
            .sqrt();
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:762-768,789-797
        // ```c
        //         if ( !(rnorm > tol) ) {
        //             /* We converged at the current iterate */
        //             break;
        //         } else {
        //             /* we did not converge, so increment the iteration counter
        //                and start a new iteration */
        //             if (++its <= maxits) {
        // ...
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
        // ---- Stopping test per lane; a finished lane writes its result.
        let mut step = [false; MAX_LANES];
        for l in 0..V::L {
            let s = slot[l];
            if s == NO_PROBLEM {
                continue;
            }
            let lane_rnorm = rnorm[l];
            let finished = if !(lane_rnorm > K_COMPO_ADJUST_ERR_TOLERANCE) {
                true
            } else {
                its[l] += 1;
                its[l] > K_COMPO_ADJUST_ITERATION_LIMIT
            };
            if finished {
                let converged = its[l] <= K_COMPO_ADJUST_ITERATION_LIMIT
                    && lane_rnorm <= K_COMPO_ADJUST_ERR_TOLERANCE
                    && st.z[M - 1][l] < 1.0;
                let output = &mut outputs[s];
                output.status = if converged { 0 } else { 1 };
                for k in 0..NN {
                    output.x[k] = st.x[k][l];
                }
                // No NCBI counterpart: iteration counter; it does not change any value NCBI computes.
                crate::utils::xstats::add(&crate::utils::xstats::NEWTON_ITERS, its[l] as u64);
                slot[l] = NO_PROBLEM;
            } else {
                step[l] = true;
            }
        }
        if step.iter().any(|&s| s) {
            newton_step(st, &step, round);
        }
    }
    (rounds, busy)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:749-752
// ```c
//     its = 0;        /* Initialize the iteration count. Note that we may
//                        converge in zero iterations if the initial x is
//                        optimal. */
//     while (its <= maxits) {
// ```
// The four-lane loop over a whole batch (no lane holds a problem, the queue at its start).
/// Returns (rounds, rounds × lanes that held a problem).
#[inline(always)]
fn solve_lanes4_body(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st: &mut LaneState<F4>,
) -> (u64, u64) {
    let mut lanes = LaneSlots::empty();
    let mut next = 0usize;
    solve_lanes_body(inputs, outputs, st, &mut lanes, &mut next, 0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:749-752
// ```c
//     its = 0;        /* Initialize the iteration count. Note that we may
//                        converge in zero iterations if the initial x is
//                        optimal. */
//     while (its <= maxits) {
// ```
// The eight-lane loop over a whole batch until the queue is empty and at most `hand_over_at` lanes
// hold a problem (0: to the end); then those problems move to four lanes in lane order (`move_lane`,
// with their `its` and first-round flag) and finish there. Each problem keeps its operations.
/// Returns [eight-lane rounds, their busy lane rounds, four-lane rounds, their busy lane rounds].
#[inline(always)]
fn solve_lanes8_body(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st8: &mut LaneState<F8>,
    st4: &mut LaneState<F4>,
    hand_over_at: usize,
) -> [u64; 4] {
    assert!(hand_over_at <= LANES);
    let mut lanes = LaneSlots::empty();
    let mut next = 0usize;
    let (rounds8, busy8) =
        solve_lanes_body(inputs, outputs, st8, &mut lanes, &mut next, hand_over_at);
    let mut tail = LaneSlots::empty();
    let mut d = 0usize;
    for l in 0..LANES8 {
        if lanes.slot[l] != NO_PROBLEM {
            move_lane(st4, d, st8, l);
            tail.slot[d] = lanes.slot[l];
            tail.its[d] = lanes.its[l];
            tail.shortcut[d] = lanes.shortcut[l];
            d += 1;
        }
    }
    let (rounds4, busy4) = solve_lanes_body(inputs, outputs, st4, &mut tail, &mut next, 0);
    [rounds8, busy8, rounds4, busy4]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the same lane body compiled with AVX2 enabled so that the compiler may use
// 256-bit vectors. Rust does not fuse `a * b + c` or reorder sums, so every lane computes the same
// IEEE-754 operations as the body without the feature; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn solve_lanes_avx2(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st: &mut LaneState<F4>,
) -> (u64, u64) {
    solve_lanes4_body(inputs, outputs, st)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the lane body without target features (Wasm, other targets,
// LOSAT_X_NEWTONEXACT_LEVEL=0); it does not change any value NCBI computes.
fn solve_lanes_plain(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st: &mut LaneState<F4>,
) -> (u64, u64) {
    solve_lanes4_body(inputs, outputs, st)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the eight-lane loop compiled with AVX2 enabled (each `F8` operation is two
// 256-bit operations on independent values); x86-64 runs it only in the unit tests
// (the eight lanes run on aarch64). It does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn solve_lanes8_avx2(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st8: &mut LaneState<F8>,
    st4: &mut LaneState<F4>,
    hand_over_at: usize,
) -> [u64; 4] {
    solve_lanes8_body(inputs, outputs, st8, st4, hand_over_at)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-752
// ```c
//     while (its <= maxits) {
// ```
// No NCBI counterpart: the eight-lane loop without target features (aarch64 NEON through the
// baseline target; other targets in the unit tests); it does not change any value NCBI computes.
fn solve_lanes8_plain(
    inputs: &[NewtonInput],
    outputs: &mut [NewtonOutput],
    st8: &mut LaneState<F8>,
    st4: &mut LaneState<F4>,
    hand_over_at: usize,
) -> [u64; 4] {
    solve_lanes8_body(inputs, outputs, st8, st4, hand_over_at)
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
// Solves every problem of the batch, each with the result of this function (`x` and the status).
// Batches of fewer than `MIN_LANE_BATCH` problems use `x_newton_exact` (same results); with
// aarch64, batches of `MIN_LANE8_BATCH` or more use eight lanes.
/// Solves `inputs[i]` into `outputs[i]`.
pub(crate) fn solve_batch(inputs: &[NewtonInput], outputs: &mut [NewtonOutput]) {
    let lanes8 = if eight_lanes() {
        Some(HAND_OVER_LANES)
    } else {
        None
    };
    solve_batch_with(inputs, outputs, lanes8);
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
// `solve_batch` with the width given: `lanes8` is `None` for four lanes, `Some(hand_over_at)` for
// eight lanes on batches of `MIN_LANE8_BATCH` or more (see `solve_lanes8_body`).
/// Solves `inputs[i]` into `outputs[i]`.
fn solve_batch_with(inputs: &[NewtonInput], outputs: &mut [NewtonOutput], lanes8: Option<usize>) {
    assert_eq!(inputs.len(), outputs.len());
    if let Some(hand_over_at) = lanes8 {
        if inputs.len() >= MIN_LANE8_BATCH {
            solve_batch8(inputs, outputs, hand_over_at);
            return;
        }
    }
    if inputs.len() < MIN_LANE_BATCH {
        for (input, output) in inputs.iter().zip(outputs.iter_mut()) {
            output.status = optimize_target_frequencies_exact(
                &mut output.x,
                &input.q,
                &input.row,
                &input.col,
                input.relative_entropy,
            );
        }
        return;
    }
    let mut st = new_lane_state();
    #[cfg(target_arch = "x86_64")]
    let (rounds, busy) = if super::x_newton_exact::x86_level() >= 1 {
        // SAFETY: AVX2 was detected at run time (x86_level).
        unsafe { solve_lanes_avx2(inputs, outputs, &mut st) }
    } else {
        solve_lanes_plain(inputs, outputs, &mut st)
    };
    #[cfg(not(target_arch = "x86_64"))]
    let (rounds, busy) = solve_lanes_plain(inputs, outputs, &mut st);
    BATCHES.fetch_add(1, Ordering::Relaxed);
    LANE_ROUNDS.fetch_add(rounds * LANES as u64, Ordering::Relaxed);
    BUSY_LANE_ROUNDS.fetch_add(busy, Ordering::Relaxed);
}

// No NCBI counterpart: the working states of the eight-lane batches of this thread, allocated on the
// first such batch and kept (a fresh zeroed pair per batch would cost about 5 % of the Newton time).
// Every lane is loaded, moved or mirrored before its values are read, so values left by an earlier
// batch are never used. It does not change any value NCBI computes.
thread_local! {
    static LANES8_STATES: RefCell<Option<(Box<LaneState<F8>>, Box<LaneState<F4>>)>> =
        const { RefCell::new(None) };
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
// An eight-lane batch (`solve_lanes8_body`), AVX2 or plain as the four-lane batches.
/// Solves `inputs[i]` into `outputs[i]` on eight lanes.
fn solve_batch8(inputs: &[NewtonInput], outputs: &mut [NewtonOutput], hand_over_at: usize) {
    LANES8_STATES.with(|states| {
        let mut states = states.borrow_mut();
        let (st8, st4) = states.get_or_insert_with(|| (new_lane_state(), new_lane_state()));
        #[cfg(target_arch = "x86_64")]
        let counts = if super::x_newton_exact::x86_level() >= 1 {
            // SAFETY: AVX2 was detected at run time (x86_level).
            unsafe { solve_lanes8_avx2(inputs, outputs, st8, st4, hand_over_at) }
        } else {
            solve_lanes8_plain(inputs, outputs, st8, st4, hand_over_at)
        };
        #[cfg(not(target_arch = "x86_64"))]
        let counts = solve_lanes8_plain(inputs, outputs, st8, st4, hand_over_at);
        let [rounds8, busy8, rounds4, busy4] = counts;
        BATCHES.fetch_add(1, Ordering::Relaxed);
        LANE_ROUNDS.fetch_add(
            rounds8 * LANES8 as u64 + rounds4 * LANES as u64,
            Ordering::Relaxed,
        );
        BUSY_LANE_ROUNDS.fetch_add(busy8 + busy4, Ordering::Relaxed);
        BATCHES8.fetch_add(1, Ordering::Relaxed);
        ROUNDS8.fetch_add(rounds8, Ordering::Relaxed);
        TAIL_ROUNDS.fetch_add(rounds4, Ordering::Relaxed);
    });
}

// No NCBI counterpart: one solved problem kept for a later call with the same input; it does not
// change any value NCBI computes.
struct StoredResult {
    input: NewtonInput,
    x: [f64; NN],
    status: i32,
}

// No NCBI counterpart: per-thread store of results computed ahead (bounded; the oldest entry is
// dropped first); it does not change any value NCBI computes.
thread_local! {
    static STORE: RefCell<VecDeque<Box<StoredResult>>> = const { RefCell::new(VecDeque::new()) };
}

// No NCBI counterpart: solves a batch of inputs gathered ahead of their calls and keeps the results
// for `take`; it does not change any value NCBI computes.
/// Solves `inputs` and keeps the results in this thread's store.
pub(crate) fn prefetch(inputs: Vec<NewtonInput>) {
    // Dispatch point (LOSAT_X_NEWTONDROP, in the strict set; `=0` keeps the results not taken).
    x_prefetch(inputs, x_newton_drop());
}

// No NCBI counterpart: `prefetch`, with `drop` for LOSAT_X_NEWTONDROP; it does not change any value
// NCBI computes.
/// A call comes after the windows (TBLASTN) or the matches (BLASTP: the next matches of a query, or
/// a job's matches in the windowed redo) of this thread's earlier calls were walked, so their
/// results not taken are not asked for again; with `drop` they are dropped before this batch is
/// solved, without it they stay until the store holds `STORE_CAPACITY` results. A dropped result
/// asked for after all would be solved by the caller with the same bits (`take` returns `None`).
#[inline]
fn x_prefetch(inputs: Vec<NewtonInput>, drop: bool) {
    if drop {
        STORE.with(|store| {
            let mut store = store.borrow_mut();
            DROPPED.fetch_add(store.len() as u64, Ordering::Relaxed);
            store.clear();
        });
        XBATCH.with(|batch| {
            if let Some(old) = batch.borrow_mut().take() {
                DROPPED.fetch_add(old.left as u64, Ordering::Relaxed);
            }
        });
    }
    if inputs.is_empty() {
        return;
    }
    if drop && x_probe_nlbatch() {
        return x_prefetch_batch(inputs);
    }
    let mut outputs = vec![NewtonOutput::new(); inputs.len()];
    solve_batch(&inputs, &mut outputs);
    PREFETCHED.fetch_add(inputs.len() as u64, Ordering::Relaxed);
    STORE.with(|store| {
        let mut store = store.borrow_mut();
        for (input, output) in inputs.into_iter().zip(outputs) {
            if store.len() >= STORE_CAPACITY {
                store.pop_front();
            }
            store.push_back(Box::new(StoredResult {
                input,
                x: output.x,
                status: output.status,
            }));
        }
    });
}

// No NCBI counterpart (probe D14 s6, exp branch only): LOSAT_PROBE_NLBATCH=1 (with NEWTONDROP) keeps
// the one batch a thread solved ahead as its input and output vectors (no box per result, no copy);
// it does not change any value NCBI computes.
fn x_probe_nlbatch() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_PROBE_NLBATCH").is_some_and(|v| v == "1"))
}

// No NCBI counterpart (probe): the batch of this thread's last prefetch; `taken` marks the results
// taken, `next` is the first one not taken, `left` how many are not taken.
struct XBatch {
    inputs: Vec<NewtonInput>,
    outputs: Vec<NewtonOutput>,
    taken: Vec<bool>,
    next: usize,
    left: usize,
}

// No NCBI counterpart (probe): per-thread batch computed ahead (NEWTONDROP drops it at the next
// prefetch); it does not change any value NCBI computes.
thread_local! {
    static XBATCH: RefCell<Option<XBatch>> = const { RefCell::new(None) };
}

// No NCBI counterpart (probe): `prefetch` with the batch kept as it was solved.
fn x_prefetch_batch(inputs: Vec<NewtonInput>) {
    let n = inputs.len();
    let mut outputs = vec![NewtonOutput::new(); n];
    solve_batch(&inputs, &mut outputs);
    PREFETCHED.fetch_add(n as u64, Ordering::Relaxed);
    XBATCH.with(|batch| {
        *batch.borrow_mut() = Some(XBatch {
            inputs,
            outputs,
            taken: vec![false; n],
            next: 0,
            left: n,
        });
    });
}

// No NCBI counterpart (probe): `take` from the kept batch: the first result not taken, from `next`
// on, whose input equals this one by bits; the batch goes when its last result is taken.
fn x_take_batch(
    x: &mut [f64],
    q: &[f64],
    row: &[f64; N],
    col: &[f64; N],
    relative_entropy: f64,
) -> Option<i32> {
    XBATCH.with(|cell| {
        let mut cell = cell.borrow_mut();
        let found = cell.as_mut().and_then(|batch| {
            let i = (batch.next..batch.inputs.len()).find(|&i| {
                let input = &batch.inputs[i];
                !batch.taken[i]
                    && input.relative_entropy.to_bits() == relative_entropy.to_bits()
                    && same_bits(&input.col, col)
                    && same_bits(&input.row, row)
                    && same_bits(&input.q, q)
            })?;
            x.copy_from_slice(&batch.outputs[i].x);
            batch.taken[i] = true;
            batch.left -= 1;
            while batch.next < batch.taken.len() && batch.taken[batch.next] {
                batch.next += 1;
            }
            Some((batch.outputs[i].status, batch.left == 0))
        });
        match found {
            Some((status, last)) => {
                if last {
                    *cell = None;
                }
                HITS.fetch_add(1, Ordering::Relaxed);
                Some(status)
            }
            None => {
                MISSES.fetch_add(1, Ordering::Relaxed);
                None
            }
        }
    })
}

// No NCBI counterpart (memory only; the results kept are the ones `prefetch` solved):
// LOSAT_X_NEWTONDROP, in the strict set (on unless `LOSAT_X_NEWTONDROP=0`).
/// Whether `prefetch` drops the results of this thread's earlier calls that were not taken.
fn x_newton_drop() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_NEWTONDROP"))
}

// No NCBI counterpart: bitwise equality of two f64 slices.
#[inline]
fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(u, v)| u.to_bits() == v.to_bits())
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
// No NCBI counterpart: looks up a result computed ahead for exactly this input (relative entropy,
// column and row probabilities, q compared by bits); on a hit `x` receives the stored final iterate
// and the stored status is returned. It does not change any value NCBI computes.
/// The stored result for this input, if any (removed from the store).
pub(crate) fn take(
    x: &mut [f64],
    q: &[f64],
    row: &[f64; N],
    col: &[f64; N],
    relative_entropy: f64,
) -> Option<i32> {
    if x_probe_nlbatch() && x_newton_drop() {
        return x_take_batch(x, q, row, col, relative_entropy);
    }
    STORE.with(|store| {
        let mut store = store.borrow_mut();
        let position = store.iter().position(|entry| {
            entry.input.relative_entropy.to_bits() == relative_entropy.to_bits()
                && same_bits(&entry.input.col, col)
                && same_bits(&entry.input.row, row)
                && same_bits(&entry.input.q, q)
        });
        match position {
            Some(position) => {
                let entry = store.remove(position).expect("stored entry");
                x.copy_from_slice(&entry.x);
                HITS.fetch_add(1, Ordering::Relaxed);
                Some(entry.status)
            }
            None => {
                MISSES.fetch_add(1, Ordering::Relaxed);
                None
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::super::adjust_scores::{
        build_matrix_info, read_aa_composition, x_newton_input, x_optimize_reference_for_test,
        BlastCompositionWorkspace,
    };
    use super::super::redo_alignment::BlastCompoAdjustMode;
    use super::*;
    use crate::config::ScoringMatrix;

    // No NCBI counterpart: random number generator of the test; it does not change any value NCBI computes.
    fn xorshift(state: &mut u64) -> f64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1385-1399
    // ```c
    //     Blast_ApplyPseudocounts(row_probs, length1,
    //                             NRrecord->first_standard_freq, pseudocounts);
    //     Blast_ApplyPseudocounts(col_probs, length2,
    //                             NRrecord->second_standard_freq, pseudocounts);
    //
    //     status =
    //         Blast_OptimizeTargetFrequencies(&NRrecord->mat_final[0][0],
    // ```
    // Test problems: realistic ones through `x_newton_input` (BLOSUM62 joint probabilities as q,
    // compositions of random sequences with biased letter frequencies, pseudocounts), and harder
    // ones with random q, compositions and relative entropies.
    fn test_problems() -> Vec<NewtonInput> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
        let workspace = BlastCompositionWorkspace::new_blosum62();
        let mut problems = Vec::new();
        let mut attempts = 0;
        while problems.len() < 40 && attempts < 400 {
            attempts += 1;
            let mut weights = [0.0f64; 20];
            for w in weights.iter_mut() {
                *w = 0.05 + xorshift(&mut state).powi(3) * 4.0;
            }
            let total: f64 = weights.iter().sum();
            let mut sequences = Vec::new();
            for _ in 0..2 {
                let len = 30 + (xorshift(&mut state) * 900.0) as usize;
                let mut seq = Vec::with_capacity(len);
                for _ in 0..len {
                    let mut u = xorshift(&mut state) * total;
                    let mut letter = 0usize;
                    while letter < 19 && u >= weights[letter] {
                        u -= weights[letter];
                        letter += 1;
                    }
                    // NCBIstdaa codes of the 20 true amino acids (A..Y without B, X, Z, U, *).
                    const TRUE_AA: [u8; 20] = [
                        1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
                    ];
                    seq.push(TRUE_AA[letter]);
                }
                sequences.push(seq);
                // the second sequence uses other weights
                for w in weights.iter_mut() {
                    *w = 0.05 + xorshift(&mut state).powi(2) * 2.0;
                }
            }
            let query = read_aa_composition(&sequences[0]);
            let subject = read_aa_composition(&sequences[1]);
            let mode = if attempts % 3 == 0 {
                BlastCompoAdjustMode::ForceFullMatrixAdjust
            } else {
                BlastCompoAdjustMode::CompositionMatrixAdjust
            };
            if let Some(input) = x_newton_input(
                &matrix_info,
                &query,
                sequences[0].len() as i32,
                &subject,
                sequences[1].len() as i32,
                mode,
                &workspace,
            ) {
                problems.push(input);
            }
        }
        assert!(
            problems.len() >= 30,
            "realistic problems: {}",
            problems.len()
        );
        let blosum_q = problems[0].q;
        // Random joint distributions, compositions and relative entropies.
        for case in 0..24 {
            let mut input = problems[case % problems.len()].clone();
            if case % 2 == 0 {
                let mut total = 0.0;
                for v in input.q.iter_mut() {
                    *v = 1.0e-4 + xorshift(&mut state).powi(2);
                    total += *v;
                }
                for v in input.q.iter_mut() {
                    *v /= total;
                }
            } else {
                input.q = blosum_q;
            }
            let (mut sr, mut sc) = (0.0, 0.0);
            for k in 0..N {
                input.row[k] = 0.002 + xorshift(&mut state).powi(3);
                input.col[k] = 0.002 + xorshift(&mut state).powi(3);
                sr += input.row[k];
                sc += input.col[k];
            }
            for k in 0..N {
                input.row[k] /= sr;
                input.col[k] /= sc;
            }
            input.relative_entropy = match case % 4 {
                0 => 0.05 + xorshift(&mut state),
                1 => 0.44,
                2 => 1.5 + 2.0 * xorshift(&mut state),
                _ => 1.0e-3,
            };
            problems.push(input);
        }
        problems
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
    // Expected results of each test problem: the exact module and the reference port (which must
    // agree with each other).
    fn expected(problems: &[NewtonInput]) -> Vec<NewtonOutput> {
        problems
            .iter()
            .map(|input| {
                let mut exact = NewtonOutput::new();
                exact.status = optimize_target_frequencies_exact(
                    &mut exact.x,
                    &input.q,
                    &input.row,
                    &input.col,
                    input.relative_entropy,
                );
                let mut reference = NewtonOutput::new();
                reference.status = x_optimize_reference_for_test(
                    &mut reference.x,
                    &input.q,
                    &input.row,
                    &input.col,
                    input.relative_entropy,
                );
                assert_eq!(exact.status, reference.status);
                assert!(same_bits(&exact.x, &reference.x));
                reference
            })
            .collect()
    }

    // No NCBI counterpart: bitwise comparison of a batch result with the expected results.
    fn check_batch(label: &str, batch: &[usize], outputs: &[NewtonOutput], want: &[NewtonOutput]) {
        for (position, &problem) in batch.iter().enumerate() {
            assert_eq!(
                outputs[position].status, want[problem].status,
                "{label}: problem {problem} status"
            );
            for k in 0..NN {
                assert_eq!(
                    outputs[position].x[k].to_bits(),
                    want[problem].x[k].to_bits(),
                    "{label}: problem {problem} x[{k}]: {} vs {}",
                    outputs[position].x[k],
                    want[problem].x[k]
                );
            }
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
    // ```
    // The test solves batches of every size 1..=13 and one of 64 (mixed iteration counts, lanes
    // refilled as problems end) and compares status and all 400 values by bits with the exact
    // module and the reference port, then repeats batches with the literal factorisation forced on
    // one lane in one round.
    #[test]
    fn x_newton_lanes_match_exact_and_reference_bitwise() {
        let problems = test_problems();
        let want = expected(&problems);
        let failed = want.iter().filter(|o| o.status != 0).count();
        eprintln!(
            "[x_newton_lanes test] problems={} status1={}",
            problems.len(),
            failed
        );
        let mut start = 0usize;
        let mut sizes: Vec<usize> = (1..=13).collect();
        sizes.push(64);
        for &size in &sizes {
            let batch: Vec<usize> = (0..size)
                .map(|i| (start + 7 * i) % problems.len())
                .collect();
            start += 5;
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); size];
            solve_batch(&inputs, &mut outputs);
            check_batch(&format!("batch {size}"), &batch, &outputs, &want);
        }
        // The plain (non-AVX2) compilation of the lane body.
        {
            let batch: Vec<usize> = (0..9).map(|i| (3 * i + 1) % problems.len()).collect();
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); batch.len()];
            let mut st = new_lane_state();
            solve_lanes_plain(&inputs, &mut outputs, &mut st);
            check_batch("plain", &batch, &outputs, &want);
        }
        // Literal factorisation forced on one lane in one round.
        for (round, lane) in [(0usize, 0usize), (1, 3), (2, 1), (5, 2), (9, 0)] {
            let batch: Vec<usize> = (0..10).map(|i| (round + 11 * i) % problems.len()).collect();
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); batch.len()];
            FORCE_LITERAL.with(|f| f.set(Some((round, lane))));
            solve_batch(&inputs, &mut outputs);
            FORCE_LITERAL.with(|f| f.set(None));
            check_batch(
                &format!("literal round {round} lane {lane}"),
                &batch,
                &outputs,
                &want,
            );
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
    // ```
    // The eight-lane batches (aarch64): every size 1..=21 and one of 64, with the move
    // to four lanes and without it, the plain compilation, and the literal factorisation forced on one
    // lane (any of the eight) in one round; status and all 400 values compared by bits with the exact
    // module and the reference port.
    #[test]
    fn x_newton_lanes8_match_exact_and_reference_bitwise() {
        let problems = test_problems();
        let want = expected(&problems);
        let mut start = 0usize;
        let mut sizes: Vec<usize> = (1..=21).collect();
        sizes.push(64);
        for hand_over_at in [HAND_OVER_LANES, 0] {
            for &size in &sizes {
                let batch: Vec<usize> = (0..size)
                    .map(|i| (start + 7 * i) % problems.len())
                    .collect();
                start += 5;
                let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
                let mut outputs = vec![NewtonOutput::new(); size];
                solve_batch_with(&inputs, &mut outputs, Some(hand_over_at));
                check_batch(
                    &format!("lanes8 hand-over {hand_over_at} batch {size}"),
                    &batch,
                    &outputs,
                    &want,
                );
            }
        }
        // The plain (non-AVX2) compilation of the eight-lane body; the batch ends on four lanes.
        {
            let batch: Vec<usize> = (0..17).map(|i| (3 * i + 2) % problems.len()).collect();
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); batch.len()];
            let mut st8 = new_lane_state::<F8>();
            let mut st4 = new_lane_state::<F4>();
            let [rounds8, busy8, rounds4, busy4] =
                solve_lanes8_plain(&inputs, &mut outputs, &mut st8, &mut st4, HAND_OVER_LANES);
            check_batch("lanes8 plain", &batch, &outputs, &want);
            assert!(rounds8 > 0 && rounds4 > 0, "rounds {rounds8} + {rounds4}");
            assert!(busy8 <= rounds8 * LANES8 as u64 && busy4 <= rounds4 * LANES as u64);
        }
        // Literal factorisation forced on one lane in one round.
        for (round, lane) in [
            (0usize, 0usize),
            (1, 7),
            (2, 5),
            (4, 3),
            (6, 6),
            (9, 1),
            (12, 4),
        ] {
            let batch: Vec<usize> = (0..13).map(|i| (round + 11 * i) % problems.len()).collect();
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); batch.len()];
            FORCE_LITERAL.with(|f| f.set(Some((round, lane))));
            solve_batch_with(&inputs, &mut outputs, Some(HAND_OVER_LANES));
            FORCE_LITERAL.with(|f| f.set(None));
            check_batch(
                &format!("lanes8 literal round {round} lane {lane}"),
                &batch,
                &outputs,
                &want,
            );
        }
    }

    // No NCBI counterpart: the store returns a result only for a bitwise-equal input; it does not
    // change any value NCBI computes.
    #[test]
    fn x_newton_lanes_store_returns_only_exact_input_matches() {
        let problems = test_problems();
        let want = expected(&problems[..5]);
        prefetch(problems[..5].to_vec());
        let mut x = [0.0f64; NN];
        // A different relative entropy (one ulp) misses.
        let p = &problems[2];
        let re = f64::from_bits(p.relative_entropy.to_bits() + 1);
        assert_eq!(take(&mut x, &p.q, &p.row, &p.col, re), None);
        // Exact inputs hit once each, in any order.
        for i in [3usize, 0, 4, 1, 2] {
            let p = &problems[i];
            let status = take(&mut x, &p.q, &p.row, &p.col, p.relative_entropy);
            assert_eq!(status, Some(want[i].status));
            assert!(same_bits(&x, &want[i].x));
            assert_eq!(take(&mut x, &p.q, &p.row, &p.col, p.relative_entropy), None);
        }
    }

    // No NCBI counterpart: LOSAT_X_NEWTONDROP; a prefetch drops this thread's results of earlier
    // prefetches not taken (with `=0` they stay); it does not change any value NCBI computes.
    #[test]
    fn x_newton_drop_empties_the_store_at_the_next_prefetch() {
        let problems = test_problems();
        let want = expected(&problems[..8]);
        let mut x = [0.0f64; NN];
        let take_at = |x: &mut [f64; NN], i: usize| {
            let p = &problems[i];
            take(x, &p.q, &p.row, &p.col, p.relative_entropy)
        };
        // the same calls with the switch (`drop`) and with `=0`, in this process
        for drop in [true, false] {
            x_prefetch(Vec::new(), true);
            x_prefetch(problems[..5].to_vec(), drop);
            for i in [0usize, 2, 4] {
                assert_eq!(take_at(&mut x, i), Some(want[i].status), "drop {drop}");
                assert!(same_bits(&x, &want[i].x));
            }
            x_prefetch(problems[5..8].to_vec(), drop);
            // problems 1 and 3 were not taken before the second prefetch
            for i in [1usize, 3] {
                let kept = take_at(&mut x, i);
                if drop {
                    assert_eq!(kept, None, "problem {i} dropped");
                } else {
                    assert_eq!(kept, Some(want[i].status), "problem {i} kept with =0");
                    assert!(same_bits(&x, &want[i].x));
                }
            }
            for i in 5..8 {
                assert_eq!(take_at(&mut x, i), Some(want[i].status), "drop {drop}");
                assert!(same_bits(&x, &want[i].x));
            }
            // an empty prefetch drops what is left as well
            x_prefetch(problems[..2].to_vec(), drop);
            assert_eq!(take_at(&mut x, 0), Some(want[0].status), "drop {drop}");
            x_prefetch(Vec::new(), drop);
            assert_eq!(
                take_at(&mut x, 1).is_none(),
                drop,
                "problem 1 after an empty prefetch, drop {drop}"
            );
        }
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
    // The ignored test times `optimize_target_frequencies_exact` (one problem per call) against
    // four lanes and eight lanes (with and without the move to four lanes, at every batch size) on
    // the realistic test problems and checks the results by bits.
    /// Timing (run with --release --ignored --nocapture).
    #[test]
    #[ignore]
    fn x_bench_newton_lanes() {
        let problems: Vec<NewtonInput> = test_problems().into_iter().take(36).collect();
        let reps = 10;
        // rounds of the four forms in turn; the fastest round of each form is printed
        let rounds = 9;
        for batch in [5usize, 6, 7, 8, 9, 10, 11, 12, 16, 21, 32, 36] {
            let inputs: Vec<NewtonInput> = problems[..batch].to_vec();
            let mut exact = vec![NewtonOutput::new(); batch];
            let mut lanes = vec![NewtonOutput::new(); batch];
            let mut best = [f64::INFINITY; 4];
            for _ in 0..rounds {
                let t0 = std::time::Instant::now();
                for _ in 0..reps {
                    for (input, output) in inputs.iter().zip(exact.iter_mut()) {
                        output.status = optimize_target_frequencies_exact(
                            &mut output.x,
                            &input.q,
                            &input.row,
                            &input.col,
                            input.relative_entropy,
                        );
                    }
                }
                best[0] = best[0].min(t0.elapsed().as_secs_f64() / (reps * batch) as f64 * 1e6);
                for (form, t) in best[1..].iter_mut().enumerate() {
                    let t0 = std::time::Instant::now();
                    for _ in 0..reps {
                        match form {
                            0 => solve_batch_with(&inputs, &mut lanes, None),
                            1 => solve_batch8(&inputs, &mut lanes, HAND_OVER_LANES),
                            _ => solve_batch8(&inputs, &mut lanes, 0),
                        }
                    }
                    *t = t.min(t0.elapsed().as_secs_f64() / (reps * batch) as f64 * 1e6);
                    for (a, b) in exact.iter().zip(&lanes) {
                        assert_eq!(a.status, b.status);
                        assert!(same_bits(&a.x, &b.x));
                    }
                }
            }
            eprintln!(
                "[BENCH] batch {batch}: exact {:.1} us/problem, lanes {:.1}, lanes8 {:.1}, lanes8 without move {:.1} us/problem",
                best[0], best[1], best[2], best[3]
            );
        }
    }
}
