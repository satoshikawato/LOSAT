//! NCBI composition-adjustment support.
//!
//! Reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/

pub mod adjust_scores;
pub mod redo_alignment;
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
// `x_newton_exact` reproduces `Blast_OptimizeTargetFrequencies` (optimize_target_freq.c:686-815) with
// the same floating-point operations on every value, in the same order. It is used only when
// LOSAT_X_NEWTONEXACT or LOSAT_X_NEWTONEXACTSHADOW is set; see the module comment.
pub(crate) mod x_newton_exact;
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
// `x_newton_lanes` solves several independent calls of this function at once (one per vector lane,
// the operations of `x_newton_exact` per lane) ahead of the calls, and hands a result to a call whose
// input has the same bits. Used only when LOSAT_X_NEWTONLANES or LOSAT_X_NEWTONLANESSHADOW is set; see
// the module comment.
pub(crate) mod x_newton_lanes;
