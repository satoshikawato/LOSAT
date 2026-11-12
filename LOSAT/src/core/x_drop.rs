//! The X-drops of a protein search in raw scores, from the X-drops in bits that `-xdrop_ungap`,
//! `-xdrop_gap` and `-xdrop_gap_final` (or the program's defaults) give: the ungapped X-drop of a
//! context (NCBI's `BlastInitialWordParametersNew`), the preliminary and final gapped X-drops
//! (`BlastExtensionParametersNew`) and the X-drop of the composition-based redo
//! (`s_GappingParamsNew`), and the X-drops that NCBI's gapped kernels can size
//! (`kernel_x_drop_fits`). BLASTP, BLASTX, TBLASTN and TBLASTX convert their X-drops here, once
//! per search or query batch, never in an extension loop.
//!
//! NCBI converts each double to `Int4` once. A value past `Int4` is converted as the x86-64 build
//! of NCBI BLAST+ 2.17.0 (the oracle) converts it, to `INT_MIN`
//! (`blast_util::ncbi_int4_from_double`): every gapped kernel then raises it to the gap costs
//! (`x_dropoff < gap_open_extend`), and every ungapped extension stops at the first letter
//! (`(maxscore - score) >= dropoff`), with no further arithmetic on it. NCBI 2.17.0 gives these
//! outputs deterministically (line RG X1 probe: `-xdrop_ungap 1e9` to `1e300`, `-xdrop_gap 1e9`
//! and `1e300`, the redo of `-xdrop_gap_final 1e8`), so they are reproduced. Only the sizes of
//! the gapped kernels' work areas, `Int4` arithmetic that overflows (undefined behaviour), are
//! refused (`kernel_x_drop_fits`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:426-429 (the same at
//! 783-786 in `Blast_SemiGappedAlign`, 1032-1035, 1387-1390 and 2003-2006)
//! ```c
//!     x_dropoff = gap_align->gap_x_dropoff;
//!
//!     if (x_dropoff < gap_open_extend)
//!         x_dropoff = gap_open_extend;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:855-860 (the same at 911-916,
//! 960-965 and 1007-1012)
//! ```c
//!     for (i = q_off; i < query->length; i++) {
//!         score += matrix[q[i]][s[i]];
//!
//!         if (score > maxscore) {
//! ...
//!         if (score <= 0 || (maxscore - score) >= dropoff)
//!             break;
//! ```

use crate::core::blast_util::ncbi_int4_from_double;

/// NCBI's `NCBIMATH_LN2` (`0.69314718055994530941723212145818`, the same double).
const NCBIMATH_LN2: f64 = std::f64::consts::LN_2;

/// A context's ungapped X-drop (`x_dropoff_init`) from the X-drop in bits, the context's ungapped
/// Lambda (`sbp->kbp[context]`) and the score block's scale factor: one conversion of the double
/// product, `INT_MIN` past `Int4`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:218-221
/// ```c
///         p->cutoffs[context].x_dropoff_init =
///             (Int4)(sbp->scale_factor *
///                    ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
/// ```
#[inline]
pub(crate) fn ungapped_x_dropoff(bits: f64, lambda: f64, scale_factor: f64) -> i32 {
    ncbi_int4_from_double(scale_factor * (bits * NCBIMATH_LN2 / lambda).ceil())
}

/// The preliminary and final gapped X-drops (`gap_x_dropoff`, `gap_x_dropoff_final`) from the
/// X-drops in bits and the smallest gapped Lambda of the query batch: the final X-drop is the
/// larger of its own double and the preliminary raw X-drop, then converted.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:454-463
/// ```c
///    if (sbp->kbp_gap) {
///       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
///       params->gap_x_dropoff = (Int4)
///           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
///       /* Note that this conversion from bits to raw score is done prematurely
///          when rescaling and composition based statistics is applied, as we
///          lose precision. Therefore this is redone in Kappa_RedoAlignmentCore */
///       params->gap_x_dropoff_final = (Int4)
///           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
///    }
/// ```
#[inline]
pub(crate) fn gapped_x_dropoffs(gap_bits: f64, final_bits: f64, min_lambda: f64) -> (i32, i32) {
    let gapped = ncbi_int4_from_double(gap_bits * NCBIMATH_LN2 / min_lambda);
    let finale =
        ncbi_int4_from_double((final_bits * NCBIMATH_LN2 / min_lambda).max(f64::from(gapped)));
    (gapped, finale)
}

/// The X-drop of the composition-based redo (`gapping_params->x_dropoff`) from the final X-drop in
/// bits, the smallest gapped Lambda of the scaled score block and the search's final raw X-drop.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2374-2382
/// ```c
///     for (i = 0;  i < num_queries;  i++) {
///         if (context->sbp->kbp_gap[i] != NULL &&
///             context->sbp->kbp_gap[i]->Lambda < min_lambda) {
///             min_lambda = context->sbp->kbp_gap[i]->Lambda;
///         }
///     }
///     gapping_params->x_dropoff = (Int4)
///         MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
///             extendParams->gap_x_dropoff_final);
/// ```
#[inline]
pub(crate) fn redo_x_dropoff(final_bits: f64, scaled_min_lambda: f64, extension_final: i32) -> i32 {
    ncbi_int4_from_double(
        (final_bits * NCBIMATH_LN2 / scaled_min_lambda).max(f64::from(extension_final)),
    )
}

/// The largest quotient of a raw X-drop by the gap extension that NCBI's gapped kernels keep in
/// `Int4`, plus one: `num_extra_cells = x_dropoff / gap_extend + 3`, the DP window
/// `dp_mem_alloc = num_extra_cells + 100`, and its growth `2 * dp_mem_alloc`, which leaves `Int4`
/// once `dp_mem_alloc` is 2^30 (at the first row that passes the window's first 97 cells).
pub(crate) const KERNEL_X_DROP_QUOTIENT_LIMIT: i64 = (1 << 30) - 103;

/// Whether NCBI's gapped kernels (`Blast_SemiGappedAlign`, the preliminary stage's, and
/// `ALIGN_EX`, the traceback's) keep their sizes in `Int4` for the raw X-drop `x_dropoff` with the
/// gap costs they run with: the X-drop raised to `gap_open + gap_extend`, over `gap_extend`
/// (`Int4` division), below `KERNEL_X_DROP_QUOTIENT_LIMIT`. LOSAT's kernels size their memory by
/// the sequences, but NCBI's overflow is undefined behaviour, so such an X-drop is refused. A
/// later growth of the window (`last_b_index + num_extra_cells + 100`, doubled) overflows only
/// for alignments of more than 2^29 letters of one sequence, the limit of NCBI's kernels at any
/// X-drop.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:457-464
/// ```c
///     if (gap_extend > 0)
///         num_extra_cells = x_dropoff / gap_extend + 3;
///     else
///         num_extra_cells = N + 3;
///
///     if (num_extra_cells > gap_align->dp_mem_alloc) {
///         gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
///                                       2 * gap_align->dp_mem_alloc);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:641-644 (the same at
/// 923-926 in `Blast_SemiGappedAlign`)
/// ```c
///         if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {
///
///             gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
///                                           2 * gap_align->dp_mem_alloc);
/// ```
pub(crate) fn kernel_x_drop_fits(x_dropoff: i32, gap_open: i32, gap_extend: i32) -> bool {
    let gap_open_extend = i64::from(gap_open) + i64::from(gap_extend);
    let x_dropoff = i64::from(x_dropoff).max(gap_open_extend);
    gap_extend <= 0 || x_dropoff / i64::from(gap_extend) < KERNEL_X_DROP_QUOTIENT_LIMIT
}

/// The refusal of an X-drop in bits (`option`, as the command line names it) whose raw X-drop
/// NCBI's gapped kernels cannot size (`kernel_x_drop_fits`) in `program`'s search.
pub(crate) fn kernel_x_drop_error(program: &str, option: &str, bits: f64) -> anyhow::Error {
    anyhow::anyhow!(
        "an {option} of {bits:?} bits, whose raw X-drop makes the sizes of NCBI BLAST+'s gapped alignment kernels leave 32-bit integers (undefined behaviour), is not supported by LOSAT's {program}"
    )
}

/// The check of the preliminary and final gapped X-drops (raw, `gapped_x_dropoffs`) against
/// NCBI's kernels with the search's gap costs, before the search: the preliminary stage runs
/// `Blast_SemiGappedAlign` with the preliminary X-drop, and the traceback without
/// composition-based statistics `ALIGN_EX` with the final one (the composition-based redo runs
/// with its scaled gap costs, under which every `Int4` X-drop fits).
pub(crate) fn check_gapped_x_dropoffs(
    program: &str,
    (gap_bits, final_bits): (f64, f64),
    (gapped, finale): (i32, i32),
    gap_open: i32,
    gap_extend: i32,
) -> anyhow::Result<()> {
    if !kernel_x_drop_fits(gapped, gap_open, gap_extend) {
        return Err(kernel_x_drop_error(program, "-xdrop_gap", gap_bits));
    }
    if !kernel_x_drop_fits(finale, gap_open, gap_extend) {
        return Err(kernel_x_drop_error(program, "-xdrop_gap_final", final_bits));
    }
    Ok(())
}

#[cfg(test)]
#[path = "x_drop_tests.rs"]
mod x_drop_tests;
