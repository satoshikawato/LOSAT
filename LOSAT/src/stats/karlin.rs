use super::tables::KarlinParams;
use crate::utils::glibc_fma_libm::GlibcLibm;

/// Calculate bit score from raw score using Karlin-Altschul statistics
///
/// Formula: S' = (lambda * S - ln(K)) / ln(2)
/// where S is the raw score, S' is the bit score
pub fn bit_score(raw_score: i32, params: &KarlinParams) -> f64 {
    let ln2 = 2.0_f64.glibc_ln();
    (params.lambda * (raw_score as f64) - params.k.glibc_ln()) / ln2
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
/// `BLAST_KarlinStoE_simple`: the e-value of one alignment of raw score `raw_score` in the
/// effective search space `search_space`, from the block's Lambda and `log_k`, its `logK`; -1
/// for a block with a negative Lambda, K or H. The one owner of this e-value: every caller of
/// NCBI's function (`Blast_HSPListGetEvalues` without a Gumbel block, the first e-values of the
/// linking) calls this one or `evalue_from_raw_score`.
///
/// A negative value comes from a PSSM's file (`-in_pssm`), an RPS database's K, or the failed
/// Karlin calculation of a PSSM's round (-1 in all three); NCBI's other searches leave no hit in
/// a context without valid statistics (its block is freed and the context masked).
#[inline]
pub fn evalue_from_raw_score_with_log_k(
    raw_score: i32,
    params: &KarlinParams,
    log_k: f64,
    search_space: f64,
) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 || params.h < 0.0 {
        return -1.0;
    }
    search_space * ((-params.lambda * (raw_score as f64)) + log_k).glibc_exp()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2722-2725
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2729-2733
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3619-3622
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3883-3899
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2876-2878
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2123-2128
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2750-2752
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1192-1193
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:97-121
/// `BLAST_KarlinStoE_simple` (`evalue_from_raw_score_with_log_k`) with `logK` taken as the log
/// of the block's K. NCBI sets `logK` to `log(K)` wherever it sets K (the ungapped calculation,
/// the tables, a copy, the rescaling of composition-based statistics, a PSSM's blocks, an RPS
/// database's K); its one other value, `HUGE_VAL` beside the -1s of a failed calculation, is
/// behind the check.
///
/// No NCBI counterpart (code placement only): out of line, so that its check does not change the
/// register allocation of its callers' loops (inlined, the linking loop of TBLASTN's
/// `x_link_values` reloaded two values from the stack, +3 % of its instructions in a search
/// with a Gumbel block, which never calls this function).
#[inline(never)]
pub fn evalue_from_raw_score(raw_score: i32, params: &KarlinParams, search_space: f64) -> f64 {
    evalue_from_raw_score_with_log_k(raw_score, params, params.k.glibc_ln(), search_space)
}

#[cfg(test)]
#[path = "karlin_tests.rs"]
mod tests;
