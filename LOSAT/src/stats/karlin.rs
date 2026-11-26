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

/// Calculate E-value from raw score, Karlin parameters, precomputed `log(K)`, and search space
/// (NCBI BLAST_KarlinStoE_simple compatible).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
///
/// Formula: E = searchsp * exp(-Lambda * S + logK)
/// where S is the raw score, `log_k` is the natural logarithm of K, and searchsp is the effective
/// search space. Returns -1.0 when Lambda, K, or H is negative.
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

/// Calculate E-value from raw score and search space (NCBI BLAST_KarlinStoE_simple compatible)
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
///
/// Formula: E = searchsp * exp(-Lambda * S + logK)
/// where S is the raw score and searchsp is the effective search space. Returns -1.0 when
/// Lambda, K, or H is negative. Every valid NCBI Karlin block stores the natural logarithm of K
/// in `logK`, so `params.k.glibc_ln()` matches `kbp->logK` bit-for-bit.
#[inline]
pub fn evalue_from_raw_score(raw_score: i32, params: &KarlinParams, search_space: f64) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 || params.h < 0.0 {
        return -1.0;
    }
    evalue_from_raw_score_with_log_k(raw_score, params, params.k.glibc_ln(), search_space)
}

#[cfg(test)]
#[path = "karlin_tests.rs"]
mod tests;
