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

/// Calculate E-value from raw score and search space (NCBI BLAST_KarlinStoE_simple compatible)
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
///
/// Formula: E = searchsp * exp(-Lambda * S + logK)
/// where S is the raw score, searchsp is the effective search space
pub fn evalue_from_raw_score(raw_score: i32, params: &KarlinParams, search_space: f64) -> f64 {
    if params.lambda < 0.0 || params.k < 0.0 {
        return -1.0;
    }
    search_space * ((-params.lambda * (raw_score as f64)) + params.k.glibc_ln()).glibc_exp()
}

#[cfg(test)]
#[path = "karlin_tests.rs"]
mod tests;
