//! NCBI-style sum-statistics (even-gap) HSP linking for TBLASTX.
//!
//! Currently simplified - returns hits with their individual E-values unchanged.
//! TODO: Implement proper NCBI-style iterative chain extraction with set-level E-values.

use crate::stats::KarlinParams;

use super::chaining::ExtendedHit;

/// Apply sum-statistics linking.
/// 
/// Current implementation: Pass through with no modification (individual E-values preserved).
/// This is faster and avoids the complexity of NCBI's iterative chain extraction,
/// but may result in slightly different hit counts compared to NCBI BLAST.
pub fn apply_sum_stats_even_gap_linking(
    hits: Vec<ExtendedHit>,
    _params: &KarlinParams,
) -> Vec<ExtendedHit> {
    // Simply return hits unchanged - individual E-values are already computed
    hits
}
