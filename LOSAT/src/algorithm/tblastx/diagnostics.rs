//! Diagnostic counters for TBLASTX pipeline
//!
//! This module re-exports diagnostic functionality from the common module
//! for backward compatibility.

pub use crate::algorithm::common::diagnostics::{
    diagnostics_enabled, DiagnosticCounters, ProteinDiagnosticCounters,
};

/// Print diagnostic summary for TBLASTX
/// 
/// Note: NCBI BLAST does not use a fixed MIN_UNGAPPED_SCORE threshold for TBLASTX.
/// Instead, it uses cutoffs->cutoff_score which is dynamically calculated from E-value.
/// Since we removed the fixed threshold to match NCBI BLAST behavior, the diagnostic
/// will show 0 for "Filtered (low score)" since no such filtering occurs.
pub fn print_summary(counters: &DiagnosticCounters) {
    // Pass 0 since we don't use a fixed threshold (matches NCBI BLAST behavior)
    counters.print_summary(0);
}

