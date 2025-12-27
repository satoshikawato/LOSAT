//! Diagnostic counters for TBLASTX pipeline
//!
//! This module re-exports diagnostic functionality from the common module
//! for backward compatibility.

pub use crate::algorithm::common::diagnostics::{
    diagnostics_enabled, DiagnosticCounters, ProteinDiagnosticCounters,
};

use super::constants::MIN_UNGAPPED_SCORE;

/// Print diagnostic summary for TBLASTX (with TBLASTX-specific constants)
pub fn print_summary(counters: &DiagnosticCounters) {
    counters.print_summary(MIN_UNGAPPED_SCORE);
}

