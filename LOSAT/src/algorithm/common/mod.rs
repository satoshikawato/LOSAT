//! Common modules shared between BLASTN and TBLASTX algorithms
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - E-value calculation functions
//! - HSP chaining and filtering logic

pub mod diagnostics;
pub mod evalue;
pub mod chaining;


