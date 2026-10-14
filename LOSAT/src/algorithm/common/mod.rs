//! Common modules shared between BLASTN and TBLASTX algorithms
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - E-value calculation functions
//! - HSP chaining and filtering logic
//! - The search of an NCBI BLAST database (`-db`): opening, reader errors, stored titles

pub mod chaining;
pub(crate) mod database;
pub mod diagnostics;
pub mod evalue;
pub mod protein_query_split;
