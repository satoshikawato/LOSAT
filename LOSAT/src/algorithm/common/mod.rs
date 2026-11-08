//! Common modules shared between the programs (no program module is used from here)
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - E-value calculation functions
//! - The search of an NCBI BLAST database (`-db`): opening, reader errors, stored titles
//! - The subjects of the tabular reports with a field list (`report::tabular`)
//! - A `-db` search's database read where it is searched (`x_db_stream`)
//! - Mask helpers (`masks`)
//! - The chunks of a translated subject frame and their scan ranges (`translated_chunks`)

pub(crate) mod database;
pub mod diagnostics;
pub mod evalue;
pub(crate) mod masks;
pub(crate) mod tabular;
pub(crate) mod translated_chunks;
pub(crate) mod x_db_stream;

#[cfg(test)]
mod helpers_tests;
