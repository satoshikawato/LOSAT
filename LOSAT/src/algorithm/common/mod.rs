//! Common modules shared between BLASTN and TBLASTX algorithms
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - E-value calculation functions
//! - HSP chaining and filtering logic
//! - The search of an NCBI BLAST database (`-db`): opening, reader errors, stored titles
//! - The subjects of the tabular reports with a field list (`report::tabular`)
//! - A `-db` search's database read where it is searched (`x_db_stream`)

pub mod chaining;
pub(crate) mod database;
pub mod diagnostics;
pub mod evalue;
pub mod protein_query_split;
pub(crate) mod tabular;
pub(crate) mod x_db_stream;
