//! Common modules shared between the programs (no program module is used from here)
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - The search of an NCBI BLAST database (`-db`): opening, reader errors, stored titles
//! - The subjects of the tabular reports with a field list (`report::tabular`)
//! - A `-db` search's database read where it is searched (`x_db_stream`)
//! - Mask helpers (`masks`)
//! - The chunks of a translated subject frame and their scan ranges (`translated_chunks`)
//! - The context of an offset of a concatenated query without the binary search (`x_context_table`)
//! - The two-hit walks on keys of score and step, of blastp, BLASTX and TBLASTX (`x_key_walk`)

pub(crate) mod database;
pub(crate) mod db_mask_search;
pub mod diagnostics;
pub(crate) mod masks;
pub(crate) mod tabular;
pub(crate) mod translated_chunks;
pub(crate) mod unknown_matrix;
pub(crate) mod x_context_table;
pub(crate) mod x_db_stream;
pub(crate) mod x_key_walk;
