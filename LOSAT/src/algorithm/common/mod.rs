//! Common modules shared between the programs (no program module is used from here)
//!
//! This module contains shared functionality that is used by both nucleotide
//! and protein alignment algorithms, including:
//! - Diagnostic counters for tracking pipeline statistics
//! - The search of an NCBI BLAST database (`-db`): opening, reader errors, stored titles
//! - The subjects of the tabular reports with a field list (`report::tabular`)
//! - A `-db` search's database read where it is searched (`x_db_stream`)
//! - Mask helpers (`masks`)
//! - The diagonal table's entry of the protein word finders and its two-hit test (`diag_entry`)
//! - The walk of a protein sequence's diagonal table over its `-db_hard_mask` chunks (`hard_chunk_walk`)
//!   and the one word of a chunk shorter than the word, with LOSAT's refusals there
//! (`hard_chunk_scan`)
//! - The chunks of a translated subject frame and their scan ranges (`translated_chunks`)
//! - The context of an offset of a concatenated query without the binary search (`x_context_table`)

pub(crate) mod database;
pub(crate) mod db_mask_search;
pub(crate) mod diag_entry;
pub mod diagnostics;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276
// The hard chunks of a protein database sequence with `-db_hard_mask`, whose word finders BLASTP
// and BLASTX share (series NP, M4b-3b).
pub(crate) mod hard_chunk_scan;
pub(crate) mod hard_chunk_walk;
pub(crate) mod masks;
pub(crate) mod tabular;
pub(crate) mod translated_chunks;
pub(crate) mod unknown_matrix;
pub(crate) mod x_context_table;
pub(crate) mod x_db_stream;
