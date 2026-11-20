//! The reader of a stored Blast4 request for `-import_search_strategy` (`CImportStrategy`): the
//! request as the typed model of `request`, read by the generic reader of the type table
//! (`report::serial::asn_types`).
//!
//! What NCBI makes of the parts of a stored request (a search strategy of
//! `-import_search_strategy`):
//!
//! - `options`: the request's options lists as the program's options (`CBlastOptionsBuilder`);
//! - `queries`: the request's query Bioseqs as the FASTA text NCBI writes and reads back
//!   (`s_ImportQueries`, `CFastaOstream`).

pub(crate) mod import;
pub mod options;
pub mod queries;
pub(crate) mod request;

#[cfg(test)]
mod strategy_tests;
