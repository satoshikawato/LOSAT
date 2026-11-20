//! Reading a stored BLAST request (a search strategy of `-import_search_strategy`, the request of
//! a BLAST archive): what NCBI makes of its parts.
//!
//! - `options`: the request's options lists as the program's options (`CBlastOptionsBuilder`);
//! - `queries`: the request's query Bioseqs as the FASTA text NCBI writes and reads back
//!   (`s_ImportQueries`, `CFastaOstream`).

pub mod options;
pub mod queries;
