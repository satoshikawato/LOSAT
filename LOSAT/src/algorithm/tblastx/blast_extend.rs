//! Diagonal table management for TBLASTX
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-180
//!
//! The entry of the diagonal table (`DiagStruct`), which the protein programs share, lives in
//! `algorithm::common::diag_entry`.

pub use crate::algorithm::common::diag_entry::DiagStruct;
