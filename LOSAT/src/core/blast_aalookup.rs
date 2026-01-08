//! Amino Acid Lookup Table
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//!
//! This module re-exports the amino acid lookup table implementation
//! from the tblastx algorithm module during migration.

pub use crate::algorithm::tblastx::lookup::*;
