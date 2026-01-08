//! Nucleotide Lookup Table
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c
//!
//! This module re-exports the nucleotide lookup table implementation
//! from the blastn algorithm module during migration.

pub use crate::algorithm::blastn::lookup::*;
