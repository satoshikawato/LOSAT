//! Sequence representation and manipulation utilities
//!
//! This module provides efficient sequence representations for BLAST operations,
//! including the 2-bit packed nucleotide format used by NCBI BLAST.

pub mod packed_nucleotide;

pub use packed_nucleotide::*;
