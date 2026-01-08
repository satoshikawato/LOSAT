//! Amino Acid Sequence Scanning
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c
//!
//! This module provides amino acid sequence scanning algorithms
//! for finding seed matches in TBLASTX.
//!
//! During migration, this re-exports from the original location.

pub use crate::algorithm::tblastx::blast_aascan::*;
