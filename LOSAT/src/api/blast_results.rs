//! BLAST Results
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp
//!
//! This module provides result structures for BLAST searches
//! including Hit and HSP representations.
//! During migration, this re-exports from the common module.

pub use crate::common::Hit;
