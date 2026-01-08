//! BLAST Extension Infrastructure
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c
//!
//! This module provides common extension infrastructure including
//! diagonal tracking structures (DiagStruct) used for two-hit detection.
//!
//! During migration, this re-exports from the original locations.

// Common extension infrastructure will be consolidated here
// For now, re-export from both AA and NA extension modules

pub use crate::algorithm::tblastx::blast_extend::*;
