//! BLAST Gapped Alignment
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c
//!
//! This module provides gapped alignment algorithms including
//! dynamic programming and traceback generation.
//!
//! During migration, this re-exports from the original locations.

pub use crate::algorithm::tblastx::blast_gapalign::*;
