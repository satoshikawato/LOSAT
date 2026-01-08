//! BLAST Search Engine
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c
//!
//! This module provides the main BLAST search engine that coordinates
//! seeding, extension, and HSP processing.
//!
//! During migration, this re-exports from the original locations.

// Re-export TBLASTX engine
pub use crate::algorithm::tblastx::blast_engine::*;
