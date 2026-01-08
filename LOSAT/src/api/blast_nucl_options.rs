//! BLASTN Options
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp
//!
//! This module provides BLASTN-specific options and configuration.
//! During migration, this re-exports from the algorithm modules.

pub use crate::algorithm::blastn::args::*;
pub use crate::algorithm::blastn::constants::*;
