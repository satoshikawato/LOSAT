//! BLASTN Arguments
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp
//!
//! This module provides BLASTN-specific argument handling.
//! During migration, this re-exports from the algorithm modules.

pub use crate::algorithm::blastn::args::*;
