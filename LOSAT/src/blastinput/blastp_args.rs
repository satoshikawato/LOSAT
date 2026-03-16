//! BLASTP Arguments
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp
//!
//! This module provides BLASTP-specific argument handling.
//! During migration, this re-exports from the algorithm modules.

pub use crate::algorithm::blastp::args::*;
