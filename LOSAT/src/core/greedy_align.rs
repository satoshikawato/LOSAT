//! Greedy Alignment Algorithm
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/greedy_align.c
//!
//! This module provides the greedy alignment algorithm used in
//! megablast for efficient nucleotide alignment.
//!
//! During migration, this re-exports from the original location.

pub use crate::algorithm::blastn::alignment::greedy::*;
