//! Common BLAST Arguments
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp
//!
//! This module provides common argument handling shared between
//! different BLAST programs.

// Common constants and utilities for argument handling
pub const DEFAULT_EVALUE: f64 = 10.0;
pub const DEFAULT_NUM_THREADS: usize = 1;
pub const DEFAULT_OUTFMT: u32 = 6;
