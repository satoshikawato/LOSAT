//! BLAST Input Handling
//!
//! Reference: ncbi-blast/c++/src/algo/blast/blastinput/
//!
//! This module provides input handling for BLAST searches
//! including argument parsing and FASTA input processing.
//!
//! # Structure
//!
//! - `blast_args` - Common argument handling
//! - `blastn_args` - BLASTN-specific arguments
//! - `tblastx_args` - TBLASTX-specific arguments

pub mod blast_args;
pub mod blastn_args;
pub mod tblastx_args;
