//! BLAST API Layer
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/
//!
//! This module provides high-level API for running BLAST searches.
//! It wraps the core algorithms and provides a user-friendly interface.
//!
//! # Structure
//!
//! - `local_blast` - Main entry point for local BLAST searches
//! - `blast_nucl_options` - BLASTN-specific options
//! - `tblastx_options` - TBLASTX-specific options
//! - `blast_results` - Result structures (Hit, HSP)
//! - `blast_options_handle` - Options handling and validation

pub mod local_blast;
pub mod blast_nucl_options;
pub mod tblastx_options;
pub mod blast_results;
pub mod blast_options_handle;
