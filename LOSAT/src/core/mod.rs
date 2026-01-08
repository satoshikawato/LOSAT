//! BLAST Core Algorithms
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/
//!
//! This module contains all core BLAST algorithm implementations,
//! independent of program type (BLASTN, TBLASTX, etc).
//!
//! # Structure
//!
//! The module is organized to match NCBI BLAST's core/ directory:
//!
//! - **Statistics** (`blast_stat`, `blast_parameters`)
//!   - Karlin-Altschul statistics, E-value calculation
//!   - Length adjustment, search space
//!   - Cutoff score calculations
//!
//! - **Lookup Tables** (`blast_lookup`, `blast_aalookup`, `blast_nalookup`)
//!   - Lookup table construction and management
//!   - Word finding for seeds
//!
//! - **Extension** (`blast_extend`, `aa_ungapped`, `na_ungapped`)
//!   - Ungapped extension algorithms
//!   - Two-hit extension logic
//!
//! - **Scanning** (`blast_aascan`, `blast_nascan`)
//!   - Sequence scanning for seeds
//!
//! - **Gapped Alignment** (`blast_gapalign`, `blast_sw`, `blast_traceback`, `greedy_align`)
//!   - Dynamic programming alignment
//!   - Smith-Waterman implementation
//!   - Traceback generation
//!
//! - **HSP Management** (`blast_hits`, `link_hsps`, `blast_itree`, `hspfilter_*`)
//!   - HSP storage and sorting
//!   - Sum-statistics linking
//!   - HSP filtering and culling
//!
//! - **Utilities** (`blast_util`, `blast_encoding`, `gencode_singleton`, etc.)
//!   - General utility functions
//!   - Sequence encoding
//!   - Genetic code tables

// Statistics
pub mod blast_stat;

// Utilities
pub mod blast_util;
pub mod blast_encoding;
pub mod gencode_singleton;

// Filtering/Masking
pub mod blast_filter;
mod blast_seg_lnfact;  // Internal - ln(n!) lookup table for SEG
pub mod blast_seg;

// Lookup Tables
pub mod blast_lookup;
pub mod blast_aalookup;
pub mod blast_nalookup;

// Extension
pub mod blast_extend;
pub mod aa_ungapped;
pub mod na_ungapped;

// Scanning
pub mod blast_aascan;

// Gapped Alignment
pub mod blast_gapalign;
pub mod greedy_align;

// HSP Management
pub mod blast_hits;
pub mod link_hsps;
pub mod hspfilter_culling;

// Parameters and Options
pub mod blast_parameters;
pub mod blast_options;

// Engine
pub mod blast_engine;

// Diagnostics
pub mod blast_diagnostics;

// TODO: Additional modules to be added as migration progresses:
// pub mod blast_nascan;
// pub mod blast_sw;
// pub mod blast_traceback;
// pub mod blast_itree;
// pub mod hspfilter_collector;
// pub mod blast_program;
// pub mod ncbi_math;
// pub mod blast_setup;
