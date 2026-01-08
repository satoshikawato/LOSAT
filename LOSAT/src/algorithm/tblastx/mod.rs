//! TBLASTX algorithm module
//!
//! This module implements TBLASTX (translated DNA vs translated DNA search),
//! which translates both query and subject sequences in all six reading frames
//! and performs protein-protein alignment.

pub mod args;
pub mod constants;
pub mod translation;
pub mod lookup;
pub mod extension;
pub mod chaining;
pub mod sum_stats_linking;
pub mod reevaluate;
pub mod diagnostics;
pub mod ncbi_cutoffs;
pub mod hsp_culling;
pub mod tracing;
pub mod filtering;
mod scan;
pub mod blast_extend;
pub mod blast_gapalign;
pub mod blast_aascan;
pub mod blast_engine;

pub use args::TblastxArgs;
pub use blast_engine::run;


