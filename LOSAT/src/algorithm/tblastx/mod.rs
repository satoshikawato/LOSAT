//! TBLASTX algorithm module
//!
//! This module implements TBLASTX (translated DNA vs translated DNA search),
//! which translates both query and subject sequences in all six reading frames
//! and performs protein-protein alignment.

pub mod args;
pub mod blast_aascan;
pub mod blast_engine;
pub mod blast_extend;
pub mod blast_gapalign;
pub mod chaining;
pub mod constants;
pub mod diagnostics;
pub mod extension;
pub mod filtering;
pub mod hsp_culling;
pub mod lookup;
pub mod ncbi_cutoffs;
pub mod reevaluate;
mod scan;
pub mod sum_stats_linking;
pub mod tracing;
pub mod translation;

pub use args::TblastxArgs;
pub use blast_engine::run;
#[cfg(target_arch = "wasm32")]
pub use blast_engine::run_web_pair;
