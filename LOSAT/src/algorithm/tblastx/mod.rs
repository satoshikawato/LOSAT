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
mod database;
pub mod diagnostics;
pub mod extension;
pub mod filtering;
pub mod hsp_culling;
pub mod lookup;
pub mod ncbi_cutoffs;
pub mod reevaluate;
pub(crate) mod report;
mod scan;
pub mod stage_dump;
pub(crate) mod subjects;
pub mod sum_stats_linking;
pub mod tracing;
pub mod translation;
mod ungapped_hit_sort;
// Diagnostic (LOSAT_TIMING): the stage times of the two linking calls of a subject. No NCBI
// counterpart; it computes no value the search reads.
mod x_link_timing;
// EXPERIMENT (LOSAT_X_TBXSCANBF): the lookup scan without a branch on the presence vector.
mod x_scan_bf;
pub mod x_seed_bucket;

pub use args::TblastxArgs;
pub use blast_engine::run;
