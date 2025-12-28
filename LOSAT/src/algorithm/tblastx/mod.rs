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
pub mod diagnostics;
pub mod utils;

pub use args::TblastxArgs;
pub use utils::run;


