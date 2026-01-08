//! HSP Linking (Sum Statistics)
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c
//!
//! This module provides HSP linking algorithms for computing
//! sum statistics (small-gap and large-gap E-values).
//!
//! During migration, this re-exports from the original location.

pub use crate::algorithm::tblastx::sum_stats_linking::*;
