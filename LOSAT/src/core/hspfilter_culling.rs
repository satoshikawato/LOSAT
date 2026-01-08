//! HSP Culling Filter
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
//!
//! This module provides HSP culling algorithms using interval trees
//! to remove redundant HSPs.
//!
//! During migration, this re-exports from the original location.

pub use crate::algorithm::tblastx::hsp_culling::*;
