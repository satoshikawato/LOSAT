//! BLASTP algorithm module
//!
//! This module implements BLASTP (protein vs protein search) in pure Rust,
//! using NCBI BLAST as the behavioral reference.

mod alignment;
pub mod args;
pub mod blast_engine;
mod encoding;
mod extension;
mod gapalign;
mod hsp;
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:102-3691
// ```c
// static void s_HSPListNormalizeScores(...);
// static void s_HitlistReapContained(...);
// static int s_HitlistEvaluateAndPurge(...);
// static void s_ComputeNumIdentities(...);
// ```
mod kappa;

pub use args::BlastpArgs;
pub use blast_engine::run;
