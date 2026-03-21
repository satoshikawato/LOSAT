//! Local BLAST Search
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp
//!
//! This module provides the main entry point for running local BLAST searches.
//! During migration, this re-exports from the algorithm modules.

// Re-export TBLASTX run function
pub use crate::algorithm::tblastx::blast_engine::run as run_tblastx;

// Re-export BLASTN run function
pub use crate::algorithm::blastn::blast_engine::run as run_blastn;

// Re-export BLASTP run function
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:1196-1218
// ```c
// CLocalBlast::Run()
// {
//     x_SetupSearch();
//     x_RunSearch();
// }
// ```
pub use crate::algorithm::blastp::blast_engine::run as run_blastp;
