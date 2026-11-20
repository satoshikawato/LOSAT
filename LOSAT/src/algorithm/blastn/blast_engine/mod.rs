//! BLASTN search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c
//!            ncbi-blast/c++/src/algo/blast/core/blast_engine.c
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Seed finding with two-hit extension
//! - HSP filtering and output generation
//!
//! # Module Structure
//!
//! - `run` - Main run() function for BLASTN search

mod run;

// Re-export the main run function
pub use run::run;
pub use run::run_local;
pub use run::run_local_database;
pub(crate) use run::QueryBatchSizes;
// LOSAT_X_BNLUTDROP (graduated), EXPERIMENT (LOSAT_X_QGROUPPAR, LOSAT_X_BNQPACK, LOSAT_X_QPREPPAR,
// LOSAT_X_BNTBRANGELOG): the summary lines `main` prints at exit (no NCBI counterpart).
pub use run::x_lut_drop_print_summary;
pub use run::x_qgroup_par_print_summary;
pub use run::x_query_drop_print_summary;
pub use run::x_query_prep_par_print_summary;
pub use run::x_tb_range_print_summary;
// ABI v1 (`web_api::v1_blastn`) checks a search without subjects as the CLI does.
#[cfg(target_arch = "wasm32")]
pub(crate) use run::check_subjects_not_empty;
// The end of a query reading and NCBI's message for a sequence without letters, shared with
// TBLASTX (port plan step S5).
pub(crate) use crate::blastinput::query_batch::{QueryReading, NO_DATA_MESSAGE};
#[cfg(test)]
mod megablast_xdrop_tests;
