//! BLAST API Layer
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/
//!
//! The local search entry shared by the command line and both web ABIs: each program's run
//! function and the routing of its report outputs.
//!
//! # Structure
//!
//! - `local_blast` - Main entry point for local BLAST searches

pub mod local_blast;
