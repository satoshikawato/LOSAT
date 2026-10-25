//! BLAST Core Algorithms
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/
//!
//! Ports of NCBI core/ (and composition_adjustment/) code that several programs share. The
//! lookup, scan, extension, gapped alignment and linking ports are per program, under
//! `algorithm::<program>`.
//!
//! # Structure
//!
//! - **Statistics** (`blast_stat`): Karlin-Altschul statistics, E-values, length
//!   adjustment, search space
//! - **Utilities** (`blast_util`, `blast_encoding`): sequence encoding and helpers
//! - **Filtering/Masking** (`blast_seg`, a re-export of `utils::seg`)
//! - **Composition adjustment** (`composition_adjustment`, NCBI algo/blast/composition_adjustment/)
//! - **HSP filters** (`blast_hits`: `-max_hsps`, `-qcov_hsp_perc`; `hspfilter_culling`:
//!   `-culling_limit`)

// Statistics
pub mod blast_stat;

// Utilities
pub mod blast_encoding;
pub mod blast_util;

// Filtering/Masking
pub mod blast_seg;

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-190
// ```c
// typedef struct BlastCompo_Alignment { ... } BlastCompo_Alignment;
// typedef struct BlastCompo_GappingParams { ... } BlastCompo_GappingParams;
// typedef struct BlastCompo_SequenceRange { ... } BlastCompo_SequenceRange;
// ```
pub mod composition_adjustment;

// HSP Management
pub mod blast_hits;
pub mod hspfilter_culling;
