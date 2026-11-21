//! BLAST Core Algorithms
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/
//!
//! Ports of NCBI core/ (and composition_adjustment/) code that several programs share. The
//! lookup, scan, extension, gapped alignment and linking ports are per program, under
//! `algorithm::<program>`.
//!
//! # Structure
//!
//! - **Utilities** (`blast_util`, `blast_encoding`): sequence encoding and helpers
//! - **Filtering/Masking** (`blast_seg`, a re-export of `utils::seg`)
//! - **Composition adjustment** (`composition_adjustment`, NCBI algo/blast/composition_adjustment/)
//! - **HSP filters** (`blast_hits`: `-max_hsps`, `-qcov_hsp_perc`; `hspfilter_besthit`:
//!   `-best_hit_overhang`, `-best_hit_score_edge`; `hspfilter_culling`:
//!   `-culling_limit`)
//! - **HSP stream** (`blast_hspstream`): the preliminary stage's close with `ADAPTIVE_CBS`
//! - **Subject chunks with hard ranges** (`subject_split`, NCBI core/blast_engine.c's
//!   `s_GetNextSubjectChunk`): BLASTN's `-db_hard_mask` search and its refusals
//! - **PSI-BLAST's PSSM engine** (`psi`, NCBI core/blast_psi.c and blast_psi_priv.c): psiblast's
//!   `-out_ascii_pssm` and `-out_pssm`

// Utilities
pub mod blast_encoding;
pub mod blast_util;
/// The even-gap linking of a gapped search's HSPs with sum statistics (psiblast's `-sum_stats true`).
pub(crate) mod link_hsps;
/// The X-drops in raw scores of the protein programs (`-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`).
pub(crate) mod x_drop;

// Filtering/Masking
pub mod blast_seg;

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/redo_alignment.h:49-190
pub mod composition_adjustment;

// HSP Management
pub mod blast_hits;
// The close of the preliminary stage's HSP stream with `ADAPTIVE_CBS` (NCBI core/blast_hspstream.c)
pub(crate) mod blast_hspstream;
pub(crate) mod hspfilter_besthit;
pub mod hspfilter_culling;

// The chunks of a database sequence with hard ranges (NCBI core/blast_engine.c's
// `s_GetNextSubjectChunk`): BLASTN's `-db_hard_mask` search and its refusals
pub mod subject_split;

// PSI-BLAST's PSSM engine (NCBI core/blast_psi.c, blast_psi_priv.c)
pub(crate) mod psi;
