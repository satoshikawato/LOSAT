//! BLAST Statistical Functions
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c
//!
//! This is the old statistics layer. The live, NCBI-faithful layer is
//! `crate::stats`; what remains here is what blastp, blastx and tblastn
//! composition-based statistics still call.
//!
//! # Organization
//!
//! - `karlin_params` - Parameter definitions (KarlinParams)
//! - `sum_statistics` - `p_to_e` / `e_to_p` (libm versions, see the module doc)
//! - `composition` - Karlin parameter calculation from composition

pub mod composition;
pub mod karlin_params;
pub mod sum_statistics;

// Re-export main types and functions for convenience
pub use composition::{
    apply_check_ideal, compute_aa_composition, compute_blosum62_ideal_karlin_params,
    compute_karlin_params_ungapped, compute_score_freq_profile, compute_std_aa_composition,
    ScoreFreqProfile,
};
pub use karlin_params::KarlinParams;
pub use sum_statistics::{e_to_p, p_to_e};
