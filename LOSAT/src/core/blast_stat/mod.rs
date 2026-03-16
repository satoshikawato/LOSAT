//! BLAST Statistical Functions
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c
//!
//! This module provides Karlin-Altschul statistical calculations for BLAST:
//! - Bit score and E-value calculation
//! - Length adjustment for effective search space
//! - Sum statistics for HSP linking
//! - Statistical parameter lookup tables
//!
//! # Organization
//!
//! The module is organized into submodules matching NCBI BLAST's blast_stat.c:
//! - `karlin_params` - Parameter definitions (KarlinParams)
//! - `search_space` - Effective search space (SearchSpace)
//! - `score_calc` - Score/E-value calculations
//! - `length_adjustment` - Length adjustment algorithm
//! - `sum_statistics` - Sum statistics for multi-HSP E-values
//! - `composition` - Karlin parameter calculation from composition
//! - `lookup_tables` - Parameter lookup tables

pub mod composition;
pub mod karlin_params;
pub mod length_adjustment;
pub mod lookup_tables;
pub mod score_calc;
pub mod search_space;
pub mod sum_statistics;

// Re-export main types and functions for convenience
pub use composition::{
    apply_check_ideal, compute_aa_composition, compute_blosum62_ideal_karlin_params,
    compute_karlin_params_ungapped, compute_score_freq_profile, compute_std_aa_composition,
    ScoreFreqProfile,
};
pub use karlin_params::KarlinParams;
pub use length_adjustment::{
    compute_length_adjustment_ncbi, compute_length_adjustment_simple, LengthAdjustmentResult,
};
pub use lookup_tables::{
    lookup_nucl_params, lookup_protein_params, lookup_protein_params_gapped,
    lookup_protein_params_ungapped, requires_even_scores,
};
pub use score_calc::{
    bit_score, calculate_statistics, evalue, evalue_from_raw_score, raw_score_from_bit_score,
    raw_score_from_evalue, raw_score_from_evalue_with_decay, simple_evalue,
};
pub use search_space::{compute_tblastx_search_space, SearchSpace};
pub use sum_statistics::{
    defaults, e_to_p, gap_decay_divisor, large_gap_sum_e, ln_factorial, ln_gamma_int,
    normalize_score, p_to_e, small_gap_sum_e, uneven_gap_sum_e,
};
