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

pub mod karlin_params;
pub mod search_space;
pub mod score_calc;
pub mod length_adjustment;
pub mod sum_statistics;
pub mod composition;
pub mod lookup_tables;

// Re-export main types and functions for convenience
pub use karlin_params::KarlinParams;
pub use search_space::{SearchSpace, compute_tblastx_search_space};
pub use score_calc::{
    bit_score, evalue, calculate_statistics,
    raw_score_from_evalue, raw_score_from_evalue_with_decay,
    raw_score_from_bit_score, evalue_from_raw_score, simple_evalue,
};
pub use length_adjustment::{
    LengthAdjustmentResult,
    compute_length_adjustment_ncbi, compute_length_adjustment_simple,
};
pub use sum_statistics::{
    gap_decay_divisor, ln_factorial, ln_gamma_int, p_to_e, e_to_p,
    small_gap_sum_e, uneven_gap_sum_e, large_gap_sum_e, normalize_score,
    defaults,
};
pub use composition::{
    compute_aa_composition, compute_std_aa_composition,
    ScoreFreqProfile, compute_score_freq_profile,
    compute_karlin_params_ungapped, apply_check_ideal,
};
pub use lookup_tables::{
    lookup_nucl_params, lookup_protein_params,
    requires_even_scores, lookup_protein_params_ungapped, lookup_protein_params_gapped,
};
